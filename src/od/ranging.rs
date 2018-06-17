extern crate hifitime;
extern crate nalgebra as na;
extern crate rand;

use self::hifitime::instant::Instant;
use self::hifitime::julian::*;
use self::na::{Matrix2x6, U2, U6, Vector2};
use self::rand::distributions::{Distribution, Normal};
use self::rand::thread_rng;
use super::serde::ser::SerializeSeq;
use super::serde::{Serialize, Serializer};
use super::Measurement;
use celestia::{CoordinateFrame, State, ECEF};
use utils::{r2, r3};

/// GroundStation defines a Two Way ranging equipment.
#[derive(Debug, Copy, Clone)]
pub struct GroundStation {
    pub name: &'static str,
    pub elevation_mask: f64,
    pub latitude: f64,
    pub longitude: f64,
    pub height: f64,
    range_noise: Normal,
    range_rate_noise: Normal,
}

impl GroundStation {
    /// Initializes a new Two Way ranging equipment from the noise values.
    pub fn from_noise_values(
        name: &'static str,
        elevation_mask: f64,
        latitude: f64,
        longitude: f64,
        height: f64,
        range_noise: f64,
        range_rate_noise: f64,
    ) -> GroundStation {
        GroundStation {
            name,
            elevation_mask,
            latitude,
            longitude,
            height,
            range_noise: Normal::new(0.0, range_noise),
            range_rate_noise: Normal::new(0.0, range_rate_noise),
        }
    }

    pub fn dss65_madrid(elevation_mask: f64, range_noise: f64, range_rate_noise: f64) -> GroundStation {
        GroundStation::from_noise_values(
            "Madrid",
            elevation_mask,
            40.427222,
            4.250556,
            0.834939,
            range_noise,
            range_rate_noise,
        )
    }

    pub fn dss34_canberra(elevation_mask: f64, range_noise: f64, range_rate_noise: f64) -> GroundStation {
        GroundStation::from_noise_values(
            "Canberra",
            elevation_mask,
            -35.398333,
            148.981944,
            0.691750,
            range_noise,
            range_rate_noise,
        )
    }

    pub fn dss13_goldstone(elevation_mask: f64, range_noise: f64, range_rate_noise: f64) -> GroundStation {
        GroundStation::from_noise_values(
            "Goldstone",
            elevation_mask,
            35.247164,
            243.205,
            1.07114904,
            range_noise,
            range_rate_noise,
        )
    }

    /// Perform a measurement from the ground station to the receiver (rx).
    pub fn measure<F: CoordinateFrame>(self, rx: State<F>, dt: Instant) -> StdMeasurement {
        use std::f64::consts::PI;
        let tx_ecef =
            State::<ECEF>::from_geodesic(self.latitude, self.longitude, self.height, ModifiedJulian::from_instant(dt));
        let rx_ecef = rx.in_frame(ECEF {});
        // Let's start by computing the range and range rate
        let rho_ecef = rx_ecef.radius() - tx_ecef.radius();
        let delta_v_ecef = (rx_ecef.velocity() - tx_ecef.velocity()) / rho_ecef.norm();
        let rho_dot = rho_ecef.dot(&delta_v_ecef);

        // Convert to SEZ to compute elevation
        let rho_sez = r2(PI / 2.0 - self.latitude.to_radians()) * r3(self.longitude.to_radians()) * rho_ecef;
        let elevation = (rho_sez[(2, 0)] / rho_ecef.norm()).asin().to_degrees();

        StdMeasurement::new(
            dt,
            tx_ecef,
            rx_ecef,
            Vector2::new(
                rho_ecef.norm() + self.range_noise.sample(&mut thread_rng()),
                rho_dot + self.range_rate_noise.sample(&mut thread_rng()),
            ),
            elevation >= self.elevation_mask,
        )
    }
}

/// Stores a standard measurement of range (km) and range rate (km/s)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StdMeasurement {
    dt: Instant,
    obs: Vector2<f64>,
    visible: bool,
    h_tilde: Matrix2x6<f64>,
}

impl StdMeasurement {
    pub fn range(&self) -> f64 {
        self.obs[(0, 0)]
    }
    pub fn range_rate(&self) -> f64 {
        self.obs[(1, 0)]
    }
}

impl Measurement for StdMeasurement {
    type StateSize = U6;
    type MeasurementSize = U2;

    fn new<F: CoordinateFrame>(dt: Instant, tx: State<F>, rx: State<F>, obs: Vector2<f64>, visible: bool) -> StdMeasurement {
        // Compute the H tilde matrix
        // Let's define some helper variables.
        let range = obs[(0, 0)];
        let range_rate = obs[(1, 0)];
        let mut h_tilde = Matrix2x6::zeros();
        // \partial \rho / \partial {x,y,z}
        h_tilde[(0, 0)] = (rx.x - tx.x) / range;
        h_tilde[(0, 1)] = (rx.y - tx.y) / range;
        h_tilde[(0, 2)] = (rx.z - tx.z) / range;
        // \partial \dot\rho / \partial {x,y,z}
        h_tilde[(1, 0)] = (rx.vx - tx.vx) / range + (range_rate / range.powi(2)) * (rx.x - tx.x);
        h_tilde[(1, 1)] = (rx.vy - tx.vy) / range + (range_rate / range.powi(2)) * (rx.y - tx.y);
        h_tilde[(1, 2)] = (rx.vz - tx.vz) / range + (range_rate / range.powi(2)) * (rx.z - tx.z);
        h_tilde[(1, 3)] = (rx.x - tx.x) / range;
        h_tilde[(1, 4)] = (rx.y - tx.y) / range;
        h_tilde[(1, 5)] = (rx.z - tx.z) / range;

        StdMeasurement {
            dt,
            obs,
            visible,
            h_tilde,
        }
    }

    /// Returns this measurement as a vector of Range and Range Rate
    ///
    /// **Units:** km, km/s
    fn observation(&self) -> &Vector2<f64> {
        &self.obs
    }

    fn sensitivity(&self) -> &Matrix2x6<f64> {
        &self.h_tilde
    }

    fn visible(&self) -> bool {
        self.visible
    }

    fn at(&self) -> Instant {
        self.dt
    }
}

impl Serialize for StdMeasurement {
    /// Serializes the observation vector at the given time.
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(3))?;
        seq.serialize_element(&ModifiedJulian::from_instant(self.dt).days)?;
        let obs = self.observation();
        seq.serialize_element(&obs[(0, 0)])?;
        seq.serialize_element(&obs[(1, 0)])?;
        seq.end()
    }
}
