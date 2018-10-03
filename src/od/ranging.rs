extern crate dual_num;
extern crate hifitime;
extern crate nalgebra as na;
extern crate rand;

use self::dual_num::linalg::norm;
use self::dual_num::{partials, Dual};
use self::hifitime::instant::Instant;
use self::hifitime::julian::*;
use self::na::{Matrix2x6, Matrix6, U2, U3, U6, Vector2, Vector3, Vector6};
use self::rand::distributions::Normal;
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

        // Convert to SEZ to compute elevation
        let rho_sez = r2(PI / 2.0 - self.latitude.to_radians()) * r3(self.longitude.to_radians()) * rho_ecef;
        let elevation = (rho_sez[(2, 0)] / rho_ecef.norm()).asin().to_degrees();

        StdMeasurement::new(dt, tx_ecef, rx_ecef, elevation >= self.elevation_mask)
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

    // This is an example of the sensitivity matrix (H tilde) of a ranging method.
    pub fn compute_dual_sensitivity(state: &Matrix6<Dual<f64>>) -> Matrix2x6<Dual<f64>> {
        let range_mat = state.fixed_slice::<U3, U6>(0, 0).into_owned();
        let velocity_mat = state.fixed_slice::<U3, U6>(3, 0).into_owned();
        let mut range_slice = Vec::with_capacity(6);
        let mut range_rate_slice = Vec::with_capacity(6);

        for j in 0..6 {
            let rho_vec = Vector3::new(range_mat[(0, j)], range_mat[(1, j)], range_mat[(2, j)]);
            let range = norm(&rho_vec);
            let delta_v_vec = (Vector3::new(velocity_mat[(0, j)], velocity_mat[(1, j)], velocity_mat[(2, j)])) / range;
            let rho_dot = rho_vec.dot(&delta_v_vec);

            range_slice.push(range);
            range_rate_slice.push(rho_dot);
        }

        let mut rtn = Matrix2x6::zeros();

        rtn.set_row(0, &Vector6::from_row_slice(&range_slice).transpose());
        rtn.set_row(1, &Vector6::from_row_slice(&range_rate_slice).transpose());
        rtn
    }
}

impl Measurement for StdMeasurement {
    type StateSize = U6;
    type MeasurementSize = U2;

    fn new<F: CoordinateFrame>(dt: Instant, tx: State<F>, rx: State<F>, visible: bool) -> StdMeasurement {
        let (obs, h_tilde) = partials((rx - tx).to_cartesian_vec(), StdMeasurement::compute_dual_sensitivity);

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
