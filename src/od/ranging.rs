extern crate dual_num;
extern crate hifitime;
extern crate nalgebra as na;
extern crate rand;

use self::dual_num::linalg::norm;
use self::dual_num::{hyperspace_from_vector, Hyperdual};
use self::hifitime::instant::Instant;
use self::hifitime::julian::*;
use self::na::{DimName, Matrix1x6, Matrix2x6, Vector1, Vector2, VectorN, U1, U2, U3, U6, U7};
use self::rand::distributions::Normal;
use super::serde::ser::SerializeSeq;
use super::serde::{Serialize, Serializer};
use super::Measurement;
use celestia::{Frame, Geoid, State};
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
            40.427_222,
            4.250_556,
            0.834_939,
            range_noise,
            range_rate_noise,
        )
    }

    pub fn dss34_canberra(elevation_mask: f64, range_noise: f64, range_rate_noise: f64) -> GroundStation {
        GroundStation::from_noise_values(
            "Canberra",
            elevation_mask,
            -35.398_333,
            148.981_944,
            0.691_750,
            range_noise,
            range_rate_noise,
        )
    }

    pub fn dss13_goldstone(elevation_mask: f64, range_noise: f64, range_rate_noise: f64) -> GroundStation {
        GroundStation::from_noise_values(
            "Goldstone",
            elevation_mask,
            35.247_164,
            243.205,
            1.071_149_04,
            range_noise,
            range_rate_noise,
        )
    }

    /// Perform a measurement from the ground station to the receiver (rx).
    pub fn measure(self, rx: State<Geoid>, dt: Instant) -> StdMeasurement {
        use std::f64::consts::PI;
        // TODO: Get the frame from cosm instead of using the one from Rx!
        // TODO: Also change the frame number based on the axes, right now, ECI frame == ECEF!
        if rx.frame.id() != 3 {
            unimplemented!("the receiver is not around the Earth");
        }
        let mjd_dt = ModifiedJulian::from_instant(dt);
        let tx = State::from_geodesic(self.latitude, self.longitude, self.height, mjd_dt, rx.frame);
        /*
        // Convert the station to "ECEF"
        let theta = gast(dt);
        let tx_ecef_r = r3(-theta) * tx.radius();
        let tx_ecef_v = r3(-theta) * tx.velocity(); // XXX: Shouldn't we be using the transport theorem?
                                                    // HACK: change after Cosm supports rotations
        let tx_ecef = State::from_cartesian_vec(
            &Vector6::from_iterator(tx_ecef_r.iter().chain(tx_ecef_v.iter()).cloned()),
            mjd_dt,
            rx.frame,
        );
        */
        // Let's start by computing the range and range rate
        let rho_ecef = rx.radius() - tx.radius();

        // Convert to SEZ to compute elevation
        let rho_sez = r2(PI / 2.0 - self.latitude.to_radians()) * r3(self.longitude.to_radians()) * rho_ecef;
        let elevation = (rho_sez[(2, 0)] / rho_ecef.norm()).asin().to_degrees();

        StdMeasurement::new(dt, tx, rx, elevation >= self.elevation_mask)
    }
}

/// Computes the (approximate) Greenwich Apparent Sideral Time as per IAU2000.
///
/// NOTE: This is an approximation valid to within 0.9 seconds in absolute value.
/// In fact, hifitime does not support UT1, but according to the [IERS](https://www.iers.org/IERS/EN/Science/EarthRotation/UTC.html;jsessionid=A6E88EB4CF0FC2E1A3C10D807F51B829.live2?nn=12932)
/// UTC with leap seconds is always within 0.9 seconds to UT1, and hifitime inherently supports leap seconds.
/// Reference: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016
fn gast(at: Instant) -> f64 {
    use std::f64::consts::PI;
    let tu = ModifiedJulian::from_instant(at).days - 51_544.5;
    2.0 * PI * (0.779_057_273_264_0 + 1.002_737_811_911_354_6 * tu)
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

    fn compute_sensitivity(state: &VectorN<Hyperdual<f64, U7>, U6>) -> (Vector2<f64>, Matrix2x6<f64>) {
        // Extract data from hyperspace
        let range_vec = state.fixed_rows::<U3>(0).into_owned();
        let velocity_vec = state.fixed_rows::<U3>(3).into_owned();

        // Code up math as usual
        let delta_v_vec = velocity_vec / norm(&range_vec);
        let range = norm(&range_vec);
        let range_rate = range_vec.dot(&delta_v_vec);

        // Extract result into Vector2 and Matrix2x6
        let mut fx = Vector2::zeros();
        let mut pmat = Matrix2x6::zeros();
        for i in 0..U2::dim() {
            fx[i] = if i == 0 { range.real() } else { range_rate.real() };
            for j in 1..U7::dim() {
                pmat[(i, j - 1)] = if i == 0 { range[j] } else { range_rate[j] };
            }
        }
        (fx, pmat)
    }
}

impl Measurement for StdMeasurement {
    type StateSize = U6;
    type MeasurementSize = U2;

    fn new<F: Frame>(dt: Instant, tx: State<F>, rx: State<F>, visible: bool) -> StdMeasurement {
        let hyperstate = hyperspace_from_vector(&(rx - tx).to_cartesian_vec());
        let (obs, h_tilde) = Self::compute_sensitivity(&hyperstate);

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

/// Stores a standard measurement of range (km)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RangeMsr {
    dt: Instant,
    obs: Vector1<f64>,
    visible: bool,
    h_tilde: Matrix1x6<f64>,
}

impl RangeMsr {
    pub fn range(&self) -> f64 {
        self.obs[(0, 0)]
    }

    fn compute_sensitivity(state: &VectorN<Hyperdual<f64, U7>, U6>) -> (Vector1<f64>, Matrix1x6<f64>) {
        // Extract data from hyperspace
        let range_vec = state.fixed_rows::<U3>(0).into_owned();

        // Code up math as usual
        let range = norm(&range_vec);

        // Extract result into Vector2 and Matrix2x6
        let fx = Vector1::new(range.real());
        let mut pmat = Matrix1x6::zeros();

        for j in 1..U7::dim() {
            pmat[(j - 1)] = range[j];
        }

        (fx, pmat)
    }
}

impl Measurement for RangeMsr {
    type StateSize = U6;
    type MeasurementSize = U1;

    fn new<F: Frame>(_: Instant, tx: State<F>, rx: State<F>, visible: bool) -> RangeMsr {
        assert_eq!(tx.frame.id(), rx.frame.id(), "tx and rx in different frames");
        assert_eq!(tx.dt, rx.dt, "tx and rx states have different times");

        let dt = tx.dt;
        let hyperstate = hyperspace_from_vector(&(rx - tx).to_cartesian_vec());
        let (obs, h_tilde) = Self::compute_sensitivity(&hyperstate);

        RangeMsr {
            dt,
            obs: Vector1::new(obs[0]),
            visible,
            h_tilde: h_tilde.fixed_columns::<U6>(0).into_owned(),
        }
    }

    /// Returns this measurement as a vector of Range and Range Rate
    ///
    /// **Units:** km, km/s
    fn observation(&self) -> &Vector1<f64> {
        &self.obs
    }

    fn sensitivity(&self) -> &Matrix1x6<f64> {
        &self.h_tilde
    }

    fn visible(&self) -> bool {
        self.visible
    }

    fn at(&self) -> Instant {
        self.dt
    }
}

impl Serialize for RangeMsr {
    /// Serializes the observation vector at the given time.
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(3))?;
        seq.serialize_element(&ModifiedJulian::from_instant(self.dt).days)?;
        let obs = self.observation();
        seq.serialize_element(&obs[(0, 0)])?;
        seq.end()
    }
}

/// Stores a standard measurement of range (km)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DopplerMsr {
    dt: Instant,
    obs: Vector1<f64>,
    visible: bool,
    h_tilde: Matrix1x6<f64>,
}

impl DopplerMsr {
    pub fn range_rate(&self) -> f64 {
        self.obs[(0, 0)]
    }

    fn compute_sensitivity(state: &VectorN<Hyperdual<f64, U7>, U6>) -> (Vector1<f64>, Matrix1x6<f64>) {
        // Extract data from hyperspace
        let range_vec = state.fixed_rows::<U3>(0).into_owned();
        let velocity_vec = state.fixed_rows::<U3>(3).into_owned();

        // Code up math as usual
        let delta_v_vec = velocity_vec / norm(&range_vec);
        let range_rate = range_vec.dot(&delta_v_vec);

        // Extract result into Vector2 and Matrix2x6
        let fx = Vector1::new(range_rate.real());
        let mut pmat = Matrix1x6::zeros();

        for j in 1..U7::dim() {
            pmat[(0, j - 1)] = range_rate[j];
        }
        (fx, pmat)
    }
}

impl Measurement for DopplerMsr {
    type StateSize = U6;
    type MeasurementSize = U1;

    fn new<F: Frame>(_: Instant, tx: State<F>, rx: State<F>, visible: bool) -> DopplerMsr {
        assert_eq!(tx.frame.id(), rx.frame.id(), "tx and rx in different frames");
        assert_eq!(tx.dt, rx.dt, "tx and rx states have different times");

        let dt = tx.dt;
        let hyperstate = hyperspace_from_vector(&(rx - tx).to_cartesian_vec());
        let (obs, h_tilde) = Self::compute_sensitivity(&hyperstate);

        DopplerMsr {
            dt,
            obs: Vector1::new(obs[1]),
            visible,
            h_tilde: h_tilde.fixed_columns::<U6>(1).into_owned(),
        }
    }

    /// Returns this measurement as a vector of Range and Range Rate
    ///
    /// **Units:** km, km/s
    fn observation(&self) -> &Vector1<f64> {
        &self.obs
    }

    fn sensitivity(&self) -> &Matrix1x6<f64> {
        &self.h_tilde
    }

    fn visible(&self) -> bool {
        self.visible
    }

    fn at(&self) -> Instant {
        self.dt
    }
}

impl Serialize for DopplerMsr {
    /// Serializes the observation vector at the given time.
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(3))?;
        seq.serialize_element(&ModifiedJulian::from_instant(self.dt).days)?;
        let obs = self.observation();
        seq.serialize_element(&obs[(0, 0)])?;
        seq.end()
    }
}
