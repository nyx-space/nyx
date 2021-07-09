/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

use super::hyperdual::linalg::norm;
use super::hyperdual::{hyperspace_from_vector, Hyperdual};
use super::rand::thread_rng;
use super::rand_distr::{Distribution, Normal};
use super::serde::ser::SerializeSeq;
use super::serde::{Serialize, Serializer};
use super::{Measurement, MeasurementDevice, TimeTagged};
use crate::cosmic::{Cosm, Frame, Orbit};
use crate::dimensions::{DimName, Matrix1x6, Matrix2x6, OVector, Vector1, Vector2, U1, U2, U6, U7};
use crate::time::Epoch;
use crate::utils::dcm_assemble;
use crate::Spacecraft;
use std::fmt;
use std::sync::Arc;

/// GroundStation defines a Two Way ranging equipment.
#[derive(Debug, Clone)]
pub struct GroundStation {
    pub name: String,
    /// in degrees
    pub elevation_mask: f64,
    /// in degrees
    pub latitude: f64,
    /// in degrees
    pub longitude: f64,
    /// in km
    pub height: f64,
    /// Frame in which this station is defined
    pub frame: Frame,
    pub cosm: Arc<Cosm>,
    range_noise: Normal<f64>,
    range_rate_noise: Normal<f64>,
}

impl GroundStation {
    /// Initializes a new Two Way ranging equipment from the noise values.
    pub fn from_noise_values(
        name: String,
        elevation_mask: f64,
        latitude: f64,
        longitude: f64,
        height: f64,
        range_noise: f64,
        range_rate_noise: f64,
        frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Self {
        Self {
            name,
            elevation_mask,
            latitude,
            longitude,
            height,
            frame,
            cosm,
            range_noise: Normal::new(0.0, range_noise).unwrap(),
            range_rate_noise: Normal::new(0.0, range_rate_noise).unwrap(),
        }
    }

    /// Initializes a point on the surface of a celestial object.
    /// This is meant for analysis, not for spacecraft navigation.
    pub fn from_point(
        name: String,
        latitude: f64,
        longitude: f64,
        height: f64,
        frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Self {
        Self::from_noise_values(
            name, 0.0, latitude, longitude, height, 0.0, 0.0, frame, cosm,
        )
    }

    pub fn dss65_madrid(
        elevation_mask: f64,
        range_noise: f64,
        range_rate_noise: f64,
        cosm: Arc<Cosm>,
    ) -> Self {
        Self::from_noise_values(
            "Madrid".to_string(),
            elevation_mask,
            40.427_222,
            4.250_556,
            0.834_939,
            range_noise,
            range_rate_noise,
            cosm.frame("IAU Earth"),
            cosm,
        )
    }

    pub fn dss34_canberra(
        elevation_mask: f64,
        range_noise: f64,
        range_rate_noise: f64,
        cosm: Arc<Cosm>,
    ) -> Self {
        Self::from_noise_values(
            "Canberra".to_string(),
            elevation_mask,
            -35.398_333,
            148.981_944,
            0.691_750,
            range_noise,
            range_rate_noise,
            cosm.frame("IAU Earth"),
            cosm,
        )
    }

    pub fn dss13_goldstone(
        elevation_mask: f64,
        range_noise: f64,
        range_rate_noise: f64,
        cosm: Arc<Cosm>,
    ) -> Self {
        Self::from_noise_values(
            "Goldstone".to_string(),
            elevation_mask,
            35.247_164,
            243.205,
            1.071_149_04,
            range_noise,
            range_rate_noise,
            cosm.frame("IAU Earth"),
            cosm,
        )
    }

    /// Computes the elevation of the provided object seen from this ground station.
    /// Also returns the ground station's orbit in the Solar System Barycenter frame
    pub fn elevation_of(&self, rx: &Orbit) -> (f64, Orbit, Orbit) {
        // Start by computing the rotation matrix from the inertial frame to the topocentric frame of this ground station.
        let dt = rx.dt;
        let station_state =
            Orbit::from_geodesic(self.latitude, self.longitude, self.height, dt, self.frame);

        // This is called the RFT in GMAT's TopocentricAxes.cpp.
        let dcm_topo2fixed = station_state.dcm_from_traj_frame(Frame::SEZ).unwrap();

        let ssb_frame = self.cosm.frame("SSB");

        // These are called RIF and RIFDot in GMAT.
        let (dcm_inertial2fixed, dcm_inertial2fixed_dt) = self
            .cosm
            .try_dcm_from_to_in_parts(&self.frame, &ssb_frame, dt)
            .unwrap();

        let dcm_inertial2topo = dcm_inertial2fixed * dcm_topo2fixed.transpose();
        let dcm_inertial2topo_dt = dcm_inertial2fixed_dt * dcm_topo2fixed.transpose();

        let dcm_topo2inertial = dcm_assemble(dcm_inertial2topo, dcm_inertial2topo_dt);

        // Now, compute the range vector in the solar system barycenter frame for both objects.
        // As per GMAT MathSpec, these ranges are different.
        let tx_ssb = self.cosm.frame_chg(&station_state, ssb_frame);
        let rx_ssb = self.cosm.frame_chg(rx, ssb_frame);

        let mut rho_sez_unit = (rx_ssb - tx_ssb)
            .with_rotation_by(dcm_topo2inertial.transpose())
            .radius();
        rho_sez_unit /= rho_sez_unit.norm();

        // NOTE: We don't need to normalize it here because we've already normalized rho_sez!
        let elevation = rho_sez_unit[2].asin().to_degrees();

        // Return elevation in degrees and tx
        (elevation, rx_ssb, tx_ssb)
    }
}
impl MeasurementDevice<Orbit, StdMeasurement> for GroundStation {
    /// Perform a measurement from the ground station to the receiver (rx).
    fn measure(&self, rx: &Orbit) -> Option<StdMeasurement> {
        let (elevation, rx_ssb, tx_ssb) = self.elevation_of(rx);

        Some(StdMeasurement::new(
            rx.dt,
            tx_ssb,
            rx_ssb,
            elevation >= self.elevation_mask,
            &self.range_noise,
            &self.range_rate_noise,
        ))
    }
}

impl MeasurementDevice<Spacecraft, StdMeasurement> for GroundStation {
    /// Perform a measurement from the ground station to the receiver (rx).
    fn measure(&self, sc_rx: &Spacecraft) -> Option<StdMeasurement> {
        let (elevation, rx_ssb, tx_ssb) = self.elevation_of(&sc_rx.orbit);

        Some(StdMeasurement::new(
            rx_ssb.dt,
            tx_ssb,
            rx_ssb,
            elevation >= self.elevation_mask,
            &self.range_noise,
            &self.range_rate_noise,
        ))
    }
}

impl fmt::Display for GroundStation {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "[{}] {} (lat.: {:.2} deg    long.: {:.2} deg    alt.: {:.2} m)",
            self.frame,
            self.name,
            self.latitude,
            self.longitude,
            self.height * 1e3,
        )
    }
}

/// Stores a standard measurement of range (km) and range rate (km/s)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StdMeasurement {
    pub dt: Epoch,
    pub obs: Vector2<f64>,
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

    fn compute_sensitivity(
        state: &OVector<Hyperdual<f64, U7>, U6>,
        range_noise: f64,
        range_rate_noise: f64,
    ) -> (Vector2<f64>, Matrix2x6<f64>) {
        // Extract data from hyperspace
        let range_vec = state.fixed_rows::<3>(0).into_owned();
        let velocity_vec = state.fixed_rows::<3>(3).into_owned();

        // Code up math as usual
        let delta_v_vec = velocity_vec / norm(&range_vec);
        let range = norm(&range_vec) + Hyperdual::from(range_noise);
        let range_rate = range_vec.dot(&delta_v_vec) + Hyperdual::from(range_rate_noise);

        // Extract result into Vector2 and Matrix2x6
        let mut fx = Vector2::zeros();
        let mut pmat = Matrix2x6::zeros();
        for i in 0..U2::dim() {
            fx[i] = if i == 0 {
                range.real()
            } else {
                range_rate.real()
            };
            for j in 1..U7::dim() {
                pmat[(i, j - 1)] = if i == 0 { range[j] } else { range_rate[j] };
            }
        }
        (fx, pmat)
    }

    /// Generate noiseless measurement
    pub fn noiseless(dt: Epoch, tx: Orbit, rx: Orbit, visible: bool) -> StdMeasurement {
        Self::raw(dt, tx, rx, visible, 0.0, 0.0)
    }

    /// Generate a new measurement with the provided noise distribution.
    pub fn new<D: Distribution<f64>>(
        dt: Epoch,
        tx: Orbit,
        rx: Orbit,
        visible: bool,
        range_dist: &D,
        range_rate_dist: &D,
    ) -> StdMeasurement {
        Self::raw(
            dt,
            tx,
            rx,
            visible,
            range_dist.sample(&mut thread_rng()),
            range_rate_dist.sample(&mut thread_rng()),
        )
    }

    /// Generate a new measurement with the provided noise values.
    pub fn raw(
        dt: Epoch,
        tx: Orbit,
        rx: Orbit,
        visible: bool,
        range_noise: f64,
        range_rate_noise: f64,
    ) -> StdMeasurement {
        assert_eq!(tx.frame, rx.frame, "tx & rx in different frames");
        assert_eq!(tx.dt, rx.dt, "tx & rx states have different times");

        let hyperstate = hyperspace_from_vector(&(rx - tx).to_cartesian_vec());
        let (obs, h_tilde) = Self::compute_sensitivity(&hyperstate, range_noise, range_rate_noise);

        StdMeasurement {
            dt,
            obs,
            visible,
            h_tilde,
        }
    }

    /// Initializes a StdMeasurement from real tracking data (sensitivity is zero)
    pub fn real(dt: Epoch, range: f64, range_rate: f64) -> Self {
        Self {
            dt,
            obs: Vector2::new(range, range_rate),
            visible: true,
            h_tilde: Matrix2x6::zeros(),
        }
    }
}

impl Measurement for StdMeasurement {
    type StateSize = U6;
    type MeasurementSize = U2;

    /// Returns this measurement as a vector of Range and Range Rate
    ///
    /// **Units:** km, km/s
    fn observation(&self) -> Vector2<f64> {
        self.obs
    }

    fn sensitivity(&self) -> Matrix2x6<f64> {
        self.h_tilde
    }

    fn visible(&self) -> bool {
        self.visible
    }
}

impl TimeTagged for StdMeasurement {
    fn epoch(&self) -> Epoch {
        self.dt
    }

    fn set_epoch(&mut self, dt: Epoch) {
        self.dt = dt
    }
}

impl Serialize for StdMeasurement {
    /// Serializes the observation vector at the given time.
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(3))?;
        seq.serialize_element(&self.dt.as_mjd_tai_days())?;
        let obs = self.observation();
        seq.serialize_element(&obs[(0, 0)])?;
        seq.serialize_element(&obs[(1, 0)])?;
        seq.end()
    }
}

/// Stores a standard measurement of range (km)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RangeMsr {
    pub dt: Epoch,
    pub obs: Vector1<f64>,
    visible: bool,
    h_tilde: Matrix1x6<f64>,
}

impl RangeMsr {
    pub fn range(&self) -> f64 {
        self.obs[(0, 0)]
    }

    fn compute_sensitivity(
        state: &OVector<Hyperdual<f64, U7>, U6>,
    ) -> (Vector1<f64>, Matrix1x6<f64>) {
        // Extract data from hyperspace
        let range_vec = state.fixed_rows::<3>(0).into_owned();

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

    pub fn new(tx: Orbit, rx: Orbit, visible: bool) -> RangeMsr {
        assert_eq!(tx.frame, rx.frame, "tx and rx in different frames");
        assert_eq!(tx.dt, rx.dt, "tx and rx states have different times");

        let dt = tx.dt;
        let hyperstate = hyperspace_from_vector(&(rx - tx).to_cartesian_vec());
        let (obs, h_tilde) = Self::compute_sensitivity(&hyperstate);

        RangeMsr {
            dt,
            obs,
            visible,
            h_tilde,
        }
    }
}

impl Measurement for RangeMsr {
    type StateSize = U6;
    type MeasurementSize = U1;

    /// Returns this measurement as a vector of Range and Range Rate
    ///
    /// **Units:** km
    fn observation(&self) -> Vector1<f64> {
        self.obs
    }

    fn sensitivity(&self) -> Matrix1x6<f64> {
        self.h_tilde
    }

    fn visible(&self) -> bool {
        self.visible
    }
}

impl TimeTagged for RangeMsr {
    fn epoch(&self) -> Epoch {
        self.dt
    }

    fn set_epoch(&mut self, dt: Epoch) {
        self.dt = dt
    }
}

impl Serialize for RangeMsr {
    /// Serializes the observation vector at the given time.
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(3))?;
        seq.serialize_element(&self.dt.as_mjd_tai_days())?;
        let obs = self.observation();
        seq.serialize_element(&obs[(0, 0)])?;
        seq.end()
    }
}

/// Stores a standard measurement of range (km)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DopplerMsr {
    pub dt: Epoch,
    pub obs: Vector1<f64>,
    visible: bool,
    h_tilde: Matrix1x6<f64>,
}

impl DopplerMsr {
    pub fn range_rate(&self) -> f64 {
        self.obs[(0, 0)]
    }

    fn compute_sensitivity(
        state: &OVector<Hyperdual<f64, U7>, U6>,
    ) -> (Vector1<f64>, Matrix1x6<f64>) {
        // Extract data from hyperspace
        let range_vec = state.fixed_rows::<3>(0).into_owned();
        let velocity_vec = state.fixed_rows::<3>(3).into_owned();

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

    pub fn new(tx: Orbit, rx: Orbit, visible: bool) -> DopplerMsr {
        assert_eq!(tx.frame, rx.frame, "tx and rx in different frames");
        assert_eq!(tx.dt, rx.dt, "tx and rx states have different times");

        let dt = tx.dt;
        let hyperstate = hyperspace_from_vector(&(rx - tx).to_cartesian_vec());
        let (obs, h_tilde) = Self::compute_sensitivity(&hyperstate);

        DopplerMsr {
            dt,
            obs,
            visible,
            h_tilde,
        }
    }
}

impl Measurement for DopplerMsr {
    type StateSize = U6;
    type MeasurementSize = U1;

    /// Returns this measurement as a vector of Range and Range Rate
    ///
    /// **Units:** km/s
    fn observation(&self) -> Vector1<f64> {
        self.obs
    }

    fn sensitivity(&self) -> Matrix1x6<f64> {
        self.h_tilde
    }

    fn visible(&self) -> bool {
        self.visible
    }
}

impl TimeTagged for DopplerMsr {
    fn epoch(&self) -> Epoch {
        self.dt
    }

    fn set_epoch(&mut self, dt: Epoch) {
        self.dt = dt
    }
}

impl Serialize for DopplerMsr {
    /// Serializes the observation vector at the given time.
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(3))?;
        seq.serialize_element(&self.dt.as_mjd_tai_days())?;
        let obs = self.observation();
        seq.serialize_element(&obs[(0, 0)])?;
        seq.end()
    }
}
