/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::msr::RangeDoppler;
use super::noise::GaussMarkov;
use super::TrackingDeviceSim;
use crate::cosmic::{Cosm, Frame, Orbit};
use crate::io::{frame_from_str, frame_to_str, ConfigRepr, Configurable};
use crate::md::prelude::{Interpolatable, Traj};
use crate::md::EventEvaluator;
use crate::time::Epoch;
use crate::utils::between_0_360;
use crate::{NyxError, Spacecraft};
use hifitime::{Duration, Unit};
use nalgebra::{allocator::Allocator, DefaultAllocator};
use rand_pcg::Pcg64Mcg;
use serde_derive::{Deserialize, Serialize};
use std::fmt;
use std::sync::Arc;

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// GroundStation defines a two-way ranging and doppler station.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(module = "nyx_space.orbit_determination"))]
pub struct GroundStation {
    pub name: String,
    /// in degrees
    pub elevation_mask_deg: f64,
    /// in degrees
    pub latitude_deg: f64,
    /// in degrees
    pub longitude_deg: f64,
    /// in km
    pub height_km: f64,
    /// Frame in which this station is defined
    #[serde(serialize_with = "frame_to_str", deserialize_with = "frame_from_str")]
    pub frame: Frame,
    /// Duration needed to generate a measurement (if unset, it is assumed to be instantaneous)
    #[serde(skip)]
    pub integration_time: Option<Duration>,
    /// Whether to correct for light travel time
    pub light_time_correction: bool,
    /// Noise on the timestamp of the measurement
    pub timestamp_noise_s: Option<GaussMarkov>,
    /// Noise on the range data of the measurement
    pub range_noise_km: Option<GaussMarkov>,
    /// Noise on the Doppler data of the measurement
    pub doppler_noise_km_s: Option<GaussMarkov>,
}

impl GroundStation {
    /// Initializes a point on the surface of a celestial object.
    /// This is meant for analysis, not for spacecraft navigation.
    pub fn from_point(
        name: String,
        latitude_deg: f64,
        longitude_deg: f64,
        height_km: f64,
        frame: Frame,
    ) -> Self {
        Self {
            name,
            elevation_mask_deg: 0.0,
            latitude_deg,
            longitude_deg,
            height_km,
            frame,
            integration_time: None,
            light_time_correction: false,
            timestamp_noise_s: None,
            range_noise_km: None,
            doppler_noise_km_s: None,
        }
    }

    pub fn dss65_madrid(
        elevation_mask: f64,
        range_noise_km: GaussMarkov,
        doppler_noise_km_s: GaussMarkov,
        iau_earth: Frame,
    ) -> Self {
        Self {
            name: "Madrid".to_string(),
            elevation_mask_deg: elevation_mask,
            latitude_deg: 40.427_222,
            longitude_deg: 4.250_556,
            height_km: 0.834_939,
            frame: iau_earth,
            integration_time: None,
            light_time_correction: false,
            timestamp_noise_s: None,
            range_noise_km: Some(range_noise_km),
            doppler_noise_km_s: Some(doppler_noise_km_s),
        }
    }

    pub fn dss34_canberra(
        elevation_mask: f64,
        range_noise_km: GaussMarkov,
        doppler_noise_km_s: GaussMarkov,
        iau_earth: Frame,
    ) -> Self {
        Self {
            name: "Canberra".to_string(),
            elevation_mask_deg: elevation_mask,
            latitude_deg: -35.398_333,
            longitude_deg: 148.981_944,
            height_km: 0.691_750,
            frame: iau_earth,
            integration_time: None,
            light_time_correction: false,
            timestamp_noise_s: None,
            range_noise_km: Some(range_noise_km),
            doppler_noise_km_s: Some(doppler_noise_km_s),
        }
    }

    pub fn dss13_goldstone(
        elevation_mask: f64,
        range_noise_km: GaussMarkov,
        doppler_noise_km_s: GaussMarkov,
        iau_earth: Frame,
    ) -> Self {
        Self {
            name: "Goldstone".to_string(),
            elevation_mask_deg: elevation_mask,
            latitude_deg: 35.247_164,
            longitude_deg: 243.205,
            height_km: 1.071_149_04,
            frame: iau_earth,
            integration_time: None,
            light_time_correction: false,
            timestamp_noise_s: None,
            range_noise_km: Some(range_noise_km),
            doppler_noise_km_s: Some(doppler_noise_km_s),
        }
    }

    /// Computes the azimuth and elevation of the provided object seen from this ground station, both in degrees.
    /// Also returns the ground station's orbit in the frame of the receiver
    pub fn azimuth_elevation_of(&self, rx: Orbit, cosm: &Cosm) -> (f64, f64, Orbit, Orbit) {
        // Start by converting the receiver spacecraft into the ground station frame.
        let rx_gs_frame = cosm.frame_chg(&rx, self.frame);

        let dt = rx.epoch;
        // Then, compute the rotation matrix from the body fixed frame of the ground station to its topocentric frame SEZ.
        let tx_gs_frame = self.to_orbit(dt);
        // Note: we're only looking at the radii so we don't need to apply the transport theorem here.
        let dcm_topo2fixed = tx_gs_frame.dcm_from_traj_frame(Frame::SEZ).unwrap();

        // Now, rotate the spacecraft in the SEZ frame to compute its elevation as seen from the ground station.
        // We transpose the DCM so that it's the fixed to topocentric rotation.
        let rx_sez = rx_gs_frame.with_position_rotated_by(dcm_topo2fixed.transpose());
        let tx_sez = tx_gs_frame.with_position_rotated_by(dcm_topo2fixed.transpose());
        // Now, let's compute the range ρ.
        let rho_sez = rx_sez - tx_sez;

        // Finally, compute the elevation (math is the same as declination)
        // Source: Vallado, section 4.4.3
        // Only the sine is needed as per Vallado, and the formula is the same as the declination
        // because we're in the SEZ frame.
        let elevation_deg = rho_sez.declination_deg();
        if (elevation_deg - 90.0).abs() < 1e-6 {
            warn!("object nearly overhead (el = {elevation_deg} deg), azimuth may be incorrect");
        }
        // For the elevation, we need to perform a quadrant check because it's measured from 0 to 360 degrees.
        let azimuth_deg = between_0_360((-rho_sez.y_km.atan2(rho_sez.x_km)).to_degrees());

        // Return elevation in degrees and rx/tx in the inertial frame of the spacecraft
        (
            azimuth_deg,
            elevation_deg,
            rx,
            cosm.frame_chg(&tx_gs_frame, rx.frame),
        )
    }

    /// Return this ground station as an orbit in its current frame
    pub fn to_orbit(&self, epoch: Epoch) -> Orbit {
        Orbit::from_geodesic(
            self.latitude_deg,
            self.longitude_deg,
            self.height_km,
            epoch,
            self.frame,
        )
    }

    /// Returns the timestamp noise, range noise, and doppler noise for this ground station at the provided epoch.
    fn noises(
        &mut self,
        epoch: Epoch,
        rng: Option<&mut Pcg64Mcg>,
    ) -> Result<(f64, f64, f64), NyxError> {
        let timestamp_noise_s;
        let range_noise_km;
        let doppler_noise_km_s;

        match rng {
            Some(rng) => {
                // Add the range noise, or return an error if it's not configured.
                range_noise_km = self
                    .range_noise_km
                    .ok_or_else(|| NyxError::CustomError("Range noise not configured".to_string()))?
                    .next_bias(epoch, rng);

                // Add the Doppler noise, or return an error if it's not configured.
                doppler_noise_km_s = self
                    .doppler_noise_km_s
                    .ok_or_else(|| {
                        NyxError::CustomError("Doppler noise not configured".to_string())
                    })?
                    .next_bias(epoch, rng);

                // Only add the epoch noise if it's configured, it's valid to not have any noise on the clock.
                if let Some(mut timestamp_noise) = self.timestamp_noise_s {
                    timestamp_noise_s = timestamp_noise.next_bias(epoch, rng);
                } else {
                    timestamp_noise_s = 0.0;
                }
            }
            None => {
                timestamp_noise_s = 0.0;
                range_noise_km = 0.0;
                doppler_noise_km_s = 0.0;
            }
        };

        Ok((timestamp_noise_s, range_noise_km, doppler_noise_km_s))
    }
}

impl ConfigRepr for GroundStation {}

impl Configurable for GroundStation {
    type IntermediateRepr = GroundStation;

    fn from_config(
        cfg: Self::IntermediateRepr,
        _cosm: Arc<Cosm>,
    ) -> Result<Self, crate::io::ConfigError>
    where
        Self: Sized,
    {
        Ok(cfg)
    }

    fn to_config(&self) -> Result<Self::IntermediateRepr, crate::io::ConfigError> {
        Ok(self.clone())
    }
}

impl TrackingDeviceSim<Orbit, RangeDoppler> for GroundStation {
    /// Perform a measurement from the ground station to the receiver. If there is no integration time of the measurement, then this is assumed to be an instantaneous measurement instead of a two way measurement.
    ///
    /// WARNING: For StdMeasurement, we just call the instantaneous measurement function.
    fn measure(
        &mut self,
        epoch: Epoch,
        traj: &Traj<Orbit>,
        rng: Option<&mut Pcg64Mcg>,
        cosm: Arc<Cosm>,
    ) -> Result<Option<RangeDoppler>, NyxError> {
        match self.integration_time {
            Some(integration_time) => {
                let rx_0 = traj.at(epoch - integration_time)?;
                let rx_1 = traj.at(epoch)?;

                let (_, elevation_0, rx_0, tx_0) = self.azimuth_elevation_of(rx_0, &cosm);
                let (_, elevation_1, rx_1, tx_1) = self.azimuth_elevation_of(rx_1, &cosm);

                if elevation_0 < self.elevation_mask_deg || elevation_1 < self.elevation_mask_deg {
                    debug!(
                        "{} (el. mask {:.3} deg) but object moves from {elevation_0:.3} to {elevation_1:.3} deg -- no measurement",
                        self.name, self.elevation_mask_deg
                    );
                    return Ok(None);
                }

                // Noises are computed at the midpoint of the integration time.
                let (timestamp_noise_s, range_noise_km, doppler_noise_km_s) =
                    self.noises(epoch - integration_time * 0.5, rng)?;

                Ok(Some(RangeDoppler::two_way(
                    (tx_0, tx_1),
                    (rx_0, rx_1),
                    timestamp_noise_s,
                    range_noise_km,
                    doppler_noise_km_s,
                )))
            }
            None => self.measure_instantaneous(traj.at(epoch)?, rng, cosm),
        }
    }

    fn name(&self) -> String {
        self.name.clone()
    }

    fn location(&self, epoch: Epoch, frame: Frame, cosm: &Cosm) -> Orbit {
        cosm.frame_chg(&self.to_orbit(epoch), frame)
    }

    fn measure_instantaneous(
        &mut self,
        rx: Orbit,
        rng: Option<&mut Pcg64Mcg>,
        cosm: Arc<Cosm>,
    ) -> Result<Option<RangeDoppler>, NyxError> {
        let (_, elevation, rx, tx) = self.azimuth_elevation_of(rx, &cosm);

        if elevation >= self.elevation_mask_deg {
            // Only update the noises if the measurement is valid.
            let (timestamp_noise_s, range_noise_km, doppler_noise_km_s) =
                self.noises(rx.epoch, rng)?;

            Ok(Some(RangeDoppler::one_way(
                tx,
                rx,
                timestamp_noise_s,
                range_noise_km,
                doppler_noise_km_s,
            )))
        } else {
            debug!(
                "{} (el. mask {:.3} deg), object at {elevation:.3} deg -- no measurement",
                self.name, self.elevation_mask_deg
            );
            Ok(None)
        }
    }
}

impl TrackingDeviceSim<Spacecraft, RangeDoppler> for GroundStation {
    /// Perform a measurement from the ground station to the receiver (rx).
    fn measure(
        &mut self,
        epoch: Epoch,
        traj: &Traj<Spacecraft>,
        rng: Option<&mut Pcg64Mcg>,
        cosm: Arc<Cosm>,
    ) -> Result<Option<RangeDoppler>, NyxError> {
        let rx = traj.at(epoch)?;
        self.measure_instantaneous(rx, rng, cosm)
    }

    fn name(&self) -> String {
        self.name.clone()
    }

    fn location(&self, epoch: Epoch, frame: Frame, cosm: &Cosm) -> Orbit {
        cosm.frame_chg(&self.to_orbit(epoch), frame)
    }

    fn measure_instantaneous(
        &mut self,
        rx: Spacecraft,
        rng: Option<&mut Pcg64Mcg>,
        cosm: Arc<Cosm>,
    ) -> Result<Option<RangeDoppler>, NyxError> {
        let (_, elevation, rx, tx) = self.azimuth_elevation_of(rx.orbit, &cosm);

        if elevation >= self.elevation_mask_deg {
            // Only update the noises if the measurement is valid.
            let (timestamp_noise_s, range_noise_km, doppler_noise_km_s) =
                self.noises(rx.epoch, rng)?;

            Ok(Some(RangeDoppler::one_way(
                tx,
                rx,
                timestamp_noise_s,
                range_noise_km,
                doppler_noise_km_s,
            )))
        } else {
            debug!(
                "{} (el. mask {:.3} deg), object at {elevation:.3} deg -- no measurement",
                self.name, self.elevation_mask_deg
            );
            Ok(None)
        }
    }
}

impl fmt::Display for GroundStation {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{} (lat.: {:.4} deg    long.: {:.4} deg    alt.: {:.3} m) [{}]",
            self.name,
            self.latitude_deg,
            self.longitude_deg,
            self.height_km * 1e3,
            self.frame,
        )
    }
}

impl<S: Interpolatable> EventEvaluator<S> for &GroundStation
where
    DefaultAllocator:
        Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size> + Allocator<f64, S::VecLength>,
{
    /// Compute the elevation in the SEZ frame. This call will panic if the frame of the input state does not match that of the ground station.
    fn eval(&self, rx_gs_frame: &S) -> f64 {
        let dt = rx_gs_frame.epoch();
        // Then, compute the rotation matrix from the body fixed frame of the ground station to its topocentric frame SEZ.
        let tx_gs_frame = self.to_orbit(dt);
        // Note: we're only looking at the radii so we don't need to apply the transport theorem here.
        let dcm_topo2fixed = tx_gs_frame.dcm_from_traj_frame(Frame::SEZ).unwrap();

        // Now, rotate the spacecraft in the SEZ frame to compute its elevation as seen from the ground station.
        // We transpose the DCM so that it's the fixed to topocentric rotation.
        let rx_sez = rx_gs_frame
            .orbit()
            .with_position_rotated_by(dcm_topo2fixed.transpose());
        let tx_sez = tx_gs_frame.with_position_rotated_by(dcm_topo2fixed.transpose());
        // Now, let's compute the range ρ.
        let rho_sez = rx_sez - tx_sez;

        // Finally, compute the elevation (math is the same as declination)
        // Source: Vallado, section 4.4.3
        // Only the sine is needed as per Vallado, and the formula is the same as the declination
        // because we're in the SEZ frame.
        rho_sez.declination_deg() - self.elevation_mask_deg
    }

    fn eval_string(&self, state: &S) -> String {
        format!(
            "Elevation from {} is {:.6} deg on {}",
            self.name,
            self.eval(state) + self.elevation_mask_deg,
            state.epoch()
        )
    }

    /// Epoch precision of the election evaluator is 1 ms
    fn epoch_precision(&self) -> Duration {
        1 * Unit::Second
    }

    /// Angle precision of the elevation evaluator is 1 millidegree.
    fn value_precision(&self) -> f64 {
        1e-3
    }
}

#[cfg(test)]
mod gs_ut {
    use crate::io::ConfigRepr;
    use crate::od::prelude::*;

    #[test]
    fn test_load_single() {
        use std::env;
        use std::path::PathBuf;

        use hifitime::TimeUnits;

        let cosm = Cosm::de438();

        // Get the path to the root directory of the current Cargo project
        let test_data: PathBuf = [
            env::var("CARGO_MANIFEST_DIR").unwrap(),
            "data".to_string(),
            "tests".to_string(),
            "config".to_string(),
            "one_ground_station.yaml".to_string(),
        ]
        .iter()
        .collect();

        assert!(test_data.exists(), "Could not find the test data");

        let gs = GroundStation::load(test_data).unwrap();

        dbg!(&gs);

        let expected_gs = GroundStation {
            name: "Demo ground station".to_string(),
            frame: cosm.frame("IAU Earth"),
            elevation_mask_deg: 5.0,
            range_noise_km: Some(GaussMarkov::new(1.days(), 5e-3, 1e-4).unwrap()),
            doppler_noise_km_s: Some(GaussMarkov::new(1.days(), 5e-5, 1.5e-6).unwrap()),
            latitude_deg: 2.3522,
            longitude_deg: 48.8566,
            height_km: 0.4,
            light_time_correction: false,
            timestamp_noise_s: None,
            integration_time: None,
        };

        assert_eq!(expected_gs, gs);
    }

    #[test]
    fn test_load_many() {
        use crate::cosmic::Cosm;
        use hifitime::TimeUnits;
        use std::env;
        use std::path::PathBuf;

        // Get the path to the root directory of the current Cargo project

        let test_file: PathBuf = [
            env::var("CARGO_MANIFEST_DIR").unwrap(),
            "data".to_string(),
            "tests".to_string(),
            "config".to_string(),
            "many_ground_stations.yaml".to_string(),
        ]
        .iter()
        .collect();

        let stations = GroundStation::load_many(test_file).unwrap();

        dbg!(&stations);

        let cosm = Cosm::de438();

        let expected = vec![
            GroundStation {
                name: "Demo ground station".to_string(),
                frame: cosm.frame("IAU Earth"),
                elevation_mask_deg: 5.0,
                range_noise_km: Some(GaussMarkov::new(1.days(), 5e-3, 1e-4).unwrap()),
                doppler_noise_km_s: Some(GaussMarkov::new(1.days(), 5e-5, 1.5e-6).unwrap()),
                latitude_deg: 2.3522,
                longitude_deg: 48.8566,
                height_km: 0.4,
                light_time_correction: false,
                timestamp_noise_s: None,
                integration_time: None,
            },
            GroundStation {
                name: "Canberra".to_string(),
                frame: cosm.frame("IAU Earth"),
                elevation_mask_deg: 5.0,
                range_noise_km: Some(GaussMarkov::new(1.days(), 5e-3, 1e-4).unwrap()),
                doppler_noise_km_s: Some(GaussMarkov::new(1.days(), 5e-5, 1.5e-6).unwrap()),
                latitude_deg: -35.398333,
                longitude_deg: 148.981944,
                height_km: 0.691750,
                light_time_correction: false,
                timestamp_noise_s: None,
                integration_time: None,
            },
        ];

        assert_eq!(expected, stations);
    }
}
