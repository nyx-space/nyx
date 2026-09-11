/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::{ODAlmanacSnafu, ODError, ODTrajSnafu, TrackingDevice};
use crate::Spacecraft;
use crate::io::ConfigError;
use crate::md::prelude::Traj;
use crate::od::msr::measurement::Measurement;
use crate::od::msr::two_way::solve_two_way_picard;
use crate::od::msr::{IntegrationRef, MeasurementType};
use crate::time::Epoch;
use anise::astro::Aberration;
use anise::errors::AlmanacResult;
use anise::frames::Frame;
use anise::prelude::{Almanac, Orbit};
use hifitime::TimeUnits;
use indexmap::IndexSet;
use log::debug;
use rand_pcg::Pcg64Mcg;
use snafu::ResultExt;

use super::GroundStation;

impl TrackingDevice<Spacecraft> for GroundStation {
    fn measurement_types(&self) -> &IndexSet<MeasurementType> {
        &self.measurement_types
    }

    /// Perform a measurement from the ground station to the receiver (rx).
    /// The epoch MUST be the station reception epoch.
    fn measure(
        &mut self,
        epoch: Epoch,
        traj: &Traj<Spacecraft>,
        rng: Option<&mut Pcg64Mcg>,
        almanac: &Almanac,
    ) -> Result<Option<Measurement>, ODError> {
        let mut msr = Measurement::new(self.name.clone(), epoch);
        msr.doppler_config = self.doppler_config;

        if self.light_time_correction {
            // Solve for Relativistic Picard Light-Time

            // Step 2a: Solve Downlink Leg (t3 -> t2) for Range & Angles
            let two_way_sol = match solve_two_way_picard(epoch, self, traj, almanac) {
                Ok(sol) => sol,
                Err(_) => return Ok(None),
            };

            let rx = traj.at(two_way_sol.t2_bounce).context(ODTrajSnafu {
                details: "fetching state for bounce epoch".to_string(),
            })?;

            if let Some(obstruction_body) = self.obstruction_body {
                // Check if there was an obstruction at the bounce epoch with LT
                let observer =
                    Spacecraft::from(self.to_orbit(epoch, almanac).context(ODAlmanacSnafu {
                        action: "building ground station orbit",
                    })?);
                let is_obstructed = almanac
                    .line_of_sight_obstructed(
                        observer.orbit,
                        rx.orbit,
                        obstruction_body.into(),
                        Aberration::LT,
                    )
                    .context(ODAlmanacSnafu {
                        action: "computing line of sight",
                    })?;

                if is_obstructed {
                    return Ok(None);
                }
            }

            // Evaluate Azimuth/Elevation from Downlink Look Direction at t3
            // Construct the apparent target state using the solved bounce state r_sc(t2)
            let aer_downlink = almanac
                .azimuth_elevation_range_sez_from_location(
                    rx.orbit,
                    self.location.clone(),
                    None,
                    None, // Position r_sc(t2) is already retarded; do not apply LT twice
                )
                .context(ODAlmanacSnafu {
                    action: "computing downlink AER",
                })?;

            if aer_downlink.elevation_above_mask_deg() < 0.0 || aer_downlink.is_obstructed() {
                return Ok(None);
            }

            let noises = self.noises(epoch, rng)?;

            for (ii, msr_type) in self.measurement_types.iter().enumerate() {
                let noise = noises[ii + 1];
                let val = match msr_type {
                    MeasurementType::Range => two_way_sol.range_km() + noise,

                    MeasurementType::Azimuth => aer_downlink.azimuth_deg + noise,
                    MeasurementType::Elevation => aer_downlink.elevation_deg + noise,

                    MeasurementType::Doppler => {
                        let doppler_cfg = msr
                            .doppler_config
                            .ok_or_else(|| ODError::ODConfigError { source: ConfigError::InvalidConfig {
                                msg: "Doppler measurement requires doppler_config on GroundStation".to_string()
                            }})?;

                        let integr_time = doppler_cfg.integration_time;

                        // Compute integration window boundaries from integration_ref
                        let (t_start, t_end) = match doppler_cfg.integration_ref {
                            IntegrationRef::Start => (epoch, epoch + integr_time),
                            IntegrationRef::Middle => {
                                (epoch - integr_time * 0.5, epoch + integr_time * 0.5)
                            }
                            IntegrationRef::End => (epoch - integr_time, epoch),
                        };

                        // Evaluate two-way ranges at window boundaries
                        let r_start =
                            solve_two_way_picard(t_start, self, traj, almanac)?.range_km();
                        let r_end = solve_two_way_picard(t_end, self, traj, almanac)?.range_km();

                        // Differenced range rate + Doppler noise
                        ((r_end - r_start) / integr_time.to_seconds()) + noise
                    }

                    _ => {
                        return Err(ODError::ODLimitation {
                            action: format!("MeasurementType::{msr_type:?} is unsupported"),
                        });
                    }
                };

                msr.push(*msr_type, val);
            }

            Ok(Some(msr))
        } else {
            let rx = traj.at(epoch).context(ODTrajSnafu {
                details: "fetching state for instantaneous measurement".to_string(),
            })?;

            if let Some(obstruction_body) = self.obstruction_body {
                let observer =
                    Spacecraft::from(self.to_orbit(epoch, almanac).context(ODAlmanacSnafu {
                        action: "building ground station orbit",
                    })?);
                let ab_corr = Aberration::NONE;
                let is_obstructed = almanac
                    .line_of_sight_obstructed(
                        observer.orbit,
                        rx.orbit,
                        obstruction_body.into(),
                        ab_corr,
                    )
                    .context(ODAlmanacSnafu {
                        action: "computing line of sight",
                    })?;

                if is_obstructed {
                    return Ok(None);
                }
            }

            let aer = almanac
                .azimuth_elevation_range_sez_from_location(
                    rx.orbit,
                    self.location.clone(),
                    None,
                    None,
                )
                .context(ODAlmanacSnafu {
                    action: "computing AER",
                })?;

            if aer.elevation_above_mask_deg() >= 0.0 && !aer.is_obstructed() {
                // Only update the noises if the measurement is valid.
                let noises = self.noises(rx.orbit.epoch, rng)?;

                let mut msr =
                    Measurement::new(self.name.clone(), rx.orbit.epoch + noises[0].seconds());
                msr.doppler_config = self.doppler_config;

                for (ii, msr_type) in self.measurement_types.iter().enumerate() {
                    let msr_value = if msr_type == &MeasurementType::Doppler {
                        if let Some(doppler_cfg) = msr.doppler_config {
                            let integr_time = doppler_cfg.integration_time;

                            // Compute integration window boundaries from integration_ref
                            let (t_start, t_end) = match doppler_cfg.integration_ref {
                                IntegrationRef::Start => (epoch, epoch + integr_time),
                                IntegrationRef::Middle => {
                                    (epoch - integr_time * 0.5, epoch + integr_time * 0.5)
                                }
                                IntegrationRef::End => (epoch - integr_time, epoch),
                            };

                            // Evaluate two-way ranges at window boundaries
                            let sc_start = traj.at(t_start).context(ODTrajSnafu {
                                details: "fetching state for start of integration".to_string(),
                            })?;
                            let sc_end = traj.at(t_end).context(ODTrajSnafu {
                                details: "fetching state for end of integration".to_string(),
                            })?;
                            let aer_start = almanac
                                .azimuth_elevation_range_sez_from_location(
                                    sc_start.orbit,
                                    self.location.clone(),
                                    None,
                                    None,
                                )
                                .context(ODAlmanacSnafu {
                                    action: "computing AER at start of integration time",
                                })?;

                            let aer_end = almanac
                                .azimuth_elevation_range_sez_from_location(
                                    sc_end.orbit,
                                    self.location.clone(),
                                    None,
                                    None,
                                )
                                .context(ODAlmanacSnafu {
                                    action: "computing AER at end of integration time",
                                })?;

                            // Differenced range rate + Doppler noise
                            ((aer_end.range_km - aer_start.range_km) / integr_time.to_seconds())
                                + noises[ii + 1]
                        } else {
                            msr_type.compute_one_way(aer, noises[ii + 1])?
                        }
                    } else {
                        msr_type.compute_one_way(aer, noises[ii + 1])?
                    };
                    msr.push(*msr_type, msr_value);
                }

                Ok(Some(msr))
            } else {
                debug!(
                    "{} {} object at {:.3} deg -- no measurement",
                    self.name,
                    rx.orbit.epoch,
                    aer.elevation_above_mask_deg(),
                );
                Ok(None)
            }
        }
    }

    fn name(&self) -> String {
        self.name.clone()
    }

    fn location(&self, epoch: Epoch, frame: Frame, almanac: &Almanac) -> AlmanacResult<Orbit> {
        almanac.transform_to(self.to_orbit(epoch, almanac).unwrap(), frame, None)
    }

    fn measure_instantaneous(
        &mut self,
        rx: Spacecraft,
        rng: Option<&mut Pcg64Mcg>,
        almanac: &Almanac,
    ) -> Result<Option<Measurement>, ODError> {
        // HACK This function should be avoided. A future version will remove the instantaneous measurement
        // because it isn't physically adequate.
        if let Some(obstruction_body) = self.obstruction_body {
            let observer = Spacecraft::from(self.to_orbit(rx.orbit.epoch, almanac).context(
                ODAlmanacSnafu {
                    action: "building ground station orbit",
                },
            )?);
            let ab_corr = if self.light_time_correction {
                Aberration::LT
            } else {
                Aberration::NONE
            };
            let is_obstructed = almanac
                .line_of_sight_obstructed(
                    observer.orbit,
                    rx.orbit,
                    obstruction_body.into(),
                    ab_corr,
                )
                .context(ODAlmanacSnafu {
                    action: "computing line of sight",
                })?;

            if is_obstructed {
                return Ok(None);
            }
        }

        let aer = almanac
            .azimuth_elevation_range_sez_from_location(rx.orbit, self.location.clone(), None, None)
            .context(ODAlmanacSnafu {
                action: "computing AER",
            })?;

        if aer.elevation_above_mask_deg() >= 0.0 && !aer.is_obstructed() {
            // Only update the noises if the measurement is valid.
            let noises = self.noises(rx.orbit.epoch, rng)?;

            let mut msr = Measurement::new(self.name.clone(), rx.orbit.epoch + noises[0].seconds());
            msr.doppler_config = self.doppler_config;

            for (ii, msr_type) in self.measurement_types.iter().enumerate() {
                let msr_value = msr_type.compute_one_way(aer, noises[ii + 1])?;
                msr.push(*msr_type, msr_value);
            }

            Ok(Some(msr))
        } else {
            debug!(
                "{} {} object at {:.3} deg -- no measurement",
                self.name,
                rx.orbit.epoch,
                aer.elevation_above_mask_deg(),
            );
            Ok(None)
        }
    }

    /// Returns the measurement noise of this ground station.
    ///
    /// # Methodology
    /// Noises are modeled using a [StochasticNoise] process, defined by the sigma on the turn-on bias and on the steady state noise.
    /// The measurement noise is computed assuming that all measurements are independent variables, i.e. the measurement matrix is
    /// a diagonal matrix. The first item in the diagonal is the range noise (in km), set to the square of the steady state sigma. The
    /// second item is the Doppler noise (in km/s), set to the square of the steady state sigma of that Gauss Markov process.
    fn measurement_covar(&self, msr_type: MeasurementType, epoch: Epoch) -> Result<f64, ODError> {
        let stochastics = self.stochastic_noises.as_ref().unwrap();

        Ok(stochastics
            .get(&msr_type)
            .ok_or(ODError::NoiseNotConfigured {
                kind: format!("{msr_type:?}"),
            })?
            .covariance(epoch))
    }

    fn measurement_bias(&self, msr_type: MeasurementType, _epoch: Epoch) -> Result<f64, ODError> {
        let stochastics = self.stochastic_noises.as_ref().unwrap();

        if let Some(gm) = stochastics
            .get(&msr_type)
            .ok_or(ODError::NoiseNotConfigured {
                kind: format!("{msr_type:?}"),
            })?
            .bias
        {
            Ok(gm.constant.unwrap_or(0.0))
        } else {
            Ok(0.0)
        }
    }
}
