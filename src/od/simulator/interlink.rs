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
use anise::almanac::Almanac;
use anise::astro::Aberration;
use anise::errors::AlmanacResult;
use anise::prelude::{Frame, Orbit};
use hifitime::{Duration, Epoch, TimeSeries, TimeUnits};
use indexmap::{IndexMap, IndexSet};
use rand_pcg::Pcg64Mcg;
use serde::{Deserialize, Serialize};
use snafu::{ensure, ResultExt};

use crate::io::ConfigRepr;
use crate::md::prelude::Traj;
use crate::md::Trajectory;
use crate::od::msr::MeasurementType;
use crate::od::noise::StochasticNoise;
use crate::od::prelude::{Measurement, NoiseNotConfiguredSnafu, ODError};
use crate::od::TrackingDevice;
use crate::od::{ODAlmanacSnafu, ODTrajSnafu};
use crate::Spacecraft;
use crate::State;

use std::collections::BTreeMap;
use std::sync::Arc;

use super::TrkConfig;

#[derive(Clone)]
pub struct InterlinkArcSim {
    /// Receiver spacecraft in this link
    pub rx_spacecraft: Vec<Trajectory>,
    /// Transmitter spacercaft
    pub trajectory: InterlinkSpacecraft,
    /// Configuration of each device
    pub configs: BTreeMap<String, TrkConfig>,
    /// Random number generator used for this tracking arc, ensures repeatability
    rng: Pcg64Mcg,
    /// Greatest common denominator time series that allows this arc to meet all of the conditions.
    time_series: TimeSeries,
}

// Defines a (transmitter) spacecraft capable of inter-satellite links
#[derive(Clone, Debug)]
pub struct InterlinkSpacecraft {
    /// Trajectory of the spacercaft orbit
    pub traj: Trajectory,
    /// Measurement types supported by the link
    pub measurement_types: IndexSet<MeasurementType>,
    /// Integration time used to generate the measurement.
    pub integration_time: Option<Duration>,
    pub timestamp_noise_s: Option<StochasticNoise>,
    pub stochastic_noises: Option<IndexMap<MeasurementType, StochasticNoise>>,
    /// Aberration correction used in the interlink
    pub ab_corr: Option<Aberration>,
}

impl InterlinkSpacecraft {
    /// Returns the noises for all measurement types configured for this ground station at the provided epoch, timestamp noise is the first entry.
    fn noises(&mut self, epoch: Epoch, rng: Option<&mut Pcg64Mcg>) -> Result<Vec<f64>, ODError> {
        let mut noises = vec![0.0; self.measurement_types.len() + 1];

        if let Some(rng) = rng {
            ensure!(
                self.stochastic_noises.is_some(),
                NoiseNotConfiguredSnafu {
                    kind: "ground station stochastics".to_string(),
                }
            );
            // Add the timestamp noise first

            if let Some(mut timestamp_noise) = self.timestamp_noise_s {
                noises[0] = timestamp_noise.sample(epoch, rng);
            }

            let stochastics = self.stochastic_noises.as_mut().unwrap();

            for (ii, msr_type) in self.measurement_types.iter().enumerate() {
                noises[ii + 1] = stochastics
                    .get_mut(msr_type)
                    .ok_or(ODError::NoiseNotConfigured {
                        kind: format!("{msr_type:?}"),
                    })?
                    .sample(epoch, rng);
            }
        }

        Ok(noises)
    }
}

impl TrackingDevice<Spacecraft> for InterlinkSpacecraft {
    fn name(&self) -> String {
        self.traj.name.clone().unwrap_or("unnamed".to_string())
    }

    fn measurement_types(&self) -> &IndexSet<MeasurementType> {
        &self.measurement_types
    }

    fn location(&self, epoch: Epoch, frame: Frame, almanac: Arc<Almanac>) -> AlmanacResult<Orbit> {
        almanac.transform_to(self.traj.at(epoch).unwrap().orbit, frame, self.ab_corr)
    }

    fn measure(
        &mut self,
        epoch: Epoch,
        traj: &Traj<Spacecraft>,
        rng: Option<&mut Pcg64Mcg>,
        almanac: Arc<Almanac>,
    ) -> Result<Option<Measurement>, ODError> {
        match self.integration_time {
            Some(integration_time) => {
                // TODO: This should support measurement alignment
                // If out of traj bounds, return None, else the whole strand is rejected.
                let rx_0 = match traj.at(epoch - integration_time).context(ODTrajSnafu {
                    details: format!(
                        "fetching state {epoch} at start of ground station integration time {integration_time}"
                    ),
                }) {
                    Ok(rx) => rx,
                    Err(_) => return Ok(None),
                };

                let rx_1 = match traj.at(epoch).context(ODTrajSnafu {
                    details: format!(
                        "fetching state {epoch} at end of ground station integration time"
                    ),
                }) {
                    Ok(rx) => rx,
                    Err(_) => return Ok(None),
                };

                // Start of integration time
                let msr_t0 = self.measure_instantaneous(rx_0, None, almanac.clone())?;

                // End of integration time
                let msr_t1 = self.measure_instantaneous(rx_1, None, almanac.clone())?;

                if msr_t0.is_some() && msr_t1.is_some() {
                    // Line of sight in both cases

                    // Noises are computed at the midpoint of the integration time.
                    let noises = self.noises(epoch - integration_time * 0.5, rng)?;

                    let mut msr = Measurement::new(self.name(), epoch + noises[0].seconds());

                    for (ii, msr_type) in self.measurement_types.iter().enumerate() {
                        let msr_value_0 = msr_t0.as_ref().unwrap().data[msr_type];
                        let msr_value_1 = msr_t1.as_ref().unwrap().data[msr_type];

                        let msr_value =
                            (msr_value_1 - msr_value_0) * 0.5 + noises[ii + 1] / 2.0_f64.sqrt();
                        msr.push(*msr_type, msr_value);
                    }

                    Ok(Some(msr))
                } else {
                    Ok(None)
                }
            }
            None => self.measure_instantaneous(
                traj.at(epoch).context(ODTrajSnafu {
                    details: "fetching state for instantaneous measurement".to_string(),
                })?,
                rng,
                almanac,
            ),
        }
    }

    fn measure_instantaneous(
        &mut self,
        rx: Spacecraft,
        rng: Option<&mut Pcg64Mcg>,
        almanac: Arc<Almanac>,
    ) -> Result<Option<Measurement>, ODError> {
        let observer = self.traj.at(rx.epoch()).context(ODTrajSnafu {
            details: format!("fetching state {} for interlink", rx.epoch()),
        })?;

        let is_obstructed = almanac
            .line_of_sight_obstructed(observer.orbit, rx.orbit, observer.orbit.frame, self.ab_corr)
            .context(ODAlmanacSnafu {
                action: "computing line of sight",
            })?;

        if is_obstructed {
            Ok(None)
        } else {
            // Convert the receiver into the body fixed transmitter frame.
            let rx_in_tx_frame = almanac
                .transform_to(rx.orbit, observer.orbit.frame, self.ab_corr)
                .context(ODAlmanacSnafu {
                    action: "transforming receiver to transmitter frame",
                })?;

            let rho_tx_frame = rx_in_tx_frame.radius_km - observer.orbit.radius_km;

            // Compute the range-rate \dot ρ. Note that rx_in_tx_frame is already the relative velocity of rx wrt tx!
            let range_rate_km_s =
                rho_tx_frame.dot(&rx_in_tx_frame.velocity_km_s) / rho_tx_frame.norm();

            let noises = self.noises(observer.epoch(), rng)?;

            let mut msr = Measurement::new(self.name(), rx.orbit.epoch + noises[0].seconds());

            for (ii, msr_type) in self.measurement_types.iter().enumerate() {
                let msr_value = match ii {
                    0 => rho_tx_frame.norm(),
                    1 => range_rate_km_s,
                    _ => unreachable!(),
                } + noises[ii + 1];
                msr.push(*msr_type, msr_value);
            }

            Ok(Some(msr))
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

impl Serialize for InterlinkSpacecraft {
    fn serialize<S>(&self, _serializer: S) -> Result<S::Ok, S::Error>
    where
        S: serde::Serializer,
    {
        unimplemented!("interlink spacecraft cannot be serialized")
    }
}

impl<'de> Deserialize<'de> for InterlinkSpacecraft {
    fn deserialize<D>(_deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        unimplemented!("interlink spacecraft cannot be deserialized")
    }
}

impl ConfigRepr for InterlinkSpacecraft {}
