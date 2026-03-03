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
use anise::errors::AlmanacResult;
use anise::prelude::{Frame, Orbit};
use hifitime::{Epoch, TimeUnits};
use indexmap::IndexSet;
use rand_pcg::Pcg64Mcg;
use snafu::ResultExt;

use crate::md::prelude::Traj;
use crate::od::groundpnt::GroundAsset;
use crate::od::interlink::InterlinkTxSpacecraft;
use crate::od::msr::MeasurementType;
use crate::od::prelude::{Measurement, ODError};
use crate::od::TrackingDevice;
use crate::od::{ODAlmanacSnafu, ODTrajSnafu};
use crate::State;

use std::sync::Arc;

impl TrackingDevice<GroundAsset> for InterlinkTxSpacecraft {
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
        traj: &Traj<GroundAsset>,
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
                let msr_t0_opt = self.measure_instantaneous(rx_0, None, almanac.clone())?;

                // End of integration time
                let msr_t1_opt = self.measure_instantaneous(rx_1, None, almanac.clone())?;

                if let Some(msr_t0) = msr_t0_opt {
                    if let Some(msr_t1) = msr_t1_opt {
                        // Line of sight in both cases

                        // Noises are computed at the midpoint of the integration time.
                        let noises = self.noises(epoch - integration_time * 0.5, rng)?;

                        let mut msr = Measurement::new(
                            "GroundAsset".to_string(),
                            epoch + noises[0].seconds(),
                        );

                        for (ii, msr_type) in self.measurement_types.iter().enumerate() {
                            let msr_value_0 = msr_t0.data[msr_type];
                            let msr_value_1 = msr_t1.data[msr_type];

                            let msr_value =
                                (msr_value_1 + msr_value_0) * 0.5 + noises[ii + 1] / 2.0_f64.sqrt();
                            msr.push(*msr_type, msr_value);
                        }

                        Ok(Some(msr))
                    } else {
                        Ok(None)
                    }
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

    /// Returns the Range and Doppler from pnt vehicle to ground asset (i.e. the opposed AER range and doppler data.)
    fn measure_instantaneous(
        &mut self,
        rx: GroundAsset,
        rng: Option<&mut Pcg64Mcg>,
        almanac: Arc<Almanac>,
    ) -> Result<Option<Measurement>, ODError> {
        let pnt_veh = self.traj.at(rx.epoch()).context(ODTrajSnafu {
            details: format!("fetching state {} for interlink", rx.epoch()),
        })?;

        let asset_loc = rx.to_location();

        let aer = almanac
            .azimuth_elevation_range_sez_from_location(pnt_veh.orbit, asset_loc, None, None)
            .context(ODAlmanacSnafu {
                action: "transforming receiver to transmitter frame",
            })?;

        if aer.elevation_above_mask_deg() < 0.0 {
            Ok(None)
        } else {
            let noises = self.noises(rx.epoch, rng)?;

            let mut msr =
                Measurement::new("GroundAsset".to_string(), rx.epoch + noises[0].seconds());

            for (ii, msr_type) in self.measurement_types.iter().enumerate() {
                let msr_value = match *msr_type {
                    MeasurementType::Range => -aer.range_km,
                    MeasurementType::Doppler => -aer.range_rate_km_s,
                    // Or return an error for unsupported types
                    _ => unreachable!("unsupported measurement type for interlink: {:?}", msr_type),
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
