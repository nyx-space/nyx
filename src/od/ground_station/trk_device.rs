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
use crate::md::prelude::{Interpolatable, Traj};
use crate::od::msr::measurement::Measurement;
use crate::od::msr::MeasurementType;
use crate::time::Epoch;
use crate::Spacecraft;
use anise::errors::AlmanacResult;
use anise::frames::Frame;
use anise::prelude::{Almanac, Orbit};
use hifitime::TimeUnits;
use indexmap::IndexSet;
use rand_pcg::Pcg64Mcg;
use snafu::ResultExt;
use std::sync::Arc;

use super::GroundStation;

impl TrackingDevice<Spacecraft> for GroundStation {
    fn measurement_types(&self) -> &IndexSet<MeasurementType> {
        &self.measurement_types
    }

    /// Perform a measurement from the ground station to the receiver (rx).
    fn measure(
        &mut self,
        epoch: Epoch,
        traj: &Traj<Spacecraft>,
        rng: Option<&mut Pcg64Mcg>,
        almanac: Arc<Almanac>,
    ) -> Result<Option<Measurement>, ODError> {
        match self.integration_time {
            Some(integration_time) => {
                // If out of traj bounds, return None.
                let rx_0 = traj.at(epoch - integration_time).context(ODTrajSnafu)?;
                let rx_1 = traj.at(epoch).context(ODTrajSnafu)?;

                let obstructing_body = if !self.frame.ephem_origin_match(rx_0.frame()) {
                    Some(rx_0.frame())
                } else {
                    None
                };

                let aer_t0 = self
                    .azimuth_elevation_of(rx_0.orbit, obstructing_body, &almanac)
                    .context(ODAlmanacSnafu {
                        action: "computing AER",
                    })?;
                let aer_t1 = self
                    .azimuth_elevation_of(rx_1.orbit, obstructing_body, &almanac)
                    .context(ODAlmanacSnafu {
                        action: "computing AER",
                    })?;

                if aer_t0.elevation_deg < self.elevation_mask_deg
                    || aer_t1.elevation_deg < self.elevation_mask_deg
                {
                    debug!(
                        "{} (el. mask {:.3} deg) but object moves from {:.3} to {:.3} deg -- no measurement",
                        self.name, self.elevation_mask_deg, aer_t0.elevation_deg, aer_t1.elevation_deg
                    );
                    return Ok(None);
                }

                // Noises are computed at the midpoint of the integration time.
                let noises = self.noises(epoch - integration_time * 0.5, rng)?;

                let mut msr = Measurement::new(self.name.clone(), epoch + noises[0].seconds());

                for (ii, msr_type) in self.measurement_types.iter().enumerate() {
                    let msr_value = msr_type.compute_two_way(aer_t0, aer_t1, noises[ii + 1])?;
                    msr.push(*msr_type, msr_value);
                }

                Ok(Some(msr))
            }
            None => self.measure_instantaneous(traj.at(epoch).context(ODTrajSnafu)?, rng, almanac),
        }
    }

    fn name(&self) -> String {
        self.name.clone()
    }

    fn location(&self, epoch: Epoch, frame: Frame, almanac: Arc<Almanac>) -> AlmanacResult<Orbit> {
        almanac.transform_to(self.to_orbit(epoch, &almanac).unwrap(), frame, None)
    }

    fn measure_instantaneous(
        &mut self,
        rx: Spacecraft,
        rng: Option<&mut Pcg64Mcg>,
        almanac: Arc<Almanac>,
    ) -> Result<Option<Measurement>, ODError> {
        let obstructing_body = if !self.frame.ephem_origin_match(rx.frame()) {
            Some(rx.frame())
        } else {
            None
        };

        let aer = self
            .azimuth_elevation_of(rx.orbit, obstructing_body, &almanac)
            .context(ODAlmanacSnafu {
                action: "computing AER",
            })?;

        if aer.elevation_deg >= self.elevation_mask_deg {
            // Only update the noises if the measurement is valid.
            let noises = self.noises(rx.orbit.epoch, rng)?;

            let mut msr = Measurement::new(self.name.clone(), rx.orbit.epoch + noises[0].seconds());

            for (ii, msr_type) in self.measurement_types.iter().enumerate() {
                let msr_value = msr_type.compute_one_way(aer, noises[ii + 1])?;
                msr.push(*msr_type, msr_value);
            }

            Ok(Some(msr))
        } else {
            debug!(
                "{} {} (el. mask {:.3} deg), object at {:.3} deg -- no measurement",
                self.name, rx.orbit.epoch, self.elevation_mask_deg, aer.elevation_deg
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
        Ok(match msr_type {
            MeasurementType::Range => self
                .range_noise_km
                .ok_or(ODError::NoiseNotConfigured { kind: "Range" })?
                .covariance(epoch),
            MeasurementType::Doppler => self
                .doppler_noise_km_s
                .ok_or(ODError::NoiseNotConfigured { kind: "Doppler" })?
                .covariance(epoch),
        })
    }
}
