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

use anise::prelude::Almanac;
use hifitime::Duration;
use snafu::ResultExt;
use std::sync::Arc;

use super::{Event, EventEvaluator};
use crate::errors::{EventAlmanacSnafu, EventError, EventPhysicsSnafu, EventStateSnafu};
use crate::md::StateParameter;
use crate::utils::between_pm_x;
use crate::{Spacecraft, State};

pub(crate) fn angled_value(cur_angle: f64, desired_angle: f64) -> f64 {
    if between_pm_x(cur_angle, desired_angle) > 0.0 {
        cur_angle - desired_angle
    } else {
        cur_angle + 2.0 * desired_angle
    }
}

impl EventEvaluator<Spacecraft> for Event {
    fn eval(&self, state: &Spacecraft, almanac: Arc<Almanac>) -> Result<f64, EventError> {
        let state = if let Some(frame) = self.obs_frame {
            if state.orbit.frame == frame {
                *state
            } else {
                state.with_orbit(
                    almanac
                        .transform_to(state.orbit, frame, None)
                        .context(EventAlmanacSnafu)?,
                )
            }
        } else {
            *state
        };

        // Return the parameter centered around the desired value
        match self.parameter {
            StateParameter::Apoapsis => Ok(angled_value(
                state.orbit.ta_deg().context(EventPhysicsSnafu)?,
                180.0,
            )),
            StateParameter::Periapsis => Ok(between_pm_x(
                state.orbit.ta_deg().context(EventPhysicsSnafu)?,
                180.0,
            )),
            StateParameter::FuelMass => Ok(state.fuel_mass_kg - self.desired_value),
            _ => Ok(state.value(self.parameter).context(EventStateSnafu {
                param: self.parameter,
            })? - self.desired_value),
        }
    }

    #[allow(clippy::identity_op)]
    fn epoch_precision(&self) -> Duration {
        1 * self.epoch_precision
    }

    fn value_precision(&self) -> f64 {
        self.value_precision
    }

    fn eval_string(
        &self,
        state: &Spacecraft,
        _almanac: Arc<Almanac>,
    ) -> Result<String, EventError> {
        match self.parameter {
            StateParameter::Apoapsis | StateParameter::Periapsis => {
                Ok(format!("{}", self.parameter))
            }
            _ => {
                let unit = if self.parameter.unit().is_empty() {
                    String::new()
                } else {
                    format!(" ({})", self.parameter.unit())
                };
                let val = state.value(self.parameter).context(EventStateSnafu {
                    param: self.parameter,
                })?;

                Ok(format!("{}{} = {:.3}{}", self.parameter, unit, val, unit))
            }
        }
    }
}
