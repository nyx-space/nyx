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

use hifitime::Duration;

use super::{Event, EventEvaluator};
use crate::cosmic::Orbit;
use crate::md::StateParameter;
use crate::utils::between_pm_x;
use crate::{Spacecraft, State};

fn angled_value(cur_angle: f64, desired_angle: f64) -> f64 {
    if between_pm_x(cur_angle, desired_angle) > 0.0 {
        cur_angle - desired_angle
    } else {
        cur_angle + 2.0 * desired_angle
    }
}

impl EventEvaluator<Orbit> for Event {
    #[allow(clippy::identity_op)]
    fn epoch_precision(&self) -> Duration {
        1 * self.epoch_precision
    }

    fn value_precision(&self) -> f64 {
        self.value_precision
    }

    fn eval(&self, state: &Orbit) -> f64 {
        // Transform the state if needed
        let state = if let Some((frame, cosm)) = &self.in_frame {
            if state.frame == *frame {
                *state
            } else {
                cosm.frame_chg(state, *frame)
            }
        } else {
            *state
        };

        // Return the parameter centered around the desired value
        match self.parameter {
            StateParameter::Apoapsis => angled_value(state.ta_deg(), 180.0),
            StateParameter::Periapsis => between_pm_x(state.ta_deg(), 180.0),
            _ => state.value(self.parameter).unwrap() - self.desired_value,
        }
    }

    fn eval_string(&self, state: &Orbit) -> String {
        match self.parameter {
            StateParameter::Apoapsis | StateParameter::Periapsis => {
                format!("{}", self.parameter)
            }
            _ => {
                let unit = if self.parameter.unit().is_empty() {
                    String::new()
                } else {
                    format!(" ({})", self.parameter.unit())
                };

                let val = state.value(self.parameter).unwrap();
                format!("{}{} = {:.3}{}", self.parameter, unit, val, unit)
            }
        }
    }
}

impl EventEvaluator<Spacecraft> for Event {
    fn eval(&self, state: &Spacecraft) -> f64 {
        match self.parameter {
            StateParameter::FuelMass => state.fuel_mass_kg - self.desired_value,
            _ => self.eval(&state.orbit),
        }
    }

    #[allow(clippy::identity_op)]
    fn epoch_precision(&self) -> Duration {
        1 * self.epoch_precision
    }

    fn value_precision(&self) -> f64 {
        self.value_precision
    }

    fn eval_string(&self, state: &Spacecraft) -> String {
        match self.parameter {
            StateParameter::Apoapsis | StateParameter::Periapsis => {
                format!("{}", self.parameter)
            }
            _ => {
                let unit = if self.parameter.unit().is_empty() {
                    String::new()
                } else {
                    format!(" ({})", self.parameter.unit())
                };
                let val = state.value(self.parameter).unwrap();
                format!("{}{} = {:.3}{}", self.parameter, unit, val, unit)
            }
        }
    }
}
