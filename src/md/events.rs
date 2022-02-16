/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::StateParameter;
use crate::cosmic::{Cosm, Frame, Orbit};
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, Vector3};
use crate::time::{Duration, Unit};
use crate::utils::between_pm_x;
use crate::{Spacecraft, State};
use std::default::Default;
use std::fmt;
use std::sync::Arc;

fn angled_value(cur_angle: f64, desired_angle: f64) -> f64 {
    if between_pm_x(cur_angle, desired_angle) > 0.0 {
        cur_angle - desired_angle
    } else {
        cur_angle + 2.0 * desired_angle
    }
}

/// A trait to specify how a specific event must be evaluated.
pub trait EventEvaluator<S: State>: fmt::Display + Send + Sync
where
    DefaultAllocator:
        Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size> + Allocator<f64, S::VecLength>,
{
    // Evaluation of event crossing, must return whether the condition happened between between both states.
    fn eval_crossing(&self, prev_state: &S, next_state: &S) -> bool {
        self.eval(prev_state) * self.eval(next_state) < 0.0
    }

    // Evaluation of the event, must return a value corresponding to whether the state is before or after the event
    fn eval(&self, state: &S) -> f64;

    fn epoch_precision(&self) -> Duration;
    fn value_precision(&self) -> f64;
}

/// Defines a state parameter event finder
#[derive(Clone, Debug)]
pub struct Event {
    /// The state parameter
    pub parameter: StateParameter,
    /// The desired self.desired_value, must be in the same units as the state parameter
    pub desired_value: f64,
    /// The time precision after which the solver will report that it cannot find any more precise
    pub epoch_precision: Unit,
    /// The precision on the desired value
    pub value_precision: f64,
    /// An optional frame in which to search this -- it IS recommended to convert the whole trajectory instead of searching in a given frame!
    pub in_frame: Option<(Frame, Arc<Cosm>)>,
}

impl fmt::Display for Event {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.parameter)?;
        if self.parameter != StateParameter::Apoapsis && self.parameter != StateParameter::Periapsis
        {
            if self.desired_value.abs() > 1e3 {
                write!(f, " = {:e}", self.desired_value)?;
            } else {
                write!(f, " = {}", self.desired_value,)?;
            }
        }
        if let Some((frame, _)) = self.in_frame {
            write!(f, "in frame {}", frame)?;
        }
        fmt::Result::Ok(())
    }
}

impl Event {
    /// Match a specific event for the parameter to hit the specified value.
    /// By default, the time precision is 1 millisecond and the value precision is 1e-3 of whatever
    /// unit is the default for that parameter. For example, a radius event will seek the requested
    /// value at the meter level, and an angle event will seek it at the thousands of a degree.
    pub fn new(parameter: StateParameter, desired_value: f64) -> Self {
        Self::within_tolerance(
            parameter,
            desired_value,
            parameter.default_event_precision(),
        )
    }

    /// Match a specific event for the parameter to hit the specified value with the provided tolerance on the value
    pub fn within_tolerance(
        parameter: StateParameter,
        desired_value: f64,
        value_precision: f64,
    ) -> Self {
        Self::specific(parameter, desired_value, value_precision, Unit::Millisecond)
    }

    /// Match a specific event for the parameter to hit the specified value with the provided tolerance on the value and time
    pub fn specific(
        parameter: StateParameter,
        desired_value: f64,
        value_precision: f64,
        epoch_precision: Unit,
    ) -> Self {
        Self {
            parameter,
            desired_value,
            epoch_precision,
            value_precision,
            in_frame: None,
        }
    }

    /// Match the periapasis i.e. True Anomaly == 0
    pub fn periapsis() -> Self {
        Self::new(StateParameter::Periapsis, 0.0)
    }

    /// Match the apoapasis i.e. True Anomaly == 180
    pub fn apoapsis() -> Self {
        Self::new(StateParameter::Apoapsis, 180.0)
    }

    /// Match a specific event in another frame, using the default epoch precision and value.
    pub fn in_frame(
        parameter: StateParameter,
        desired_value: f64,
        target_frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Self {
        warn!("Searching for an event in another frame is slow: you should instead convert the trajectory into that other frame");
        Self {
            parameter,
            desired_value,
            epoch_precision: Unit::Millisecond,
            value_precision: 1e-3,
            in_frame: Some((target_frame, cosm)),
        }
    }
}

impl Default for Event {
    fn default() -> Self {
        Self {
            parameter: StateParameter::Custom { mapping: 0 },
            desired_value: 0.0,
            value_precision: 1e-3,
            epoch_precision: Unit::Second,
            in_frame: None,
        }
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
            StateParameter::AoL
            | StateParameter::AoP
            | StateParameter::Declination
            | StateParameter::EccentricAnomaly
            | StateParameter::FlightPathAngle
            | StateParameter::HyperbolicAnomaly
            | StateParameter::Inclination
            | StateParameter::MeanAnomaly
            | StateParameter::RightAscension
            | StateParameter::RAAN
            | StateParameter::TrueAnomaly
            | StateParameter::TrueLongitude
            | StateParameter::VelocityDeclination => {
                angled_value(state.value(&self.parameter).unwrap(), self.desired_value)
            }
            StateParameter::Apoapsis => angled_value(state.ta(), 180.0),
            StateParameter::Periapsis => between_pm_x(state.ta(), 180.0),
            StateParameter::SlantAngle { x, y, z } => {
                let mut tgt = Vector3::new(x, y, z);
                tgt /= tgt.norm();
                angled_value(
                    state.r_hat().dot(&tgt).acos().to_degrees(),
                    self.desired_value,
                )
            }
            _ => state.value(&self.parameter).unwrap() - self.desired_value,
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
}
