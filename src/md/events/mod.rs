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

pub mod details;
use log::warn;
pub mod evaluators;
pub mod search;
use super::StateParameter;
use crate::errors::EventError;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::time::{Duration, Unit};
use crate::State;
use anise::prelude::{Almanac, Frame};
use anise::structure::planetocentric::ellipsoid::Ellipsoid;
use serde::{Deserialize, Serialize};

use std::default::Default;
use std::fmt;
use std::sync::Arc;

/// A trait to specify how a specific event must be evaluated.
pub trait EventEvaluator<S: State>: fmt::Display + Send + Sync
where
    DefaultAllocator: Allocator<S::Size> + Allocator<S::Size, S::Size> + Allocator<S::VecLength>,
{
    // Evaluation of event crossing, must return whether the condition happened between between both states.
    fn eval_crossing(
        &self,
        prev_state: &S,
        next_state: &S,
        almanac: Arc<Almanac>,
    ) -> Result<bool, EventError> {
        let prev = self.eval(prev_state, almanac.clone())?;
        let next = self.eval(next_state, almanac)?;

        Ok(prev * next < 0.0)
    }

    /// Evaluation of the event, must return a value corresponding to whether the state is before or after the event
    fn eval(&self, state: &S, almanac: Arc<Almanac>) -> Result<f64, EventError>;
    /// Returns a string representation of the event evaluation for the given state
    fn eval_string(&self, state: &S, almanac: Arc<Almanac>) -> Result<String, EventError>;
    fn epoch_precision(&self) -> Duration;
    fn value_precision(&self) -> f64;
}

/// Defines a state parameter event finder
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub struct Event {
    /// The state parameter
    pub parameter: StateParameter,
    /// The desired self.desired_value, must be in the same units as the state parameter
    pub desired_value: f64,
    /// The duration precision after which the solver will report that it cannot find any more precise
    pub epoch_precision: Duration,
    /// The precision on the desired value
    pub value_precision: f64,
    /// An optional frame in which to search this -- it IS recommended to convert the whole trajectory instead of searching in a given frame!
    pub obs_frame: Option<Frame>,
}

impl fmt::Display for Event {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.parameter)?;
        if self.parameter != StateParameter::Apoapsis && self.parameter != StateParameter::Periapsis
        {
            if self.desired_value.abs() > 1e3 {
                write!(
                    f,
                    " = {:e} {} (± {:e} {})",
                    self.desired_value,
                    self.parameter.unit(),
                    self.value_precision,
                    self.parameter.unit()
                )?;
            } else {
                write!(
                    f,
                    " = {} {} (± {} {})",
                    self.desired_value,
                    self.parameter.unit(),
                    self.value_precision,
                    self.parameter.unit()
                )?;
            }
        }
        if let Some(frame) = self.obs_frame {
            write!(f, "in frame {frame}")?;
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
        unit_precision: Unit,
    ) -> Self {
        Self {
            parameter,
            desired_value,
            epoch_precision: 1 * unit_precision,
            value_precision,
            obs_frame: None,
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

    /// Match the central body's mean equatorial radius.
    /// This is useful for detecting when an object might impact the central body.
    pub fn mean_surface(body: &Ellipsoid) -> Self {
        Self::new(StateParameter::Rmag, body.mean_equatorial_radius_km())
    }

    /// Match a specific event in another frame, using the default epoch precision and value.
    pub fn in_frame(parameter: StateParameter, desired_value: f64, target_frame: Frame) -> Self {
        warn!("Searching for an event in another frame is slow: you should instead convert the trajectory into that other frame");
        Self {
            parameter,
            desired_value,
            epoch_precision: Unit::Millisecond * 1,
            value_precision: 1e-3,
            obs_frame: Some(target_frame),
        }
    }
}

impl Default for Event {
    fn default() -> Self {
        Self {
            parameter: StateParameter::Periapsis,
            desired_value: 0.0,
            value_precision: 1e-3,
            epoch_precision: Unit::Second * 1,
            obs_frame: None,
        }
    }
}
