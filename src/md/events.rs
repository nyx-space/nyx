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

use super::StateParameter;
use crate::celestia::{Cosm, Frame, Orbit};
use crate::dimensions::allocator::Allocator;
use crate::dimensions::DefaultAllocator;
use crate::time::{Duration, TimeUnit};
use crate::utils::between_pm_x;
use crate::{SpacecraftState, State};
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
    DefaultAllocator: Allocator<f64, S::Size>,
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
    pub epoch_precision: TimeUnit,
    /// The precision on the desired value
    pub value_precision: f64,
    /// An optional frame in which to search this -- it IS recommended to convert the whole trajectory instead of searching in a given frame!
    pub in_frame: Option<(Frame, Arc<Cosm>)>,
}

impl fmt::Display for Event {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.parameter)?;
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
        Self {
            parameter,
            desired_value,
            epoch_precision: TimeUnit::Millisecond,
            value_precision: 1e-3,
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
            epoch_precision: TimeUnit::Millisecond,
            value_precision: 1e-3,
            in_frame: Some((target_frame, cosm)),
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
            StateParameter::AoL => angled_value(state.aol(), self.desired_value),
            StateParameter::AoP => angled_value(state.aop(), self.desired_value),
            StateParameter::Apoapsis => angled_value(state.ta(), 180.0),
            StateParameter::Declination => angled_value(state.declination(), self.desired_value),
            StateParameter::ApoapsisRadius => state.apoapsis() - self.desired_value,
            StateParameter::EccentricAnomaly => angled_value(state.ea(), self.desired_value),
            StateParameter::Eccentricity => state.ecc() - self.desired_value,
            StateParameter::Energy => state.energy() - self.desired_value,
            StateParameter::GeodeticHeight => state.geodetic_height() - self.desired_value,
            StateParameter::GeodeticLatitude => state.geodetic_latitude() - self.desired_value,
            StateParameter::GeodeticLongitude => state.geodetic_longitude() - self.desired_value,
            StateParameter::Hmag => state.hmag() - self.desired_value,
            StateParameter::HX => state.hx() - self.desired_value,
            StateParameter::HY => state.hy() - self.desired_value,
            StateParameter::HZ => state.hz() - self.desired_value,
            StateParameter::Inclination => angled_value(state.inc(), self.desired_value),
            StateParameter::MeanAnomaly => angled_value(state.ma(), self.desired_value),
            StateParameter::Periapsis => between_pm_x(state.ta(), 180.0),
            StateParameter::PeriapsisRadius => state.periapsis() - self.desired_value,
            StateParameter::Period => state.period().in_seconds() - self.desired_value,
            StateParameter::RightAscension => {
                angled_value(state.right_ascension(), self.desired_value)
            }
            StateParameter::RAAN => angled_value(state.raan(), self.desired_value),
            StateParameter::Rmag => state.rmag() - self.desired_value,
            StateParameter::SemiParameter => state.semi_parameter() - self.desired_value,
            StateParameter::SemiMinorAxis => state.semi_minor_axis() - self.desired_value,
            StateParameter::SMA => state.sma() - self.desired_value,
            StateParameter::TrueAnomaly => angled_value(state.ta(), self.desired_value),
            StateParameter::TrueLongitude => angled_value(state.tlong(), self.desired_value),
            StateParameter::Vmag => state.vmag() - self.desired_value,
            StateParameter::X => state.x - self.desired_value,
            StateParameter::Y => state.y - self.desired_value,
            StateParameter::Z => state.z - self.desired_value,
            StateParameter::VX => state.vx - self.desired_value,
            StateParameter::VY => state.vy - self.desired_value,
            StateParameter::VZ => state.vz - self.desired_value,
            _ => unimplemented!(),
        }
    }
}

impl EventEvaluator<SpacecraftState> for Event {
    fn eval(&self, state: &SpacecraftState) -> f64 {
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

/// Computes the elevation between a body fixed point and something else.
#[derive(Clone, Debug)]
pub struct ElevationEvent {}
