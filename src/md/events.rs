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

#[derive(Clone, Debug)]
pub struct Event {
    pub parameter: StateParameter,
    pub epoch_precision: TimeUnit,
    pub value_precision: f64,
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
    pub fn new(parameter: StateParameter) -> Self {
        Self {
            parameter,
            epoch_precision: TimeUnit::Millisecond,
            value_precision: 1e-3,
            in_frame: None,
        }
    }

    /// Allows only reporting the current value
    pub fn status(parameter: StateParameter) -> Self {
        Self {
            parameter,
            epoch_precision: TimeUnit::Millisecond,
            value_precision: 0.0,
            in_frame: None,
        }
    }

    /// Match a specific event in another frame, using the default epoch precision and value.
    pub fn in_frame(parameter: StateParameter, target_frame: Frame, cosm: Arc<Cosm>) -> Self {
        Self {
            parameter,
            epoch_precision: TimeUnit::Millisecond,
            value_precision: 1e-3,
            in_frame: Some((target_frame, cosm)),
        }
    }

    /// Report the status in a specific frame
    pub fn status_in_frame(
        parameter: StateParameter,
        target_frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Self {
        Self {
            parameter,
            epoch_precision: TimeUnit::Millisecond,
            value_precision: 0.0,
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
            StateParameter::AoL(value) => angled_value(state.aol(), value),
            StateParameter::AoP(value) => angled_value(state.aop(), value),
            StateParameter::Apoapsis => angled_value(state.ta(), 180.0),
            StateParameter::Declination(value) => angled_value(state.declination(), value),
            StateParameter::ApoapsisRadius(value) => state.apoapsis() - value,
            StateParameter::EccentricAnomaly(value) => angled_value(state.ea(), value),
            StateParameter::Eccentricity(value) => state.ecc() - value,
            StateParameter::Energy(value) => state.energy() - value,
            StateParameter::GeodeticHeight(value) => state.geodetic_height() - value,
            StateParameter::GeodeticLatitude(value) => state.geodetic_latitude() - value,
            StateParameter::GeodeticLongitude(value) => state.geodetic_longitude() - value,
            StateParameter::Hmag(value) => state.hmag() - value,
            StateParameter::HX(value) => state.hx() - value,
            StateParameter::HY(value) => state.hy() - value,
            StateParameter::HZ(value) => state.hz() - value,
            StateParameter::Inclination(value) => angled_value(state.inc(), value),
            StateParameter::MeanAnomaly(value) => angled_value(state.ma(), value),
            StateParameter::Periapsis => between_pm_x(state.ta(), 180.0),
            StateParameter::PeriapsisRadius(value) => state.periapsis() - value,
            StateParameter::Period(value) => state.period().in_seconds() - value,
            StateParameter::RightAscension(value) => angled_value(state.right_ascension(), value),
            StateParameter::RAAN(value) => angled_value(state.raan(), value),
            StateParameter::Rmag(value) => state.rmag() - value,
            StateParameter::SemiParameter(value) => state.semi_parameter() - value,
            StateParameter::SemiMinorAxis(value) => state.semi_minor_axis() - value,
            StateParameter::SMA(value) => state.sma() - value,
            StateParameter::TrueAnomaly(value) => angled_value(state.ta(), value),
            StateParameter::TrueLongitude(value) => angled_value(state.tlong(), value),
            StateParameter::Vmag(value) => state.vmag() - value,
            StateParameter::X(value) => state.x - value,
            StateParameter::Y(value) => state.y - value,
            StateParameter::Z(value) => state.z - value,
            StateParameter::VX(value) => state.vx - value,
            StateParameter::VY(value) => state.vy - value,
            StateParameter::VZ(value) => state.vz - value,
            _ => unimplemented!(),
        }
    }
}

impl EventEvaluator<SpacecraftState> for Event {
    fn eval(&self, state: &SpacecraftState) -> f64 {
        match self.parameter {
            StateParameter::FuelMass(value) => state.fuel_mass_kg - value,
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

#[derive(Clone, Debug)]
pub struct EclipseEvent {}
