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
use crate::cosmic::OrbitPartial;
use std::fmt;

/// Defines a state parameter event finder
#[derive(Copy, Clone, Debug)]
pub struct Objective {
    /// The state parameter to target
    pub parameter: StateParameter,
    /// The desired self.desired_value, must be in the same units as the state parameter
    pub desired_value: f64,
    /// The precision on the desired value
    pub tolerance: f64,
    /// A multiplicative factor this parameter's error in the targeting (defaults to 1.0)
    pub multiplicative_factor: f64,
    /// An additive factor to this parameters's error in the targeting (defaults to 0.0)
    pub additive_factor: f64,
}

impl Objective {
    /// Match a specific value for the parameter.
    /// By default, the tolerance on the parameter is 0.1 times whatever unit is the default for that parameter.
    /// For example, a radius event will seek the requested value at the decimeter level, and an angle event will seek it at the tenth of a degree.
    pub fn new(parameter: StateParameter, desired_value: f64) -> Self {
        Self::within_tolerance(
            parameter,
            desired_value,
            parameter.default_event_precision(),
        )
    }

    /// Match a specific value for the parameter to hit the specified value with the provided tolerance on the value
    pub fn within_tolerance(parameter: StateParameter, desired_value: f64, tolerance: f64) -> Self {
        Self {
            parameter,
            desired_value,
            tolerance,
            multiplicative_factor: 1.0,
            additive_factor: 0.0,
        }
    }

    /// Returns whether this objective has been achieved, and the associated parameter error.
    pub fn assess(&self, achieved: OrbitPartial) -> (bool, f64) {
        self.assess_raw(achieved.real())
    }

    /// Returns whether this objective has been achieved, and the associated parameter error.
    /// Warning: the parameter `achieved` must be in the same unit as the objective.
    pub fn assess_raw(&self, achieved: f64) -> (bool, f64) {
        let param_err =
            self.multiplicative_factor * (self.desired_value - achieved) + self.additive_factor;

        (param_err.abs() <= self.tolerance, param_err)
    }
}

impl fmt::Display for Objective {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "\t{:x}", self)
    }
}

impl fmt::LowerHex for Objective {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        let max_obj_tol = self.tolerance.log10().abs().ceil() as usize;

        write!(
            f,
            "{:?} → {:.prec$} ",
            self.parameter,
            self.desired_value,
            prec = max_obj_tol
        )?;

        if self.tolerance.abs() < 1e-1 {
            write!(f, "(± {:.1e})", self.tolerance)
        } else {
            write!(f, " (± {:.2})", self.tolerance)
        }
    }
}
