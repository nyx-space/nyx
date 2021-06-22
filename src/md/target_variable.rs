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

use super::NyxError;
use crate::cosmic::Frame;
use std::convert::TryFrom;
use std::default::Default;

/// Defines the kind of correction to apply in the targeter
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Vary {
    /// Vary position component X
    PositionX,
    /// Vary position component Y
    PositionY,
    /// Vary position component Z
    PositionZ,
    /// Vary velocity component X
    VelocityX,
    /// Vary velocity component Y
    VelocityY,
    /// Vary velocity component Z
    VelocityZ,
}

impl Vary {
    pub fn vec_index(&self) -> usize {
        match self {
            Vary::PositionX => 0,
            Vary::PositionY => 1,
            Vary::PositionZ => 2,
            Vary::VelocityX => 3,
            Vary::VelocityY => 4,
            Vary::VelocityZ => 5,
        }
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Variable {
    /// The component that will be varied in the targeter
    pub component: Vary,
    /// The perturbation for the finite differencing algorithm
    pub perturbation: f64,
    /// The initial guess of this variable
    pub init_guess: f64,
    /// The maximum step this variable may have between each iteration
    pub max_step: f64,
    /// The absolute maximum value this parameter can ever have
    pub max_value: f64,
    /// The absolute minimum value this parameter can ever have
    pub min_value: f64,
    /// The frame in which this variable should be applied, must be either a local frame or inertial
    pub frame: Option<Frame>,
}

impl Variable {
    /// Returns whether the configuration of this variable is valid
    pub fn valid(&self) -> Result<(), NyxError> {
        if self.max_step < 0.0 {
            let msg = format!(
                "{:?}: max step is negative: {}",
                self.component, self.max_step
            );
            error!("{}", msg);
            return Err(NyxError::TargetError(msg));
        }
        if self.max_value < 0.0 {
            let msg = format!(
                "{:?}: max value is negative: {}",
                self.component, self.max_value
            );
            error!("{}", msg);
            return Err(NyxError::TargetError(msg));
        }
        if self.min_value > self.max_value {
            let msg = format!(
                "{:?}: min value is greater than max value: {} > {}",
                self.component, self.min_value, self.max_value
            );
            error!("{}", msg);
            return Err(NyxError::TargetError(msg));
        }
        Ok(())
    }
}

impl Default for Variable {
    fn default() -> Self {
        Self {
            component: Vary::VelocityX,
            perturbation: 0.0001,
            init_guess: 0.0,
            max_step: 0.2,
            max_value: 5.0,
            min_value: -5.0,
            frame: None,
        }
    }
}

impl TryFrom<Vary> for Variable {
    type Error = NyxError;

    fn try_from(vary: Vary) -> Result<Self, Self::Error> {
        match vary {
            Vary::PositionX
            | Vary::PositionY
            | Vary::PositionZ
            | Vary::VelocityX
            | Vary::VelocityY
            | Vary::VelocityZ => Ok(Self {
                component: vary,
                ..Default::default()
            }),
        }
    }
}
