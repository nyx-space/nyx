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
use std::convert::TryFrom;
use std::default::Default;

/// Defines the kind of correction to apply in the targeter
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Vary {
    /// Vary position component X in the integration frame
    PositionX,
    /// Vary position component Y in the integration frame
    PositionY,
    /// Vary position component Z in the integration frame
    PositionZ,
    /// Vary velocity component X in the integration frame
    VelocityX,
    /// Vary velocity component Y in the integration frame
    VelocityY,
    /// Vary velocity component Z in the integration frame
    VelocityZ,
    /// Vary velocity component V in the VNC frame
    VelocityV,
    /// Vary velocity component N in the VNC frame
    VelocityN,
    /// Vary velocity component C in the VNC frame
    VelocityC,
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
}

impl Variable {
    /// Returns whether the configuration of this variable is valid
    pub fn valid(&self) -> bool {
        if self.max_step < 0.0 {
            error!(
                "{:?}: max step is negative: {}",
                self.component, self.max_step
            );
            return false;
        }
        if self.max_value < 0.0 {
            error!(
                "{:?}: max value is negative: {}",
                self.component, self.max_value
            );
            return false;
        }
        if self.min_value > self.max_value {
            error!(
                "{:?}: min value is greater than max value: {} > {}",
                self.component, self.min_value, self.max_value
            );
            return false;
        }
        true
    }
}

impl Default for Variable {
    fn default() -> Self {
        Self {
            component: Vary::VelocityX,
            perturbation: 0.0001,
            init_guess: 0.0,
            max_step: 0.5,
            max_value: std::f64::INFINITY,
            min_value: std::f64::NEG_INFINITY,
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
            | Vary::VelocityZ
            | Vary::VelocityV
            | Vary::VelocityN
            | Vary::VelocityC => Ok(Self {
                component: vary,
                ..Default::default()
            }),
        }
    }
}
