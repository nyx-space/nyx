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

use crate::cosmic::Frame;
use crate::errors::TargetingError;
use std::default::Default;

/// Defines the kind of correction to apply in the targeter
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
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
    /// Maneuver's in-plane alpha
    MnvrAlpha,
    /// Maneuver's in-plane alpha-dot
    MnvrAlphaDot,
    /// Maneuver's in-plane alpha double-dot
    MnvrAlphaDDot,
    /// Maneuver's out-of-plane beta
    MnvrBeta,
    /// Maneuver's out-of-plane beta dot
    MnvrBetaDot,
    /// Maneuver's out-of-plane beta double-dot
    MnvrBetaDDot,
    /// Start epoch difference in seconds
    StartEpoch,
    /// Burn duration difference in seconds
    Duration,
    /// End epoch difference in seconds
    EndEpoch,
    /// Thrust direction in X
    Tx,
    /// Thrust direction in Y
    Ty,
    /// Thrust direction in Z
    Tz,
    /// Thrust level during the burn
    ThrustLevel,
}

impl Vary {
    pub fn is_finite_burn(&self) -> bool {
        *self == Self::MnvrAlpha
            || *self == Self::MnvrAlphaDDot
            || *self == Self::MnvrAlphaDDot
            || *self == Self::MnvrBeta
            || *self == Self::MnvrBetaDDot
            || *self == Self::MnvrBetaDDot
            || *self == Self::StartEpoch
            || *self == Self::Duration
            || *self == Self::EndEpoch
            || *self == Self::Tx
            || *self == Self::Ty
            || *self == Self::Tz
            || *self == Self::ThrustLevel
    }

    pub fn vec_index(&self) -> usize {
        match self {
            Self::PositionX | Self::Tx | Self::MnvrAlphaDDot | Self::MnvrBetaDDot => 0,
            Self::PositionY | Self::Ty | Self::MnvrAlphaDot | Self::MnvrBetaDot => 1,
            Self::PositionZ | Self::Tz | Self::MnvrAlpha | Self::MnvrBeta => 2,
            Self::VelocityX | Self::ThrustLevel => 3,
            Self::VelocityY => 4,
            Self::VelocityZ => 5,
            Self::StartEpoch => 6,
            Self::Duration | Self::EndEpoch => 7,
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
    pub fn valid(&self) -> Result<(), TargetingError> {
        if self.max_step < 0.0 {
            let msg = format!(
                "{:?}: max step is negative: {}",
                self.component, self.max_step
            );
            error!("{}", msg);
            return Err(TargetingError::VariableError(msg));
        }
        if self.max_value < 0.0 {
            let msg = format!(
                "{:?}: max value is negative: {}",
                self.component, self.max_value
            );
            error!("{}", msg);
            return Err(TargetingError::VariableError(msg));
        }
        if self.min_value > self.max_value {
            let msg = format!(
                "{:?}: min value is greater than max value: {} > {}",
                self.component, self.min_value, self.max_value
            );
            error!("{}", msg);
            return Err(TargetingError::VariableError(msg));
        }
        Ok(())
    }

    pub fn with_initial_guess(self, guess: f64) -> Self {
        let mut me = self;
        me.init_guess = guess;
        me
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

impl From<Vary> for Variable {
    fn from(vary: Vary) -> Self {
        match vary {
            Vary::PositionX
            | Vary::PositionY
            | Vary::PositionZ
            | Vary::VelocityX
            | Vary::VelocityY
            | Vary::VelocityZ => Self {
                component: vary,
                ..Default::default()
            },
            Vary::MnvrAlpha | Vary::MnvrBeta => Self {
                component: vary,
                perturbation: 0.01_f64.to_radians(),
                max_step: 45.0_f64.to_radians(),
                max_value: 360.0_f64.to_radians(),
                min_value: -360.0_f64.to_radians(),
                ..Default::default()
            },
            Vary::MnvrAlphaDot | Vary::MnvrBetaDot => Self {
                component: vary,
                perturbation: 0.01_f64.to_radians(),
                max_step: 45.0_f64.to_radians(),
                max_value: 360.0_f64.to_radians(),
                min_value: -360.0_f64.to_radians(),
                ..Default::default()
            },
            Vary::MnvrAlphaDDot | Vary::MnvrBetaDDot => Self {
                component: vary,
                perturbation: 0.01_f64.to_radians(),
                max_step: 45.0_f64.to_radians(),
                max_value: 360.0_f64.to_radians(),
                min_value: -360.0_f64.to_radians(),
                ..Default::default()
            },
            Vary::StartEpoch | Vary::EndEpoch => Self {
                component: vary,
                perturbation: 0.5,
                max_step: 60.0,
                max_value: 600.0,
                min_value: -600.0,
                ..Default::default()
            },
            Vary::Duration => Self {
                component: vary,
                perturbation: 1.0,
                max_step: 60.0,
                max_value: 600.0,
                min_value: 0.0,
                ..Default::default()
            },
            Vary::Tx | Vary::Ty | Vary::Tz => Self {
                component: vary,
                max_value: 1.0,
                min_value: -1.0,
                ..Default::default()
            },
            Vary::ThrustLevel => Self {
                component: vary,
                min_value: 0.0,
                max_value: 1.0,
                init_guess: 1.0,
                ..Default::default()
            },
        }
    }
}
