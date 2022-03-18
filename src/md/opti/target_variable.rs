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

use crate::cosmic::Frame;
use crate::errors::TargetingError;
use std::default::Default;
use std::f64::consts::{FRAC_PI_2, FRAC_PI_8, PI};
use std::fmt;

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
    MnvrDelta,
    /// Maneuver's out-of-plane beta dot
    MnvrDeltaDot,
    /// Maneuver's out-of-plane beta double-dot
    MnvrDeltaDDot,
    /// Start epoch difference in seconds
    StartEpoch,
    /// Burn duration difference in seconds
    Duration,
    /// End epoch difference in seconds
    EndEpoch,
    /// Thrust direction in X
    ThrustX,
    /// Thrust direction in Y
    ThrustY,
    /// Thrust direction in Z
    ThrustZ,
    /// Thrust level during the burn.
    ThrustLevel,
    /// Thrust direction rate in X
    ThrustRateX,
    /// Thrust direction rate in Y
    ThrustRateY,
    /// Thrust direction rate in Z
    ThrustRateZ,
    /// Thrust direction acceleration in X
    ThrustAccelX,
    /// Thrust direction acceleration in Y
    ThrustAccelY,
    /// Thrust direction acceleration in Z
    ThrustAccelZ,
}

impl Vary {
    pub fn is_finite_burn(&self) -> bool {
        *self == Self::MnvrAlpha
            || *self == Self::MnvrAlphaDDot
            || *self == Self::MnvrAlphaDDot
            || *self == Self::MnvrDelta
            || *self == Self::MnvrDeltaDDot
            || *self == Self::MnvrDeltaDDot
            || *self == Self::StartEpoch
            || *self == Self::Duration
            || *self == Self::EndEpoch
            || *self == Self::ThrustX
            || *self == Self::ThrustY
            || *self == Self::ThrustZ
            || *self == Self::ThrustLevel
            || *self == Self::ThrustRateX
            || *self == Self::ThrustRateY
            || *self == Self::ThrustRateZ
            || *self == Self::ThrustAccelX
            || *self == Self::ThrustAccelY
            || *self == Self::ThrustAccelZ
    }

    pub fn vec_index(&self) -> usize {
        match self {
            Self::PositionX | Self::ThrustX | Self::MnvrAlphaDDot | Self::MnvrDeltaDDot => 0,
            Self::PositionY | Self::ThrustY | Self::MnvrAlphaDot | Self::MnvrDeltaDot => 1,
            Self::PositionZ | Self::ThrustZ | Self::MnvrAlpha | Self::MnvrDelta => 2,
            Self::VelocityX | Self::ThrustRateX => 3,
            Self::VelocityY | Self::ThrustRateY => 4,
            Self::VelocityZ | Self::ThrustRateZ => 5,
            Self::StartEpoch | Self::ThrustAccelX => 6,
            Self::Duration | Self::EndEpoch | Self::ThrustAccelY => 7,
            Self::ThrustAccelZ => 8,
            _ => unreachable!(),
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

    pub fn with_pert(self, pert: f64) -> Self {
        let mut me = self;
        me.perturbation = pert;
        me
    }

    pub fn with_min(self, min_val: f64) -> Self {
        let mut me = self;
        me.min_value = min_val;
        me
    }

    pub fn with_max(self, max_val: f64) -> Self {
        let mut me = self;
        me.max_value = max_val;
        me
    }

    /// Ensure that `val` is within the variable bounds
    pub fn apply_bounds(&self, val: f64) -> f64 {
        if val > self.max_value {
            self.max_value
        } else if val < self.min_value {
            self.min_value
        } else {
            val
        }
    }

    /// Ensure that `val` is within the variable bounds
    pub fn ensure_bounds(&self, val: &mut f64) {
        *val = self.check_bounds(*val).0;
    }

    /// Returns the input value unless it is out of bounds, then it returns the bound, and whether the input value was OK
    pub fn check_bounds(&self, val: f64) -> (f64, bool) {
        if val > self.max_value {
            (self.max_value, false)
        } else if val < self.min_value {
            (self.min_value, false)
        } else {
            (val, true)
        }
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
            Vary::MnvrAlpha | Vary::MnvrAlphaDot | Vary::MnvrAlphaDDot => Self {
                component: vary,
                perturbation: 0.1 * PI,
                max_step: FRAC_PI_8,
                max_value: PI,
                min_value: 0.0,
                ..Default::default()
            },
            Vary::MnvrDelta | Vary::MnvrDeltaDot | Vary::MnvrDeltaDDot => Self {
                component: vary,
                perturbation: 0.1 * PI,
                max_step: FRAC_PI_8,
                max_value: FRAC_PI_2,
                min_value: -FRAC_PI_2,
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
            Vary::ThrustX | Vary::ThrustY | Vary::ThrustZ => Self {
                component: vary,
                max_value: 1.0,
                min_value: -1.0,
                ..Default::default()
            },
            Vary::ThrustRateX | Vary::ThrustRateY | Vary::ThrustRateZ => Self {
                component: vary,
                perturbation: 1e-10,
                max_value: 1.0,
                min_value: -1.0,
                ..Default::default()
            },
            Vary::ThrustAccelX | Vary::ThrustAccelY | Vary::ThrustAccelZ => Self {
                component: vary,
                perturbation: 1e-15,
                max_value: 1.0,
                min_value: -1.0,
                ..Default::default()
            },
            Vary::ThrustLevel => Self {
                component: vary,
                perturbation: -0.0001, // Perturb the thrust by -1%
                min_value: 0.0001,
                max_value: 1.0,
                init_guess: 1.0,
                ..Default::default()
            },
        }
    }
}

impl fmt::Display for Variable {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}{:?} = {} ± {:} ∈ [{}; {}]",
            match self.frame {
                Some(f) => format!("{}", f),
                None => "".to_string(),
            },
            self.component,
            self.init_guess,
            self.perturbation,
            self.min_value,
            self.max_value
        )
    }
}
