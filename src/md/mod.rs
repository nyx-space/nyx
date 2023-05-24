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

use crate::errors::NyxError;
use crate::{Orbit, Spacecraft};
use std::error::Error;
use std::fmt;

pub mod prelude {
    pub use super::{
        optimizer::*,
        trajectory::{ExportCfg, Interpolatable, Traj},
        Ephemeris, Event, ScTraj, StateParameter,
    };
    pub use crate::cosmic::{
        try_achieve_b_plane, BPlane, BPlaneTarget, Bodies, Cosm, Frame, GuidanceMode,
        LightTimeCalc, Orbit, OrbitDual,
    };
    pub use crate::dynamics::{
        Drag, Harmonics, OrbitalDynamics, PointMasses, SolarPressure, SpacecraftDynamics,
    };
    pub use crate::dynamics::{Dynamics, NyxError};
    pub use crate::io::gravity::HarmonicsMem;
    pub use crate::md::objective::Objective;
    pub use crate::propagators::{PropOpts, Propagator};
    pub use crate::time::{Duration, Epoch, TimeUnits, Unit};
    pub use crate::Spacecraft;
    pub use crate::{State, TimeTagged};
    pub use std::sync::Arc;
}

pub mod trajectory;

mod events;
pub use events::{Event, EventEvaluator};

pub mod objective;
pub mod opti;
pub use opti::optimizer;
pub type ScTraj = trajectory::Traj<Spacecraft>;
pub type Ephemeris = trajectory::Traj<Orbit>;

mod param;
pub use param::StateParameter;

pub use opti::target_variable::{Variable, Vary};

#[allow(clippy::result_large_err)]
#[derive(Clone, PartialEq, Debug)]
pub enum TargetingError {
    /// Raised if the variables to be adjusted lead to an under-determined of the problem for the targeter
    UnderdeterminedProblem,
    /// Raised if the variables of the problem are incorrectly configured
    VariableError(String),
    /// Raised in case of a targeting frame error
    FrameError(String),
    /// Raised if the targeted was configured with a variable that isn't supported (e.g. a maneuver alpha variable in the multiple shooting)
    UnsupportedVariable(Variable),
    /// Raised when the verification of a solution has failed
    Verification(String),
}

impl fmt::Display for TargetingError {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::UnderdeterminedProblem => write!(f, "The variables to be adjusted lead to an under-determined of the problem for the targeter"),
            Self::VariableError(e) => write!(f, "Incorrectly configured variable: {e}"),
            Self::FrameError(e) => write!(f, "Frame error in targeter: {e}"),
            Self::UnsupportedVariable(v) => write!(f, "Unsupported variable in problem: {v:?}"),
            Self::Verification(e) => write!(f, "Verification of targeting solution failed: {e}")
        }
    }
}

impl Error for TargetingError {}

impl From<TargetingError> for NyxError {
    fn from(e: TargetingError) -> Self {
        NyxError::Targeter(Box::new(e))
    }
}
