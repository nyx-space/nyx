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

use crate::cosmic::AstroError;
use crate::dynamics::guidance::GuidanceErrors;
use crate::errors::NyxError;
use crate::propagators::PropagationError;
use crate::{Orbit, Spacecraft};
use snafu::prelude::*;

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

pub(crate) mod events;
pub use events::{Event, EventEvaluator};

pub mod objective;
pub mod opti;
pub use opti::optimizer;
pub type ScTraj = trajectory::Traj<Spacecraft>;
pub type Ephemeris = trajectory::Traj<Orbit>;

mod param;
pub use param::StateParameter;

pub use opti::target_variable::{Variable, Vary};

use self::trajectory::TrajError;

#[allow(clippy::result_large_err)]
#[derive(Clone, PartialEq, Debug, Snafu)]
#[snafu(visibility(pub(crate)))]
pub enum TargetingError {
    #[snafu(display(
        "The variables to be adjusted lead to an under-determined of the problem for the targeter"
    ))]
    UnderdeterminedProblem,
    /// Raised if the variables of the problem are incorrectly configured
    #[snafu(display("Incorrectly configured variable: {msg}"))]
    VariableError { msg: String },
    /// Raised in case of a targeting frame error
    #[snafu(display("Frame error in targeter: {msg}"))]
    FrameError { msg: String },
    /// Raised if the targeted was configured with a variable that isn't supported (e.g. a maneuver alpha variable in the multiple shooting)
    #[snafu(display("Unsupported variable in problem: {var:?}"))]
    UnsupportedVariable { var: Variable },
    #[snafu(display("Verification of targeting solution failed: {msg}"))]
    Verification { msg: String },
    #[snafu(display("astro error during targeting: {source}"))]
    Astro { source: AstroError },
    #[snafu(display("targeting aborted, too many iterations"))]
    TooManyIterations,
    #[snafu(display("correction is ineffective at {action}: value at previous iteration {prev_val}, current value: {cur_val}"))]
    CorrectionIneffective {
        prev_val: f64,
        cur_val: f64,
        action: &'static str,
    },
    #[snafu(display("encountered a guidance error: {source}"))]
    GuidanceError { source: GuidanceErrors },
    #[snafu(display("not a finite burn"))]
    NotFinite,
    #[snafu(display("Jacobian is signular"))]
    SingularJacobian,
    #[snafu(display("propagation error during targeting: {source}"))]
    PropError { source: PropagationError },
    #[snafu(display("during an optimization,  encountered {source}"))]
    TargetingTrajError { source: TrajError },
}

impl From<TargetingError> for NyxError {
    fn from(e: TargetingError) -> Self {
        NyxError::Targeter {
            source: Box::new(e),
        }
    }
}
