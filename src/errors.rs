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

pub use crate::time::Errors as TimeErrors;
use crate::Spacecraft;
use std::convert::From;
use std::error::Error;
use std::fmt;

#[derive(Clone, PartialEq, Debug)]
pub enum NyxError {
    /// STM is singular, check the automatic differentiation function and file a bug
    SingularStateTransitionMatrix,
    /// Fuel exhausted error, try running without fuel depletion, and then adding it.
    /// Parameter is the state at which the fuel has exhaused
    FuelExhausted(Box<Spacecraft>),
    /// Propagation event not triggered withinin the max propagation time
    ConditionNeverTriggered,
    /// Propagation event not hit enough times (requested, found).
    UnsufficientTriggers(usize, usize),
    /// Maximum iterations reached, value corresponds to the number of iterations used
    MaxIterReached(String),
    /// Event not in braket
    EventNotInEpochBraket(String, String),
    /// The operation was expecting the state to have an STM, but it isn't present.
    StateTransitionMatrixUnset,
    /// The sensitivity matrix was not updated prior to requesting a filter measurement update
    SensitivityNotUpdated,
    /// Kalman gain is singular, file a bug if this is encountered
    SingularKalmanGain,
    /// Covariance is singular
    SingularCovarianceMatrix,
    /// Jacobian of some optimization problem is singular
    SingularJacobian,
    TargetsTooClose,
    LambertNotReasonablePhi,
    LambertMultiRevNotSupported,
    /// Returns this error if the partials for this model are not defined, thereby preventing the computation of the STM
    PartialsUndefined,
    /// Returned if trying to set a parameter for something which does not have that parameter.
    ParameterUnavailableForType,
    LoadingError(String),
    FileUnreadable(String),
    ObjectNotFound(String),
    NoInterpolationData(String),
    InvalidInterpolationData(String),
    OutOfInterpolationWindow(String),
    NoStateData(String),
    DisjointFrameOrientations(String, String),
    /// When there is a controller but there isn't any thruster available
    CtrlExistsButNoThrusterAvail,
    /// The control vector returned by a controller must be a unit vector. Use the throttle() function to specify the amount.
    CtrlNotAUnitVector(f64),
    /// The control throttle range must be between 0.0 and 1.0 (both included) as it represents a percentage.
    CtrlThrottleRangeErr(f64),
    /// An objective based analysis or control was attempted, but no objective was defined.
    NoObjectiveDefined,
    /// Error when exporting data
    ExportError(String),
    /// This computation requires the orbit to be hyperbolic
    NotHyperbolic(String),
    /// Raised if a differential corrector is not decreasing the error
    CorrectionIneffective(String),
    /// When there is an error during a Monte Carlo or in the conditions starting a Monte Carlo run
    MonteCarlo(String),
    /// Raised if the variables to be adjusted lead to an over-determined of the problem for the targeter
    TargetError(String),
    /// Raised if the variables to be adjusted lead to an under-determined of the problem for the targeter
    UnderdeterminedProblem,
    /// Returned if CCSDS encountered an error
    CCSDS(String),
    /// Returned if the targeter for `node_no` has failed
    MultipleShootingTargeter(usize, Box<NyxError>),
    /// Returned when the trajectory could not be created
    TrajectoryCreationError,
    /// Some custom error for new dynamics
    CustomError(String),
    /// Hifitime errors that rose upward
    TimeError(TimeErrors),
}

impl fmt::Display for NyxError {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::SensitivityNotUpdated => write!(
                f,
                "The measurement matrix H_tilde was not updated prior to measurement update"
            ),
            Self::SingularKalmanGain => write!(
                f,
                "Gain could not be computed because H*P_bar*H + R is singular"
            ),
            Self::SingularStateTransitionMatrix => write!(
                f,
                "STM is singular, propagation or smoothing cannot proceed"
            ),
            Self::SingularCovarianceMatrix => {
                write!(f, "Covariance is singular, smoothing cannot proceed")
            }
            Self::FuelExhausted(sc) => write!(
                f,
                "Spacecraft fuel exhausted, disable fuel depletion and place maneuvers\n{}",
                sc
            ),
            Self::ConditionNeverTriggered => write!(
                f,
                "Try increasing the search space, i.e. increase the maximum propagation time"
            ),
            Self::TargetsTooClose => write!(f, "Lambert too close: Δν ~=0 and A ~=0"),
            Self::LambertMultiRevNotSupported => {
                write!(f, "Use the Izzo algorithm for multi-rev transfers")
            }
            Self::LambertNotReasonablePhi => {
                write!(f, "No reasonable phi found to connect both radii")
            }
            Self::MultipleShootingTargeter(n, e) => {
                write!(f, "Multiple shooting failed on node {} with {}", n, e)
            }
            _ => write!(f, "{:?}", self),
        }
    }
}

impl Error for NyxError {}

impl From<TimeErrors> for NyxError {
    fn from(e: TimeErrors) -> Self {
        NyxError::TimeError(e)
    }
}
