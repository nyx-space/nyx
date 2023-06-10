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

use super::thiserror::Error;
use crate::io::ConfigError;
use crate::md::trajectory::TrajError;
use crate::md::StateParameter;
pub use crate::md::TargetingError;
pub use crate::time::Errors as TimeErrors;
use crate::Spacecraft;
use std::convert::From;

/// Represents all possible errors that can occur in the Nyx library.
#[derive(Error, Debug, PartialEq)]
pub enum NyxError {
    /// Occurs when there is a CCSDS error.
    #[error("CCSDS error: {0}")]
    CCSDS(String),

    /// Occurs when a propagation event does not trigger within the maximum propagation time.
    #[error("Propagation event not triggered within the max propagation time")]
    ConditionNeverTriggered,

    /// Occurs when there is a configuration file error.
    #[error("Config error: {0}")]
    ConfigError(ConfigError),

    /// Occurs when control variables do not decrease targeting error in differential corrector.
    #[error("Control variables to not decrease targeting error in differential corrector: {0}")]
    CorrectionIneffective(String),

    /// Occurs when the covariance matrix is not positive semi-definite.
    #[error("Covariance is not positive semi definite")]
    CovarianceMatrixNotPsd,

    /// Occurs when the control vector is not a unit vector.
    #[error("Control vector is not a unit vector: {0}")]
    CtrlNotAUnitVector(f64),

    /// Occurs when the throttle is not between 0.0 and 1.0.
    #[error("Throttle is not between 0.0 and 1.0: {0}")]
    CtrlThrottleRangeErr(f64),

    /// Occurs when a custom error is raised.
    #[error("Custom error: {0}")]
    CustomError(String),

    /// Occurs when it is not possible to convert the state to another frame as the frames are disjoint.
    #[error("Cannot convert between disjoint frames: {0} <-> {1}")]
    DisjointFrameOrientations(String, String),

    /// Occurs when there is an error when exporting data.
    #[error("Error when exporting data: {0}")]
    ExportError(String),

    /// Occurs when a file could not be read.
    #[error("Could not read file: {0}")]
    FileUnreadable(String),

    /// Occurs when the fuel of the spacecraft is exhausted.
    #[error("Fuel exhausted at {0}")]
    FuelExhausted(Box<Spacecraft>),

    /// Occurs when there is a guidance law configuration error.
    #[error("Guidance law config error: {0}")]
    GuidanceConfigError(String),

    /// Occurs when the interpolation data is invalid.
    #[error("Invalid interpolation data: {0}")]
    InvalidInterpolationData(String),

    /// Occurs when a file could not be loaded.
    #[error("Could not load file: {0}")]
    LoadingError(String),

    /// Occurs when multi-revolution Lambert is not supported. The Izzo algorithm should be used for multi-revolution transfers.
    #[error("Use the Izzo algorithm for multi-rev transfers")]
    LambertMultiRevNotSupported,

    /// Occurs when no reasonable phi is found to connect both radii.
    #[error("No reasonable phi found to connect both radii")]
    LambertNotReasonablePhi,

    /// Occurs when there is a math domain error.
    #[error("Math domain error: {0}")]
    MathDomain(String),

    /// Occurs when there is a Monte Carlo error.
    #[error("Monte Carlo error: {0}")]
    MonteCarlo(String),

    /// Occurs when the maximum number of iterations is reached for a particular operation.
    #[error("Maximum iterations of {0} reached")]
    MaxIterReached(String),

    /// Occurs when there is no interpolation data.
    #[error("No interpolation data: {0}")]
    NoInterpolationData(String),

    /// Occurs when an objective based analysis or control was attempted, but no objective was defined.
    #[error("An objective based analysis or control was attempted, but no objective was defined")]
    NoObjectiveDefined,

    /// Occurs when there is no state data.
    #[error("No state data: {0}")]
    NoStateData(String),

    /// Occurs when no thruster is attached to the spacecraft.
    #[error("No thruster attached to spacecraft")]
    NoThrusterAvail,

    /// Occurs when a celestial object or spacecraft is not found.
    #[error("Cosm object not found: `{0}` (available: {1:?})")]
    ObjectNotFound(String, Vec<String>),

    /// Occurs when partials for a particular dynamical model are not defined.
    #[error("Partials for this model are not defined")]
    PartialsUndefined,

    /// Occurs when trying to modify a polynomial's (error)-th error but the polynomial has less orders than that.
    #[error("Happens when trying to modify a polynomial's (error)-th error but the polynomial has less orders than that")]
    PolynomialOrderError(usize),

    /// Occurs when the sensitivity matrix must be updated prior to a filter measurement update, but it is not.
    #[error("The sensitivity matrix must be updated prior to a filter measurement update")]
    SensitivityNotUpdated,

    /// Occurs when the State Transition Matrix (STM) is singular, which means that propagation or smoothing cannot proceed.
    #[error("STM is singular, propagation or smoothing cannot proceed")]
    SingularStateTransitionMatrix,

    /// Occurs when the covariance matrix is singular.
    #[error("Singular Covariance")]
    SingularCovarianceMatrix,

    /// Occurs when the Kalman Gain could not be computed because the matrix H*P_bar*H + R is singular.
    #[error("Gain could not be computed because H*P_bar*H + R is singular")]
    SingularKalmanGain,

    /// Occurs when the Jacobian matrix is singular.
    #[error("Singular Jacobian")]
    SingularJacobian,

    /// Occurs when a state parameter cannot be used in a particular function.
    #[error("Unavailable parameter {0:?}: {1}")]
    StateParameterUnavailable(StateParameter, String),

    /// Occurs when there is a targeting error.
    #[error("Targeting error: {0}")]
    Targeter(Box<TargetingError>),

    /// Occurs when the targets in the Lambert solver are too close, i.e., Δν ~=0 and A ~=0.
    #[error("Lambert too close: Δν ~=0 and A ~=0")]
    TargetsTooClose,

    /// Occurs when there is a time related error.
    #[error("Time related error: {0}")]
    TimeError(TimeErrors),

    /// Occurs when there is a trajectory error.
    #[error("Trajectory error: {0}")]
    Trajectory(TrajError),

    /// Occurs when a propagation event does not hit enough times. The error includes the number of requested and found triggers.
    #[error("Propagation event not hit enough times (requested, found).")]
    UnsufficientTriggers(usize, usize),
}

/// Converts a `TimeErrors` into a `NyxError`.
impl From<TimeErrors> for NyxError {
    fn from(e: TimeErrors) -> Self {
        NyxError::TimeError(e)
    }
}

/// Converts a `TrajError` into a `NyxError`.
impl From<TrajError> for NyxError {
    fn from(e: TrajError) -> Self {
        NyxError::Trajectory(e)
    }
}

/// Converts a `ConfigError` into a `NyxError`.
impl From<ConfigError> for NyxError {
    fn from(e: ConfigError) -> Self {
        NyxError::ConfigError(e)
    }
}
