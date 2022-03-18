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

use super::thiserror::Error;
use crate::md::trajectory::TrajError;
pub use crate::md::TargetingError;
pub use crate::time::Errors as TimeErrors;
use crate::Spacecraft;
use std::convert::From;

#[derive(Clone, PartialEq, Error, Debug)]
pub enum NyxError {
    #[error("STM is singular, propagation or smoothing cannot proceed")]
    SingularStateTransitionMatrix,
    #[error("Fuel exhausted at {0}")]
    FuelExhausted(Box<Spacecraft>),
    #[error("Propagation event not triggered within the max propagation time")]
    ConditionNeverTriggered,
    #[error("Propagation event not hit enough times (requested, found).")]
    UnsufficientTriggers(usize, usize),
    #[error("Maximum iterations of {0} reached")]
    MaxIterReached(String),
    #[error("Event not in braket: {0} <=> {1}")]
    EventNotInEpochBraket(String, String),
    #[error("The operation was expecting the state to have an STM, but it isn't present.")]
    StateTransitionMatrixUnset,
    #[error("The sensitivity matrix must be updated prior to a filter measurement update")]
    SensitivityNotUpdated,
    #[error("Gain could not be computed because H*P_bar*H + R is singular")]
    SingularKalmanGain,
    #[error("Singular Covariance")]
    SingularCovarianceMatrix,
    #[error("Singular Jacobian")]
    SingularJacobian,
    #[error("Lambert too close: Δν ~=0 and A ~=0")]
    TargetsTooClose,
    #[error("No reasonable phi found to connect both radii")]
    LambertNotReasonablePhi,
    #[error("Use the Izzo algorithm for multi-rev transfers")]
    LambertMultiRevNotSupported,
    #[error("Partials for this model are not defined")]
    PartialsUndefined,
    #[error("State parameter cannot be used in this function")]
    StateParameterUnavailable,
    #[error("Could not load file: {0}")]
    LoadingError(String),
    #[error("Could not read file: {0}")]
    FileUnreadable(String),
    #[error("Cosm object not found")]
    ObjectNotFound(String),
    #[error("No interpolation data: {0}")]
    NoInterpolationData(String),
    #[error("Invalid interpolation data: {0}")]
    InvalidInterpolationData(String),
    #[error("No state data: {0}")]
    NoStateData(String),
    #[error("Cannot convert between disjoint frames: {0} <-> {1}")]
    DisjointFrameOrientations(String, String),
    #[error("No thruster attached to spacecraft")]
    NoThrusterAvail,
    #[error("Control vector is not a unit vector: {0}")]
    CtrlNotAUnitVector(f64),
    #[error("Throttle is not between 0.0 and 1.0: {0}")]
    CtrlThrottleRangeErr(f64),
    #[error("Happens when trying to modify a polynomial's (error)-th error but the polynomial has less orders than that")]
    PolynomialOrderError(usize),
    #[error("An objective based analysis or control was attempted, but no objective was defined")]
    NoObjectiveDefined,
    #[error("Error when exporting data: {0}")]
    ExportError(String),
    #[error("This computation requires the orbit to be hyperbolic: {0}")]
    NotHyperbolic(String),
    #[error("Control variables to not decrease targeting error in differential corrector: {0}")]
    CorrectionIneffective(String),
    #[error("Monte Carlo error: {0}")]
    MonteCarlo(String),
    #[error("CCSDS error: {0}")]
    CCSDS(String),
    #[error("Multiple shooting failed on node {0} with {1}")]
    MultipleShootingTargeter(usize, Box<NyxError>),
    #[error("Custom error: {0}")]
    CustomError(String),
    #[error("Time related error: {0}")]
    TimeError(TimeErrors),
    #[error("Targeting error: {0}")]
    Targeter(TargetingError),
    #[error("Trajectory error: {0}")]
    Trajectory(TrajError),
    #[error("Math domain error: {0}")]
    MathDomain(String),
    #[error("Guidance law config error: {0}")]
    GuidanceConfigError(String),
}

impl From<TimeErrors> for NyxError {
    fn from(e: TimeErrors) -> Self {
        NyxError::TimeError(e)
    }
}

impl From<TrajError> for NyxError {
    fn from(e: TrajError) -> Self {
        NyxError::Trajectory(e)
    }
}
