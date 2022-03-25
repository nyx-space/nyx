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
    /// STM is singular, propagation or smoothing cannot proceed
    #[error("STM is singular, propagation or smoothing cannot proceed")]
    SingularStateTransitionMatrix,
    /// Fuel exhausted at the provided spacecraft state
    #[error("Fuel exhausted at {0}")]
    FuelExhausted(Box<Spacecraft>),
    /// Propagation event not triggered within the max propagation time
    #[error("Propagation event not triggered within the max propagation time")]
    ConditionNeverTriggered,
    /// Propagation event not hit enough times (requested, found).
    #[error("Propagation event not hit enough times (requested, found).")]
    UnsufficientTriggers(usize, usize),
    /// Maximum iterations reached
    #[error("Maximum iterations of {0} reached")]
    MaxIterReached(String),
    /// Event not found within the provided epochs
    #[error("Event not in braket: {0} <=> {1}")]
    EventNotInEpochBraket(String, String),
    /// The operation was expecting the state to have an STM, but it isn't present.
    #[error("The operation was expecting the state to have an STM, but it isn't present.")]
    StateTransitionMatrixUnset,
    /// The sensitivity matrix must be updated prior to a filter measurement update
    #[error("The sensitivity matrix must be updated prior to a filter measurement update")]
    SensitivityNotUpdated,
    /// Kalman Gain could not be computed because H*P_bar*H + R is singular
    #[error("Gain could not be computed because H*P_bar*H + R is singular")]
    SingularKalmanGain,
    /// Singular Covariance
    #[error("Singular Covariance")]
    SingularCovarianceMatrix,
    /// Covariance is not positive semi definite
    #[error("Covariance is not positive semi definite")]
    CovarianceMatrixNotPsd,
    /// Singular Jacobian
    #[error("Singular Jacobian")]
    SingularJacobian,
    /// Targets in Lambert solver too close: Δν ~=0 and A ~=0
    #[error("Lambert too close: Δν ~=0 and A ~=0")]
    TargetsTooClose,
    /// No reasonable phi found to connect both radii
    #[error("No reasonable phi found to connect both radii")]
    LambertNotReasonablePhi,
    /// Multi revolution Lambert not supported, use the Izzo algorithm for multi-rev transfers
    #[error("Use the Izzo algorithm for multi-rev transfers")]
    LambertMultiRevNotSupported,
    /// Partials for this dynamical model are not defined
    #[error("Partials for this model are not defined")]
    PartialsUndefined,
    /// State parameter cannot be used in this function
    #[error("State parameter cannot be used in this function")]
    StateParameterUnavailable,
    /// Could not load file
    #[error("Could not load file: {0}")]
    LoadingError(String),
    /// Could not read file
    #[error("Could not read file: {0}")]
    FileUnreadable(String),
    /// Celestial object or spacecraft not found
    #[error("Cosm object not found: {0}")]
    ObjectNotFound(String),
    /// No interpolation data
    #[error("No interpolation data: {0}")]
    NoInterpolationData(String),
    /// Invalid interpolation data
    #[error("Invalid interpolation data: {0}")]
    InvalidInterpolationData(String),
    /// No state data
    #[error("No state data: {0}")]
    NoStateData(String),
    /// Cannot convert the state to another frame as the frames are disjoint
    #[error("Cannot convert between disjoint frames: {0} <-> {1}")]
    DisjointFrameOrientations(String, String),
    /// No thruster attached to spacecraft
    #[error("No thruster attached to spacecraft")]
    NoThrusterAvail,
    /// Control vector is not a unit vector
    #[error("Control vector is not a unit vector: {0}")]
    CtrlNotAUnitVector(f64),
    /// Throttle is not between 0.0 and 1.0
    #[error("Throttle is not between 0.0 and 1.0: {0}")]
    CtrlThrottleRangeErr(f64),
    /// Happens when trying to modify a polynomial's (error)-th error but the polynomial has less orders than that
    #[error("Happens when trying to modify a polynomial's (error)-th error but the polynomial has less orders than that")]
    PolynomialOrderError(usize),
    /// An objective based analysis or control was attempted, but no objective was defined
    #[error("An objective based analysis or control was attempted, but no objective was defined")]
    NoObjectiveDefined,
    /// Error when exporting data
    #[error("Error when exporting data: {0}")]
    ExportError(String),
    /// This computation requires the orbit to be hyperbolic
    #[error("This computation requires the orbit to be hyperbolic: {0}")]
    NotHyperbolic(String),
    /// Control variables to not decrease targeting error in differential corrector
    #[error("Control variables to not decrease targeting error in differential corrector: {0}")]
    CorrectionIneffective(String),
    /// Monte Carlo error
    #[error("Monte Carlo error: {0}")]
    MonteCarlo(String),
    /// CCSDS error
    #[error("CCSDS error: {0}")]
    CCSDS(String),
    /// Multiple shooting failed with the provided error at the provided node computation
    #[error("Multiple shooting failed on node {0} with {1}")]
    MultipleShootingTargeter(usize, Box<NyxError>),
    #[error("Custom error: {0}")]
    CustomError(String),
    /// Time related error
    #[error("Time related error: {0}")]
    TimeError(TimeErrors),
    /// Targeting error
    #[error("Targeting error: {0}")]
    Targeter(TargetingError),
    /// Trajectory error
    #[error("Trajectory error: {0}")]
    Trajectory(TrajError),
    /// Math domain
    #[error("Math domain error: {0}")]
    MathDomain(String),
    /// Guidance law config error
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
