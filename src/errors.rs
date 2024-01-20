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
use snafu::prelude::*;
use std::convert::From;

// #[derive(Error, Debug, PartialEq, Snafu)]
#[derive(Error, Debug, PartialEq)]
pub enum NyxError {
    /// Maximum iterations reached
    #[error("Maximum iterations of {msg} reached")]
    MaxIterReached { msg: String },
    /// Covariance is not positive semi definite
    #[error("Covariance is not positive semi definite")]
    CovarianceMatrixNotPsd,
    /// Targets in Lambert solver too close: Δν ~=0 and A ~=0
    #[error("Lambert too close: Δν ~=0 and A ~=0")]
    TargetsTooClose,
    /// No reasonable phi found to connect both radii
    #[error("No reasonable phi found to connect both radii")]
    LambertNotReasonablePhi,
    /// Multi revolution Lambert not supported, use the Izzo algorithm for multi-rev transfers
    #[error("Use the Izzo algorithm for multi-rev transfers")]
    LambertMultiRevNotSupported,
    /// State parameter cannot be used in this function
    #[error("Unavailable parameter {param:?}: {msg}")]
    StateParameterUnavailable { param: StateParameter, msg: String },
    /// Could not load file
    #[error("Could not load file: {msg}")]
    LoadingError { msg: String },
    /// Could not read file
    #[error("Could not read file: {msg}")]
    FileUnreadable { msg: String },
    /// Celestial object or spacecraft not found
    #[error("Cosm object not found: `{needle}` (available: {haystack:?})")]
    ObjectNotFound {
        needle: String,
        haystack: Vec<String>,
    },
    /// No interpolation data
    #[error("No interpolation data: {msg}")]
    NoInterpolationData { msg: String },
    /// Invalid interpolation data
    #[error("Invalid interpolation data: {msg}")]
    InvalidInterpolationData { msg: String },
    /// No state data
    #[error("No state data: {msg}")]
    NoStateData { msg: String },
    /// No thruster attached to spacecraft
    #[error("No thruster attached to spacecraft")]
    NoThrusterAvail,
    /// Happens when trying to modify a polynomial's (error)-th error but the polynomial has less orders than that
    #[error("Happens when trying to modify a polynomial's (error)-th error but the polynomial has less orders than that")]
    PolynomialOrderError { order: usize },
    /// An objective based analysis or control was attempted, but no objective was defined
    #[error("An objective based analysis or control was attempted, but no objective was defined")]
    NoObjectiveDefined,
    /// This computation requires the orbit to be hyperbolic
    #[error("This computation requires the orbit to be hyperbolic: {msg}")]
    NotHyperbolic { msg: String },
    /// Monte Carlo error
    #[error("Monte Carlo error: {msg}")]
    MonteCarlo { msg: String },
    /// CCSDS error
    #[error("CCSDS error: {msg}")]
    CCSDS { msg: String },
    #[error("Custom error: {0}")]
    CustomError(String),
    /// Trajectory error
    #[error("Trajectory error: {source}")]
    Trajectory { source: TrajError },
    /// Math domain
    #[error("Math domain error: {msg}")]
    MathDomain { msg: String },
    /// Guidance law config error
    #[error("Guidance law config error: {msg}")]
    GuidanceConfigError { msg: String },
    /// Configuration file error
    #[error("Config error: {source}")]
    ConfigError { source: ConfigError },
}

impl From<TrajError> for NyxError {
    fn from(source: TrajError) -> Self {
        NyxError::Trajectory { source }
    }
}

impl From<ConfigError> for NyxError {
    fn from(source: ConfigError) -> Self {
        NyxError::ConfigError { source }
    }
}
