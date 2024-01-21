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

use crate::io::ConfigError;
use crate::md::trajectory::TrajError;
use crate::md::StateParameter;
pub use crate::md::TargetingError;
use snafu::prelude::*;
use std::convert::From;

#[derive(Debug, PartialEq, Snafu)]
// #[derive(Error, Debug, PartialEq)]
pub enum NyxError {
    /// Maximum iterations reached
    #[snafu(display("Maximum iterations of {msg} reached"))]
    MaxIterReached { msg: String },
    /// Covariance is not positive semi definite
    #[snafu(display("Covariance is not positive semi definite"))]
    CovarianceMatrixNotPsd,
    /// Targets in Lambert solver too close: Δν ~=0 and A ~=0
    #[snafu(display("Lambert too close: Δν ~=0 and A ~=0"))]
    TargetsTooClose,
    /// No reasonable phi found to connect both radii
    #[snafu(display("No reasonable phi found to connect both radii"))]
    LambertNotReasonablePhi,
    /// Multi revolution Lambert not supported, use the Izzo algorithm for multi-rev transfers
    #[snafu(display("Use the Izzo algorithm for multi-rev transfers"))]
    LambertMultiRevNotSupported,
    /// State parameter cannot be used in this function
    #[snafu(display("Unavailable parameter {param:?}: {msg}"))]
    StateParameterUnavailable { param: StateParameter, msg: String },
    /// Could not load file
    #[snafu(display("Could not load file: {msg}"))]
    LoadingError { msg: String },
    /// Could not read file
    #[snafu(display("Could not read file: {msg}"))]
    FileUnreadable { msg: String },
    /// Celestial object or spacecraft not found
    #[snafu(display("Cosm object not found: `{needle}` (available: {haystack:?})"))]
    ObjectNotFound {
        needle: String,
        haystack: Vec<String>,
    },
    /// No interpolation data
    #[snafu(display("No interpolation data: {msg}"))]
    NoInterpolationData { msg: String },
    /// Invalid interpolation data
    #[snafu(display("Invalid interpolation data: {msg}"))]
    InvalidInterpolationData { msg: String },
    /// No state data
    #[snafu(display("No state data: {msg}"))]
    NoStateData { msg: String },
    /// No thruster attached to spacecraft
    #[snafu(display("No thruster attached to spacecraft"))]
    NoThrusterAvail,
    /// Happens when trying to modify a polynomial's (error)-th error but the polynomial has less orders than that
    #[snafu(display("Happens when trying to modify a polynomial's (error)-th error but the polynomial has less orders than that"))]
    PolynomialOrderError { order: usize },
    /// An objective based analysis or control was attempted, but no objective was defined
    #[snafu(display(
        "An objective based analysis or control was attempted, but no objective was defined"
    ))]
    NoObjectiveDefined,
    /// This computation requires the orbit to be hyperbolic
    #[snafu(display("This computation requires the orbit to be hyperbolic: {msg}"))]
    NotHyperbolic { msg: String },
    /// Monte Carlo error
    #[snafu(display("Monte Carlo error: {msg}"))]
    MonteCarlo { msg: String },
    /// CCSDS error
    #[snafu(display("CCSDS error: {msg}"))]
    CCSDS { msg: String },
    #[snafu(display("Error: {msg}"))]
    CustomError { msg: String },
    /// Trajectory error
    #[snafu(display("Trajectory error: {source}"))]
    Trajectory { source: TrajError },
    /// Math domain
    #[snafu(display("Math domain error: {msg}"))]
    MathDomain { msg: String },
    /// Guidance law config error
    #[snafu(display("Guidance law config error: {msg}"))]
    GuidanceConfigError { msg: String },
    /// Configuration file error
    #[snafu(display("Config error: {source}"))]
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
