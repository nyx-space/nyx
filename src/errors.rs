/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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

use crate::md::trajectory::TrajError;
use crate::md::StateParameter;
pub use crate::md::TargetingError;
use crate::{cosmic::AstroError, io::ConfigError};
use anise::errors::{AlmanacError, PhysicsError};
use hifitime::Epoch;
use snafu::prelude::*;
use std::convert::From;

#[derive(Debug, Snafu)]
#[snafu(visibility(pub(crate)))]
pub enum NyxError {
    #[snafu(display("Maximum iterations of {msg} reached"))]
    MaxIterReached { msg: String },
    #[snafu(display("Covariance is not positive semi definite"))]
    CovarianceMatrixNotPsd,
    #[snafu(display("Lambert too close: Δν ~=0 and A ~=0"))]
    TargetsTooClose,
    #[snafu(display("No reasonable phi found to connect both radii"))]
    LambertNotReasonablePhi,
    #[snafu(display("Use the Izzo algorithm for multi-rev transfers"))]
    LambertMultiRevNotSupported,
    #[snafu(display("Unavailable parameter {param:?}: {msg}"))]
    StateParameterUnavailable { param: StateParameter, msg: String },
    #[snafu(display("Could not load file: {msg}"))]
    LoadingError { msg: String },
    #[snafu(display("Could not read file: {msg}"))]
    FileUnreadable { msg: String },
    #[snafu(display("Cosm object not found: `{needle}` (available: {haystack:?})"))]
    ObjectNotFound {
        needle: String,
        haystack: Vec<String>,
    },
    #[snafu(display("No interpolation data: {msg}"))]
    NoInterpolationData { msg: String },
    #[snafu(display("Invalid interpolation data: {msg}"))]
    InvalidInterpolationData { msg: String },
    #[snafu(display("No state data: {msg}"))]
    NoStateData { msg: String },
    #[snafu(display("Happens when trying to modify a polynomial's (error)-th error but the polynomial has less orders than that"))]
    PolynomialOrderError { order: usize },
    #[snafu(display(
        "An objective based analysis or control was attempted, but no objective was defined"
    ))]
    NoObjectiveDefined,
    #[snafu(display("This computation requires the orbit to be hyperbolic: {msg}"))]
    NotHyperbolic { msg: String },
    #[snafu(display("Monte Carlo error: {msg}"))]
    MonteCarlo { msg: String },
    #[snafu(display("CCSDS error: {msg}"))]
    CCSDS { msg: String },
    #[snafu(display("Error: {msg}"))]
    CustomError { msg: String },
    #[snafu(display("Trajectory error: {source}"))]
    Trajectory { source: TrajError },
    #[snafu(display("Math domain error: {msg}"))]
    MathDomain { msg: String },
    #[snafu(display("Guidance law config error: {msg}"))]
    GuidanceConfigError { msg: String },
    #[snafu(display("Config error: {source}"))]
    ConfigError { source: ConfigError },
    #[snafu(display("issue due to Almanac: {action} {source}"))]
    FromAlmanacError {
        #[snafu(source(from(AlmanacError, Box::new)))]
        source: Box<AlmanacError>,
        action: &'static str,
    },
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

#[derive(Debug, PartialEq, Snafu)]
#[snafu(visibility(pub(crate)))]
pub enum StateError {
    #[snafu(display("{param} is unavailable for this kind of state"))]
    Unavailable { param: StateParameter },
    #[snafu(display("{param} is read only for this kind of state"))]
    ReadOnly { param: StateParameter },
    #[snafu(display("{param} computation caused {source}"))]
    StateAstroError {
        param: StateParameter,
        source: AstroError,
    },
    #[snafu(display("No thruster attached to spacecraft"))]
    NoThrusterAvail,
}

#[derive(Debug, PartialEq, Snafu)]
#[snafu(visibility(pub(crate)))]
pub enum EventError {
    #[snafu(display("during event computation: {source}"))]
    EventAlmanacError {
        #[snafu(source(from(AlmanacError, Box::new)))]
        source: Box<AlmanacError>,
    },
    #[snafu(display("during event computation: {source}"))]
    EventStateError {
        param: StateParameter,
        source: StateError,
    },
    #[snafu(display("during event computation: {source}"))]
    EventPhysicsError { source: PhysicsError },
    #[snafu(display("when computing an event in a trajectory {source}"))]
    EventTrajError { source: TrajError },
    #[snafu(display("Event {event} not found between {start} and {end}"))]
    NotFound {
        start: Epoch,
        end: Epoch,
        event: String,
    },
}

#[derive(Debug, Snafu)]
#[snafu(visibility(pub(crate)))]
pub enum MonteCarloError {
    #[snafu(display("Monte Carlo caused {source}"))]
    StateError { source: StateError },
    #[snafu(display("for {param}, expected percentage between 0.0 and 1.0 but got {prct}"))]
    ParamPercentage { param: StateParameter, prct: f64 },
}
