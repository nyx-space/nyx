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

use crate::dynamics::DynamicsError;
pub use crate::dynamics::{Dynamics, NyxError};
use crate::errors::StateError;
use crate::io::{ConfigError, InputOutputError};
use crate::linalg::OVector;
use crate::md::trajectory::TrajError;
use crate::propagators::PropagationError;
use crate::time::Epoch;
pub use crate::{State, TimeTagged};
use anise::almanac::planetary::PlanetaryDataError;
use anise::errors::AlmanacError;
use hifitime::Duration;
use snafu::prelude::Snafu;
use std::sync::Arc;

pub mod filter;
pub use filter::Filter;

/// Provides a range and range rate measuring models.
mod ground_station;
pub use ground_station::GroundStation;

/// Provides Estimate handling functionalities.
pub mod estimate;

/// Provides noise modeling
pub mod noise;

/// Provides all of the support measurement models
pub mod msr;

/// Provides all of the functionality to simulate measurements from ground stations
pub mod simulator;

/// Provides the interfaces to the orbit determination process
pub mod process;

pub use simulator::TrackingDevice;

/// Provides all state noise compensation functionality
pub mod snc;

/// A helper type for spacecraft orbit determination.
pub type SpacecraftODProcess<'a> = self::process::ODProcess<
    'a,
    crate::md::prelude::SpacecraftDynamics,
    nalgebra::Const<2>,
    nalgebra::Const<3>,
    filter::kalman::KF<crate::Spacecraft, nalgebra::Const<3>, nalgebra::Const<2>>,
    GroundStation,
>;

/// A helper type for spacecraft orbit determination sequentially processing measurements
pub type SpacecraftODProcessSeq<'a> = self::process::ODProcess<
    'a,
    crate::md::prelude::SpacecraftDynamics,
    nalgebra::Const<1>,
    nalgebra::Const<3>,
    filter::kalman::KF<crate::Spacecraft, nalgebra::Const<3>, nalgebra::Const<1>>,
    GroundStation,
>;

#[allow(unused_imports)]
pub mod prelude {
    pub use super::estimate::*;
    pub use super::filter::kalman::*;
    pub use super::ground_station::*;
    pub use super::msr::*;
    pub use super::noise::{GaussMarkov, StochasticNoise, WhiteNoise};
    pub use super::process::*;
    pub use super::simulator::TrackingArcSim;
    pub use super::simulator::*;
    pub use super::snc::*;
    pub use super::*;

    pub use crate::time::{Duration, Epoch, TimeUnits, Unit};
}

#[derive(Debug, PartialEq, Snafu)]
#[snafu(visibility(pub(crate)))]
pub enum ODError {
    #[snafu(display("during an orbit determination, encountered {source}"))]
    ODPropError { source: PropagationError },
    #[snafu(display("during an orbit determination, encountered {source}"))]
    ODDynamicsError { source: DynamicsError },
    #[snafu(display("at least {need} measurements required for {action}"))]
    TooFewMeasurements { need: usize, action: &'static str },
    #[snafu(display("invalid step size: {step}"))]
    StepSizeError { step: Duration },
    #[snafu(display("filter iteration did not converge in {loops} iterations"))]
    Diverged { loops: usize },
    #[snafu(display("STM is singular"))]
    SingularStateTransitionMatrix,
    #[snafu(display("invalid measurement @ {epoch} = {val}"))]
    InvalidMeasurement { epoch: Epoch, val: f64 },
    #[snafu(display("sensitivity matrix must be updated before this call"))]
    SensitivityNotUpdated,
    #[snafu(display("Kalman gain is singular"))]
    SingularKalmanGain,
    #[snafu(display("noise matrix is singular"))]
    SingularNoiseRk,
    #[snafu(display("{kind} noise not configured"))]
    NoiseNotConfigured { kind: String },
    #[snafu(display("measurement sim error: {details}"))]
    MeasurementSimError { details: String },
    #[snafu(display("during an OD encountered {source}: {details}"))]
    ODTrajError { source: TrajError, details: String },
    #[snafu(display("OD failed because {source}"))]
    ODConfigError { source: ConfigError },
    #[snafu(display("OD failed because of an I/O error: {source}"))]
    ODIOError { source: InputOutputError },
    #[snafu(display("OD failed due to Almanac: {action} {source}"))]
    ODAlmanac {
        #[snafu(source(from(AlmanacError, Box::new)))]
        source: Box<AlmanacError>,
        action: &'static str,
    },
    #[snafu(display("OD failed due to planetary data in Almanac: {action} {source}"))]
    ODPlanetaryData {
        #[snafu(source(from(PlanetaryDataError, Box::new)))]
        source: Box<PlanetaryDataError>,
        action: &'static str,
    },
    #[snafu(display("not enough residuals to {action}"))]
    ODNoResiduals { action: &'static str },
    #[snafu(display("Could not {action} OD results: {source}"))]
    ODStateError {
        #[snafu(source(from(StateError, Box::new)))]
        source: Box<StateError>,
        action: &'static str,
    },
}
