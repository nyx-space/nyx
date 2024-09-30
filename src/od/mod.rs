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
use crate::io::{ConfigError, InputOutputError};
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector};
use crate::md::trajectory::TrajError;
use crate::propagators::PropagationError;
use crate::time::Epoch;
use crate::Orbit;
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

use arrow::datatypes::Field;
pub use simulator::TrackingDeviceSim;

/// Provides all state noise compensation functionality
pub mod snc;

/// A helper type for spacecraft orbit determination.
pub type SpacecraftODProcess<'a> = self::process::ODProcess<
    'a,
    crate::md::prelude::SpacecraftDynamics,
    msr::RangeDoppler,
    nalgebra::Const<3>,
    crate::Spacecraft,
    filter::kalman::KF<crate::Spacecraft, nalgebra::Const<3>, nalgebra::Const<2>>,
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

/// A trait defining a measurement that can be used in the orbit determination process.
pub trait Measurement: Copy + TimeTagged {
    /// Defines how much data is measured. For example, if measuring range and range rate, this should be of size 2 (nalgebra::U2).
    type MeasurementSize: DimName;

    /// Returns the fields for this kind of measurement.
    /// The metadata must include a `unit` field with the unit.
    fn fields() -> Vec<Field>;

    /// Initializes a new measurement from the provided data.
    fn from_observation(epoch: Epoch, obs: OVector<f64, Self::MeasurementSize>) -> Self
    where
        DefaultAllocator: Allocator<Self::MeasurementSize>;

    /// Returns the measurement/observation as a vector.
    fn observation(&self) -> OVector<f64, Self::MeasurementSize>
    where
        DefaultAllocator: Allocator<Self::MeasurementSize>;
}

/// The Estimate trait defines the interface that is the opposite of a `SolveFor`.
/// For example, `impl EstimateFrom<Spacecraft> for Orbit` means that the `Orbit` can be estimated (i.e. "solved for") from a `Spacecraft`.
///
/// In the future, there will be a way to estimate ground station biases, for example. This will need a new State that includes both the Spacecraft and
/// the ground station bias information. Then, the `impl EstimateFrom<SpacecraftAndBias> for OrbitAndBias` will be added, where `OrbitAndBias` is the
/// new State that includes the orbit and the bias of one ground station.
pub trait EstimateFrom<O: State, M: Measurement>
where
    Self: State,
    DefaultAllocator: Allocator<<O as State>::Size>
        + Allocator<<O as State>::VecLength>
        + Allocator<<O as State>::Size, <O as State>::Size>
        + Allocator<Self::Size>
        + Allocator<Self::VecLength>
        + Allocator<Self::Size, Self::Size>,
{
    /// From the state extract the state to be estimated
    fn extract(from: O) -> Self;

    /// Returns the measurement sensitivity (often referred to as H tilde).
    ///
    /// # Limitations
    /// The transmitter is necessarily an Orbit. This implies that any non-orbit parameter in the estimation vector must be a zero-bias estimator, i.e. it must be assumed that the parameter should be zero.
    /// This is a limitation of the current implementation. It could be fixed by specifying another State like trait in the EstimateFrom trait, significantly adding complexity with little practical use.
    /// To solve for non zero bias parameters, one ought to be able to estimate the _delta_ of that parameter and want that delta to return to zero, thereby becoming a zero-bias estimator.
    fn sensitivity(
        msr: &M,
        receiver: Self,
        transmitter: Orbit,
    ) -> OMatrix<f64, M::MeasurementSize, Self::Size>
    where
        DefaultAllocator: Allocator<M::MeasurementSize, Self::Size>;
}

/// A generic implementation of EstimateFrom for any State that is also a Measurement, e.g. if there is a direct observation of the full state.
/// WARNING: The frame of the full state measurement is _not_ checked to match that of `Self` or of the filtering frame.
impl<O> EstimateFrom<O, O> for O
where
    O: State + Measurement,
    Self: State,
    DefaultAllocator: Allocator<<O as State>::Size>
        + Allocator<<O as State>::VecLength>
        + Allocator<<O as State>::Size, <O as State>::Size>
        + Allocator<Self::Size>
        + Allocator<Self::VecLength>
        + Allocator<Self::Size, Self::Size>,
{
    fn extract(from: O) -> Self {
        from
    }

    fn sensitivity(
        _full_state_msr: &O,
        _receiver: Self,
        _transmitter: Orbit,
    ) -> OMatrix<f64, <O as Measurement>::MeasurementSize, Self::Size>
    where
        DefaultAllocator: Allocator<<O as Measurement>::MeasurementSize, Self::Size>,
    {
        OMatrix::<f64, O::MeasurementSize, Self::Size>::identity()
    }
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
    #[snafu(display("Noise matrix is singular"))]
    SingularNoiseRk,
    #[snafu(display("{kind} noise not configured"))]
    NoiseNotConfigured { kind: &'static str },
    #[snafu(display("during an OD encountered {source}"))]
    ODTrajError { source: TrajError },
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
}
