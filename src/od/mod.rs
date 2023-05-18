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

pub use crate::dynamics::{Dynamics, NyxError};
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector};
use crate::time::Epoch;
use crate::Orbit;
pub use crate::{cosmic::Cosm, State, TimeTagged};
use std::sync::Arc;

use crate::io::{CovarFormat, EpochFormat};

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
pub use simulator::trackdata::TrackingDeviceSim;

/// Provides all state noise compensation functionality
pub mod snc;

pub mod prelude {
    pub use super::estimate::*;
    pub use super::filter::kalman::*;
    pub use super::ground_station::*;
    pub use super::msr::*;
    pub use super::process::*;
    pub use super::simulator::arc::TrackingArcSim;
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
        DefaultAllocator: Allocator<f64, Self::MeasurementSize>;

    /// Returns the measurement/observation as a vector.
    fn observation(&self) -> OVector<f64, Self::MeasurementSize>
    where
        DefaultAllocator: Allocator<f64, Self::MeasurementSize>;
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
    DefaultAllocator: Allocator<f64, <O as State>::Size>
        + Allocator<f64, <O as State>::VecLength>
        + Allocator<f64, <O as State>::Size, <O as State>::Size>
        + Allocator<f64, Self::Size>
        + Allocator<f64, Self::VecLength>
        + Allocator<f64, Self::Size, Self::Size>,
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
        DefaultAllocator: Allocator<f64, M::MeasurementSize, Self::Size>;
}

/// A generic implementation of EstimateFrom for any State that is also a Measurement, e.g. if there is a direct observation of the full state.
/// WARNING: The frame of the full state measurement is _not_ checked to match that of `Self` or of the filtering frame.
impl<O> EstimateFrom<O, O> for O
where
    O: State + Measurement,
    Self: State,
    DefaultAllocator: Allocator<f64, <O as State>::Size>
        + Allocator<f64, <O as State>::VecLength>
        + Allocator<f64, <O as State>::Size, <O as State>::Size>
        + Allocator<f64, Self::Size>
        + Allocator<f64, Self::VecLength>
        + Allocator<f64, Self::Size, Self::Size>,
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
        DefaultAllocator: Allocator<f64, <O as Measurement>::MeasurementSize, Self::Size>,
    {
        OMatrix::<f64, O::MeasurementSize, Self::Size>::identity()
    }
}
