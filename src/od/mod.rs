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

extern crate hyperdual;
extern crate rand;
extern crate rand_distr;
extern crate serde;

pub use crate::dynamics::{Dynamics, NyxError};
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector};
use crate::time::Epoch;
pub use crate::{State, TimeTagged};

use crate::io::{CovarFormat, EpochFormat};

/// Provides the Kalman filters. The [examples](https://github.com/ChristopherRabotin/nyx/tree/master/examples) folder may help in the setup.
pub mod kalman;

/// Provides a range and range rate measuring models.
pub mod measurement;

/// Provides Estimate handling functionalities.
pub mod estimate;

/// Provide Residual handling functionalities.
pub mod residual;

/// Provides some helper for filtering.
pub mod ui;

/// Provides all of the measurements functionality, including measurement simulation
// pub mod msr;

/// Provides all state noise compensation functionality
pub mod snc;

/// Defines a Filter trait where S is the size of the estimated state, A the number of acceleration components of the EOMs (used for process noise matrix size), M the size of the measurements.
pub trait Filter<T, A, M>
where
    A: DimName,
    M: DimName,
    T: State,
    DefaultAllocator: Allocator<f64, M>
        + Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::VecLength>
        + Allocator<f64, A>
        + Allocator<f64, M, M>
        + Allocator<f64, M, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<f64, A, A>
        + Allocator<f64, <T as State>::Size, A>
        + Allocator<f64, A, <T as State>::Size>,
{
    type Estimate: estimate::Estimate<T>;

    /// Returns the previous estimate
    fn previous_estimate(&self) -> &Self::Estimate;

    /// Set the previous estimate
    fn set_previous_estimate(&mut self, est: &Self::Estimate);

    /// Update the sensitivity matrix (or "H tilde"). This function **must** be called prior to each
    /// call to `measurement_update`.
    fn update_h_tilde(&mut self, h_tilde: OMatrix<f64, M, <T as State>::Size>);

    /// Computes a time update/prediction at the provided nominal state (i.e. advances the filter estimate with the updated STM).
    ///
    /// Returns an error if the STM was not updated.
    fn time_update(&mut self, nominal_state: T) -> Result<Self::Estimate, NyxError>;

    /// Computes the measurement update with a provided real observation and computed observation.
    ///
    ///Returns an error if the STM or sensitivity matrices were not updated.
    fn measurement_update(
        &mut self,
        nominal_state: T,
        real_obs: &OVector<f64, M>,
        computed_obs: &OVector<f64, M>,
    ) -> Result<(Self::Estimate, residual::Residual<M>), NyxError>;

    /// Returns whether the filter is an extended filter (e.g. EKF)
    fn is_extended(&self) -> bool;

    /// Sets the filter to be extended or not depending on the value of status
    fn set_extended(&mut self, status: bool);

    /// Sets the process noise matrix of the estimated state
    fn set_process_noise(&mut self, snc: snc::SNC<A>);

    /// Returns the measurement noise used at this given epoch
    fn measurement_noise(&self, epoch: Epoch) -> &OMatrix<f64, M, M>;
}

/// A trait defining a measurement of size `MeasurementSize`
pub trait Measurement: TimeTagged
where
    Self: Sized,
    DefaultAllocator: Allocator<f64, Self::MeasurementSize>
        + Allocator<f64, <Self::State as State>::Size>
        + Allocator<f64, <Self::State as State>::Size, <Self::State as State>::Size>
        + Allocator<f64, <Self::State as State>::VecLength>
        + Allocator<f64, Self::MeasurementSize, <Self::State as State>::Size>,
{
    /// Defines the state size of the estimated state
    type State: State;
    /// Defines how much data is measured. For example, if measuring range and range rate, this should be of size 2 (nalgebra::U2).
    type MeasurementSize: DimName;

    /// Returns the measurement/observation as a vector.
    fn observation(&self) -> OVector<f64, Self::MeasurementSize>
    where
        DefaultAllocator: Allocator<f64, Self::MeasurementSize>;

    /// Returns the measurement sensitivity (often referred to as H tilde).
    fn sensitivity(&self) -> OMatrix<f64, Self::MeasurementSize, <Self::State as State>::Size>
    where
        DefaultAllocator: Allocator<f64, <Self::State as State>::Size, Self::MeasurementSize>;

    /// Returns whether the transmitter and receiver where in line of sight.
    fn visible(&self) -> bool;
}

/// A trait to generalize measurement devices such as a ground station
pub trait MeasurementDevice<MsrIn, Msr>
where
    Self: Sized,
    Msr: Measurement,
    DefaultAllocator: Allocator<f64, <Msr::State as State>::Size>
        + Allocator<f64, <Msr::State as State>::Size, <Msr::State as State>::Size>
        + Allocator<f64, <Msr::State as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, <Msr::State as State>::Size>,
{
    /// Returns the measurement if the device and generate one, else returns None
    fn measure(&self, input: &MsrIn) -> Option<Msr>;
}

pub trait EstimateFrom<O: State>
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
}
