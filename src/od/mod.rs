extern crate serde;

use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, DimName, MatrixMN, VectorN};
use crate::time::Epoch;
use crate::{State, TimeTagged};
pub use dynamics::{Dynamics, NyxError};

use crate::io::{CovarFormat, EpochFormat};

/// Provides the Kalman filters. The [examples](https://github.com/ChristopherRabotin/nyx/tree/master/examples) folder may help in the setup.
pub mod kalman;

/// Provides a range and range rate measuring models.
pub mod ranging;

/// Provides Estimate handling functionalities.
pub mod estimate;

/// Provide Residual handling functionalities.
pub mod residual;

/// Provides some helper for filtering.
pub mod ui;

/// Provides the Square Root Information Filter
pub mod srif;

/// Provides all state noise compensation functionality
pub mod snc;

// / A trait container to specify that given dynamics support linearization, and can be used for state transition matrix computation.
// /
// / This trait will likely be made obsolete after the implementation of [#32](https://github.com/ChristopherRabotin/nyx/issues/32).
// pub trait Estimable<N>
// where
//     Self: Dynamics + Sized,
// {
//     /// Defines the state size of the estimated state
//     type LinStateSize: DimName;
//     /// Returns the estimated state
//     fn extract_estimated_state(
//         &self,
//         prop_state: &Self::StateType,
//     ) -> VectorN<f64, Self::LinStateSize>
//     where
//         DefaultAllocator: Allocator<f64, Self::LinStateSize>;

//     /// Returns the estimated state
//     fn estimated_state(&self) -> VectorN<f64, Self::LinStateSize>
//     where
//         DefaultAllocator: Allocator<f64, Self::LinStateSize>,
//     {
//         self.extract_estimated_state(&self.state())
//     }

//     /// Sets the estimated state
//     fn set_estimated_state(&mut self, new_state: VectorN<f64, Self::LinStateSize>)
//     where
//         DefaultAllocator: Allocator<f64, Self::LinStateSize>;

//     /// Defines the gradient of the equations of motion for these dynamics.
//     fn stm(&self) -> MatrixMN<f64, Self::LinStateSize, Self::LinStateSize>
//     where
//         DefaultAllocator: Allocator<f64, Self::LinStateSize>
//             + Allocator<f64, Self::LinStateSize, Self::LinStateSize>,
//     {
//         self.extract_stm(&self.state())
//     }

//     /// Converts the Dynamics' state type to a measurement to be ingested in a filter
//     fn to_measurement(&self, prop_state: &Self::StateType) -> N;

//     /// Extracts the STM from the dynamics state
//     fn extract_stm(
//         &self,
//         prop_state: &Self::StateType,
//     ) -> MatrixMN<f64, Self::LinStateSize, Self::LinStateSize>
//     where
//         DefaultAllocator: Allocator<f64, Self::LinStateSize>
//             + Allocator<f64, Self::LinStateSize, Self::LinStateSize>;
// }

/// Defines a Filter trait where S is the size of the estimated state, A the number of acceleration components of the EOMs (used for process noise matrix size), M the size of the measurements.
pub trait Filter<T, A, M>
where
    A: DimName,
    M: DimName,
    T: State,
    DefaultAllocator: Allocator<f64, M>
        + Allocator<f64, <T as State>::Size>
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

    /// Update the State Transition Matrix (STM). This function **must** be called in between each
    /// call to `time_update` or `measurement_update`.
    fn update_stm(&mut self, new_stm: MatrixMN<f64, <T as State>::Size, <T as State>::Size>);

    /// Update the sensitivity matrix (or "H tilde"). This function **must** be called prior to each
    /// call to `measurement_update`.
    fn update_h_tilde(&mut self, h_tilde: MatrixMN<f64, M, <T as State>::Size>);

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
        real_obs: VectorN<f64, M>,
        computed_obs: VectorN<f64, M>,
    ) -> Result<(Self::Estimate, residual::Residual<M>), NyxError>;

    /// Returns whether the filter is an extended filter (e.g. EKF)
    fn is_extended(&self) -> bool;

    /// Sets the filter to be extended or not depending on the value of status
    fn set_extended(&mut self, status: bool);

    /// Sets the process noise matrix of the estimated state
    fn set_process_noise(&mut self, snc: snc::SNC<A>);
}

/// A trait defining a measurement of size `MeasurementSize`
pub trait Measurement: TimeTagged
where
    Self: Sized,
    DefaultAllocator: Allocator<f64, Self::MeasurementSize>
        + Allocator<f64, Self::MeasurementSize, Self::StateSize>,
{
    /// Defines the state size of the estimated state
    type StateSize: DimName;
    /// Defines how much data is measured. For example, if measuring range and range rate, this should be of size 2 (nalgebra::U2).
    type MeasurementSize: DimName;

    /// Returns the measurement/observation as a vector.
    fn observation(&self) -> VectorN<f64, Self::MeasurementSize>
    where
        DefaultAllocator: Allocator<f64, Self::MeasurementSize>;

    /// Returns the measurement sensitivity (often referred to as H tilde).
    fn sensitivity(&self) -> MatrixMN<f64, Self::MeasurementSize, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize, Self::MeasurementSize>;

    /// Returns whether the transmitter and receiver where in line of sight.
    fn visible(&self) -> bool;
}

/// A trait to generalize measurement devices such as a ground station
pub trait MeasurementDevice<MsrIn, Msr>
where
    Self: Sized,
    Msr: Measurement,
    DefaultAllocator: Allocator<f64, Msr::StateSize>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::StateSize>,
{
    /// Returns the measurement if the device and generate one, else returns None
    fn measure(&self, input: &MsrIn) -> Option<Msr>;
}

pub trait EstimateFrom<O: State>
where
    Self: State,
{
    // TODO: I need a `from` as well because somehow I need to update the state of the propagator... ugh
    fn extract(from: &O) -> Self;
}

// impl EstimateAs<Orbit> for SpacecraftState;
// impl EstimateFrom<SpacecraftState> for Orbit;
