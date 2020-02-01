extern crate hyperdual;
extern crate nalgebra as na;
extern crate serde;

use self::hyperdual::{hyperspace_from_vector, Hyperdual, Owned};
use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, DimName, MatrixMN, VectorN};
use crate::hifitime::Epoch;
use celestia::{Frame, State};
use dynamics::Dynamics;
use std::fmt;

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

/// A trait container to specify that given dynamics support linearization, and can be used for state transition matrix computation.
///
/// This trait will likely be made obsolete after the implementation of [#32](https://github.com/ChristopherRabotin/nyx/issues/32).
pub trait Estimable<N>
where
    Self: Dynamics + Sized,
{
    /// Defines the state size of the estimated state
    type LinStateSize: DimName;
    /// Returns the estimated state
    fn estimated_state(&self) -> VectorN<f64, Self::LinStateSize>
    where
        DefaultAllocator: Allocator<f64, Self::LinStateSize>;
    /// Sets the estimated state
    fn set_estimated_state(&mut self, new_state: VectorN<f64, Self::LinStateSize>)
    where
        DefaultAllocator: Allocator<f64, Self::LinStateSize>;
    /// Defines the gradient of the equations of motion for these dynamics.
    fn stm(&self) -> MatrixMN<f64, Self::LinStateSize, Self::LinStateSize>
    where
        DefaultAllocator: Allocator<f64, Self::LinStateSize>
            + Allocator<f64, Self::LinStateSize, Self::LinStateSize>;

    fn to_measurement(&self, prop_state: &Self::StateType) -> (Epoch, N);
}

/// A trait defining a measurement of size `MeasurementSize`
pub trait Measurement
where
    Self: Sized,
    DefaultAllocator: Allocator<f64, Self::MeasurementSize>
        + Allocator<f64, Self::MeasurementSize, Self::StateSize>,
{
    /// Defines the state size of the estimated state
    type StateSize: DimName;
    /// Defines how much data is measured. For example, if measuring range and range rate, this should be of size 2 (nalgebra::U2).
    type MeasurementSize: DimName;

    /// Computes a new measurement from the provided information.
    fn new<F: Frame>(dt: Epoch, tx: State<F>, rx: State<F>, visible: bool) -> Self;

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

    /// Returns the time at which the measurement was performed.
    fn at(&self) -> Epoch;
}

/// A trait to generalize measurement devices such as a ground station
pub trait MeasurementDevice<N>
where
    Self: Sized,
    N: Measurement,
    DefaultAllocator: Allocator<f64, N::StateSize>
        + Allocator<f64, N::StateSize, N::MeasurementSize>
        + Allocator<f64, N::MeasurementSize>
        + Allocator<f64, N::MeasurementSize, N::StateSize>,
{
    type MeasurementInput: Copy;
    fn measure(&self, state: &Self::MeasurementInput) -> N;
}

/// A trait container to specify that given dynamics support linearization, and can be used for state transition matrix computation.
///
/// This trait will likely be made obsolete after the implementation of [#32](https://github.com/ChristopherRabotin/nyx/issues/32).
pub trait AutoDiffDynamics: Dynamics
where
    Self: Sized,
{
    /// Defines the state size of the estimated state
    type HyperStateSize: DimName;
    type STMSize: DimName;

    /// Defines the equations of motion for Dual numbers for these dynamics.
    fn dual_eom(
        &self,
        t: f64,
        state: &VectorN<Hyperdual<f64, Self::HyperStateSize>, Self::STMSize>,
    ) -> (
        VectorN<f64, Self::STMSize>,
        MatrixMN<f64, Self::STMSize, Self::STMSize>,
    )
    where
        DefaultAllocator: Allocator<f64, Self::HyperStateSize>
            + Allocator<f64, Self::STMSize>
            + Allocator<f64, Self::STMSize, Self::STMSize>
            + Allocator<Hyperdual<f64, Self::HyperStateSize>, Self::STMSize>,
        Owned<f64, Self::HyperStateSize>: Copy;

    /// Computes both the state and the gradient of the dynamics. These may be accessed by the
    /// related getters.
    fn compute(
        &self,
        t: f64,
        state: &VectorN<f64, Self::STMSize>,
    ) -> (
        VectorN<f64, Self::STMSize>,
        MatrixMN<f64, Self::STMSize, Self::STMSize>,
    )
    where
        DefaultAllocator: Allocator<f64, Self::HyperStateSize>
            + Allocator<f64, Self::STMSize>
            + Allocator<f64, Self::STMSize, Self::STMSize>
            + Allocator<Hyperdual<f64, Self::HyperStateSize>, Self::STMSize>,
        Owned<f64, Self::HyperStateSize>: Copy,
    {
        let hyperstate = hyperspace_from_vector(&state);

        let (state, grad) = self.dual_eom(t, &hyperstate);

        (state, grad)
    }
}
/*
impl<T: AutoDiffDynamics> Linearization for T
where
    DefaultAllocator:
        Allocator<Dual<f64>, T::StateSize> + Allocator<Dual<f64>, T::StateSize, T::StateSize>,
{
    type StateSize = T::HyperStateSize;

    /// Returns the gradient of the dynamics at the given state.
    ///
    /// **WARNING:** Requires a prior call to self.compute() ! This is where the auto-differentiation happens.
    fn gradient(
        &self,
        _t: f64,
        _state: &VectorN<f64, Self::StateSize>,
    ) -> MatrixMN<f64, Self::StateSize, Self::StateSize>
    where
        DefaultAllocator:
            Allocator<f64, Self::StateSize> + Allocator<f64, Self::StateSize, Self::StateSize>,
    {
        panic!("retrieve the gradient by calling self.compute(...)");
    }
}*/

/// Specifies the format of the Epoch during serialization
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EpochFormat {
    /// Default is MJD TAI, as defined in [hifitime](https://docs.rs/hifitime/).
    MjdTai,
    MjdTt,
    MjdUtc,
    JdeEt,
    JdeTai,
    JdeTt,
    JdeUtc,
    /// Seconds past TAI Epoch
    TaiSecs,
    /// Days past TAI Epoch
    TaiDays,
}

impl fmt::Display for EpochFormat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            EpochFormat::MjdTai => write!(f, "MJD TAI"),
            EpochFormat::MjdTt => write!(f, "MJD TT"),
            EpochFormat::MjdUtc => write!(f, "MJD UTC"),
            EpochFormat::JdeEt => write!(f, "JDE ET"),
            EpochFormat::JdeTai => write!(f, "JDE TAI"),
            EpochFormat::JdeTt => write!(f, "JDE TT"),
            EpochFormat::JdeUtc => write!(f, "JDE UTC"),
            EpochFormat::TaiSecs => write!(f, "TAI+ s"),
            EpochFormat::TaiDays => write!(f, "TAI+ days"),
        }
    }
}

/// Specifies the format of the covariance during serialization
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CovarFormat {
    /// Default: allows plotting the variance of the elements instead of the covariance
    Sqrt,
    /// Keeps the covariance as computed, i.e. one sigma (~68%), causes e.g. positional elements in km^2.
    Sigma1,
    /// Three sigma covers about 99.7% of the distribution
    Sigma3,
    /// Allows specifying a custom multiplication factor of each element of the covariance.
    MulSigma(f64),
}

impl fmt::Display for CovarFormat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            CovarFormat::Sqrt => write!(f, "exptd_val_"),
            CovarFormat::Sigma1 => write!(f, "covar_"),
            CovarFormat::Sigma3 => write!(f, "3sig_covar"),
            CovarFormat::MulSigma(x) => write!(f, "{}sig_covar", x),
        }
    }
}
