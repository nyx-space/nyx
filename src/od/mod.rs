extern crate dual_num;
extern crate hifitime;
extern crate nalgebra as na;
extern crate serde;

use self::dual_num::{hyperspace_from_vector, Dual, DualN, Float, FloatConst, Hyperdual, Owned};
use self::hifitime::instant::Instant;
use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, DimName, MatrixMN, VectorN};
use celestia::{CoordinateFrame, State};
use dynamics::Dynamics;

/// Provides the Kalman filters. The [examples](https://github.com/ChristopherRabotin/nyx/tree/master/examples) folder may help in the setup.
pub mod kalman;

/// Provides a range and range rate measuring models.
pub mod ranging;

/// A trait container to specify that given dynamics support linearization, and can be used for state transition matrix computation.
///
/// This trait will likely be made obsolete after the implementation of [#32](https://github.com/ChristopherRabotin/nyx/issues/32).
pub trait Linearization
where
    Self: Sized,
{
    /// Defines the state size of the estimated state
    type StateSize: DimName;
    /// Defines the gradient of the equations of motion for these dynamics.
    fn gradient(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> MatrixMN<f64, Self::StateSize, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize> + Allocator<f64, Self::StateSize, Self::StateSize>;
}

/// A trait defining a measurement of size `MeasurementSize`
pub trait Measurement
where
    Self: Sized,
    DefaultAllocator: Allocator<f64, Self::MeasurementSize> + Allocator<f64, Self::MeasurementSize, Self::StateSize>,
{
    /// Defines the state size of the estimated state
    type StateSize: DimName;
    /// Defines how much data is measured. For example, if measuring range and range rate, this should be of size 2 (nalgebra::U2).
    type MeasurementSize: DimName;

    /// Computes a new measurement from the provided information.
    fn new<F: CoordinateFrame>(dt: Instant, tx: State<F>, rx: State<F>, visible: bool) -> Self;

    /// Returns the measurement/observation as a vector.
    fn observation(&self) -> &VectorN<f64, Self::MeasurementSize>
    where
        DefaultAllocator: Allocator<f64, Self::MeasurementSize>;

    /// Returns the measurement sensitivity (often referred to as H tilde).
    fn sensitivity(&self) -> &MatrixMN<f64, Self::MeasurementSize, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize, Self::MeasurementSize>;

    /// Returns whether the transmitter and receiver where in line of sight.
    fn visible(&self) -> bool;

    /// Returns the time at which the measurement was performed.
    fn at(&self) -> Instant;
}

/// A trait container to specify that given dynamics support linearization, and can be used for state transition matrix computation.
///
/// This trait will likely be made obsolete after the implementation of [#32](https://github.com/ChristopherRabotin/nyx/issues/32).
pub trait AutoDiffDynamics: Dynamics
where
    Self: Copy + Sized,
{
    /// Defines the state size of the estimated state
    type HyperStateSize: DimName;
    type STMSize: DimName;

    /// Defines the equations of motion for Dual numbers for these dynamics.
    fn dual_eom(
        &self,
        t: f64,
        state: &VectorN<Hyperdual<f64, Self::HyperStateSize>, Self::STMSize>,
    ) -> (VectorN<f64, Self::STMSize>, MatrixMN<f64, Self::STMSize, Self::STMSize>)
    where
        DefaultAllocator: Allocator<f64, Self::HyperStateSize>
            + Allocator<f64, Self::STMSize>
            + Allocator<f64, Self::STMSize, Self::STMSize>
            + Allocator<Hyperdual<f64, Self::HyperStateSize>, Self::STMSize>,
        Owned<f64, Self::HyperStateSize>: Copy;

    /// Computes both the state and the gradient of the dynamics. These may be accessed by the related
    /// getters.
    fn compute(
        &self,
        t: f64,
        state: &VectorN<f64, Self::STMSize>,
    ) -> (VectorN<f64, Self::STMSize>, MatrixMN<f64, Self::STMSize, Self::STMSize>)
    where
        DefaultAllocator: Allocator<f64, Self::STMSize> + Allocator<f64, Self::STMSize, Self::STMSize>,
    {
        let hyperstate = hyperspace_from_vector(&state);

        let (state, grad) = self.dual_eom(t, &hyperstate);

        (state, grad)
    }
}

impl<T: AutoDiffDynamics> Linearization for T
where
    DefaultAllocator: Allocator<Dual<f64>, T::StateSize> + Allocator<Dual<f64>, T::StateSize, T::StateSize>,
{
    type StateSize = T::HyperStateSize;

    /// Returns the gradient of the dynamics at the given state.
    ///
    /// **WARNING:** Requires a prior call to self.compute() ! This is where the auto-differentiation happens.
    fn gradient(&self, _t: f64, _state: &VectorN<f64, Self::StateSize>) -> MatrixMN<f64, Self::StateSize, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize> + Allocator<f64, Self::StateSize, Self::StateSize>,
    {
        panic!("retrieve the gradient by calling self.compute(...)");
    }
}
