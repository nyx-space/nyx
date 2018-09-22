extern crate dual_num;
extern crate hifitime;
extern crate nalgebra as na;
extern crate serde;

use self::dual_num::Dual;
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
    fn new<F: CoordinateFrame>(
        dt: Instant,
        tx: State<F>,
        rx: State<F>,
        obs: VectorN<f64, Self::MeasurementSize>,
        visible: bool,
    ) -> Self;

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
    Self: Sized,
{
    /// Defines the state size of the estimated state
    type HyperStateSize: DimName;

    /// Defines the equations of motion for Dual numbers for these dynamics.
    fn dual_eom(
        &self,
        t: f64,
        state: &MatrixMN<Dual<f64>, Self::HyperStateSize, Self::HyperStateSize>,
    ) -> MatrixMN<Dual<f64>, Self::HyperStateSize, Self::HyperStateSize>
    where
        DefaultAllocator: Allocator<Dual<f64>, Self::HyperStateSize>
            + Allocator<Dual<f64>, Self::HyperStateSize, Self::HyperStateSize>
            + Allocator<f64, Self::HyperStateSize>
            + Allocator<f64, Self::HyperStateSize, Self::HyperStateSize>;

    /// Returns the state of the dynamics (does **not** include the STM, use Dynamics::state() for that)
    fn dynamics_state(&self) -> VectorN<f64, Self::HyperStateSize>
    where
        DefaultAllocator: Allocator<f64, Self::HyperStateSize>;

    /// Set the state of the dynamics (does **not** include the STM, use Dynamics::set_state() for that)
    fn set_dynamics_state(&mut self, state: VectorN<f64, Self::HyperStateSize>)
    where
        DefaultAllocator: Allocator<f64, Self::HyperStateSize>;

    /// Returns the gradient of the dynamics.
    fn dynamics_gradient(&self) -> MatrixMN<f64, Self::HyperStateSize, Self::HyperStateSize>
    where
        DefaultAllocator: Allocator<f64, Self::HyperStateSize, Self::HyperStateSize>;

    /// Set the gradient of the dynamics
    fn set_dynamics_gradient(&mut self, gradient: MatrixMN<f64, Self::HyperStateSize, Self::HyperStateSize>)
    where
        DefaultAllocator: Allocator<f64, Self::HyperStateSize, Self::HyperStateSize>;

    /// Computes both the state and the gradient of the dynamics. These may be accessed by the related
    /// getters.
    fn compute(&mut self, t: f64, state: &VectorN<f64, Self::HyperStateSize>)
    where
        DefaultAllocator: Allocator<Dual<f64>, Self::HyperStateSize>
            + Allocator<Dual<f64>, Self::HyperStateSize, Self::HyperStateSize>
            + Allocator<f64, Self::HyperStateSize>
            + Allocator<f64, Self::HyperStateSize, Self::HyperStateSize>,
    {
        // Create a Matrix for the hyperdual space
        let mut hyperdual_space = MatrixMN::<Dual<f64>, Self::HyperStateSize, Self::HyperStateSize>::zeros();

        for i in 0..Self::HyperStateSize::dim() {
            let mut v_i = VectorN::<Dual<f64>, Self::HyperStateSize>::zeros();

            for j in 0..Self::HyperStateSize::dim() {
                v_i[(j, 0)] = Dual::new(state[(j, 0)], if i == j { 1.0 } else { 0.0 });
            }

            hyperdual_space.set_column(i, &v_i);
        }

        let state_n_grad = self.dual_eom(t, &hyperdual_space);

        // Extract the dual part
        let mut state = VectorN::<f64, Self::HyperStateSize>::zeros();
        let mut grad = MatrixMN::<f64, Self::HyperStateSize, Self::HyperStateSize>::zeros();

        for i in 0..Self::HyperStateSize::dim() {
            for j in 0..Self::HyperStateSize::dim() {
                if j == 0 {
                    // The state is duplicated in every column
                    state[(i, 0)] = state_n_grad[(i, 0)].real();
                }

                grad[(i, j)] = state_n_grad[(i, j)].dual();
            }
        }
        self.set_dynamics_state(state);
        self.set_dynamics_gradient(grad);
    }

    /// Returns the real part of the equations of motion for these dynamics.
    fn real_eom(&mut self, t: f64, state: &VectorN<f64, Self::HyperStateSize>) -> VectorN<f64, Self::HyperStateSize>
    where
        DefaultAllocator: Allocator<Dual<f64>, Self::HyperStateSize>
            + Allocator<Dual<f64>, Self::HyperStateSize, Self::HyperStateSize>
            + Allocator<f64, Self::HyperStateSize>
            + Allocator<f64, Self::HyperStateSize, Self::HyperStateSize>,
    {
        self.compute(t, state);
        self.dynamics_state()
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
    fn gradient(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> MatrixMN<f64, Self::StateSize, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize> + Allocator<f64, Self::StateSize, Self::StateSize>,
    {
        self.dynamics_gradient()
    }
}
