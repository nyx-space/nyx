extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, Dim, DimName, MatrixNM, VectorN};

pub mod kalman;

/// The `Dynamics` trait handles and stores any equation of motion *and* the state is integrated.
///
/// Its design is such that several of the provided dynamics can be combined fairly easily. However,
/// when combining the dynamics (e.g. integrating both the attitude of a spaceraft and its orbital
///  parameters), it is up to the implementor to handle time and state organization correctly.
/// For time management, I highly recommend using `hifitime` which is thoroughly validated.
pub trait Linearization
where
    Self: Sized,
{
    /// Defines the state size for these dynamics. It must be imported from `nalgebra`.
    type StateSize: Dim + DimName;
    /// Defines the equations of motion for these dynamics, or a combination of provided dynamics.
    fn gradient(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> MatrixNM<f64, Self::StateSize, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize>;
}
