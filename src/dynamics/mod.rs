extern crate nalgebra as na;

use self::na::{DefaultAllocator, Dim, DimName, VectorN};
use self::na::allocator::Allocator;

pub mod celestial;

/// The `Dynamics` trait handles and stores any equation of motion *and* the state is integrated.
///
/// Its design is such that several of the provided dynamics can be combined fairly easily. However,
/// when combining the dynamics (e.g. integrating both the attitude of a spaceraft and its orbital
///  parameters), it is up to the implementor to handle time and state organization correctly.
/// For time management, I highly recommend using `hifitime` which is thoroughly validated.
pub trait Dynamics
where
    Self: Sized,
{
    /// Defines the state size for these dynamics. It must be imported from `nalgebra`.
    type StateSize: Dim + DimName;
    /// Returns the time of the current state
    fn time(&self) -> f64;

    /// Returns the current state of the dynamics so it can be integrated.
    fn state(&self) -> &VectorN<f64, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize>;

    /// Defines the equations of motion for these dynamics, or a combination of provided dynamics. XXX: Is this going to work in `derive`?!
    fn eom(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize>;

    /// Updates the internal state of the dynamics.
    ///
    /// NOTE: I do not think that this should be a Box<Self> because my understanding is that it
    /// would invalidate the previous state entirely. That would also mean that `fn state` would
    /// return a Box<VectorN> which would then transfer the ownership to whoever calls it. I do not
    /// think this is a good idea in this case. In fact, we can assume that when integrating the
    /// attitude of two instruments, we might need to "read" the attitude of the spacecraft, "read"
    /// that of an instrument, and compute the DCM between both. If the first attitude lost ownership
    /// of its state, it wouldn't be able to integrate it anymore.
    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>)
    where
        DefaultAllocator: Allocator<f64, Self::StateSize>;
}
