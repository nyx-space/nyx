extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, DimName, VectorN};

/// The celestial module handles all Cartesian based dynamics.
///
/// It is up to the engineer to ensure that the coordinate frames of the different dynamics borrowed
/// from this module match, or perform the appropriate coordinate transformations.
pub mod celestial;

/// The gravity module handles spherical harmonics only. It _must_ be combined with a CelestialDynamics dynamics
///
/// This module allows loading gravity models from [PDS](http://pds-geosciences.wustl.edu/), [EGM2008](http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/) and GMAT's own COF files.
pub mod gravity;

/// The drag module handles drag in a very basic fashion. Do not use for high fidelity dynamics.
pub mod drag;

/// The angular momentum module handles all angular momentum dynamics.
///
/// Note that this module does not handle attitude parameters or control. Refer to the relevant modules.
pub mod momentum;

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
    type StateSize: DimName;
    /// Defines the type which will be published on the propagator channel
    type StateType;
    /// Returns the time of the current state
    fn time(&self) -> f64;

    /// Returns the current state of the dynamics as a vector so it can be integrated.
    fn state_vector(&self) -> VectorN<f64, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize>;

    /// Defines the equations of motion for these dynamics, or a combination of provided dynamics.
    fn eom(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize>;

    /// Updates the internal state of the dynamics.
    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>)
    where
        DefaultAllocator: Allocator<f64, Self::StateSize>;

    /// Returns the state of the dynamics
    fn state(&self) -> Self::StateType;
}
