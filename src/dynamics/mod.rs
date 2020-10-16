extern crate hyperdual;

use self::hyperdual::{hyperspace_from_vector, Hyperdual, Owned};
use crate::celestia::{Frame, State};
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, DimName, MatrixMN, Vector3, VectorN, U3, U7};
use crate::time::Epoch;

pub use crate::errors::NyxError;

/// The orbital module handles all Cartesian based orbital dynamics.
///
/// It is up to the engineer to ensure that the coordinate frames of the different dynamics borrowed
/// from this module match, or perform the appropriate coordinate transformations.
pub mod orbital;

/// The gravity module handles spherical harmonics only. It _must_ be combined with a OrbitalDynamics dynamics
///
/// This module allows loading gravity models from [PDS](http://pds-geosciences.wustl.edu/), [EGM2008](http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/) and GMAT's own COF files.
// pub mod gravity;

/// The drag module handles drag in a very basic fashion. Do not use for high fidelity dynamics.
// pub mod drag;

/// The angular momentum module handles all angular momentum dynamics.
///
/// Note that this module does not handle attitude parameters or control. Refer to the relevant modules.
pub mod momentum;

/// The spacecraft module allows for simulation of spacecraft dynamics in general, including propulsion/maneuvers.
pub mod spacecraft;

/// Defines what a propulsion subsystem must implement, with some common propulsion systems.
pub mod propulsion;

/// Defines a few examples of thrust controllers.
pub mod thrustctrl;

/// Defines some velocity change controllers.
pub mod deltavctrl;

/// Defines a MissionArc, i.e. a section of a spacecraft mission. Enables maneuvers design.
pub mod missionarc;

/// Defines solar radiation pressure models
pub mod solarpressure;

/// Define drag models
pub mod drag;

/// Define the spherical harmonic models.
pub mod sph_harmonics;

/// The `Dynamics` trait handles and stores any equation of motion *and* the state is integrated.
///
/// Its design is such that several of the provided dynamics can be combined fairly easily. However,
/// when combining the dynamics (e.g. integrating both the attitude of a spaceraft and its orbital
///  parameters), it is up to the implementor to handle time and state organization correctly.
/// For time management, I highly recommend using `hifitime` which is thoroughly validated.
pub trait Dynamics {
    /// Defines the state size for these dynamics. It must be imported from `nalgebra`.
    type StateSize: DimName;
    /// Defines the type which will be published on the propagator channel
    type StateType: Copy;
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
    fn set_state(
        &mut self,
        new_t: f64,
        new_state: &VectorN<f64, Self::StateSize>,
    ) -> Result<(), NyxError>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize>;

    /// Returns the state of the dynamics
    fn state(&self) -> Self::StateType;
}

/// A trait to specify that given dynamics support linearization, and can be used for state transition matrix computation.
pub trait AutoDiff: Send + Sync {
    /// Defines the state size of the estimated state
    type HyperStateSize: DimName;
    type STMSize: DimName;

    /// Computes both the state and the gradient of the dynamics.
    fn eom_grad(
        &self,
        epoch: Epoch,
        integr_frame: Frame,
        state: &VectorN<f64, Self::STMSize>,
    ) -> (
        VectorN<f64, Self::STMSize>,
        MatrixMN<f64, Self::STMSize, Self::STMSize>,
    )
    where
        DefaultAllocator: Allocator<f64, Self::STMSize>
            + Allocator<f64, Self::STMSize, Self::STMSize>
            + Allocator<f64, Self::HyperStateSize>
            + Allocator<Hyperdual<f64, Self::HyperStateSize>, Self::STMSize>,
        Owned<f64, Self::HyperStateSize>: Copy,
    {
        let hyperstate: VectorN<Hyperdual<f64, Self::HyperStateSize>, Self::STMSize> =
            hyperspace_from_vector(&state);

        let (state, grad) = self.dual_eom(epoch, integr_frame, &hyperstate);

        (state, grad)
    }

    /// Defines the equations of motion for Dual numbers for these dynamics.
    fn dual_eom(
        &self,
        epoch: Epoch,
        integr_frame: Frame,
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
}

/// The `ForceModel` trait handles immutable dynamics which return a force. Those will be divided by the mass of the spacecraft to compute the acceleration (F = ma).
///
/// Examples include Solar Radiation Pressure, drag, etc., i.e. forces which do not need to save the current state, only act on it.
pub trait ForceModel: AutoDiff<STMSize = U3, HyperStateSize = U7> + Send + Sync {
    /// Defines the equations of motion for this force model from the provided osculating state.
    fn eom(&self, osc: &State) -> Vector3<f64>;
}

/// The `AccelModel` trait handles immutable dynamics which return an acceleration. Those can be added directly to Celestial Dynamics for example.
///
/// Examples include spherical harmonics, i.e. accelerations which do not need to save the current state, only act on it.
pub trait AccelModel: Send + Sync {
    /// Defines the equations of motion for this force model from the provided osculating state in the integration frame.
    fn eom(&self, osc: &State) -> Vector3<f64>;
}
