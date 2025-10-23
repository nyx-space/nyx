/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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

use anise::errors::{AlmanacError, PhysicsError};
pub use anise::prelude::*;

pub use crate::cosmic::{DragData, GuidanceMode, Mass, SRPData, Spacecraft};
use crate::dynamics::DynamicsError;
pub use crate::errors::NyxError;
use crate::errors::StateError;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector};
use crate::md::StateParameter;
use snafu::Snafu;
use std::fmt;

/// A trait allowing for something to have an epoch
pub trait TimeTagged {
    /// Retrieve the Epoch
    fn epoch(&self) -> Epoch;
    /// Set the Epoch
    fn set_epoch(&mut self, epoch: Epoch);
}

/// A trait for generate propagation and estimation state.
/// The first parameter is the size of the state, the second is the size of the propagated state including STM and extra items.
pub trait State: Default + Copy + PartialEq + fmt::Display + fmt::LowerExp + Send + Sync
where
    Self: Sized,
    DefaultAllocator:
        Allocator<Self::Size> + Allocator<Self::Size, Self::Size> + Allocator<Self::VecLength>,
{
    /// Size of the state and its STM
    type Size: DimName;
    type VecLength: DimName;

    /// Initialize an empty state
    /// By default, this is not implemented. This function must be implemented when filtering on this state.
    fn zeros() -> Self {
        unimplemented!()
    }

    /// Return this state as a vector for the propagation/estimation
    fn to_vector(&self) -> OVector<f64, Self::VecLength>;

    /// Returns strictly the state vector without any STM, if set.
    fn to_state_vector(&self) -> OVector<f64, Self::Size> {
        OVector::<f64, Self::Size>::from_iterator(
            self.to_vector().iter().copied().take(Self::Size::DIM),
        )
    }

    /// Return the state's _step_ state transition matrix.
    /// By default, this is not implemented. This function must be implemented when filtering on this state.
    fn stm(&self) -> Result<OMatrix<f64, Self::Size, Self::Size>, DynamicsError> {
        Err(DynamicsError::StateTransitionMatrixUnset)
    }

    /// Copies the current state but sets the STM to identity
    fn with_stm(self) -> Self;

    /// Resets the STM, unimplemented by default.
    fn reset_stm(&mut self) {
        unimplemented!()
    }

    /// Unsets the STM for this state
    fn unset_stm(&mut self);

    /// Set this state
    fn set(&mut self, epoch: Epoch, vector: &OVector<f64, Self::VecLength>);

    /// Reconstruct a new State from the provided delta time in seconds compared to the current state
    /// and with the provided vector.
    fn set_with_delta_seconds(
        mut self,
        delta_t_s: f64,
        vector: &OVector<f64, Self::VecLength>,
    ) -> Self
    where
        DefaultAllocator: Allocator<Self::VecLength>,
    {
        self.set(self.epoch() + delta_t_s, vector);
        self
    }

    /// Retrieve the Epoch
    fn epoch(&self) -> Epoch;

    /// Set the Epoch
    fn set_epoch(&mut self, epoch: Epoch);

    /// By default, this is not implemented. This function must be implemented when filtering on this state.
    fn add(self, _other: OVector<f64, Self::Size>) -> Self {
        unimplemented!()
    }

    /// Return the value of the parameter, returns an error by default
    fn value(&self, param: StateParameter) -> Result<f64, StateError> {
        Err(StateError::Unavailable { param })
    }

    /// Allows setting the value of the given parameter.
    /// NOTE: Most parameters where the `value` is available CANNOT be also set for that parameter (it's a much harder problem!)
    fn set_value(&mut self, param: StateParameter, _val: f64) -> Result<(), StateError> {
        Err(StateError::Unavailable { param })
    }

    /// Returns a copy of the orbit
    fn orbit(&self) -> Orbit;

    /// Modifies this state's orbit
    fn set_orbit(&mut self, _orbit: Orbit) {}
}

pub fn assert_orbit_eq_or_abs(left: &Orbit, right: &Orbit, epsilon: f64, msg: &str) {
    if !left.eq_within(right, epsilon, epsilon) {
        panic!(
            r#"assertion failed: {msg}
  left: `{left:?}`,
 right: `{right:?}`"#
        )
    }
}

#[derive(Debug, PartialEq, Snafu)]
#[snafu(visibility(pub(crate)))]
pub enum AstroError {
    #[snafu(display("B Plane jacobian invariant must be either VX, VY or VZ"))]
    BPlaneInvariant,
    #[snafu(display("operation requires a local frame"))]
    NotLocalFrame,
    #[snafu(display("partial derivatives not defined for this parameter"))]
    PartialsUndefined,
    #[snafu(display("Orbit is not hyperbolic so there is no hyperbolic anomaly."))]
    NotHyperbolic,
    #[snafu(display("physics error occured during astro computation: {source}"))]
    AstroPhysics { source: PhysicsError },
    #[snafu(display("ANISE Almanac error occured during astro computation: {source}"))]
    AstroAlmanac {
        #[snafu(source(from(AlmanacError, Box::new)))]
        source: Box<AlmanacError>,
    },
}

// Re-Export OrbitDual
mod orbitdual;
pub use self::orbitdual::*;

// Re-Export B Plane
mod bplane;
pub use self::bplane::*;

// Re-Export spacecraft
mod spacecraft;
pub use self::spacecraft::*;

/// The eclipse module allows finding eclipses and (conversely) visibility between a state and another one (e.g. a planet or the Sun).
pub mod eclipse;

/// Speed of light in meters per second
pub const SPEED_OF_LIGHT_M_S: f64 = SPEED_OF_LIGHT_KM_S * 1e3;
pub use anise::constants::SPEED_OF_LIGHT_KM_S;

/// Astronomical unit, in kilometers, according to the [IAU](https://www.iau.org/public/themes/measuring/).
pub const AU: f64 = 149_597_870.700;

/// From NIST special publication 330, 2008 edition, in meters per second squared
pub const STD_GRAVITY: f64 = 9.80665;
