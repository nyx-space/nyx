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

use crate::cosmic::{AstroError, Orbit};
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, Matrix3, Matrix4x3, OMatrix, OVector, Vector3};
use crate::State;
use anise::almanac::planetary::PlanetaryDataError;
use anise::almanac::Almanac;
use anise::errors::AlmanacError;
use hyperdual::Owned;
use snafu::Snafu;

use std::fmt;
use std::sync::Arc;

pub use crate::errors::NyxError;

/// Cartesian-based orbital dynamics.
///
/// Ensure coordinate frames match or perform transformations when combining dynamics.
pub mod orbital;
use self::guidance::GuidanceError;
pub use self::orbital::*;

/// Spacecraft dynamics, including propulsion and maneuvers.
pub mod spacecraft;
pub use self::spacecraft::*;

/// Guidance laws.
pub mod guidance;

/// Velocity change controllers.
pub mod deltavctrl;

/// Solar radiation pressure models.
pub mod solarpressure;
pub use self::solarpressure::*;

/// Atmospheric drag models.
pub mod drag;
pub use self::drag::*;

/// Spherical harmonic gravity models.
///
/// Supports loading models from PDS, EGM2008, and GMAT's COF files.
pub mod sph_harmonics;
pub use self::sph_harmonics::*;

/// A trait for models with equations of motion that can be integrated.
///
/// This trait is designed for composition, allowing different dynamics to be combined.
/// When combining dynamics, ensure that time and state are handled consistently.
/// `hifitime` is recommended for time management.
#[allow(clippy::type_complexity)]
pub trait Dynamics: Clone + Sync + Send
where
    DefaultAllocator: Allocator<<Self::StateType as State>::Size>
        + Allocator<<Self::StateType as State>::VecLength>
        + Allocator<<Self::StateType as State>::Size, <Self::StateType as State>::Size>,
{
    /// The state of the associated hyperdual state, almost always StateType + U1
    type HyperdualSize: DimName;
    type StateType: State;

    /// Defines the equations of motion.
    ///
    /// - `delta_t`: Time in seconds past the context epoch.
    /// - `state_vec`: The state vector, which changes at each integration step.
    /// - `state_ctx`: The state context, used to rebuild the state from the state vector.
    fn eom(
        &self,
        delta_t: f64,
        state_vec: &OVector<f64, <Self::StateType as State>::VecLength>,
        state_ctx: &Self::StateType,
        almanac: Arc<Almanac>,
    ) -> Result<OVector<f64, <Self::StateType as State>::VecLength>, DynamicsError>
    where
        DefaultAllocator: Allocator<<Self::StateType as State>::VecLength>;

    /// Defines the equations of motion for dual numbers, enabling automatic differentiation.
    ///
    /// If differentiation is not supported, this function should prevent initialization with a context that has an STM defined.
    fn dual_eom(
        &self,
        _delta_t: f64,
        _osculating_state: &Self::StateType,
        _almanac: Arc<Almanac>,
    ) -> Result<
        (
            OVector<f64, <Self::StateType as State>::Size>,
            OMatrix<f64, <Self::StateType as State>::Size, <Self::StateType as State>::Size>,
        ),
        DynamicsError,
    >
    where
        DefaultAllocator: Allocator<Self::HyperdualSize>
            + Allocator<<Self::StateType as State>::Size>
            + Allocator<<Self::StateType as State>::Size, <Self::StateType as State>::Size>,
        Owned<f64, Self::HyperdualSize>: Copy,
    {
        Err(DynamicsError::StateTransitionMatrixUnset)
    }

    /// Performs final changes after each successful integration step.
    ///
    /// Also called before the first integration step to update the initial state if needed.
    fn finally(
        &self,
        next_state: Self::StateType,
        _almanac: Arc<Almanac>,
    ) -> Result<Self::StateType, DynamicsError> {
        Ok(next_state)
    }
}

/// A trait for immutable dynamics that return a force (e.g., solar radiation pressure, drag).
///
/// The force is divided by the spacecraft's mass to compute acceleration (F=ma).
pub trait ForceModel: Send + Sync + fmt::Display {
    /// Returns the estimation index if a parameter of this force model is stored in the spacecraft state.
    fn estimation_index(&self) -> Option<usize>;

    /// Defines the equations of motion for this force model.
    fn eom(&self, ctx: &Spacecraft, almanac: Arc<Almanac>) -> Result<Vector3<f64>, DynamicsError>;

    /// Defines the partial derivatives of the equations of motion.
    ///
    /// The last row corresponds to the partials of the parameter of this force model with respect to position, which only applies to conservative forces.
    fn dual_eom(
        &self,
        osc_ctx: &Spacecraft,
        almanac: Arc<Almanac>,
    ) -> Result<(Vector3<f64>, Matrix4x3<f64>), DynamicsError>;
}

/// A trait for immutable dynamics that return an acceleration (e.g., spherical harmonics).
pub trait AccelModel: Send + Sync + fmt::Display {
    /// Defines the equations of motion for this acceleration model.
    fn eom(&self, osc: &Orbit, almanac: Arc<Almanac>) -> Result<Vector3<f64>, DynamicsError>;

    /// Defines the partial derivatives of the equations of motion.
    fn dual_eom(
        &self,
        osc_ctx: &Orbit,
        almanac: Arc<Almanac>,
    ) -> Result<(Vector3<f64>, Matrix3<f64>), DynamicsError>;
}

/// Dynamical model errors.
#[derive(Debug, PartialEq, Snafu)]
#[snafu(visibility(pub(crate)))]
pub enum DynamicsError {
    /// Fuel exhausted.
    #[snafu(display("fuel exhausted at {sc}"))]
    FuelExhausted { sc: Box<Spacecraft> },
    /// State Transition Matrix (STM) was expected but not set.
    #[snafu(display("expected STM to be set"))]
    StateTransitionMatrixUnset,
    /// Astrodynamics error.
    #[snafu(display("dynamical model encountered an astro error: {source}"))]
    DynamicsAstro { source: AstroError },
    /// Guidance error.
    #[snafu(display("dynamical model encountered an issue with the guidance: {source}"))]
    DynamicsGuidance { source: GuidanceError },
    /// Almanac error.
    #[snafu(display("dynamical model issue due to Almanac: {action} {source}"))]
    DynamicsAlmanacError {
        action: &'static str,
        #[snafu(source(from(AlmanacError, Box::new)))]
        source: Box<AlmanacError>,
    },
    /// Planetary data error.
    #[snafu(display("dynamical model issue due to planetary data: {action} {source}"))]
    DynamicsPlanetaryError {
        action: &'static str,
        source: PlanetaryDataError,
    },
}
