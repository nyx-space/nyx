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

pub use crate::errors::NyxError;

/// The orbital module handles all Cartesian based orbital dynamics.
///
/// It is up to the engineer to ensure that the coordinate frames of the different dynamics borrowed
/// from this module match, or perform the appropriate coordinate transformations.
pub mod orbital;
use self::guidance::GuidanceError;
pub use self::orbital::*;

/// The spacecraft module allows for simulation of spacecraft dynamics in general, including propulsion/maneuvers.
pub mod spacecraft;
pub use self::spacecraft::*;

/// Defines a few examples of guidance laws.
pub mod guidance;

/// Defines some velocity change controllers.
pub mod deltavctrl;

/// Defines solar radiation pressure models
pub mod solarpressure;
pub use self::solarpressure::*;

/// The drag module handles drag in a very basic fashion. Do not use for high fidelity dynamics.
pub mod drag;
pub use self::drag::*;

/// Define the gravity field models.
/// This module allows loading gravity models from [PDS](http://pds-geosciences.wustl.edu/), [EGM2008](http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/) and GMAT's own COF files.
pub mod gravity_field;
pub use self::gravity_field::*;

/// Define the solid tide models.
#[cfg(feature = "premium")]
pub mod solid_tides;
#[cfg(feature = "premium")]
pub use self::solid_tides::*;

pub mod sequence;

/// The `Dynamics` trait handles and stores any equation of motion *and* the state is integrated.
///
/// Its design is such that several of the provided dynamics can be combined fairly easily. However,
/// when combining the dynamics (e.g. integrating both the attitude of a spaceraft and its orbital
///  parameters), it is up to the implementor to handle time and state organization correctly.
/// For time management, I highly recommend using `hifitime` which is thoroughly validated.
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

    /// Defines the equations of motion for these dynamics, or a combination of provided dynamics.
    /// The time delta_t is in **seconds** PAST the context epoch. The state vector is the state which
    /// changes for every intermediate step of the integration. The state context is the state of
    /// what is being propagated, it should allow rebuilding a new state context from the
    /// provided state vector.
    fn eom(
        &self,
        delta_t: f64,
        state_vec: &OVector<f64, <Self::StateType as State>::VecLength>,
        state_ctx: &Self::StateType,
        almanac: &Almanac,
    ) -> Result<OVector<f64, <Self::StateType as State>::VecLength>, DynamicsError>
    where
        DefaultAllocator: Allocator<<Self::StateType as State>::VecLength>;

    /// Defines the equations of motion for Dual numbers for these dynamics.
    /// _All_ dynamics need to allow for automatic differentiation. However, if differentiation is not supported,
    /// then the dynamics should prevent initialization with a context which has an STM defined.
    fn dual_eom(
        &self,
        _delta_t: f64,
        _osculating_state: &Self::StateType,
        _almanac: &Almanac,
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

    /// Optionally performs some final changes after each successful integration of the equations of motion.
    /// For example, this can be used to update the Guidance mode.
    /// NOTE: This function is also called just prior to very first integration step in order to update the initial state if needed.
    fn finally(
        &self,
        next_state: Self::StateType,
        _almanac: &Almanac,
    ) -> Result<Self::StateType, DynamicsError> {
        Ok(next_state)
    }
}

/// Evaluates mass-dependent forces acting on a spacecraft.
///
/// Implementations of `ForceModel` operate on a full [`Spacecraft`] context to account for
/// physical properties such as mass, cross-sectional area, and surface coefficients (e.g.,
/// aerodynamic drag, solar radiation pressure). The evaluated force vector $\mathbf{F}$ is
/// scaled by the inverse spacecraft mass ($1/m$) and unit conversions during numerical
/// integration to yield acceleration in $\text{km/s}^2$.
pub trait ForceModel: Send + Sync + fmt::Display {
    /// Returns the state-vector index of an estimable parameter associated with this model, if configured.
    ///
    /// For example, if a drag coefficient ($C_D$) or radiation pressure coefficient ($C_R$) is
    /// actively estimated in the filter state, this returns its corresponding index in the state vector.
    fn estimation_index(&self) -> Option<usize>;

    /// Evaluates the force vector $\mathbf{F}$ at the provided state and epoch. Must be in kg*km/s^2 (or kN).
    fn eom(&self, ctx: &Spacecraft, almanac: &Almanac) -> Result<Vector3<f64>, DynamicsError>;

    /// Evaluates the nominal force vector $\mathbf{F}$ and its partial derivatives for State Transition Matrix (STM) propagation.
    ///
    /// Returns a tuple containing:
    /// 1. The nominal force vector $\mathbf{F}$.
    /// 2. The $4 \times 3$ Jacobian matrix containing spatial partial derivatives ($\partial \mathbf{F}/\partial \mathbf{r}$)
    ///    in the first three rows, and parameter partial derivatives ($\partial \mathbf{F}/\partial p$) in the fourth row.
    fn gradient(
        &self,
        osc_ctx: &Spacecraft,
        almanac: &Almanac,
    ) -> Result<(Vector3<f64>, Matrix4x3<f64>), DynamicsError>;
}

/// Evaluates mass-independent accelerations acting directly on an orbit.
///
/// Unlike [`ForceModel`], implementations of `AccelModel` operate strictly on an [`Orbit`]
/// state because the evaluated acceleration vector $\mathbf{a}$ is independent of spacecraft mass
/// or surface geometry (e.g., central-body spherical harmonics, point-mass third-body gravity).
pub trait AccelModel: Send + Sync + fmt::Display {
    /// Evaluates the acceleration vector $\mathbf{a}$ ($\text{km/s}^2$) in the integration frame at the provided orbital state and epoch.
    fn eom(&self, osc: &Orbit, almanac: &Almanac) -> Result<Vector3<f64>, DynamicsError>;

    /// Evaluates the nominal acceleration vector $\mathbf{a}$ and its spatial partial derivatives for State Transition Matrix (STM) propagation.
    ///
    /// Returns a tuple containing:
    /// 1. The nominal acceleration vector $\mathbf{a}$ ($\text{km/s}^2$).
    /// 2. The $3 \times 3$ Jacobian matrix of spatial partial derivatives ($\partial \mathbf{a}/\partial \mathbf{r}$, in $\text{s}^{-2}$).
    fn gradient(
        &self,
        osc_ctx: &Orbit,
        almanac: &Almanac,
    ) -> Result<(Vector3<f64>, Matrix3<f64>), DynamicsError>;
}

/// Stores dynamical model errors
#[derive(Debug, PartialEq, Snafu)]
#[snafu(visibility(pub(crate)))]
pub enum DynamicsError {
    #[snafu(display("spacecraft total mass is zero, cannot compute any force model"))]
    MasslessSpacecraft,
    /// Fuel exhausted at the provided spacecraft state
    #[snafu(display("fuel exhausted at {sc}"))]
    FuelExhausted { sc: Box<Spacecraft> },
    #[snafu(display("expected STM to be set"))]
    StateTransitionMatrixUnset,
    #[snafu(display("dynamical model encountered an astro error: {source}"))]
    DynamicsAstro { source: AstroError },
    #[snafu(display("dynamical model encountered an issue with the guidance: {source}"))]
    DynamicsGuidance { source: GuidanceError },
    #[snafu(display("dynamical model issue due to Almanac: {action} {source}"))]
    DynamicsAlmanacError {
        action: &'static str,
        #[snafu(source(from(AlmanacError, Box::new)))]
        source: Box<AlmanacError>,
    },
    #[snafu(display("dynamical model issue due to planetary data: {action} {source}"))]
    DynamicsPlanetaryError {
        action: &'static str,
        #[snafu(source(from(PlanetaryDataError, Box::new)))]
        source: Box<PlanetaryDataError>,
    },
}
