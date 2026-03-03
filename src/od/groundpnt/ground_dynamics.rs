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

use crate::dynamics::{Dynamics, DynamicsError};
use crate::od::groundpnt::GroundAsset;
use crate::State;
use anise::prelude::Almanac;
use hyperdual::{hyperspace_from_vector, OHyperdual};
use nalgebra::allocator::Allocator;
use nalgebra::{Const, DefaultAllocator, Matrix6, OMatrix, OVector, Vector6};
use std::sync::Arc;

#[derive(Clone)]
pub struct GroundDynamics {}

impl Dynamics for GroundDynamics {
    type StateType = GroundAsset;
    type HyperdualSize = Const<6>;

    fn eom(
        &self,
        _delta_t: f64,
        _state_vec: &OVector<f64, <Self::StateType as crate::State>::VecLength>,
        state_ctx: &Self::StateType,
        _almanac: Arc<Almanac>,
    ) -> Result<OVector<f64, <Self::StateType as crate::State>::VecLength>, DynamicsError>
    where
        DefaultAllocator: Allocator<<Self::StateType as crate::State>::VecLength>,
    {
        let d_x = Vector6::from_iterator([
            state_ctx.latitude_vel_deg_s,
            state_ctx.longitude_vel_deg_s,
            state_ctx.height_vel_km_s,
            0.0,
            0.0,
            0.0,
        ]);

        Ok(OVector::<f64, Const<42>>::from_iterator(
            d_x.iter()
                .chain(OVector::<f64, Const<36>>::zeros().iter())
                .cloned(),
        ))
    }

    fn dual_eom(
        &self,
        _delta_t: f64,
        osc: &Self::StateType,
        _almanac: Arc<Almanac>,
    ) -> Result<
        (
            OVector<f64, <Self::StateType as crate::State>::Size>,
            OMatrix<
                f64,
                <Self::StateType as crate::State>::Size,
                <Self::StateType as crate::State>::Size,
            >,
        ),
        DynamicsError,
    >
    where
        DefaultAllocator: Allocator<Self::HyperdualSize>
            + Allocator<<Self::StateType as crate::State>::Size>
            + Allocator<
                <Self::StateType as crate::State>::Size,
                <Self::StateType as crate::State>::Size,
            >,
        nalgebra::Owned<f64, Self::HyperdualSize>: Copy,
    {
        // Extract data from hyperspace
        // Build full state vector with partials in the right position (hence building with all six components)
        let state: Vector6<OHyperdual<f64, Const<7>>> =
            hyperspace_from_vector(&osc.to_state_vector());

        let velocity = state.fixed_rows::<3>(3).into_owned();

        // Extract result into Vector6 and Matrix6
        let mut dx = Vector6::zeros();
        let mut grad = Matrix6::zeros();
        for i in 0..6 {
            dx[i] = if i < 3 { velocity[i].real() } else { 0.0 };
            for j in 1..7 {
                grad[(i, j - 1)] = if i < 3 { velocity[i][j] } else { 0.0 };
            }
        }

        // This function returns the time derivative of each function. The propagator will add this to the state vector (which has the previous STM).
        // This is why we don't multiply the gradient (A matrix) with the previous STM
        Ok((dx, grad))
    }
}
