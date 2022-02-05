/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::hyperdual::linalg::norm;
use super::hyperdual::{extract_jacobian_and_result, hyperspace_from_vector, Float, OHyperdual};
use super::{Dynamics, NyxError};
use crate::linalg::{Const, Matrix3, Matrix6, OVector, SVector, Vector3, Vector6};
use crate::State;
use std::f64;
use std::fmt;
use std::sync::Arc;

pub use super::sph_harmonics::Harmonics;

/// `OrbitalDynamics` provides the equations of motion for any celestial dynamic, without state transition matrix computation.
#[derive(Clone)]
pub struct Ecrv<const N: usize> {
    pub tau: SVector<f64, N>,
}

impl<const N: usize> fmt::Display for Ecrv<N> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Ecrv {}", N)
    }
}

impl<const N: usize> Dynamics for Ecrv<N> {
    type HyperdualSize = Const<7>;
    type StateType = SVector<f64, N>;

    fn eom(
        &self,
        delta_t_s: f64,
        state: &OVector<f64, Const<42>>,
        ctx: &SVector<f64, N>,
    ) -> Result<OVector<f64, Const<42>>, NyxError> {
        let osc = ctx.set_with_delta_seconds(delta_t_s, state);
        let (new_state, new_stm) = if ctx.stm.is_some() {
            let (state, grad) = self.dual_eom(delta_t_s, &osc)?;

            let stm_dt = ctx.stm()? * grad;
            // Rebuild the STM as a vector.
            let stm_as_vec = OVector::<f64, Const<36>>::from_column_slice(stm_dt.as_slice());
            (state, stm_as_vec)
        } else {
            // Still return something of size 42, but the STM will be zeros.
            let body_acceleration = (-osc.frame.gm() / osc.rmag().powi(3)) * osc.radius();
            let mut d_x = Vector6::from_iterator(
                osc.velocity()
                    .iter()
                    .chain(body_acceleration.iter())
                    .cloned(),
            );

            // Apply the acceleration models
            for model in &self.accel_models {
                let model_acc = model.eom(&osc)?;
                for i in 0..3 {
                    d_x[i + 3] += model_acc[i];
                }
            }

            (d_x, OVector::<f64, Const<36>>::zeros())
        };
        Ok(OVector::<f64, Const<42>>::from_iterator(
            new_state.iter().chain(new_stm.iter()).cloned(),
        ))
    }

    fn dual_eom(
        &self,
        _delta_t_s: f64,
        osc: &Orbit,
    ) -> Result<(Vector6<f64>, Matrix6<f64>), NyxError> {
        // Extract data from hyperspace
        // Build full state vector with partials in the right position (hence building with all six components)
        let state: Vector6<OHyperdual<f64, Const<7>>> =
            hyperspace_from_vector(&osc.to_cartesian_vec());

        let radius = state.fixed_rows::<3>(0).into_owned();
        let velocity = state.fixed_rows::<3>(3).into_owned();

        // Code up math as usual
        let rmag = norm(&radius);
        let body_acceleration =
            radius * (OHyperdual::<f64, Const<7>>::from_real(-osc.frame.gm()) / rmag.powi(3));

        // Extract result into Vector6 and Matrix6
        let mut dx = Vector6::zeros();
        let mut grad = Matrix6::zeros();
        for i in 0..6 {
            dx[i] = if i < 3 {
                velocity[i].real()
            } else {
                body_acceleration[i - 3].real()
            };
            for j in 1..7 {
                grad[(i, j - 1)] = if i < 3 {
                    velocity[i][j]
                } else {
                    body_acceleration[i - 3][j]
                };
            }
        }

        // Apply the acceleration models
        for model in &self.accel_models {
            // let (model_acc, model_grad) = model.dual_eom(&radius, osc)?;
            let (model_acc, model_grad) = model.dual_eom(osc)?;
            for i in 0..3 {
                dx[i + 3] += model_acc[i];
                for j in 1..4 {
                    grad[(i + 3, j - 1)] += model_grad[(i, j - 1)];
                }
            }
        }

        // This function returns the time derivative of each function. The propagator will add this to the state vector (which has the previous STM).
        // This is why we don't multiply the gradient (A matrix) with the previous STM
        Ok((dx, grad))
    }
}
