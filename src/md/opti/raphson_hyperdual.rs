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

use super::solution::TargeterSolution;
use crate::errors::TargetingError;
use crate::linalg::{DMatrix, SVector};
use crate::md::ui::*;
use crate::md::StateParameter;
pub use crate::md::{Variable, Vary};
use crate::propagators::error_ctrl::ErrorCtrl;
use crate::pseudo_inverse;
use crate::utils::are_eigenvalues_stable;
use std::time::Instant;

impl<'a, E: ErrorCtrl, const V: usize, const O: usize> Optimizer<'a, E, V, O> {
    /// Differential correction using hyperdual numbers for the objectives
    #[allow(clippy::comparison_chain)]
    pub fn try_achieve_dual(
        &self,
        initial_state: Spacecraft,
        correction_epoch: Epoch,
        achievement_epoch: Epoch,
    ) -> Result<TargeterSolution<V, O>, NyxError> {
        if self.objectives.is_empty() {
            return Err(NyxError::Targeter(TargetingError::UnderdeterminedProblem));
        }

        let mut is_bplane_tgt = false;
        for obj in &self.objectives {
            if obj.parameter.is_b_plane() {
                is_bplane_tgt = true;
                break;
            }
        }

        // Now we know that the problem is correctly defined, so let's propagate as is to the epoch
        // where the correction should be applied.
        let xi_start = self
            .prop
            .with(initial_state)
            .until_epoch(correction_epoch)?;

        debug!("initial_state = {}", initial_state);
        debug!("xi_start = {}", xi_start);

        let mut xi = xi_start;

        // Store the total correction in a static vector
        let mut total_correction = SVector::<f64, V>::zeros();

        // Apply the initial guess
        for (i, var) in self.variables.iter().enumerate() {
            match var.component {
                Vary::PositionX => {
                    xi.orbit.x += var.init_guess;
                }
                Vary::PositionY => {
                    xi.orbit.y += var.init_guess;
                }
                Vary::PositionZ => {
                    xi.orbit.z += var.init_guess;
                }
                Vary::VelocityX => {
                    xi.orbit.vx += var.init_guess;
                }
                Vary::VelocityY => {
                    xi.orbit.vy += var.init_guess;
                }
                Vary::VelocityZ => {
                    xi.orbit.vz += var.init_guess;
                }
                _ => {
                    return Err(NyxError::Targeter(TargetingError::UnsupportedVariable(
                        *var,
                    )))
                }
            }
            total_correction[i] += var.init_guess;
        }

        let mut prev_err_norm = std::f64::INFINITY;

        // Determine padding in debugging info
        // For the width, we find the largest desired values and multiply it by the order of magnitude of its tolerance
        let max_obj_val = self
            .objectives
            .iter()
            .map(|obj| {
                (obj.desired_value.abs().ceil() as i32
                    * 10_i32.pow(obj.tolerance.abs().log10().ceil() as u32)) as i32
            })
            .max()
            .unwrap();

        let max_obj_tol = self
            .objectives
            .iter()
            .map(|obj| obj.tolerance.log10().abs().ceil() as usize)
            .max()
            .unwrap();

        let width = f64::from(max_obj_val).log10() as usize + 2 + max_obj_tol;

        let start_instant = Instant::now();

        for it in 0..=self.iterations {
            // Now, enable the trajectory STM for this state so we can apply the correction
            xi.enable_stm();

            // Full propagation for a half period duration is slightly more precise than a step by step one with multiplications in between.
            let xf = self.prop.with(xi).until_epoch(achievement_epoch)?.orbit;

            // Check linearization
            if !are_eigenvalues_stable(xf.stm().unwrap().complex_eigenvalues()) {
                warn!(
                    "STM linearization is broken for the requested time step of {}",
                    achievement_epoch - correction_epoch
                );
            }

            let xf_dual_obj_frame = match &self.objective_frame {
                Some((frame, cosm)) => {
                    let orbit_obj_frame = cosm.frame_chg(&xf, *frame);
                    OrbitDual::from(orbit_obj_frame)
                }
                None => OrbitDual::from(xf),
            };

            // Build the error vector
            let mut err_vector = SVector::<f64, O>::zeros();
            let mut converged = true;

            // Build the B-Plane once, if needed, and always in the objective frame
            let b_plane = if is_bplane_tgt {
                Some(BPlane::from_dual(xf_dual_obj_frame)?)
            } else {
                None
            };

            // Build debugging information
            let mut objmsg = Vec::new();

            // The Jacobian includes the sensitivity of each objective with respect to each variable for the whole trajectory.
            // As such, it includes the STM of that variable for the whole propagation arc.
            let mut jac = DMatrix::from_element(self.objectives.len(), self.variables.len(), 0.0);

            for (i, obj) in self.objectives.iter().enumerate() {
                let xf_partial = if obj.parameter.is_b_plane() {
                    match obj.parameter {
                        StateParameter::BdotR => b_plane.unwrap().b_r,
                        StateParameter::BdotT => b_plane.unwrap().b_t,
                        StateParameter::BLTOF => b_plane.unwrap().ltof_s,
                        _ => unreachable!(),
                    }
                } else {
                    xf_dual_obj_frame.partial_for(&obj.parameter)?
                };

                let achieved = xf_partial.real();

                let (ok, param_err) = obj.assess_raw(achieved);
                if !ok {
                    converged = false;
                }
                err_vector[i] = param_err;

                objmsg.push(format!(
                    "\t{:?}: achieved = {:>width$.prec$}\t desired = {:>width$.prec$}\t scaled error = {:>width$.prec$}",
                    obj.parameter,
                    achieved,
                    obj.desired_value,
                    param_err, width=width, prec=max_obj_tol
                ));

                // Build the Jacobian with the partials of the objectives with respect to all of the final state parameters
                // We localize the problem in the STM.
                // TODO: VNC (how?!)
                let mut partial_vec = DMatrix::from_element(1, 6, 0.0);
                for (i, val) in [
                    xf_partial.wtr_x(),
                    xf_partial.wtr_y(),
                    xf_partial.wtr_z(),
                    xf_partial.wtr_vx(),
                    xf_partial.wtr_vy(),
                    xf_partial.wtr_vz(),
                ]
                .iter()
                .enumerate()
                {
                    partial_vec[(0, i)] = *val;
                }

                for (j, var) in self.variables.iter().enumerate() {
                    let idx = var.component.vec_index();
                    // Compute the partial of the objective over all components wrt to all of the components in the STM of the control variable.
                    let rslt = &partial_vec * xf.stm().unwrap().fixed_columns::<1>(idx);
                    jac[(i, j)] = rslt[(0, 0)];
                }
            }

            if converged {
                let conv_dur = Instant::now() - start_instant;
                let mut state = xi_start;
                // Convert the total correction from VNC back to integration frame in case that's needed.
                for (i, var) in self.variables.iter().enumerate() {
                    match var.component {
                        Vary::PositionX => state.orbit.x += total_correction[i],
                        Vary::PositionY => state.orbit.y += total_correction[i],
                        Vary::PositionZ => state.orbit.z += total_correction[i],
                        Vary::VelocityX => state.orbit.vx += total_correction[i],
                        Vary::VelocityY => state.orbit.vy += total_correction[i],
                        Vary::VelocityZ => state.orbit.vz += total_correction[i],
                        _ => {
                            return Err(NyxError::Targeter(TargetingError::UnsupportedVariable(
                                *var,
                            )))
                        }
                    }
                }

                let sol = TargeterSolution {
                    corrected_state: state,
                    achieved_state: xi_start.with_orbit(xf),
                    correction: total_correction,
                    computation_dur: conv_dur,
                    variables: self.variables,
                    achieved_errors: err_vector,
                    achieved_objectives: self.objectives,
                    iterations: it,
                };
                info!("Targeter -- CONVERGED in {} iterations", it);
                for obj in &objmsg {
                    info!("{}", obj);
                }
                return Ok(sol);
            }

            // We haven't converged yet, so let's build the error vector
            if (err_vector.norm() - prev_err_norm).abs() < 1e-10 {
                return Err(NyxError::CorrectionIneffective(
                    "No change in objective errors".to_string(),
                ));
            }
            prev_err_norm = err_vector.norm();

            debug!("Jacobian {}", jac);

            // Perform the pseudo-inverse if needed, else just inverse
            let jac_inv = pseudo_inverse!(&jac)?;

            debug!("Inverse Jacobian {}", jac_inv);

            let mut delta = jac_inv * err_vector;

            debug!("Error vector: {}\nRaw correction: {}", err_vector, delta);

            // And finally apply it to the xi
            for (i, var) in self.variables.iter().enumerate() {
                // Choose the minimum step between the provided max step and the correction.
                if delta[i].abs() > var.max_step {
                    delta[i] = var.max_step * delta[i].signum();
                } else if delta[i] > var.max_value {
                    delta[i] = var.max_value;
                } else if delta[i] < var.min_value {
                    delta[i] = var.min_value;
                }

                info!(
                    "Correction {:?} (element {}): {}",
                    var.component, i, delta[i]
                );

                match var.component {
                    Vary::PositionX => {
                        xi.orbit.x += delta[i];
                    }
                    Vary::PositionY => {
                        xi.orbit.y += delta[i];
                    }
                    Vary::PositionZ => {
                        xi.orbit.z += delta[i];
                    }
                    Vary::VelocityX => {
                        xi.orbit.vx += delta[i];
                    }
                    Vary::VelocityY => {
                        xi.orbit.vy += delta[i];
                    }
                    Vary::VelocityZ => {
                        xi.orbit.vz += delta[i];
                    }
                    _ => {
                        return Err(NyxError::Targeter(TargetingError::UnsupportedVariable(
                            *var,
                        )))
                    }
                }
            }
            total_correction += delta;
            debug!("Total correction: {:e}", total_correction);

            // Log progress
            info!("Targeter -- Iteration #{} -- {}", it, achievement_epoch);
            for obj in &objmsg {
                info!("{}", obj);
            }
        }

        Err(NyxError::MaxIterReached(format!(
            "Failed after {} iterations:\nError: {}\n\n{}",
            self.iterations, prev_err_norm, self
        )))
    }
}
