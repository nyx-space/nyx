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

mod jacobian;
mod correction;

use super::solution::TargeterSolution;
use super::targeter::Targeter;
use crate::cosmic::{AstroAlmanacSnafu, AstroPhysicsSnafu};
use crate::dynamics::guidance::{LocalFrame, Maneuver, MnvrRepr};
use crate::errors::TargetingError;
use crate::linalg::{SMatrix, SVector, Vector6};
use crate::md::{prelude::*, AstroSnafu, UnderdeterminedProblemSnafu};
use crate::md::{PropSnafu, StateParameter};
pub use crate::md::Vary;
use crate::polyfit::CommonPolynomial;
use crate::pseudo_inverse;
use snafu::{ensure, ResultExt};
#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant;

impl<const V: usize, const O: usize> Targeter<'_, V, O> {
    /// Differential correction using finite differencing
    #[allow(clippy::comparison_chain)]
    pub fn try_achieve_fd(
        &self,
        initial_state: Spacecraft,
        correction_epoch: Epoch,
        achievement_epoch: Epoch,
        almanac: Arc<Almanac>,
    ) -> Result<TargeterSolution<V, O>, TargetingError> {
        ensure!(!self.objectives.is_empty(), UnderdeterminedProblemSnafu);

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
            .with(initial_state, almanac.clone())
            .until_epoch(correction_epoch)
            .context(PropSnafu)?;

        debug!("initial_state = {initial_state}");
        debug!("xi_start = {xi_start}");

        let mut xi = xi_start;
        // We'll store the initial state correction here.
        let _state_correction = Vector6::<f64>::zeros();

        // Store the total correction in Vector3
        let mut total_correction = SVector::<f64, V>::zeros();

        let mut mnvr = Maneuver {
            start: correction_epoch,
            end: achievement_epoch,
            thrust_prct: 1.0,
            representation: MnvrRepr::Angles {
                azimuth: CommonPolynomial::Quadratic(0.0, 0.0, 0.0),
                elevation: CommonPolynomial::Quadratic(0.0, 0.0, 0.0),
            },
            frame: LocalFrame::RCN,
        };

        let mut finite_burn_target = false;

        let mut initial_guess = SVector::<f64, V>::zeros();
        for (i, var) in self.variables.iter().enumerate() {
            initial_guess[i] = var.init_guess;
        }

        correction::apply_correction(
            self,
            &initial_guess,
            &mut xi,
            &mut mnvr,
            correction_epoch,
            achievement_epoch,
            true,
        )?;

        for (i, var) in self.variables.iter().enumerate() {
            if var.component.is_finite_burn() {
                finite_burn_target = true;
            }
            total_correction[i] += var.init_guess;
        }

        let mut prev_err_norm = f64::INFINITY;

        // Determine padding in debugging info
        // For the width, we find the largest desired values and multiply it by the order of magnitude of its tolerance
        let max_obj_val = self
            .objectives
            .iter()
            .map(|obj| {
                obj.desired_value.abs().ceil() as i32
                    * 10_i32.pow(obj.tolerance.abs().log10().ceil() as u32)
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

        #[cfg(not(target_arch = "wasm32"))]
        let start_instant = Instant::now();

        for it in 0..=self.iterations {
            // Modify each variable by the desired perturbation, propagate, compute the final parameter, and store how modifying that variable affects the final parameter
            let cur_xi = xi;

            // If we are targeting a finite burn, let's set propagate in several steps to make sure we don't miss the burn
            let xf = if finite_burn_target {
                info!("#{it} {mnvr}");
                let mut prop = self.prop.clone();
                let prop_opts = prop.opts;
                let pre_mnvr = prop
                    .with(cur_xi, almanac.clone())
                    .until_epoch(mnvr.start)
                    .context(PropSnafu)?;
                prop.dynamics = prop.dynamics.with_guidance_law(Arc::new(mnvr));
                prop.set_max_step(mnvr.duration());
                let post_mnvr = prop
                    .with(
                        pre_mnvr.with_guidance_mode(GuidanceMode::Thrust),
                        almanac.clone(),
                    )
                    .until_epoch(mnvr.end)
                    .context(PropSnafu)?;
                // Reset the propagator options to their previous configuration
                prop.opts = prop_opts;
                // And propagate until the achievement epoch
                prop.with(post_mnvr, almanac.clone())
                    .until_epoch(achievement_epoch)
                    .context(PropSnafu)?
                    .orbit
            } else {
                self.prop
                    .with(cur_xi, almanac.clone())
                    .until_epoch(achievement_epoch)
                    .context(PropSnafu)?
                    .orbit
            };

            let xf_dual_obj_frame = match &self.objective_frame {
                Some(frame) => {
                    let orbit_obj_frame = almanac
                        .transform_to(xf, *frame, None)
                        .context(AstroAlmanacSnafu)
                        .context(AstroSnafu)?;

                    OrbitDual::from(orbit_obj_frame)
                }
                None => OrbitDual::from(xf),
            };

            // Build the error vector
            let mut err_vector = SVector::<f64, O>::zeros();
            let mut converged = true;

            // Build the B-Plane once, if needed, and always in the objective frame
            let b_plane = if is_bplane_tgt {
                Some(BPlane::from_dual(xf_dual_obj_frame).context(AstroSnafu)?)
            } else {
                None
            };

            // Build debugging information
            let mut objmsg = Vec::with_capacity(self.objectives.len());

            // The Jacobian includes the sensitivity of each objective with respect to each variable for the whole trajectory.
            // As such, it includes the STM of that variable for the whole propagation arc.
            let mut jac = SMatrix::<f64, O, V>::zeros();

            for (i, obj) in self.objectives.iter().enumerate() {
                let partial = if obj.parameter.is_b_plane() {
                    match obj.parameter {
                        StateParameter::BdotR => b_plane.unwrap().b_r,
                        StateParameter::BdotT => b_plane.unwrap().b_t,
                        StateParameter::BLTOF => b_plane.unwrap().ltof_s,
                        _ => unreachable!(),
                    }
                } else {
                    xf_dual_obj_frame
                        .partial_for(obj.parameter)
                        .context(AstroSnafu)?
                };

                let achieved = partial.real();

                let (ok, param_err) = obj.assess_value(achieved);
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

                jacobian::compute_jacobian(
                    self,
                    &xi,
                    &mnvr,
                    finite_burn_target,
                    &almanac,
                    achievement_epoch,
                    &mut jac,
                    i,
                    obj,
                    achieved,
                    is_bplane_tgt,
                )?;
            }

            if converged {
                #[cfg(not(target_arch = "wasm32"))]
                let conv_dur = Instant::now() - start_instant;
                #[cfg(target_arch = "wasm32")]
                let conv_dur = Duration::ZERO.into();
                let mut corrected_state = xi_start;

                let mut state_correction = Vector6::<f64>::zeros();
                if !finite_burn_target {
                    for (i, var) in self.variables.iter().enumerate() {
                        state_correction[var.component.vec_index()] += total_correction[i];
                    }
                }
                // Now, let's apply the correction to the initial state
                if let Some(frame) = self.correction_frame {
                    let dcm_vnc2inertial = frame
                        .dcm_to_inertial(corrected_state.orbit)
                        .context(AstroPhysicsSnafu)
                        .context(AstroSnafu)?
                        .rot_mat
                        .transpose();

                    let velocity_correction =
                        dcm_vnc2inertial * state_correction.fixed_rows::<3>(3);
                    corrected_state.orbit.apply_dv_km_s(velocity_correction);
                } else {
                    corrected_state.orbit.radius_km +=
                        state_correction.fixed_rows::<3>(0).to_owned();
                    corrected_state.orbit.velocity_km_s +=
                        state_correction.fixed_rows::<3>(3).to_owned();
                }

                let sol = TargeterSolution {
                    corrected_state,
                    achieved_state: xi_start.with_orbit(xf),
                    correction: total_correction,
                    computation_dur: conv_dur,
                    variables: self.variables,
                    achieved_errors: err_vector,
                    achieved_objectives: self.objectives,
                    iterations: it,
                };
                // Log success as info
                if it == 1 {
                    info!("Targeter -- CONVERGED in 1 iteration");
                } else {
                    info!("Targeter -- CONVERGED in {it} iterations");
                }
                for obj in &objmsg {
                    info!("{obj}");
                }
                return Ok(sol);
            }

            // We haven't converged yet, so let's build t
            if (err_vector.norm() - prev_err_norm).abs() < 1e-10 {
                return Err(TargetingError::CorrectionIneffective {
                    prev_val: prev_err_norm,
                    cur_val: err_vector.norm(),
                    action: "Raphson targeter",
                });
            }
            prev_err_norm = err_vector.norm();

            debug!("Jacobian {jac}");

            // Perform the pseudo-inverse if needed, else just inverse
            let jac_inv = pseudo_inverse!(&jac)?;

            debug!("Inverse Jacobian {jac_inv}");

            let delta = jac_inv * err_vector;

            debug!(
                "Error vector (norm = {}): {}\nRaw correction: {}",
                err_vector.norm(),
                err_vector,
                delta
            );

            let applied_delta = correction::apply_correction(
                self,
                &delta,
                &mut xi,
                &mut mnvr,
                correction_epoch,
                achievement_epoch,
                false,
            )?;
            total_correction += applied_delta;
            total_correction += delta;
            debug!("Total correction: {total_correction:e}");

            // Log progress to debug
            info!("Targeter -- Iteration #{it} -- {achievement_epoch}");
            for obj in &objmsg {
                info!("{obj}");
            }
        }

        Err(TargetingError::TooManyIterations)
    }
}
