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

use super::optimizer::Optimizer;
use super::solution::TargeterSolution;
use crate::dynamics::guidance::Mnvr;
use crate::errors::TargetingError;
use crate::linalg::{SMatrix, SVector, Vector6};
use crate::md::rayon::prelude::*;
use crate::md::ui::*;
use crate::md::StateParameter;
pub use crate::md::{Variable, Vary};
use crate::polyfit::CommonPolynomial;
use crate::propagators::error_ctrl::ErrorCtrl;
use crate::pseudo_inverse;
use hifitime::TimeUnits;
use std::time::Instant;

impl<'a, E: ErrorCtrl, const V: usize, const O: usize> Optimizer<'a, E, V, O> {
    /// Differential correction using finite differencing
    #[allow(clippy::comparison_chain)]
    pub fn try_achieve_fd(
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
        // We'll store the initial state correction here.
        let mut state_correction = Vector6::<f64>::zeros();

        // Store the total correction in Vector3
        let mut total_correction = SVector::<f64, V>::zeros();

        let mut mnvr = Mnvr {
            start: correction_epoch,
            end: achievement_epoch,
            thrust_lvl: 1.0,
            alpha_inplane_radians: CommonPolynomial::Quadratic(0.0, 0.0, 0.0),
            delta_outofplane_radians: CommonPolynomial::Quadratic(0.0, 0.0, 0.0),
            frame: Frame::RCN,
        };

        let mut finite_burn_target = false;

        // Apply the initial guess
        for (i, var) in self.variables.iter().enumerate() {
            // Check the validity (this function will report to log and raise an error)
            var.valid()?;
            // Check that there is no attempt to target a position in a local frame
            if self.correction_frame.is_some() && var.component.vec_index() < 3 {
                // Then this is a position correction, which is not allowed if a frame is provided!
                let msg = format!(
                    "Variable is in frame {} but that frame cannot be used for a {:?} correction",
                    self.correction_frame.unwrap(),
                    var.component
                );
                error!("{}", msg);
                return Err(NyxError::Targeter(TargetingError::FrameError(msg)));
            }

            // Check that a thruster is provided since we'll be changing that and the burn duration
            if var.component.is_finite_burn() {
                if xi_start.thruster.is_none() {
                    // Can't do any conversion to finite burns without a thruster
                    return Err(NyxError::NoThrusterAvail);
                }
                finite_burn_target = true;
                // Modify the default maneuver
                match var.component {
                    Vary::Duration => mnvr.end = mnvr.start + var.init_guess.seconds(),
                    Vary::EndEpoch => mnvr.end += var.init_guess.seconds(),
                    Vary::StartEpoch => mnvr.start += var.init_guess.seconds(),
                    Vary::MnvrAlpha | Vary::MnvrAlphaDot | Vary::MnvrAlphaDDot => {
                        mnvr.alpha_inplane_radians = mnvr
                            .alpha_inplane_radians
                            .add_val_in_order(var.init_guess, var.component.vec_index())
                            .unwrap();
                    }
                    Vary::MnvrDelta | Vary::MnvrDeltaDot | Vary::MnvrDeltaDDot => {
                        mnvr.delta_outofplane_radians = mnvr
                            .delta_outofplane_radians
                            .add_val_in_order(var.init_guess, var.component.vec_index())
                            .unwrap();
                    }
                    Vary::ThrustX | Vary::ThrustY | Vary::ThrustZ => {
                        let mut vector = mnvr.direction();
                        vector[var.component.vec_index()] += var.perturbation;
                        mnvr.set_direction(vector)?;
                    }
                    Vary::ThrustRateX | Vary::ThrustRateY | Vary::ThrustRateZ => {
                        let mut vector = mnvr.rate();
                        vector[(var.component.vec_index() - 1) % 3] += var.perturbation;
                        mnvr.set_rate(vector)?;
                    }
                    Vary::ThrustAccelX | Vary::ThrustAccelY | Vary::ThrustAccelZ => {
                        let mut vector = mnvr.accel();
                        vector[(var.component.vec_index() - 1) % 3] += var.perturbation;
                        mnvr.set_accel(vector)?;
                    }
                    Vary::ThrustLevel => {
                        mnvr.thrust_lvl += var.perturbation;
                        if mnvr.thrust_lvl > 1.0 {
                            mnvr.thrust_lvl = 1.0
                        } else if mnvr.thrust_lvl < 0.0 {
                            mnvr.thrust_lvl = 0.0
                        }
                    }
                    _ => unreachable!(),
                }
                info!("Initial maneuver guess: {}", mnvr);
            } else {
                state_correction[var.component.vec_index()] += var.init_guess;
                // Now, let's apply the correction to the initial state
                if let Some(frame) = self.correction_frame {
                    // The following will error if the frame is not local
                    let dcm_vnc2inertial = xi.orbit.dcm_from_traj_frame(frame)?;
                    let velocity_correction =
                        dcm_vnc2inertial * state_correction.fixed_rows::<3>(3);
                    xi.orbit.apply_dv(velocity_correction);
                } else {
                    xi.orbit.x += state_correction[0];
                    xi.orbit.y += state_correction[1];
                    xi.orbit.z += state_correction[2];
                    xi.orbit.vx += state_correction[3];
                    xi.orbit.vy += state_correction[4];
                    xi.orbit.vz += state_correction[5];
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
            // Modify each variable by the desired perturbation, propagate, compute the final parameter, and store how modifying that variable affects the final parameter
            let cur_xi = xi;

            // If we are targeting a finite burn, let's set propagate in several steps to make sure we don't miss the burn
            let xf = if finite_burn_target {
                info!("#{} {}", it, mnvr);
                let mut prop = self.prop.clone();
                let prop_opts = prop.opts;
                let pre_mnvr = prop.with(cur_xi).until_epoch(mnvr.start)?;
                prop.dynamics = prop.dynamics.with_guidance_law(Arc::new(mnvr));
                prop.set_max_step(mnvr.duration());
                let post_mnvr = prop
                    .with(pre_mnvr.with_guidance_mode(GuidanceMode::Thrust))
                    .until_epoch(mnvr.end)?;
                // Reset the propagator options to their previous configuration
                prop.opts = prop_opts;
                // And propagate until the achievement epoch
                prop.with(post_mnvr).until_epoch(achievement_epoch)?.orbit
            } else {
                self.prop.with(cur_xi).until_epoch(achievement_epoch)?.orbit
            };

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
                    xf_dual_obj_frame.partial_for(&obj.parameter)?
                };

                let achieved = partial.real();

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

                let mut pert_calc: Vec<_> = self
                    .variables
                    .iter()
                    .enumerate()
                    .map(|(j, var)| (j, var, 0.0_f64))
                    .collect();

                pert_calc.par_iter_mut().for_each(|(_, var, jac_val)| {
                    let mut this_xi = xi;

                    let mut this_prop = self.prop.clone();
                    let mut this_mnvr = mnvr;

                    let mut opposed_pert = false;

                    if var.component.is_finite_burn() {
                        // Modify the burn itself
                        let pert = var.perturbation;
                        // Modify the maneuver, but do not change the epochs of the maneuver unless the change is greater than one millisecond
                        match var.component {
                            Vary::Duration => {
                                if pert.abs() > 1e-3 {
                                    this_mnvr.end = mnvr.start + pert.seconds()
                                }
                            }
                            Vary::EndEpoch => {
                                if pert.abs() > 1e-3 {
                                    this_mnvr.end = mnvr.end + pert.seconds()
                                }
                            }
                            Vary::StartEpoch => {
                                if pert.abs() > 1e-3 {
                                    this_mnvr.start = mnvr.start + pert.seconds()
                                }
                            }
                            Vary::MnvrAlpha | Vary::MnvrAlphaDot | Vary::MnvrAlphaDDot => {
                                this_mnvr.alpha_inplane_radians = mnvr
                                    .alpha_inplane_radians
                                    .add_val_in_order(pert, var.component.vec_index())
                                    .unwrap();
                            }
                            Vary::MnvrDelta | Vary::MnvrDeltaDot | Vary::MnvrDeltaDDot => {
                                this_mnvr.delta_outofplane_radians = mnvr
                                    .delta_outofplane_radians
                                    .add_val_in_order(pert, var.component.vec_index())
                                    .unwrap();
                            }
                            Vary::ThrustX | Vary::ThrustY | Vary::ThrustZ => {
                                let mut vector = this_mnvr.direction();
                                vector[var.component.vec_index()] += var.perturbation;
                                if !var.check_bounds(vector[var.component.vec_index()]).1 {
                                    // Oops, bound was hit, go the other way
                                    vector[var.component.vec_index()] -= 2.0 * var.perturbation;
                                    opposed_pert = true;
                                }
                                this_mnvr.set_direction(vector).unwrap();
                            }
                            Vary::ThrustRateX | Vary::ThrustRateY | Vary::ThrustRateZ => {
                                let mut vector = this_mnvr.rate();
                                vector[(var.component.vec_index() - 1) % 3] += var.perturbation;
                                if !var
                                    .check_bounds(vector[(var.component.vec_index() - 1) % 3])
                                    .1
                                {
                                    // Oops, bound was hit, go the other way
                                    vector[(var.component.vec_index() - 1) % 3] -=
                                        2.0 * var.perturbation;
                                    opposed_pert = true;
                                }
                                this_mnvr.set_rate(vector).unwrap();
                            }
                            Vary::ThrustAccelX | Vary::ThrustAccelY | Vary::ThrustAccelZ => {
                                let mut vector = this_mnvr.accel();
                                vector[(var.component.vec_index() - 1) % 3] += var.perturbation;
                                if !var
                                    .check_bounds(vector[(var.component.vec_index() - 1) % 3])
                                    .1
                                {
                                    // Oops, bound was hit, go the other way
                                    vector[(var.component.vec_index() - 1) % 3] -=
                                        2.0 * var.perturbation;
                                    opposed_pert = true;
                                }
                                this_mnvr.set_accel(vector).unwrap();
                            }
                            Vary::ThrustLevel => {
                                this_mnvr.thrust_lvl += var.perturbation;
                                if this_mnvr.thrust_lvl > 1.0 {
                                    this_mnvr.thrust_lvl = 1.0
                                } else if this_mnvr.thrust_lvl < 0.0 {
                                    this_mnvr.thrust_lvl = 0.0
                                }
                            }
                            _ => unreachable!(),
                        }
                    } else {
                        let mut state_correction = Vector6::<f64>::zeros();
                        state_correction[var.component.vec_index()] += var.perturbation;
                        // Now, let's apply the correction to the initial state
                        if let Some(frame) = self.correction_frame {
                            // The following will error if the frame is not local
                            let dcm_vnc2inertial =
                                this_xi.orbit.dcm_from_traj_frame(frame).unwrap();
                            let velocity_correction =
                                dcm_vnc2inertial * state_correction.fixed_rows::<3>(3);
                            this_xi.orbit.apply_dv(velocity_correction);
                        } else {
                            this_xi = xi + state_correction;
                        }
                    }

                    let this_xf = if finite_burn_target {
                        // Propagate normally until start of maneuver
                        let pre_mnvr = this_prop.with(cur_xi).until_epoch(this_mnvr.start).unwrap();
                        // Add this maneuver to the dynamics, make sure that we don't over-step this maneuver
                        let prop_opts = this_prop.opts;
                        this_prop.set_max_step(this_mnvr.duration());
                        this_prop.dynamics =
                            this_prop.dynamics.with_guidance_law(Arc::new(this_mnvr));
                        let post_mnvr = this_prop
                            .with(pre_mnvr.with_guidance_mode(GuidanceMode::Thrust))
                            .until_epoch(this_mnvr.end)
                            .unwrap();
                        // Reset the propagator options to their previous configuration
                        this_prop.opts = prop_opts;
                        // And propagate until the achievement epoch
                        this_prop
                            .with(post_mnvr)
                            .until_epoch(achievement_epoch)
                            .unwrap()
                            .orbit
                    } else {
                        this_prop
                            .with(this_xi)
                            .until_epoch(achievement_epoch)
                            .unwrap()
                            .orbit
                    };

                    let xf_dual_obj_frame = match &self.objective_frame {
                        Some((frame, cosm)) => {
                            let orbit_obj_frame = cosm.frame_chg(&this_xf, *frame);
                            OrbitDual::from(orbit_obj_frame)
                        }
                        None => OrbitDual::from(this_xf),
                    };

                    let b_plane = if is_bplane_tgt {
                        Some(BPlane::from_dual(xf_dual_obj_frame).unwrap())
                    } else {
                        None
                    };

                    let partial = if obj.parameter.is_b_plane() {
                        match obj.parameter {
                            StateParameter::BdotR => b_plane.unwrap().b_r,
                            StateParameter::BdotT => b_plane.unwrap().b_t,
                            StateParameter::BLTOF => b_plane.unwrap().ltof_s,
                            _ => unreachable!(),
                        }
                    } else {
                        xf_dual_obj_frame.partial_for(&obj.parameter).unwrap()
                    };

                    let this_achieved = partial.real();
                    *jac_val = (this_achieved - achieved) / var.perturbation;
                    if opposed_pert {
                        // We opposed the perturbation to ensure we don't over step a min/max bound
                        *jac_val = -*jac_val;
                    }
                });

                for (j, var, jac_val) in &pert_calc {
                    // If this is a thrust level, we oppose the value so that the correction can still be positive.
                    jac[(i, *j)] = if var.component == Vary::ThrustLevel {
                        -*jac_val
                    } else {
                        *jac_val
                    };
                }
            }

            if converged {
                let conv_dur = Instant::now() - start_instant;
                let mut corrected_state = xi_start;

                let mut state_correction = Vector6::<f64>::zeros();
                if !finite_burn_target {
                    for (i, var) in self.variables.iter().enumerate() {
                        state_correction[var.component.vec_index()] += total_correction[i];
                    }
                }
                // Now, let's apply the correction to the initial state
                if let Some(frame) = self.correction_frame {
                    let dcm_vnc2inertial = corrected_state
                        .orbit
                        .dcm_from_traj_frame(frame)
                        .unwrap()
                        .transpose();
                    let velocity_correction =
                        dcm_vnc2inertial * state_correction.fixed_rows::<3>(3);
                    corrected_state.orbit.apply_dv(velocity_correction);
                } else {
                    corrected_state.orbit = corrected_state.orbit + state_correction;
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
                    info!("Targeter -- CONVERGED in {} iterations", it);
                }
                for obj in &objmsg {
                    info!("{}", obj);
                }
                return Ok(sol);
            }

            // We haven't converged yet, so let's build t
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

            debug!(
                "Error vector (norm = {}): {}\nRaw correction: {}",
                err_vector.norm(),
                err_vector,
                delta
            );

            // And finally apply it to the xi
            let mut state_correction = Vector6::<f64>::zeros();
            for (i, var) in self.variables.iter().enumerate() {
                debug!(
                    "Correction {:?}{} (element {}): {}",
                    var.component,
                    match self.correction_frame {
                        Some(f) => format!(" in {:?}", f),
                        None => format!(""),
                    },
                    i,
                    delta[i]
                );

                let corr = delta[i];

                if var.component.is_finite_burn() {
                    // Modify the maneuver, but do not change the epochs of the maneuver unless the change is greater than one millisecond
                    match var.component {
                        Vary::Duration => {
                            if corr.abs() > 1e-3 {
                                // Check that we are within the bounds
                                let init_duration_s =
                                    (correction_epoch - achievement_epoch).in_seconds();
                                let acceptable_corr = var.apply_bounds(init_duration_s).seconds();
                                mnvr.end = mnvr.start + acceptable_corr;
                            }
                        }
                        Vary::EndEpoch => {
                            if corr.abs() > 1e-3 {
                                // Check that we are within the bounds
                                let total_end_corr =
                                    (mnvr.end + corr.seconds() - achievement_epoch).in_seconds();
                                let acceptable_corr = var.apply_bounds(total_end_corr).seconds();
                                mnvr.end += acceptable_corr;
                            }
                        }
                        Vary::StartEpoch => {
                            if corr.abs() > 1e-3 {
                                // Check that we are within the bounds
                                let total_start_corr =
                                    (mnvr.start + corr.seconds() - correction_epoch).in_seconds();
                                let acceptable_corr = var.apply_bounds(total_start_corr).seconds();
                                mnvr.end += acceptable_corr;

                                mnvr.start += corr.seconds()
                            }
                        }
                        Vary::MnvrAlpha | Vary::MnvrAlphaDot | Vary::MnvrAlphaDDot => {
                            mnvr.alpha_inplane_radians = mnvr
                                .alpha_inplane_radians
                                .add_val_in_order(corr, var.component.vec_index())
                                .unwrap();
                        }
                        Vary::MnvrDelta | Vary::MnvrDeltaDot | Vary::MnvrDeltaDDot => {
                            mnvr.delta_outofplane_radians = mnvr
                                .delta_outofplane_radians
                                .add_val_in_order(corr, var.component.vec_index())
                                .unwrap();
                        }
                        Vary::ThrustX | Vary::ThrustY | Vary::ThrustZ => {
                            let mut vector = mnvr.direction();
                            vector[var.component.vec_index()] += corr;
                            var.ensure_bounds(&mut vector[var.component.vec_index()]);
                            mnvr.set_direction(vector)?;
                        }
                        Vary::ThrustRateX | Vary::ThrustRateY | Vary::ThrustRateZ => {
                            let mut vector = mnvr.rate();
                            let idx = (var.component.vec_index() - 1) % 3;
                            vector[idx] += corr;
                            var.ensure_bounds(&mut vector[idx]);
                            mnvr.set_rate(vector)?;
                        }
                        Vary::ThrustAccelX | Vary::ThrustAccelY | Vary::ThrustAccelZ => {
                            let mut vector = mnvr.accel();
                            let idx = (var.component.vec_index() - 1) % 3;
                            vector[idx] += corr;
                            var.ensure_bounds(&mut vector[idx]);
                            mnvr.set_accel(vector)?;
                        }
                        Vary::ThrustLevel => {
                            mnvr.thrust_lvl += corr;
                            var.ensure_bounds(&mut mnvr.thrust_lvl);
                        }
                        _ => unreachable!(),
                    }
                } else {
                    // Choose the minimum step between the provided max step and the correction.
                    if delta[i].abs() > var.max_step.abs() {
                        delta[i] = var.max_step.abs() * delta[i].signum();
                    } else if delta[i] > var.max_value {
                        delta[i] = var.max_value;
                    } else if delta[i] < var.min_value {
                        delta[i] = var.min_value;
                    }
                    state_correction[var.component.vec_index()] += delta[i];
                }
            }

            // Now, let's apply the correction to the initial state
            if let Some(frame) = self.correction_frame {
                let dcm_vnc2inertial = xi.orbit.dcm_from_traj_frame(frame)?;
                let velocity_correction = dcm_vnc2inertial * state_correction.fixed_rows::<3>(3);
                xi.orbit.apply_dv(velocity_correction);
            } else {
                xi = xi + state_correction;
            }
            total_correction += delta;
            debug!("Total correction: {:e}", total_correction);

            // Log progress to debug
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
