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

use rayon::prelude::*;

use crate::dynamics::guidance::{ra_dec_from_unit_vector, Mnvr};
// use crate::errors::TargetingError;
use crate::linalg::{SMatrix, SVector, Vector3};
use crate::md::objective::Objective;
use crate::md::trajectory::InterpState;
use crate::md::ui::*;
use crate::md::StateParameter;
pub use crate::md::{Variable, Vary};
use crate::polyfit::CommonPolynomial;
use crate::propagators::error_ctrl::ErrorCtrl;
use crate::pseudo_inverse;
use crate::time::TimeUnits;
use core::f64::consts::TAU;
// use std::convert::TryInto;
// use std::fmt;

// use super::solution::TargeterSolution;

impl<'a, E: ErrorCtrl> Optimizer<'a, E, 3, 6> {
    /// Create a new Targeter which will apply an impulsive delta-v correction.
    /// The `spacecraft` _must_ be the spacecraft BEFORE the Δv is applied
    pub fn convert_impulsive_mnvr(
        spacecraft: Spacecraft,
        dv: Vector3<f64>,
        prop: &'a Propagator<'a, SpacecraftDynamics, E>,
    ) -> Result<Mnvr, NyxError> {
        if spacecraft.thruster.is_none() {
            // Can't do any conversion to finite burns without a thruster
            return Err(NyxError::NoThrusterAvail);
        }

        /* ************************* */
        /* Compute the initial guess */
        /* ************************* */
        // Calculate the u, dot u (=0) and ddot u from this state
        let u = dv / dv.norm();
        let r = spacecraft.orbit.radius();
        let rmag = spacecraft.orbit.rmag();
        let u_ddot = (3.0 * spacecraft.orbit.frame.gm() / rmag.powi(5))
            * (r.dot(&u) * r - (r.dot(&u).powi(2) * u));
        // Compute the control rates at the time of the impulsive maneuver (tdv)
        let (alpha_tdv, delta_tdv) = ra_dec_from_unit_vector(u);
        let (alpha_ddot_tdv, delta_ddot_tdv) = ra_dec_from_unit_vector(u_ddot);
        // Build the maneuver polynomial angles from these
        let alpha_inplane_radians = CommonPolynomial::Quadratic(0.0, 0.0, alpha_tdv);
        let beta_outofplane_radians = CommonPolynomial::Quadratic(0.0, 0.0, 0.0);

        // Compute a few thruster parameters
        let thruster = spacecraft.thruster.as_ref().unwrap();
        let v_exhaust_m_s = thruster.exhaust_velocity();

        let delta_tfb = ((v_exhaust_m_s * spacecraft.mass_kg()) / thruster.thrust_N)
            * (1.0 - (-dv.norm() * 1e3 / v_exhaust_m_s).exp());

        let impulse_epoch = spacecraft.epoch();
        // Build the estimated maneuver
        let mut mnvr = Mnvr {
            start: impulse_epoch - 0.5 * delta_tfb * Unit::Second,
            end: impulse_epoch + 0.5 * delta_tfb * Unit::Second,
            thrust_lvl: 1.0,
            alpha_inplane_radians,
            delta_outofplane_radians: beta_outofplane_radians,
            frame: Frame::Inertial,
        };

        println!("INITIAL GUESS\n{}\n\n", mnvr);

        /* ************************ */
        /* Compute the nominal traj */
        /* ************************ */
        // Pre-traj is the trajectory _before_ the impulsive maneuver
        let pre_sc = prop.with(spacecraft).for_duration(-60.0 * Unit::Minute)?;
        let (_, pre_traj) = prop.with(pre_sc).until_epoch_with_traj(impulse_epoch)?;
        // Post-traj is the trajectory _after_ the impulsive maneuver
        let (_, post_traj) = prop
            .with(spacecraft.with_dv(dv))
            .for_duration_with_traj(60.0 * Unit::Minute)?;

        println!("{}", pre_traj);
        println!("{}", post_traj);

        // Now let's setup the optimizer.
        let variables = [
            Variable::from(Vary::MnvrAlpha).with_initial_guess(alpha_tdv),
            Variable::from(Vary::MnvrAlphaDot),
            Variable::from(Vary::MnvrAlphaDDot).with_initial_guess(alpha_ddot_tdv),
            Variable::from(Vary::MnvrDelta).with_initial_guess(delta_tdv),
            Variable::from(Vary::MnvrDeltaDot),
            Variable::from(Vary::MnvrDeltaDDot).with_initial_guess(delta_ddot_tdv),
            // Variable::from(Vary::StartEpoch),
            // Variable::from(Vary::Duration),
        ];

        const NUM_VARIABLES: usize = 6;

        // The correction stores, in order, alpha_0, \dot{alpha_0}, \ddot{alpha_0}, beta_0, \dot{beta_0}, \ddot{beta_0}
        let mut prev_err_norm = std::f64::INFINITY;
        // The objectives will be updated if the duration of the maneuver is changed
        let mut sc_x0 = pre_traj.at(mnvr.start)?;
        let mut sc_xf_desired = post_traj.at(mnvr.end)?;
        let mut objectives = [
            Objective {
                parameter: StateParameter::X,
                desired_value: sc_xf_desired.orbit.x,
                tolerance: 1e-3,
                additive_factor: 0.0,
                multiplicative_factor: 1.0,
            },
            Objective {
                parameter: StateParameter::Y,
                desired_value: sc_xf_desired.orbit.y,
                tolerance: 1e-3,
                additive_factor: 0.0,
                multiplicative_factor: 1.0,
            },
            Objective {
                parameter: StateParameter::Z,
                desired_value: sc_xf_desired.orbit.z,
                tolerance: 1e-3,
                additive_factor: 0.0,
                multiplicative_factor: 1.0,
            },
            Objective {
                parameter: StateParameter::VX,
                desired_value: sc_xf_desired.orbit.vx,
                tolerance: 1e-3,
                additive_factor: 0.0,
                multiplicative_factor: 1e-3,
            },
            Objective {
                parameter: StateParameter::VY,
                desired_value: sc_xf_desired.orbit.vy,
                tolerance: 1e-3,
                additive_factor: 0.0,
                multiplicative_factor: 1e-3,
            },
            Objective {
                parameter: StateParameter::VZ,
                desired_value: sc_xf_desired.orbit.vz,
                tolerance: 1e-3,
                additive_factor: 0.0,
                multiplicative_factor: 1e-3,
            },
        ];
        const NUM_OBJ: usize = 6;

        // Determine padding in debugging info
        // For the width, we find the largest desired values and multiply it by the order of magnitude of its tolerance
        let max_obj_val = objectives
            .iter()
            .map(|obj| {
                (obj.desired_value.abs().ceil() as i32
                    * 10_i32.pow(obj.tolerance.abs().log10().ceil() as u32)) as i32
            })
            .max()
            .unwrap();

        let max_obj_tol = objectives
            .iter()
            .map(|obj| obj.tolerance.log10().abs().ceil() as usize)
            .max()
            .unwrap();

        let width = f64::from(max_obj_val).log10() as usize + 2 + max_obj_tol;

        // let start_instant = Instant::now();
        let max_iter = 10;

        for it in 0..=max_iter {
            dbg!(it);
            // Propagate with the estimated maneuver until the end of the maneuver
            let mut prop = prop.clone();
            prop.set_tolerance(1e-3);
            prop.dynamics = prop.dynamics.with_guidance_law(Arc::new(mnvr));
            let sc_xf_achieved = prop
                .with(sc_x0.with_guidance_mode(GuidanceMode::Thrust))
                .until_epoch(mnvr.end)?;

            println!(
                "#{} INIT: {}\nAchieved: {}\nDesired: {}",
                it, sc_x0, sc_xf_achieved, sc_xf_desired
            );

            // Build the error vector
            let mut err_vector = SVector::<f64, NUM_OBJ>::zeros();
            let mut converged = true;
            // Build debugging information
            let mut objmsg = Vec::with_capacity(objectives.len());

            // The Jacobian includes the sensitivity of each objective with respect to each variable for the whole trajectory.
            // As such, it includes the STM of that variable for the whole propagation arc.
            let mut jac = SMatrix::<f64, NUM_OBJ, NUM_VARIABLES>::zeros();

            // For each objective, we'll perturb the variables to compute the Jacobian with finite differencing.
            for (i, obj) in objectives.iter().enumerate() {
                let achieved = sc_xf_achieved.value_and_deriv(&obj.parameter)?.0;
                // Check if this objective has been achieved
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

                let mut pert_calc: Vec<_> = variables
                    .iter()
                    .enumerate()
                    .map(|(j, var)| (j, var, 0.0_f64))
                    .collect();

                pert_calc.par_iter_mut().for_each(|(_, var, jac_val)| {
                    let mut this_prop = prop.clone();
                    let mut this_mnvr = mnvr;

                    // Modify the burn itself
                    let pert = var.perturbation;
                    // Modify the maneuver, but do not change the epochs of the maneuver unless the change is greater than one millisecond
                    match var.component {
                        Vary::Duration => this_mnvr.end = mnvr.start + pert.seconds(),
                        Vary::EndEpoch => this_mnvr.end = mnvr.end + pert.seconds(),
                        Vary::StartEpoch => this_mnvr.start = mnvr.start + pert.seconds(),
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
                        _ => unreachable!(),
                    }

                    // Grab the nominal start time from the pre_dv trajectory
                    let this_sc_x0 = pre_traj.at(this_mnvr.start).unwrap();

                    this_prop.dynamics = this_prop.dynamics.with_guidance_law(Arc::new(this_mnvr));
                    let this_sc_xf_achieved = this_prop
                        .with(this_sc_x0.with_guidance_mode(GuidanceMode::Thrust))
                        .until_epoch(this_mnvr.end)
                        .unwrap();

                    let this_achieved = this_sc_xf_achieved
                        .value_and_deriv(&obj.parameter)
                        .unwrap()
                        .0;
                    *jac_val = (this_achieved - achieved) / var.perturbation;
                });

                for (j, _, jac_val) in &pert_calc {
                    jac[(i, *j)] = *jac_val;
                }
            }

            if converged {
                print!("\n\nI can't believe we converged\n\n");
                return Ok(mnvr);
                // let conv_dur = Instant::now() - start_instant;
                // let mut corrected_state = xi_start;

                // let mut state_correction = Vector6::<f64>::zeros();
                // for (i, var) in self.variables.iter().enumerate() {
                //     state_correction[var.component.vec_index()] += total_correction[i];
                // }
                // // Now, let's apply the correction to the initial state
                // if let Some(frame) = self.correction_frame {
                //     let dcm_vnc2inertial = corrected_state
                //         .orbit
                //         .dcm_from_traj_frame(frame)
                //         .unwrap()
                //         .transpose();
                //     let velocity_correction =
                //         dcm_vnc2inertial * state_correction.fixed_rows::<3>(3);
                //     corrected_state.orbit.apply_dv(velocity_correction);
                // } else {
                //     corrected_state.orbit = corrected_state.orbit + state_correction;
                // }

                // let sol = TargeterSolution {
                //     corrected_state,
                //     achieved_state: xi_start.with_orbit(xf),
                //     correction: total_correction,
                //     computation_dur: conv_dur,
                //     variables: self.variables.clone(),
                //     achieved_errors: err_vector,
                //     achieved_objectives: self.objectives.clone(),
                //     iterations: it,
                // };
                // // Log success as info
                // if it == 1 {
                //     info!("Targeter -- CONVERGED in 1 iteration");
                // } else {
                //     info!("Targeter -- CONVERGED in {} iterations", it);
                // }
                // for obj in &objmsg {
                //     info!("{}", obj);
                // }
                // return Ok(sol);
            }

            // We haven't converged yet, so let's build t
            if (err_vector.norm() - prev_err_norm).abs() < 1e-10 {
                return Err(NyxError::CorrectionIneffective(
                    "No change in objective errors".to_string(),
                ));
            }
            prev_err_norm = err_vector.norm();
            dbg!(prev_err_norm);

            debug!("Jacobian {}", jac);

            // Perform the pseudo-inverse if needed, else just inverse
            let jac_inv = pseudo_inverse!(&jac)?;

            debug!("Inverse Jacobian {}", jac_inv);

            let mut delta = jac_inv * err_vector;

            debug!("Error vector: {}\nRaw correction: {}", err_vector, delta);

            // And finally apply it to the maneuver
            let mut update_obj = false;
            for (i, var) in variables.iter().enumerate() {
                // Choose the minimum step between the provided max step and the correction.
                if delta[i].abs() > var.max_step.abs() {
                    delta[i] = var.max_step.abs() * delta[i].signum();
                } else if delta[i] > var.max_value {
                    delta[i] = var.max_value;
                } else if delta[i] < var.min_value {
                    delta[i] = var.min_value;
                }

                println!(
                    "Correction {:?} (element {}): {}",
                    var.component, i, delta[i]
                );

                let corr = delta[i];

                // Modify the maneuver, but do not change the epochs of the maneuver unless the change is greater than one millisecond
                match var.component {
                    Vary::Duration => {
                        mnvr.end = mnvr.start + corr.seconds();
                        update_obj = true;
                    }
                    Vary::EndEpoch => {
                        mnvr.end += corr.seconds();
                        update_obj = true;
                    }
                    Vary::StartEpoch => {
                        mnvr.start += corr.seconds();
                        update_obj = true;
                    }
                    Vary::MnvrAlpha | Vary::MnvrAlphaDot | Vary::MnvrAlphaDDot => {
                        mnvr.alpha_inplane_radians = mnvr
                            .alpha_inplane_radians
                            .add_val_in_order(corr % TAU, var.component.vec_index())
                            .unwrap();
                    }
                    Vary::MnvrDelta | Vary::MnvrDeltaDot | Vary::MnvrDeltaDDot => {
                        mnvr.delta_outofplane_radians = mnvr
                            .delta_outofplane_radians
                            .add_val_in_order(corr % TAU, var.component.vec_index())
                            .unwrap();
                    }
                    _ => unreachable!(),
                }
            }

            // Log progress to debug
            info!("Targeter -- Iteration #{}", it);
            for obj in &objmsg {
                println!("{}", obj);
            }

            println!("New mnvr {}", mnvr);
            if update_obj {
                sc_x0 = pre_traj.at(mnvr.start)?;
                sc_xf_desired = post_traj.at(mnvr.end)?;
                objectives = [
                    Objective::within_tolerance(StateParameter::X, sc_xf_desired.orbit.x, 1e-3),
                    Objective::within_tolerance(StateParameter::Y, sc_xf_desired.orbit.y, 1e-3),
                    Objective::within_tolerance(StateParameter::Z, sc_xf_desired.orbit.z, 1e-3),
                    Objective::within_tolerance(StateParameter::VX, sc_xf_desired.orbit.vx, 1e-3),
                    Objective::within_tolerance(StateParameter::VY, sc_xf_desired.orbit.vy, 1e-3),
                    Objective::within_tolerance(StateParameter::VZ, sc_xf_desired.orbit.vz, 1e-3),
                ];
            }
        }

        unreachable!();
    }
}
