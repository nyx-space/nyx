/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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
extern crate rayon;

use super::{plane_angles_from_unit_vector, unit_vector_from_plane_angles, GuidanceLaw};
use crate::cosmic::{Frame, GuidanceMode, Spacecraft};
use crate::linalg::{DMatrix, SVector, Vector3};
use crate::md::trajectory::InterpState;
use crate::md::ui::{Objective, Propagator, SpacecraftDynamics};
use crate::md::{StateParameter, Variable, Vary};
use crate::polyfit::CommonPolynomial;
use crate::propagators::ErrorCtrl;
use crate::time::{Epoch, TimeUnit};
use crate::utils::pseudo_inverse;
use crate::NyxError;
use crate::State;
use hifitime::{Duration, TimeUnitHelper};
use rayon::prelude::*;
use std::fmt;
use std::sync::Arc;
use std::time::Instant;

/// Mnvr defined a single maneuver. Direction MUST be in the VNC frame (Velocity / Normal / Cross).
/// It may be used with a maneuver scheduler.
#[derive(Copy, Clone, Debug)]
pub struct Mnvr {
    /// Start epoch of the maneuver
    pub start: Epoch,
    /// End epoch of the maneuver
    pub end: Epoch,
    /// Thrust level, if 1.0 use all thruster available at full power
    /// TODO: Convert this to a common polynomial as well to optimize throttle, throttle rate (and accel?)
    pub thrust_lvl: f64,
    /// The interpolation polynomial for the in-plane angle
    pub alpha_inplane_degrees: CommonPolynomial,
    /// The interpolation polynomial for the out-of-plane angle
    pub beta_outofplane_degrees: CommonPolynomial,
    /// The frame in which the maneuvers are defined.
    pub frame: Frame,
}

impl fmt::Display for Mnvr {
    /// Prints the polynomial with the least significant coefficients first
    #[allow(clippy::identity_op)]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.end - self.start >= 1 * TimeUnit::Millisecond {
            write!(f, "Finite burn maneuver @ {} for {}\n\tin-plane angle α: {}\n\tout-of-plane angle β: {}", self.start, self.end-self.start, self.alpha_inplane_degrees, self.beta_outofplane_degrees)
        } else {
            write!(
                f,
                "Impulsive maneuver @ {}\n\tin-plane angle α: {}\n\tout-of-plane angle β: {}",
                self.start, self.alpha_inplane_degrees, self.beta_outofplane_degrees
            )
        }
    }
}

impl Mnvr {
    /// Creates an impulsive maneuver whose vector is the deltaV.
    /// TODO: This should use William's algorithm
    pub fn from_impulsive(dt: Epoch, vector: Vector3<f64>, frame: Frame) -> Self {
        Self::from_time_invariant(dt, dt + TimeUnit::Millisecond, 1.0, vector, frame)
    }

    /// Creates a manneuver from the provided time-invariant delta-v, in km/s
    pub fn from_time_invariant(
        start: Epoch,
        end: Epoch,
        thrust_lvl: f64,
        vector: Vector3<f64>,
        frame: Frame,
    ) -> Self {
        // Convert to angles
        let (alpha, beta) = plane_angles_from_unit_vector(vector);
        Self {
            start,
            end,
            thrust_lvl,
            alpha_inplane_degrees: CommonPolynomial::Constant(alpha),
            beta_outofplane_degrees: CommonPolynomial::Constant(beta),
            frame,
        }
    }

    /// Converts the input delta-v vector in km/s at the provided Epoch to a finite burn
    /// Uses Copernicus algorithm as described in "AAS 12-236: Recent Improvements to the Copernicus Trajectory Design and Optimization System" by Williams et al.
    /// The vector is expected to be in the same frame as the spaceraft's orbit.
    /// Convergence criteria:
    ///     1. If the change in the duration of the burn is less than 1 seconds; or
    ///     2. If the magnitude of the thrust vector matches the impulsive maneuver to less than 1 mm/s
    pub fn impulsive_to_finite<'a, E: ErrorCtrl>(
        epoch: Epoch,
        dv: Vector3<f64>,
        spacecraft: Spacecraft,
        prop: &'a Propagator<'a, SpacecraftDynamics, E>,
        frame: Frame,
    ) -> Result<Self, NyxError> {
        if spacecraft.thruster.is_none() {
            // Can't do any conversion to finite burns without a thruster
            return Err(NyxError::CtrlExistsButNoThrusterAvail);
        }

        // Grab the spacecraft at the epoch of the impulsive maneuver.
        // Clone the dynamics
        let mut prop = prop.clone();
        // Propagate to the dv epoch
        prop.dynamics = prop.dynamics.without_ctrl();
        let sc_at_dv_epoch = prop.with(spacecraft).until_epoch(epoch)?;
        let xi = sc_at_dv_epoch;

        /* *** */
        /* Compute the initial guess */
        /* *** */
        // Calculate the u, dot u (=0) and ddot u from this state
        let u = dv / dv.norm();
        let r = sc_at_dv_epoch.orbit.radius();
        let rmag = sc_at_dv_epoch.orbit.rmag();
        let u_ddot = (3.0 * sc_at_dv_epoch.orbit.frame.gm() / rmag.powi(5))
            * (r.dot(&u) * r - (r.dot(&u).powi(2) * u));
        // Compute the control rates at the time of the impulsive maneuver (tdv)
        let (alpha_tdv, beta_tdv) = plane_angles_from_unit_vector(u);
        let (alpha_ddot_tdv, beta_ddot_tdv) = plane_angles_from_unit_vector(u_ddot);
        // Build the maneuver polynomial angles from these
        let alpha_inplane_degrees = CommonPolynomial::Quadratic(alpha_ddot_tdv, 0.0, alpha_tdv);
        let beta_outofplane_degrees = CommonPolynomial::Quadratic(beta_ddot_tdv, 0.0, beta_tdv);

        // Compute a few thruster parameters
        let thruster = spacecraft.thruster.as_ref().unwrap();
        let c_km_s = thruster.exhaust_velocity() * 1e-3;

        let delta_tfb = ((c_km_s * spacecraft.mass_kg()) / thruster.thrust)
            * (1.0 - (-dv.norm() / c_km_s).exp());

        // Build the estimated maneuver
        let mut mnvr = Self {
            start: epoch - 0.5 * delta_tfb * TimeUnit::Second,
            end: epoch + 0.5 * delta_tfb * TimeUnit::Second,
            thrust_lvl: 1.0,
            alpha_inplane_degrees,
            beta_outofplane_degrees,
            frame,
        };

        println!("Estimate: {}", mnvr);

        let x0_epoch = mnvr.start;
        let xf_epoch = mnvr.end;

        // Build pre-dv and post-dv trajectories so we can update the objectives at each iteration
        // HACK: Cannot yet build trajectories with a backward prop, so we go backward and generate the traj forward only
        let sc_init = prop
            .with(spacecraft)
            .until_epoch(x0_epoch - 15.0.minutes())?;
        let (_, pre_traj) = prop
            .with(sc_init)
            .until_epoch_with_traj(spacecraft.epoch())?;
        let (_, post_traj) = prop
            .with(spacecraft)
            .until_epoch_with_traj(xf_epoch + 15.0.minutes())?;

        // println!("{}\n{}", pre_traj, post_traj);

        let mut x0_constraint = pre_traj.at(x0_epoch)?;
        let mut xf_constraint = post_traj.at(xf_epoch)?;

        // The variables in this targeter
        let variables = [
            Variable::from(Vary::MnvrAlpha).with_initial_guess(alpha_tdv),
            Variable::from(Vary::MnvrAlphaDot),
            Variable::from(Vary::MnvrAlphaDDot).with_initial_guess(alpha_ddot_tdv),
            Variable::from(Vary::MnvrBeta).with_initial_guess(beta_tdv),
            Variable::from(Vary::MnvrBetaDot),
            Variable::from(Vary::MnvrBetaDDot).with_initial_guess(beta_ddot_tdv),
            Variable::from(Vary::StartEpoch),
            Variable::from(Vary::Duration),
        ];

        println!("INIT    \t{}", spacecraft);
        println!("TARGET X0\t{}", x0_constraint);
        println!("TARGET Xf\t{}", xf_constraint);

        /* *** */
        /* With the finite burn maneuver, let's target the nominal state without the maneuver enabled. */
        /* The main difference with the targeter is that this is purely statically allocated */
        /* *** */

        // The correction stores, in order, alpha_0, \dot{alpha_0}, \ddot{alpha_0}, beta_0, \dot{beta_0}, \ddot{beta_0}
        let mut prev_err_norm = std::f64::INFINITY;
        // The objectives are organized as initial state constraints first and final state last
        let mut objectives = [
            Objective::within_tolerance(StateParameter::X, x0_constraint.orbit.x, 1e-3),
            Objective::within_tolerance(StateParameter::Y, x0_constraint.orbit.y, 1e-3),
            Objective::within_tolerance(StateParameter::Z, x0_constraint.orbit.z, 1e-3),
            Objective::within_tolerance(StateParameter::VX, x0_constraint.orbit.vx, 1e-3),
            Objective::within_tolerance(StateParameter::VY, x0_constraint.orbit.vy, 1e-3),
            Objective::within_tolerance(StateParameter::VZ, x0_constraint.orbit.vz, 1e-3),
            // xf constraints
            Objective::within_tolerance(StateParameter::X, xf_constraint.orbit.x, 1e-3),
            Objective::within_tolerance(StateParameter::Y, xf_constraint.orbit.y, 1e-3),
            Objective::within_tolerance(StateParameter::Z, xf_constraint.orbit.z, 1e-3),
            Objective::within_tolerance(StateParameter::VX, xf_constraint.orbit.vx, 1e-3),
            Objective::within_tolerance(StateParameter::VY, xf_constraint.orbit.vy, 1e-3),
            Objective::within_tolerance(StateParameter::VZ, xf_constraint.orbit.vz, 1e-3),
        ];

        const N_OBJ: usize = 12;

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

        let start_instant = Instant::now();
        let max_iter = 50;
        for it in 0..=max_iter {
            // Propagate for the duration of the burn
            prop.set_max_step(mnvr.end - mnvr.start);
            prop.dynamics = prop.dynamics.with_ctrl(Arc::new(mnvr));
            let sc_x0 = prop
                .with(xi.with_guidance_mode(GuidanceMode::Thrust))
                .until_epoch(x0_epoch)?;
            let x0 = sc_x0.orbit;

            let xf = prop
                .with(sc_x0.with_guidance_mode(GuidanceMode::Thrust))
                .until_epoch(xf_epoch)?
                .orbit;

            // Build the error vector
            let mut param_errors = [0.0; N_OBJ];
            let mut converged = true;

            // Build debugging information
            let mut objmsg = Vec::with_capacity(objectives.len());

            // The Jacobian includes the sensitivity of each objective with respect to each variable for the whole trajectory.
            // As such, it includes the STM of that variable for the whole propagation arc.
            let mut jac = DMatrix::from_element(objectives.len(), variables.len(), 0.0);

            for (i, obj) in objectives.iter().enumerate() {
                let achieved = if i < objectives.len() / 2 {
                    x0.value(&obj.parameter)?
                } else {
                    xf.value(&obj.parameter)?
                };

                let (ok, param_err) = obj.assess_raw(achieved);
                if !ok {
                    converged = false;
                }
                param_errors[i] = param_err;

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

                pert_calc.par_iter_mut().for_each(|(j, var, jac_val)| {
                    let this_xi = xi;
                    let perturbation = var.perturbation;
                    let mut this_mnvr = mnvr;
                    let mut prop = prop.clone();

                    // Perturb the maneuver
                    if *j < 3 {
                        this_mnvr.alpha_inplane_degrees = this_mnvr
                            .alpha_inplane_degrees
                            .add_val_in_order(perturbation, *j)
                            .unwrap();
                    } else if *j < 6 {
                        this_mnvr.beta_outofplane_degrees = this_mnvr
                            .beta_outofplane_degrees
                            .add_val_in_order(perturbation, *j - 3)
                            .unwrap();
                    } else if var.component == Vary::StartEpoch {
                        let prev_end = this_mnvr.start;
                        // Modification of the start epoch of the burn
                        this_mnvr.start = this_mnvr.start - perturbation.seconds();
                        // println!(
                        //     "#({}, {}) START WAS {}\t\tNOW IS: {} (pert = {})",
                        //     i,
                        //     j,
                        //     prev_end,
                        //     this_mnvr.start,
                        //     perturbation.seconds(),
                        // );
                    } else if var.component == Vary::Duration {
                        let prev_end = this_mnvr.end;
                        this_mnvr.end = this_mnvr.end + perturbation.seconds();
                        // println!(
                        //     "#({}, {}) END WAS {}\t\tNOW IS: {} (pert = {})",
                        //     i,
                        //     j,
                        //     prev_end,
                        //     this_mnvr.end,
                        //     perturbation.seconds(),
                        // );
                    }

                    // Build the dynamics for this perturbation
                    prop.dynamics = prop.dynamics.with_ctrl(Arc::new(this_mnvr));

                    // Propagate for the duration of the burn: backward if we're matching the initial state, and forward otherwise
                    let this_x = if i < objectives.len() / 2 {
                        prop.with(this_xi)
                            .until_epoch(this_mnvr.start)
                            .unwrap()
                            .orbit
                    } else {
                        let xi_prime = prop.with(this_xi).until_epoch(this_mnvr.start).unwrap();
                        prop.with(xi_prime)
                            .until_epoch(this_mnvr.end)
                            .unwrap()
                            .orbit
                    };

                    let this_achieved = this_x.value(&obj.parameter).unwrap();
                    *jac_val = (this_achieved - achieved) / perturbation;
                    // if var.component == Vary::StartEpoch || var.component == Vary::Duration {
                    //     println!("#({}, {}) => {} ==> {}", i, j, this_mnvr, jac_val);
                    // }
                });

                for (j, _, jac_val) in &pert_calc {
                    jac[(i, *j)] = *jac_val;
                }
            }

            for obj in &objmsg {
                info!("{}", obj);
            }
            if converged {
                let conv_dur = Instant::now() - start_instant;
                if it == 1 {
                    info!(
                        "Targeter -- CONVERGED in 1 iteration ({:.3} seconds)",
                        conv_dur.as_secs_f64()
                    );
                } else {
                    info!(
                        "Targeter -- CONVERGED in {} iterations ({:.3} seconds)",
                        it,
                        conv_dur.as_secs_f64()
                    );
                }
                return Ok(mnvr);
            }

            // We haven't converged yet, so let's build the error vector
            let err_vector = SVector::<f64, N_OBJ>::from_row_slice(&param_errors);
            // let err_vector = DVector::from(param_errors);
            if (err_vector.norm() - prev_err_norm).abs() < 1e-10 {
                return Err(NyxError::CorrectionIneffective(
                    "No change in objective errors".to_string(),
                ));
            }
            prev_err_norm = err_vector.norm();

            debug!("Jacobian {:.6}", jac);

            // Perform the pseudo-inverse if needed, else just inverse
            let jac_inv = pseudo_inverse(&jac, NyxError::SingularJacobian)?;

            debug!("Inverse Jacobian {:.6}", jac_inv);

            let delta = jac_inv * &err_vector;

            debug!("Error vector: {}\nRaw correction: {}", err_vector, delta);

            let mut update_objs = false;
            for (i, value) in delta.iter().enumerate() {
                // Change the relevant component of the polynomial which defines this maneuver
                if i < 3 {
                    mnvr.alpha_inplane_degrees =
                        mnvr.alpha_inplane_degrees.add_val_in_order(*value, i)?;
                } else if i < 6 {
                    mnvr.beta_outofplane_degrees = mnvr
                        .beta_outofplane_degrees
                        .add_val_in_order(*value, i - 3)?;
                } else if i == 7 && value.abs() > 1e-6 {
                    // Modification of the start epoch of the burn
                    mnvr.start = mnvr.start - value.seconds();
                    update_objs = true;
                } else if i == 8 && value.abs() > 1e-5 {
                    mnvr.end = mnvr.end + value.seconds();
                    update_objs = true;
                }
            }

            // Log progress to debug
            debug!("Targeter -- Iteration #{}", it);
            for obj in &objmsg {
                debug!("{}", obj);
            }
            info!("NEW MNVR\n{}\n", mnvr);

            // Update the objectives if we've changed the start or end time of the maneuver
            if update_objs {
                x0_constraint = pre_traj.at(mnvr.start)?;
                xf_constraint = post_traj.at(mnvr.end)?;

                objectives = [
                    Objective::within_tolerance(StateParameter::X, x0_constraint.orbit.x, 1e-3),
                    Objective::within_tolerance(StateParameter::Y, x0_constraint.orbit.y, 1e-3),
                    Objective::within_tolerance(StateParameter::Z, x0_constraint.orbit.z, 1e-3),
                    Objective::within_tolerance(StateParameter::VX, x0_constraint.orbit.vx, 1e-3),
                    Objective::within_tolerance(StateParameter::VY, x0_constraint.orbit.vy, 1e-3),
                    Objective::within_tolerance(StateParameter::VZ, x0_constraint.orbit.vz, 1e-3),
                    // xf constraints
                    Objective::within_tolerance(StateParameter::X, xf_constraint.orbit.x, 1e-3),
                    Objective::within_tolerance(StateParameter::Y, xf_constraint.orbit.y, 1e-3),
                    Objective::within_tolerance(StateParameter::Z, xf_constraint.orbit.z, 1e-3),
                    Objective::within_tolerance(StateParameter::VX, xf_constraint.orbit.vx, 1e-3),
                    Objective::within_tolerance(StateParameter::VY, xf_constraint.orbit.vy, 1e-3),
                    Objective::within_tolerance(StateParameter::VZ, xf_constraint.orbit.vz, 1e-3),
                ];
            }
        }

        return Err(NyxError::MaxIterReached(format!(
            "Failed after {} iterations:\nError: {}",
            max_iter, prev_err_norm
        )));
    }

    /// Return the thrust vector computed at the provided epoch
    pub fn vector(&self, epoch: Epoch) -> Vector3<f64> {
        let t = (epoch - self.start).in_seconds();
        let alpha = self.alpha_inplane_degrees.eval(t);
        let beta = self.beta_outofplane_degrees.eval(t);
        unit_vector_from_plane_angles(alpha.to_radians(), beta.to_radians())
    }

    pub fn duration(&self) -> Duration {
        self.end - self.start
    }
}

impl GuidanceLaw for Mnvr {
    fn direction(&self, osc: &Spacecraft) -> Vector3<f64> {
        // NOTE: We do not increment the mnvr number here. The power function is called first,
        // so we let that function handle starting and stopping of the maneuver.
        match osc.mode {
            GuidanceMode::Thrust => {
                if self.start <= osc.epoch() {
                    if matches!(self.frame, Frame::Inertial) {
                        self.vector(osc.epoch())
                    } else {
                        osc.orbit.dcm_from_traj_frame(self.frame).unwrap()
                            * self.vector(osc.epoch())
                    }
                } else {
                    Vector3::zeros()
                }
            }
            _ => Vector3::zeros(),
        }
    }

    fn throttle(&self, osc: &Spacecraft) -> f64 {
        match osc.mode {
            GuidanceMode::Thrust => {
                if osc.epoch() < self.start || osc.epoch() > self.end {
                    0.0
                } else {
                    self.thrust_lvl
                }
            }
            _ => {
                // We aren't in maneuver mode, so return 0% throttle
                0.0
            }
        }
    }

    fn next(&self, sc: &Spacecraft) -> GuidanceMode {
        // Here, we're using the Custom field of the mode to store the current maneuver number we're executing
        if sc.epoch() < self.start || sc.epoch() > self.end {
            GuidanceMode::Coast
        } else {
            GuidanceMode::Thrust
        }
    }
}
