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

use super::{ra_dec_from_unit_vector, GuidanceLaw};
use crate::cosmic::{Frame, GuidanceMode, Spacecraft};
use crate::dynamics::guidance::unit_vector_from_ra_dec;
use crate::linalg::{DMatrix, SVector, Vector3};
use crate::md::trajectory::InterpState;
use crate::md::ui::{Objective, Propagator, SpacecraftDynamics};
use crate::md::{StateParameter, Variable, Vary};
use crate::polyfit::CommonPolynomial;
use crate::propagators::ErrorCtrl;
use crate::time::{Epoch, TimeUnit};
use crate::State;
use crate::{pseudo_inverse, NyxError};
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
    pub alpha_inplane_radians: CommonPolynomial,
    /// The interpolation polynomial for the out-of-plane angle
    pub delta_outofplane_radians: CommonPolynomial,
    /// The frame in which the maneuvers are defined.
    pub frame: Frame,
}

impl fmt::Display for Mnvr {
    /// Prints the polynomial with the least significant coefficients first
    #[allow(clippy::identity_op)]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.end - self.start >= 1 * TimeUnit::Millisecond {
            let start_vec = self.vector(self.start);
            let end_vec = self.vector(self.end);
            write!(
                f,
                "Finite burn maneuver @ {:.2}% on {} for {} (ending on {})",
                100.0 * self.thrust_lvl,
                self.start,
                self.end - self.start,
                self.end,
            )?;
            write!(
                f,
                "\n\tin-plane angle α: {}\n\tout-of-plane angle β: {}",
                self.alpha_inplane_radians, self.delta_outofplane_radians
            )?;
            write!(
                f,
                "\n\tinitial dir: [{:.6}, {:.6}, {:.6}]\n\tfinal dir  : [{:.6}, {:.6}, {:.6}]",
                start_vec[0], start_vec[1], start_vec[2], end_vec[0], end_vec[1], end_vec[2]
            )
        } else {
            write!(
                f,
                "Impulsive maneuver @ {}\n\tin-plane angle α: {}\n\tout-of-plane angle β: {}",
                self.start, self.alpha_inplane_radians, self.delta_outofplane_radians
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
        let (alpha, delta) = ra_dec_from_unit_vector(vector);
        Self {
            start,
            end,
            thrust_lvl,
            alpha_inplane_radians: CommonPolynomial::Constant(alpha),
            delta_outofplane_radians: CommonPolynomial::Constant(delta),
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
        let mut prop_no_ctrl = prop.clone();
        // Propagate to the dv epoch
        prop_no_ctrl.dynamics = prop.dynamics.without_ctrl();
        let sc_at_dv_epoch = prop_no_ctrl.with(spacecraft).until_epoch(epoch)?;

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
        let (alpha_tdv, delta_tdv) = ra_dec_from_unit_vector(u);
        let (alpha_ddot_tdv, delta_ddot_tdv) = ra_dec_from_unit_vector(u_ddot);
        // Build the maneuver polynomial angles from these
        let alpha_inplane_degrees = CommonPolynomial::Quadratic(alpha_ddot_tdv, 0.0, alpha_tdv);
        let delta_outofplane_degrees = CommonPolynomial::Quadratic(delta_ddot_tdv, 0.0, delta_tdv);

        // Compute a few thruster parameters
        let thruster = spacecraft.thruster.as_ref().unwrap();
        let c_km_s = thruster.exhaust_velocity();

        let delta_tfb = ((c_km_s * spacecraft.mass_kg()) / thruster.thrust)
            * (1.0 - (-dv.norm() * 1e3 / c_km_s).exp());

        // Build the estimated maneuver
        let mut mnvr = Self {
            start: epoch - 0.5 * delta_tfb * TimeUnit::Second,
            end: epoch + 0.5 * delta_tfb * TimeUnit::Second,
            thrust_lvl: 1.0,
            alpha_inplane_radians: alpha_inplane_degrees,
            delta_outofplane_radians: delta_outofplane_degrees,
            frame,
        };

        println!("Estimate: {}", mnvr);

        let spacecraft_with_dv = spacecraft.with_dv(dv);
        let mut sc_x0 = prop_no_ctrl.with(spacecraft).until_epoch(mnvr.start)?;
        let mut xf_constraint = prop_no_ctrl
            .with(spacecraft_with_dv)
            .until_epoch(mnvr.end)?;

        // The variables in this targeter
        let variables = [
            Variable::from(Vary::MnvrAlpha).with_initial_guess(alpha_tdv),
            Variable::from(Vary::MnvrAlphaDot),
            Variable::from(Vary::MnvrAlphaDDot).with_initial_guess(alpha_ddot_tdv),
            Variable::from(Vary::MnvrDelta).with_initial_guess(delta_tdv),
            Variable::from(Vary::MnvrDeltaDot),
            Variable::from(Vary::MnvrDeltaDDot).with_initial_guess(delta_ddot_tdv),
            Variable::from(Vary::StartEpoch),
            Variable::from(Vary::Duration),
        ];

        println!("INIT    \t{}", sc_x0);
        println!("TARGET Xf\t{}", xf_constraint);

        /* *** */
        /* With the finite burn maneuver, let's target the nominal state without the maneuver enabled. */
        /* The main difference with the targeter is that this is purely statically allocated */
        /* *** */

        // The correction stores, in order, alpha_0, \dot{alpha_0}, \ddot{alpha_0}, beta_0, \dot{beta_0}, \ddot{beta_0}
        let mut prev_err_norm = std::f64::INFINITY;
        // The objectives are organized as initial state constraints first and final state last
        let mut objectives = [
            Objective::within_tolerance(StateParameter::X, xf_constraint.orbit.x, 1e-3),
            Objective::within_tolerance(StateParameter::Y, xf_constraint.orbit.y, 1e-3),
            Objective::within_tolerance(StateParameter::Z, xf_constraint.orbit.z, 1e-3),
            // Objective::within_tolerance(StateParameter::VX, xf_constraint.orbit.vx, 1e-3),
            // Objective::within_tolerance(StateParameter::VY, xf_constraint.orbit.vy, 1e-3),
            // Objective::within_tolerance(StateParameter::VZ, xf_constraint.orbit.vz, 1e-3),
        ];

        const N_OBJ: usize = 3;

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
        let max_iter = 5;
        for it in 0..=max_iter {
            // Propagate for the duration of the burn
            prop.set_max_step(mnvr.end - mnvr.start);
            prop.dynamics = prop.dynamics.with_ctrl(Arc::new(mnvr));

            let xf = prop
                .with(sc_x0.with_guidance_mode(GuidanceMode::Thrust))
                .until_epoch(mnvr.end)?
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
                let achieved = xf.value(&obj.parameter)?;

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
                    let perturbation = var.perturbation;
                    let mut this_mnvr = mnvr;
                    let mut prop = prop.clone();

                    // Perturb the maneuver
                    let var = variables[*j].component;
                    // Change the relevant component of the polynomial which defines this maneuver
                    match var {
                        Vary::MnvrAlpha => {
                            this_mnvr.alpha_inplane_radians = mnvr
                                .alpha_inplane_radians
                                .add_val_in_order(perturbation, 0)
                                .unwrap();
                        }
                        Vary::MnvrAlphaDot => {
                            this_mnvr.alpha_inplane_radians = mnvr
                                .alpha_inplane_radians
                                .add_val_in_order(perturbation, 1)
                                .unwrap();
                        }
                        Vary::MnvrAlphaDDot => {
                            this_mnvr.alpha_inplane_radians = mnvr
                                .alpha_inplane_radians
                                .add_val_in_order(perturbation, 2)
                                .unwrap();
                        }

                        Vary::MnvrDelta => {
                            this_mnvr.delta_outofplane_radians = mnvr
                                .delta_outofplane_radians
                                .add_val_in_order(perturbation, 0)
                                .unwrap();
                        }
                        Vary::MnvrDeltaDot => {
                            this_mnvr.delta_outofplane_radians = mnvr
                                .delta_outofplane_radians
                                .add_val_in_order(perturbation, 1)
                                .unwrap();
                        }
                        Vary::MnvrDeltaDDot => {
                            this_mnvr.delta_outofplane_radians = mnvr
                                .delta_outofplane_radians
                                .add_val_in_order(perturbation, 2)
                                .unwrap();
                        }

                        Vary::StartEpoch => {
                            this_mnvr.start = mnvr.start + perturbation.seconds();
                        }
                        Vary::EndEpoch => {
                            this_mnvr.end = mnvr.end + perturbation.seconds();
                        }
                        Vary::Duration => {
                            this_mnvr.end = mnvr.start + perturbation.seconds();
                        }
                        _ => unreachable!(),
                    };

                    // Build the dynamics for this perturbation
                    prop.dynamics = prop.dynamics.with_ctrl(Arc::new(this_mnvr));
                    prop.set_max_step(mnvr.end - mnvr.start);

                    // Propagate for the duration of the burn
                    // Start with the spacecraft at time mvnr.start and propagate until the end of the maneuver
                    let sc_at_mnvr_start = prop_no_ctrl
                        .with(spacecraft)
                        .until_epoch(this_mnvr.start)
                        .unwrap();

                    let this_sc = prop
                        .with(sc_at_mnvr_start)
                        .until_epoch(this_mnvr.end)
                        .unwrap();
                    let this_x = this_sc.orbit;

                    let this_achieved = this_x.value(&obj.parameter).unwrap();

                    println!(
                        "{:?}\tΔ={}/{}\t{} kg used in {}\tWhere: {}",
                        var,
                        this_achieved - achieved,
                        perturbation,
                        sc_at_mnvr_start.mass_kg() - this_sc.mass_kg(),
                        this_x.epoch() - sc_at_mnvr_start.epoch(),
                        this_mnvr
                    );

                    *jac_val = (this_achieved - achieved) / perturbation;
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
            let jac_inv = pseudo_inverse!(&jac)?;

            debug!("Inverse Jacobian {:.6}", jac_inv);

            let delta = jac_inv * &err_vector;

            debug!("Error vector: {}\nRaw correction: {}", err_vector, delta);

            let mut update_objs = false;
            for (i, value) in delta.iter().enumerate() {
                let var = variables[i].component;
                // Change the relevant component of the polynomial which defines this maneuver
                match var {
                    Vary::MnvrAlpha => {
                        mnvr.alpha_inplane_radians =
                            mnvr.alpha_inplane_radians.add_val_in_order(*value, 0)?;
                    }
                    Vary::MnvrAlphaDot => {
                        mnvr.alpha_inplane_radians =
                            mnvr.alpha_inplane_radians.add_val_in_order(*value, 1)?;
                    }
                    Vary::MnvrAlphaDDot => {
                        mnvr.alpha_inplane_radians =
                            mnvr.alpha_inplane_radians.add_val_in_order(*value, 2)?;
                    }

                    Vary::MnvrDelta => {
                        mnvr.delta_outofplane_radians =
                            mnvr.delta_outofplane_radians.add_val_in_order(*value, 0)?;
                    }
                    Vary::MnvrDeltaDot => {
                        mnvr.delta_outofplane_radians =
                            mnvr.delta_outofplane_radians.add_val_in_order(*value, 1)?;
                    }
                    Vary::MnvrDeltaDDot => {
                        mnvr.delta_outofplane_radians =
                            mnvr.delta_outofplane_radians.add_val_in_order(*value, 2)?;
                    }

                    Vary::StartEpoch => {
                        mnvr.start = mnvr.start + value.seconds();
                        update_objs = true;
                    }
                    Vary::EndEpoch => {
                        mnvr.end = mnvr.end + value.seconds();
                        update_objs = true;
                    }
                    Vary::Duration => {
                        mnvr.end = mnvr.start + value.seconds();
                        update_objs = true;
                    }
                    _ => unreachable!(),
                };
            }

            // Log progress to debug
            debug!("Targeter -- Iteration #{}", it);
            for obj in &objmsg {
                debug!("{}", obj);
            }
            info!("NEW MNVR\n{}\n", mnvr);

            // Update the objectives if we've changed the start or end time of the maneuver
            if update_objs {
                sc_x0 = prop_no_ctrl
                    .with(spacecraft)
                    .until_epoch(mnvr.start)
                    .unwrap();

                xf_constraint = prop_no_ctrl
                    .with(spacecraft_with_dv)
                    .until_epoch(mnvr.end)
                    .unwrap();

                println!("INIT   X0\t{}\nTARGET Xf\t{}\n\n", sc_x0, xf_constraint);

                objectives = [
                    Objective::within_tolerance(StateParameter::X, xf_constraint.orbit.x, 1e-3),
                    Objective::within_tolerance(StateParameter::Y, xf_constraint.orbit.y, 1e-3),
                    Objective::within_tolerance(StateParameter::Z, xf_constraint.orbit.z, 1e-3),
                    // Objective::within_tolerance(StateParameter::VX, xf_constraint.orbit.vx, 1e-3),
                    // Objective::within_tolerance(StateParameter::VY, xf_constraint.orbit.vy, 1e-3),
                    // Objective::within_tolerance(StateParameter::VZ, xf_constraint.orbit.vz, 1e-3),
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
        let alpha = self.alpha_inplane_radians.eval(t);
        let delta = self.delta_outofplane_radians.eval(t);
        unit_vector_from_ra_dec(alpha, delta)
    }

    /// Return the duration of this maneuver
    pub fn duration(&self) -> Duration {
        self.end - self.start
    }

    /// Returns the direction of the burn at the start of the burn, useful for setting new angles
    pub fn direction(&self) -> Vector3<f64> {
        let alpha = self.alpha_inplane_radians.coeff_in_order(0).unwrap();
        let delta = self.delta_outofplane_radians.coeff_in_order(0).unwrap();
        unit_vector_from_ra_dec(alpha, delta)
    }

    /// Set the time-invariant direction for this finite burn while keeping the other components as they are
    pub fn set_direction(&mut self, vector: Vector3<f64>) {
        self.set_direction_and_rates(vector, self.rate(), self.accel());
    }

    /// Returns the rate of direction of the burn at the start of the burn, useful for setting new angles
    pub fn rate(&self) -> Vector3<f64> {
        match self.alpha_inplane_radians.coeff_in_order(1) {
            Ok(alpha) => {
                let delta = self.delta_outofplane_radians.coeff_in_order(1).unwrap();
                unit_vector_from_ra_dec(alpha, delta)
            }
            Err(_) => Vector3::zeros(),
        }
    }

    /// Set the rate of direction for this finite burn while keeping the other components as they are
    pub fn set_rate(&mut self, rate: Vector3<f64>) {
        self.set_direction_and_rates(self.direction(), rate, self.accel());
    }

    /// Returns the acceleration of the burn at the start of the burn, useful for setting new angles
    pub fn accel(&self) -> Vector3<f64> {
        match self.alpha_inplane_radians.coeff_in_order(2) {
            Ok(alpha) => {
                let delta = self.delta_outofplane_radians.coeff_in_order(2).unwrap();
                unit_vector_from_ra_dec(alpha, delta)
            }
            Err(_) => Vector3::zeros(),
        }
    }

    /// Set the acceleration of the direction of this finite burn while keeping the other components as they are
    pub fn set_accel(&mut self, accel: Vector3<f64>) {
        self.set_direction_and_rates(self.direction(), self.rate(), accel);
    }

    /// Set the initial direction, direction rate, and direction acceleration for this finite burn
    pub fn set_direction_and_rates(
        &mut self,
        dir: Vector3<f64>,
        rate: Vector3<f64>,
        accel: Vector3<f64>,
    ) {
        let (alpha, delta) = ra_dec_from_unit_vector(dir);
        if rate.norm() < 2e-16 && accel.norm() < 2e-16 {
            self.alpha_inplane_radians = CommonPolynomial::Constant(alpha);
            self.delta_outofplane_radians = CommonPolynomial::Constant(delta);
        } else {
            let (alpha_dt, delta_dt) = ra_dec_from_unit_vector(rate);
            if accel.norm() < 2e-16 {
                self.alpha_inplane_radians = CommonPolynomial::Linear(alpha_dt, alpha);
                self.delta_outofplane_radians = CommonPolynomial::Linear(delta_dt, delta);
            } else {
                let (alpha_ddt, delta_ddt) = ra_dec_from_unit_vector(accel);
                self.alpha_inplane_radians =
                    CommonPolynomial::Quadratic(alpha_ddt, alpha_dt, alpha);
                self.delta_outofplane_radians =
                    CommonPolynomial::Quadratic(delta_ddt, delta_dt, delta);
            }
        }
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
