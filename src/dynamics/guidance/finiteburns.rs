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
use crate::linalg::{DMatrix, Vector3, Vector6};
use crate::md::trajectory::InterpState;
use crate::md::ui::{Objective, Propagator, SpacecraftDynamics};
use crate::md::StateParameter;
use crate::polyfit::CommonPolynomial;
use crate::propagators::ErrorCtrl;
use crate::time::{Epoch, TimeUnit};
use crate::utils::pseudo_inverse;
use crate::NyxError;
use crate::State;
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
    pub thrust_lvl: f64,
    /// The interpolation polynomial for the in-plane angle
    pub alpha_inplane: CommonPolynomial,
    /// The interpolation polynomial for the out-of-plane angle
    pub beta_outofplane: CommonPolynomial,
}

impl fmt::Display for Mnvr {
    /// Prints the polynomial with the least significant coefficients first
    #[allow(clippy::identity_op)]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.end - self.start >= 1 * TimeUnit::Millisecond {
            write!(f, "Finite burn maneuver @ {} for {}\n\tin-plane angle α: {}\n\tout-of-plane angle β: {}", self.start, self.end-self.start, self.alpha_inplane, self.beta_outofplane)
        } else {
            write!(
                f,
                "Impulsive maneuver @ {}\n\tin-plane angle α: {}\n\tout-of-plane angle β: {}",
                self.start, self.alpha_inplane, self.beta_outofplane
            )
        }
    }
}

impl Mnvr {
    /// Creates an impulsive maneuver whose vector is the deltaV.
    /// TODO: This should use William's algorithm
    pub fn from_impulsive(dt: Epoch, vector: Vector3<f64>) -> Self {
        Self::from_time_invariant(dt, dt + TimeUnit::Millisecond, 1.0, vector)
    }

    /// Creates a manneuver from the provided time-invariant delta-v, in km/s
    pub fn from_time_invariant(
        start: Epoch,
        end: Epoch,
        thrust_lvl: f64,
        vector: Vector3<f64>,
    ) -> Self {
        // Convert to angles
        let (alpha, beta) = plane_angles_from_unit_vector(vector);
        Self {
            start,
            end,
            thrust_lvl,
            alpha_inplane: CommonPolynomial::Constant(alpha),
            beta_outofplane: CommonPolynomial::Constant(beta),
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
        let alpha_inplane = CommonPolynomial::Quadratic(alpha_ddot_tdv, 0.0, alpha_tdv);
        let beta_outofplane = CommonPolynomial::Quadratic(beta_ddot_tdv, 0.0, beta_tdv);

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
            alpha_inplane,
            beta_outofplane,
        };

        let target_state = prop
            .with(spacecraft)
            .until_epoch(epoch - 0.5 * delta_tfb * TimeUnit::Second)?;

        println!("INIT  \t{}", spacecraft);
        println!("TARGET\t{}", target_state);

        /* *** */
        /* With the finite burn maneuver, let's target the nominal state without the maneuver enabled. */
        /* The main difference with the targeter is that this is purely statically allocated */
        /* *** */

        // The correction stores, in order, alpha_0, \dot{alpha_0}, \ddot{alpha_0}, beta_0, \dot{beta_0}, \ddot{beta_0}
        let mut prev_err_norm = std::f64::INFINITY;
        let objectives = [
            Objective::new(StateParameter::X, target_state.orbit.x),
            Objective::new(StateParameter::Y, target_state.orbit.y),
            Objective::new(StateParameter::Z, target_state.orbit.z),
        ];

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
        let xi = sc_at_dv_epoch;
        let achievement_epoch = target_state.epoch();
        for it in 0..=25 {
            // Build the finite burn
            let fb_guess = FiniteBurns {
                mnvrs: vec![mnvr],
                frame: Frame::Inertial,
            };

            // Propagate for the duration of the burn
            prop.set_max_step(mnvr.end - mnvr.start);
            prop.dynamics = prop.dynamics.with_ctrl(Arc::new(fb_guess));
            let xf = prop
                .with(xi.with_guidance_mode(GuidanceMode::Custom(0)))
                .until_epoch(achievement_epoch)?
                .orbit;

            // Build the error vector
            let mut param_errors = [0.0; 3];
            let mut converged = true;

            // Build debugging information
            let mut objmsg = Vec::with_capacity(objectives.len());

            // The Jacobian includes the sensitivity of each objective with respect to each variable for the whole trajectory.
            // As such, it includes the STM of that variable for the whole propagation arc.
            let mut jac = DMatrix::from_element(3, 6, 0.0);
            // let mut jac = SMatrix::<f64, 3, 6>::zeros();

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

                let mut pert_calc: Vec<_> = (0..6).map(|i| (i, 0.0_f64)).collect();

                pert_calc.par_iter_mut().for_each(|(i, jac_val)| {
                    let this_xi = xi;
                    let perturbation = 0.001;
                    let mut this_mnvr = mnvr;
                    let mut prop = prop.clone();

                    // Change the relevant component of the polynomial which defines this maneuver
                    if *i < 3 {
                        this_mnvr.alpha_inplane = this_mnvr
                            .alpha_inplane
                            .add_val_in_order(perturbation, *i)
                            .unwrap();
                    } else {
                        this_mnvr.beta_outofplane = this_mnvr
                            .beta_outofplane
                            .add_val_in_order(perturbation, *i - 3)
                            .unwrap();
                    }

                    // Build the finite burn
                    let fb_guess = FiniteBurns {
                        mnvrs: vec![this_mnvr],
                        frame: Frame::Inertial,
                    };

                    // Propagate for the duration of the burn
                    prop.dynamics = prop.dynamics.with_ctrl(Arc::new(fb_guess));

                    let this_xf = prop
                        .clone()
                        .with(this_xi)
                        .until_epoch(achievement_epoch)
                        .unwrap()
                        .orbit;

                    let this_achieved = this_xf.value(&obj.parameter).unwrap();
                    *jac_val = (this_achieved - achieved) / perturbation;
                });

                for (j, jac_val) in &pert_calc {
                    jac[(i, *j)] = *jac_val;
                }
            }

            for obj in &objmsg {
                info!("{}", obj);
            }

            if converged {
                return Ok(mnvr);
            }

            // We haven't converged yet, so let's build the error vector
            let err_vector = Vector3::from_row_slice(&param_errors);
            if (err_vector.norm() - prev_err_norm).abs() < 1e-10 {
                return Err(NyxError::CorrectionIneffective(
                    "No change in objective errors".to_string(),
                ));
            }
            prev_err_norm = err_vector.norm();

            debug!("Jacobian {}", jac);

            // Perform the pseudo-inverse if needed, else just inverse
            let jac_inv = pseudo_inverse(&jac, NyxError::SingularJacobian)?;

            debug!("Inverse Jacobian {}", jac_inv);

            let delta = jac_inv * &err_vector;

            debug!("Error vector: {}\nRaw correction: {}", err_vector, delta);

            for (i, value) in delta.iter().enumerate() {
                // Change the relevant component of the polynomial which defines this maneuver
                if i < 3 {
                    mnvr.alpha_inplane = mnvr.alpha_inplane.add_val_in_order(*value, i)?;
                } else {
                    mnvr.beta_outofplane = mnvr.beta_outofplane.add_val_in_order(*value, i - 3)?;
                }
            }

            // Log progress to debug
            debug!("Targeter -- Iteration #{} -- {}", it, achievement_epoch);
            for obj in &objmsg {
                debug!("{}", obj);
            }
            debug!("{}", mnvr);
        }

        return Err(NyxError::MaxIterReached(format!(
            "Failed after {} iterations:\nError: {}",
            25, prev_err_norm
        )));
    }

    /// Return the thrust vector computed at the provided epoch
    pub fn vector(&self, epoch: Epoch) -> Vector3<f64> {
        let t = (epoch - self.start).in_seconds();
        let alpha = self.alpha_inplane.eval(t);
        let beta = self.beta_outofplane.eval(t);
        unit_vector_from_plane_angles(alpha, beta)
    }
}

/// A controller for a set of pre-determined maneuvers.
#[derive(Clone, Debug)]
pub struct FiniteBurns {
    /// Maneuvers should be provided in chronological order, first maneuver first in the list
    pub mnvrs: Vec<Mnvr>,
    /// The frame in which the maneuvers are defined.
    pub frame: Frame,
}

impl FiniteBurns {
    /// Builds a schedule from the vector of maneuvers, must be provided in chronological order.
    pub fn from_mnvrs(mnvrs: Vec<Mnvr>, frame: Frame) -> Arc<Self> {
        assert!(
            matches!(frame, Frame::Inertial | Frame::VNC),
            "Maneuvers must be either in the inertial frame or in a body frame"
        );
        Arc::new(Self { mnvrs, frame })
    }
}

impl GuidanceLaw for FiniteBurns {
    fn direction(&self, osc: &Spacecraft) -> Vector3<f64> {
        // NOTE: We do not increment the mnvr number here. The power function is called first,
        // so we let that function handle starting and stopping of the maneuver.
        match osc.mode {
            GuidanceMode::Custom(mnvr_no) => {
                let next_mnvr = self.mnvrs[mnvr_no as usize];
                if next_mnvr.start <= osc.epoch() {
                    if matches!(self.frame, Frame::Inertial) {
                        next_mnvr.vector(osc.epoch())
                    } else {
                        osc.orbit.dcm_from_traj_frame(self.frame).unwrap()
                            * next_mnvr.vector(osc.epoch())
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
            GuidanceMode::Custom(mnvr_no) => {
                let next_mnvr = self.mnvrs[mnvr_no as usize];
                if next_mnvr.start <= osc.epoch() {
                    next_mnvr.thrust_lvl
                } else {
                    0.0
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
        match sc.mode {
            GuidanceMode::Custom(mnvr_no) => {
                if (mnvr_no as usize) < self.mnvrs.len() {
                    let cur_mnvr = self.mnvrs[mnvr_no as usize];
                    if sc.epoch() >= cur_mnvr.end {
                        if mnvr_no as usize == self.mnvrs.len() - 1 {
                            // No following maneuver, so let's coast from now on.
                            GuidanceMode::Coast
                        } else {
                            GuidanceMode::Custom(mnvr_no + 1)
                        }
                    } else {
                        // Stay on the current maneuver
                        GuidanceMode::Custom(mnvr_no)
                    }
                } else {
                    // We're done with all the maneuvers, so we can coast now
                    GuidanceMode::Coast
                }
            }
            _ => {
                // If we haven't started the maneuvers yet, let's get ready to do so by switching to the mode
                // which will start the first maneuver
                GuidanceMode::Custom(0)
            }
        }
    }
}
