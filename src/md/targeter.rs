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
use super::rayon::prelude::*;
use super::StateParameter;
pub use super::{Variable, Vary};
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DMatrix, DVector, DefaultAllocator};
use crate::md::ui::*;
use crate::propagators::error_ctrl::ErrorCtrl;
use std::convert::TryInto;
use std::fmt;
use std::time::{Duration, Instant};

/// Defines a state parameter event finder
#[derive(Copy, Clone, Debug)]
pub struct Objective {
    /// The state parameter to target
    pub parameter: StateParameter,
    /// The desired self.desired_value, must be in the same units as the state parameter
    pub desired_value: f64,
    /// The precision on the desired value
    pub tolerance: f64,
    /// A multiplicative factor this parameter's error in the targeting (defaults to 1.0)
    pub multiplicative_factor: f64,
    /// An additive factor to this parameters's error in the targeting (defaults to 0.0)
    pub additive_factor: f64,
}

impl Objective {
    /// Match a specific value for the parameter.
    /// By default, the tolerance on the parameter is 0.1 times whatever unit is the default for that parameter.
    /// For example, a radius event will seek the requested value at the decimeter level, and an angle event will seek it at the tenth of a degree.
    pub fn new(parameter: StateParameter, desired_value: f64) -> Self {
        Self::within_tolerance(
            parameter,
            desired_value,
            parameter.default_event_precision(),
        )
    }

    /// Match a specific value for the parameter to hit the specified value with the provided tolerance on the value
    pub fn within_tolerance(parameter: StateParameter, desired_value: f64, tolerance: f64) -> Self {
        Self {
            parameter,
            desired_value,
            tolerance,
            multiplicative_factor: 1.0,
            additive_factor: 0.0,
        }
    }
}

/// Defines a targeter solution
#[derive(Clone, Debug)]
pub struct TargeterSolution {
    /// The corrected spacecraft state at the correction epoch
    pub state: Spacecraft,
    /// The correction vector applied
    pub correction: DVector<f64>,
    /// The kind of correction (position or velocity)
    pub variables: Vec<Variable>,
    /// The epoch at which the objectives are achieved
    pub achievement_epoch: Epoch,
    /// The errors achieved
    pub achieved_errors: Vec<f64>,
    /// The objectives set in the targeter
    pub achieved_objectives: Vec<Objective>,
    /// The number of iterations required
    pub iterations: usize,
    /// Computation duration
    pub computation_dur: Duration,
}

impl fmt::Display for TargeterSolution {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut objmsg = String::from("");
        for (i, obj) in self.achieved_objectives.iter().enumerate() {
            objmsg.push_str(&format!(
                "\n\t\t{:?} = {:.3} (wanted {:.3} +/- {:.1e})",
                obj.parameter,
                obj.desired_value + self.achieved_errors[i],
                obj.desired_value,
                obj.tolerance
            ));
        }

        let mut corrmsg = String::from("Correction:");
        let mut is_only_position = true;
        let mut is_only_velocity = true;
        for (i, var) in self.variables.iter().enumerate() {
            let unit = match var.component {
                Vary::PositionX | Vary::PositionY | Vary::PositionZ => {
                    is_only_velocity = false;
                    "m"
                }
                Vary::VelocityX
                | Vary::VelocityY
                | Vary::VelocityZ
                | Vary::VelocityV
                | Vary::VelocityN
                | Vary::VelocityC => {
                    is_only_position = false;
                    "m/s"
                }
            };
            corrmsg.push_str(&format!(
                "\n\t\t{:?} = {:.3} {}",
                var.component, self.correction[i], unit
            ));
        }

        if is_only_position {
            corrmsg.push_str(&format!(
                "\n\t\t|Δr| = {:.3} m",
                self.correction.norm() * 1e3
            ));
        } else if is_only_velocity {
            corrmsg.push_str(&format!(
                "\n\t\t|Δv| = {:.3} m/s",
                self.correction.norm() * 1e3
            ));
        }

        write!(
            f,
            "Targeter solution correcting {:?} (converged in {:.3} seconds, {} iterations):\n\t{}\n\tAchieved:{}\n\tFinal state:\n\t\t{}\n\t\t{:x}",
            self.variables.iter().map(|v| format!("{:?}", v.component)).collect::<Vec<String>>(),
            self.computation_dur.as_secs_f64(), self.iterations, corrmsg, objmsg, self.state, self.state
        )
    }
}

/// The target is a differential corrector.
#[derive(Clone, Debug)]
pub struct Targeter<'a, D: Dynamics<StateType = Spacecraft>, E: ErrorCtrl>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>,
{
    /// The propagator setup (kind, stages, etc.)
    pub prop: Arc<&'a Propagator<'a, D, E>>,
    /// The list of objectives of this targeter
    pub objectives: Vec<Objective>,
    /// An optional frame (and Cosm) to compute the objectives in.
    /// Needed if the propagation frame is separate from objectives frame (e.g. for B Plane targeting).
    pub objective_frame: Option<(Frame, Arc<Cosm>)>,
    /// The kind of correction to apply to achieve the objectives
    pub variables: Vec<Variable>,
    /// Maximum number of iterations
    pub iterations: usize,
}

impl<'a, D: Dynamics<StateType = Spacecraft>, E: ErrorCtrl> fmt::Display for Targeter<'a, D, E>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>
        + Allocator<f64, <D::StateType as State>::Size>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut objmsg = String::from("");
        for obj in &self.objectives {
            objmsg.push_str(&format!(
                "{:?} = {:.3} (+/- {:.1e}) ",
                obj.parameter, obj.desired_value, obj.tolerance
            ));
        }

        write!(
            f,
            "Targeter:\n\tObjectives: {}\n\tCorrect: {:?}",
            objmsg, self.variables
        )
    }
}

impl<'a, D: Dynamics<StateType = Spacecraft>, E: ErrorCtrl> Targeter<'a, D, E>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>
        + Allocator<f64, <D::StateType as State>::Size>,
{
    /// Create a new Targeter which will apply an impulsive delta-v correction.
    pub fn new(
        prop: Arc<&'a Propagator<'a, D, E>>,
        variables: Vec<Variable>,
        objectives: Vec<Objective>,
    ) -> Self {
        Self {
            prop,
            objectives,
            variables,
            iterations: 100,
            objective_frame: None,
        }
    }

    /// Create a new Targeter which will apply an impulsive delta-v correction.
    pub fn in_frame(
        prop: Arc<&'a Propagator<'a, D, E>>,
        variables: Vec<Variable>,
        objectives: Vec<Objective>,
        objective_frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Self {
        Self {
            prop,
            objectives,
            variables,
            iterations: 100,
            objective_frame: Some((objective_frame, cosm)),
        }
    }

    /// Create a new Targeter which will apply an impulsive delta-v correction. By default, max step is 0.5 km/s.
    pub fn delta_v_in_frame(
        prop: Arc<&'a Propagator<'a, D, E>>,
        objectives: Vec<Objective>,
        objective_frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Self {
        Self {
            prop,
            objectives,
            variables: vec![
                Vary::VelocityX.try_into().unwrap(),
                Vary::VelocityY.try_into().unwrap(),
                Vary::VelocityZ.try_into().unwrap(),
            ],
            iterations: 100,
            objective_frame: Some((objective_frame, cosm)),
        }
    }

    /// Create a new Targeter which will apply an impulsive delta-v correction.
    pub fn delta_v(prop: Arc<&'a Propagator<'a, D, E>>, objectives: Vec<Objective>) -> Self {
        Self {
            prop,
            objectives,
            variables: vec![
                Vary::VelocityX.try_into().unwrap(),
                Vary::VelocityY.try_into().unwrap(),
                Vary::VelocityZ.try_into().unwrap(),
            ],
            iterations: 100,
            objective_frame: None,
        }
    }

    /// Create a new Targeter which will MOVE the position of the spacecraft at the correction epoch
    pub fn delta_r_in_frame(
        prop: Arc<&'a Propagator<'a, D, E>>,
        objectives: Vec<Objective>,
        objective_frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Self {
        Self {
            prop,
            objectives,
            variables: vec![
                Vary::PositionX.try_into().unwrap(),
                Vary::PositionY.try_into().unwrap(),
                Vary::PositionZ.try_into().unwrap(),
            ],
            iterations: 100,
            objective_frame: Some((objective_frame, cosm)),
        }
    }

    /// Create a new Targeter which will MOVE the position of the spacecraft at the correction epoch
    pub fn delta_r(prop: Arc<&'a Propagator<'a, D, E>>, objectives: Vec<Objective>) -> Self {
        Self {
            prop,
            objectives,
            variables: vec![
                Vary::PositionX.try_into().unwrap(),
                Vary::PositionY.try_into().unwrap(),
                Vary::PositionZ.try_into().unwrap(),
            ],
            iterations: 100,
            objective_frame: None,
        }
    }

    /// Runs the targeter using finite differencing (for now).
    pub fn try_achieve_from(
        &self,
        initial_state: Spacecraft,
        correction_epoch: Epoch,
        achievement_epoch: Epoch,
    ) -> Result<TargeterSolution, NyxError> {
        self.try_achieve_fd(initial_state, correction_epoch, achievement_epoch)
    }

    /// Differential correction using finite differencing
    #[allow(clippy::comparison_chain)]
    pub fn try_achieve_fd(
        &self,
        initial_state: Spacecraft,
        correction_epoch: Epoch,
        achievement_epoch: Epoch,
    ) -> Result<TargeterSolution, NyxError> {
        if self.objectives.is_empty() {
            return Err(NyxError::UnderdeterminedProblem);
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

        // Store the total correction in Vector3
        let mut total_correction = DVector::from_element(self.variables.len(), 0.0);

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
                _ => unimplemented!(),
            }
            total_correction[i] += var.init_guess;
            // Check the validity (this function will report to log)
            var.valid();
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
            // Modify each variable by the perturbation, propagate, compute the final parameter, and store how modifying that variable affects the final parameter
            let cur_xi = xi;
            let xf = self.prop.with(cur_xi).until_epoch(achievement_epoch)?;

            let xf_dual_obj_frame = match &self.objective_frame {
                Some((frame, cosm)) => {
                    let orbit_obj_frame = cosm.frame_chg(&xf.orbit, *frame);
                    OrbitDual::from(orbit_obj_frame)
                }
                None => OrbitDual::from(xf.orbit),
            };

            // Build the error vector
            let mut param_errors = Vec::new();
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
            let mut jac = DMatrix::from_element(self.objectives.len(), self.variables.len(), 0.0);

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

                let param_err = obj.multiplicative_factor * (obj.desired_value - achieved)
                    + obj.additive_factor;

                if param_err.abs() > obj.tolerance {
                    converged = false;
                }
                param_errors.push(param_err);

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

                    match var.component {
                        Vary::PositionX => {
                            this_xi.orbit.x += var.perturbation;
                        }
                        Vary::PositionY => {
                            this_xi.orbit.y += var.perturbation;
                        }
                        Vary::PositionZ => {
                            this_xi.orbit.z += var.perturbation;
                        }
                        Vary::VelocityX | Vary::VelocityV => {
                            this_xi.orbit.vx += var.perturbation;
                        }
                        Vary::VelocityY | Vary::VelocityN => {
                            this_xi.orbit.vy += var.perturbation;
                        }
                        Vary::VelocityZ | Vary::VelocityC => {
                            this_xi.orbit.vz += var.perturbation;
                        }
                    };

                    let this_xf = self
                        .prop
                        .clone()
                        .with(this_xi)
                        .until_epoch(achievement_epoch)
                        .unwrap();

                    let xf_dual_obj_frame = match &self.objective_frame {
                        Some((frame, cosm)) => {
                            let orbit_obj_frame = cosm.frame_chg(&this_xf.orbit, *frame);
                            OrbitDual::from(orbit_obj_frame)
                        }
                        None => OrbitDual::from(this_xf.orbit),
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
                });

                for (j, _, jac_val) in &pert_calc {
                    jac[(i, *j)] = *jac_val;
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
                        _ => unimplemented!(),
                    }
                }

                let sol = TargeterSolution {
                    state,
                    correction: total_correction,
                    computation_dur: conv_dur,
                    variables: self.variables.clone(),
                    achievement_epoch,
                    achieved_errors: param_errors,
                    achieved_objectives: self.objectives.clone(),
                    iterations: it,
                };
                info!("Targeter -- CONVERGED in {} iterations", it);
                for obj in &objmsg {
                    info!("{}", obj);
                }
                return Ok(sol);
            }

            // We haven't converged yet, so let's build the error vector
            let err_vector = DVector::from(param_errors);
            if (err_vector.norm() - prev_err_norm).abs() < 1e-10 {
                return Err(NyxError::CorrectionIneffective(
                    "No change in objective errors".to_string(),
                ));
            }
            prev_err_norm = err_vector.norm();

            debug!("Jacobian {}", jac);

            // Perform the pseudo-inverse if needed, else just inverse
            let jac_inv = if self.variables.len() == self.objectives.len() {
                match jac.try_inverse() {
                    Some(inv) => inv,
                    None => return Err(NyxError::SingularStateTransitionMatrix),
                }
            } else if self.objectives.len() < self.variables.len() {
                let m1_inv = match (&jac * &jac.transpose()).try_inverse() {
                    Some(inv) => inv,
                    None => return Err(NyxError::SingularStateTransitionMatrix),
                };
                &jac.transpose() * m1_inv
            } else {
                let m2_inv = match (&jac.transpose() * &jac).try_inverse() {
                    Some(inv) => inv,
                    None => return Err(NyxError::SingularStateTransitionMatrix),
                };
                m2_inv * &jac.transpose()
            };

            debug!("Inverse Jacobian {}", jac_inv);

            let mut delta = jac_inv * &err_vector;

            debug!("Error vector: {}\nRaw correction: {}", err_vector, delta);

            // And finally apply it to the xi
            for (i, var) in self.variables.iter().enumerate() {
                // Choose the minimum step between the provided max step and the correction.
                if delta[i].abs() > var.max_step.abs() {
                    delta[i] = var.max_step.abs() * delta[i].signum();
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
                    _ => unimplemented!(),
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

    /// Differential correction using hyperdual numbers for the objectives
    #[allow(clippy::comparison_chain)]
    pub fn try_achieve_from_dual(
        &self,
        initial_state: Spacecraft,
        correction_epoch: Epoch,
        achievement_epoch: Epoch,
    ) -> Result<TargeterSolution, NyxError> {
        warn!("DO NOT USE! The dual number target is broken! https://gitlab.com/nyx-space/nyx/-/issues/207");

        if self.objectives.is_empty() {
            return Err(NyxError::UnderdeterminedProblem);
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

        // Store the total correction in Vector3
        let mut total_correction = DVector::from_element(self.variables.len(), 0.0);

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
                _ => unimplemented!(),
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
            xi.enable_traj_stm();

            let xf = self.prop.with(xi).until_epoch(achievement_epoch)?;
            // Diagonal of the STM of the trajectory (will include the frame change if needed)
            debug!("STM {}", xf.stm());

            let xf_dual_obj_frame = match &self.objective_frame {
                Some((frame, cosm)) => {
                    let orbit_obj_frame = cosm.frame_chg(&xf.orbit, *frame);
                    OrbitDual::from(orbit_obj_frame)
                }
                None => OrbitDual::from(xf.orbit),
            };

            // Build the error vector
            let mut param_errors = Vec::new();
            // let mut jac_rows = Vec::new();
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

                let param_err = obj.multiplicative_factor * (obj.desired_value - achieved)
                    + obj.additive_factor;

                if param_err.abs() > obj.tolerance {
                    converged = false;
                }
                param_errors.push(param_err);

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
                let mut partial_vec = DVector::from_element(6, 0.0);
                for (i, val) in [
                    partial.wtr_x(),
                    partial.wtr_y(),
                    partial.wtr_z(),
                    partial.wtr_vx(),
                    partial.wtr_vy(),
                    partial.wtr_vz(),
                ]
                .iter()
                .enumerate()
                {
                    partial_vec[i] = *val;
                }
                // Multiply its transpose by the STM and extract the data
                let obj_jac = -xf.stm() * partial_vec;
                println!("{}", obj_jac);
                for (j, var) in self.variables.iter().enumerate() {
                    let idx = match var.component {
                        Vary::PositionX => 0,
                        Vary::PositionY => 1,
                        Vary::PositionZ => 2,
                        Vary::VelocityX | Vary::VelocityV => 3,
                        Vary::VelocityY | Vary::VelocityN => 4,
                        Vary::VelocityZ | Vary::VelocityC => 5,
                    };
                    // We only ever need the diagonals of the STM I think?
                    jac[(i, j)] = obj_jac[idx];
                }
                // jac_rows.push();
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
                        _ => unimplemented!(),
                    }
                }

                let sol = TargeterSolution {
                    state,
                    correction: total_correction,
                    computation_dur: conv_dur,
                    variables: self.variables.clone(),
                    achievement_epoch,
                    achieved_errors: param_errors,
                    achieved_objectives: self.objectives.clone(),
                    iterations: it,
                };
                info!("Targeter -- CONVERGED in {} iterations", it);
                for obj in &objmsg {
                    info!("{}", obj);
                }
                return Ok(sol);
            }

            // We haven't converged yet, so let's build the error vector
            let err_vector = DVector::from(param_errors);
            if (err_vector.norm() - prev_err_norm).abs() < 1e-10 {
                return Err(NyxError::CorrectionIneffective(
                    "No change in objective errors".to_string(),
                ));
            }
            prev_err_norm = err_vector.norm();

            debug!("Jacobian {}", jac);

            // Perform the pseudo-inverse if needed, else just inverse
            let jac_inv = if self.variables.len() == self.objectives.len() {
                match jac.try_inverse() {
                    Some(inv) => inv,
                    None => return Err(NyxError::SingularStateTransitionMatrix),
                }
            } else if self.objectives.len() < self.variables.len() {
                let m1_inv = match (&jac * &jac.transpose()).try_inverse() {
                    Some(inv) => inv,
                    None => return Err(NyxError::SingularStateTransitionMatrix),
                };
                &jac.transpose() * m1_inv
            } else {
                let m2_inv = match (&jac.transpose() * &jac).try_inverse() {
                    Some(inv) => inv,
                    None => return Err(NyxError::SingularStateTransitionMatrix),
                };
                m2_inv * &jac.transpose()
            };

            debug!("Inverse Jacobian {:e}", jac_inv);

            let mut delta = jac_inv * err_vector;

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
                    _ => unimplemented!(),
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

    /// Apply a correction and propagate to achievement epoch. Also checks that the objectives are indeed matched
    pub fn apply(&self, solution: TargeterSolution) -> Result<Spacecraft, NyxError> {
        let (xf, _) = self.apply_with_traj(solution)?;
        Ok(xf)
    }

    /// Apply a correction and propagate to achievement epoch, return the final state and trajectory.
    /// Also checks that the objectives are indeed matched.
    /// WARNING: This checks that the final objectives are matched with TEN TIMES the initial tolerances
    /// XXX Check why that is the case.
    pub fn apply_with_traj(
        &self,
        solution: TargeterSolution,
    ) -> Result<(Spacecraft, ScTraj), NyxError> {
        // Propagate until achievement epoch
        let (xf, traj) = self
            .prop
            .with(solution.state)
            .until_epoch_with_traj(solution.achievement_epoch)?;

        // Build the partials
        let xf_dual = OrbitDual::from(xf.orbit);

        let mut is_bplane_tgt = false;
        for obj in &self.objectives {
            if obj.parameter.is_b_plane() {
                is_bplane_tgt = true;
            }
        }

        // Build the B-Plane once, if needed
        let b_plane = if is_bplane_tgt {
            Some(BPlane::from_dual(xf_dual)?)
        } else {
            None
        };

        let mut converged = true;
        let mut param_errors = Vec::new();
        for obj in &self.objectives {
            let partial = if obj.parameter.is_b_plane() {
                match obj.parameter {
                    StateParameter::BdotR => b_plane.unwrap().b_r,
                    StateParameter::BdotT => b_plane.unwrap().b_t,
                    StateParameter::BLTOF => b_plane.unwrap().ltof_s,
                    _ => unreachable!(),
                }
            } else {
                xf_dual.partial_for(&obj.parameter)?
            };

            let param_err = obj.desired_value - partial.real();

            if param_err.abs() > obj.tolerance {
                converged = false;
            }
            param_errors.push(param_err);
        }
        if converged {
            Ok((xf, traj))
        } else {
            let mut objmsg = String::from("");
            for (i, obj) in self.objectives.iter().enumerate() {
                objmsg.push_str(&format!(
                    "{:?} = {:.3} BUT should be {:.3} (+/- {:.1e}) ",
                    obj.parameter, param_errors[i], obj.desired_value, obj.tolerance
                ));
            }
            Err(NyxError::TargetError(objmsg))
        }
    }
}
