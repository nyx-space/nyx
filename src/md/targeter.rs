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
use super::StateParameter;
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DMatrix, DVector, DefaultAllocator, Vector3, U3};
use crate::md::ui::*;
use crate::propagators::error_ctrl::ErrorCtrl;
use std::fmt;

/// Defines a state parameter event finder
#[derive(Copy, Clone, Debug)]
pub struct Objective {
    /// The state parameter to target
    pub parameter: StateParameter,
    /// The desired self.desired_value, must be in the same units as the state parameter
    pub desired_value: f64,
    /// The precision on the desired value
    pub tolerance: f64,
}

impl Objective {
    /// Match a specific value for the parameter.
    /// By default, the tolerance on the parameter is 0.1 times whatever unit is the default for that parameter.
    /// For example, a radius event will seek the requested value at the decimeter level, and an angle event will seek it at the tenth of a degree.
    pub fn new(parameter: StateParameter, desired_value: f64) -> Self {
        Self::within_tolerance(
            parameter,
            desired_value,
            1e2 * parameter.default_event_precision(),
        )
    }

    /// Match a specific value for the parameter to hit the specified value with the provided tolerance on the value
    pub fn within_tolerance(parameter: StateParameter, desired_value: f64, tolerance: f64) -> Self {
        Self {
            parameter,
            desired_value,
            tolerance,
        }
    }
}

/// Defines the kind of correction to apply in the targeter
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Vary {
    PositionX,
    PositionY,
    PositionZ,
    VelocityX,
    VelocityY,
    VelocityZ,
}

/// Defines a targeter solution
#[derive(Clone, Debug)]
pub struct TargeterSolution {
    /// The corrected spacecraft state at the correction epoch
    pub state: Spacecraft,
    /// The correction vector applied
    pub correction: Vector3<f64>,
    /// The kind of correction (position or velocity)
    pub variables: Vec<Vary>,
    /// The epoch at which the objectives are achieved
    pub achievement_epoch: Epoch,
    /// The errors achieved
    pub achieved_errors: Vec<f64>,
    /// The objectives set in the targeter
    pub achieved_objectives: Vec<Objective>,
    /// The number of iterations required
    pub iterations: usize,
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
            let unit = match var {
                Vary::PositionX | Vary::PositionY | Vary::PositionZ => {
                    is_only_velocity = false;
                    "m"
                }
                Vary::VelocityX | Vary::VelocityY | Vary::VelocityZ => {
                    is_only_position = false;
                    "m/s"
                }
            };
            corrmsg.push_str(&format!(
                "\n\t\t{:?} = {:.3} {}",
                var, self.correction[i], unit
            ));
        }

        if is_only_position {
            corrmsg.push_str(&format!(
                "|Δr| = {:.3} m",
                (self.correction[0].powi(2)
                    + self.correction[1].powi(2)
                    + self.correction[2].powi(2))
                .sqrt()
                    * 1e3
            ));
        } else if is_only_velocity {
            corrmsg.push_str(&format!(
                "|Δv| = {:.3} m/s",
                (self.correction[0].powi(2)
                    + self.correction[1].powi(2)
                    + self.correction[2].powi(2))
                .sqrt()
                    * 1e3
            ));
        }

        write!(
            f,
            "Targeter solution correcting {:?} (converged in {} iterations):\n\t{}\n\tAchieved:{}\n\tInitial state:\n\t\t{}\n\t\t{:x}",
            self.variables, self.iterations, corrmsg, objmsg, self.state, self.state
        )
    }
}

/// The target is a "simple" differential corrector.
#[derive(Clone, Debug)]
pub struct Targeter<'a, D: Dynamics<StateType = Spacecraft>, E: ErrorCtrl>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::PropVecSize>,
{
    /// The propagator setup (kind, stages, etc.)
    pub prop: Arc<&'a Propagator<'a, D, E>>,
    /// The list of objectives of this targeter
    pub objectives: Vec<Objective>,
    /// The kind of correction to apply to achieve the objectives
    pub variables: Vec<Vary>,
    /// Maximum number of iterations
    pub iterations: usize,
}

impl<'a, D: Dynamics<StateType = Spacecraft>, E: ErrorCtrl> fmt::Display for Targeter<'a, D, E>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::PropVecSize>
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
        + Allocator<f64, <D::StateType as State>::PropVecSize>
        + Allocator<f64, <D::StateType as State>::Size>,
{
    /// Create a new Targeter which will apply an instantaneous delta-v correction.
    /// Defaults to 100 iterations which should be sufficient for all well-defined problems
    pub fn delta_v(prop: Arc<&'a Propagator<'a, D, E>>, objectives: Vec<Objective>) -> Self {
        Self {
            prop,
            objectives,
            variables: vec![Vary::VelocityX, Vary::VelocityY, Vary::VelocityZ],
            iterations: 300,
        }
    }

    /// Create a new Targeter which will MOVE the position of the spacecraft at the correction epoch
    /// Defaults to 100 iterations which should be sufficient for all well-defined problems
    pub fn delta_r(prop: Arc<&'a Propagator<'a, D, E>>, objectives: Vec<Objective>) -> Self {
        Self {
            prop,
            objectives,
            variables: vec![Vary::PositionX, Vary::PositionY, Vary::PositionZ],
            iterations: 100,
        }
    }

    /// Differential correction using hyperdual numbers for the objectives
    pub fn try_achieve_from(
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
            }
        }

        // Now we know that the problem is correctly defined, so let's propagate as is to the epoch
        // where the correction should be applied.
        let xi_start = self
            .prop
            .with(initial_state)
            .until_epoch(correction_epoch)?;

        let mut xi = xi_start;

        // Store the total correction in Vector3
        let mut total_correction = Vector3::zeros();

        // STM index to split on
        // let split_on = if self.variables == Vary::Velocity {
        //     3
        // } else {
        //     0
        // };

        let mut prev_err_norm = std::f64::INFINITY;

        for it in 0..=self.iterations {
            // Now, enable the trajectory STM for this state so we can apply the correction
            xi.enable_traj_stm();

            let xf = self.prop.with(xi).until_epoch(achievement_epoch)?;

            let phi_k_to_0 = xf.stm();
            println!("{}", phi_k_to_0);
            // let phi_drdv = phi_k_to_0.fixed_slice::<U3, U3>(3, split_on).to_owned();
            // let phi_inv = match phi_k_to_0
            //     .fixed_slice::<U3, U3>(0, split_on)
            //     .to_owned()
            //     .try_inverse()
            // {
            //     Some(inv) => inv,
            //     None => return Err(NyxError::SingularStateTransitionMatrix),
            // };

            // Build the partials
            let xf_dual = OrbitDual::from(xf.orbit);

            // Build the error vector
            let mut param_errors = Vec::new();
            let mut jac_rows = Vec::new();
            let mut converged = true;

            // Build the B-Plane once, if needed
            let b_plane = if is_bplane_tgt {
                Some(BPlane::from_dual(xf_dual)?)
            } else {
                None
            };

            // If we have as many objectives as variables, then we will not perform a Least Squares optimization.
            let is_full = self.objectives.len() == self.variables.len();

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

                println!("{}\n{:?}\n{}", obj.desired_value, obj.parameter, partial);

                if is_full {
                    // Only grab the partials that we need
                    let mut tmp_row = Vec::with_capacity(self.variables.len());
                    for var in &self.variables {
                        tmp_row.push(match var {
                            Vary::PositionX => partial.wtr_x(),
                            Vary::PositionY => partial.wtr_y(),
                            Vary::PositionZ => partial.wtr_z(),
                            Vary::VelocityX => partial.wtr_vx(),
                            Vary::VelocityY => partial.wtr_vy(),
                            Vary::VelocityZ => partial.wtr_vz(),
                        });
                    }
                    jac_rows.push(tmp_row);
                } else {
                    // Else, we will try a gradient descent approach on all of the state parameters.
                    jac_rows.push(vec![
                        partial.wtr_x(),
                        partial.wtr_y(),
                        partial.wtr_z(),
                        partial.wtr_vx(),
                        partial.wtr_vy(),
                        partial.wtr_vz(),
                    ]);
                }
            }

            // Build debugging information
            let mut objmsg = Vec::new();
            for (i, obj) in self.objectives.iter().enumerate() {
                objmsg.push(format!(
                    "\t{:?}: error = {:.3}\t desired = {:.3} (+/- {:.1e})",
                    obj.parameter, param_errors[i], obj.desired_value, obj.tolerance
                ));
            }

            if converged {
                let mut state = xi_start;
                // XXX: This will break is 4 or more Vary
                for (i, var) in self.variables.iter().enumerate() {
                    match var {
                        Vary::PositionX => state.orbit.x += total_correction[i],
                        Vary::PositionY => state.orbit.y += total_correction[i],
                        Vary::PositionZ => state.orbit.z += total_correction[i],
                        Vary::VelocityX => state.orbit.vx += total_correction[i],
                        Vary::VelocityY => state.orbit.vy += total_correction[i],
                        Vary::VelocityZ => state.orbit.vz += total_correction[i],
                    }
                }

                let sol = TargeterSolution {
                    state,
                    correction: total_correction,
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

            // And the Jacobian
            let mut jac = DMatrix::from_element(self.objectives.len(), 6, 0.0);
            for (i, row) in jac_rows.iter().enumerate() {
                jac[(i, 0)] = row[0];
                jac[(i, 1)] = row[1];
                jac[(i, 2)] = row[2];
                jac[(i, 3)] = row[3];
                jac[(i, 4)] = row[4];
                jac[(i, 5)] = row[5];
            }

            println!("{}", phi_k_to_0);
            let phi_obj = jac * phi_k_to_0;

            // Compute the least squares on this Phi with the objectives
            let phi_obj_inv = match (&phi_obj * &phi_obj.transpose()).try_inverse() {
                Some(inv) => inv,
                None => return Err(NyxError::SingularStateTransitionMatrix),
            };

            // Compute the correction at xf and map it to a state error in position and velocity space
            let delta = &phi_obj.transpose() * phi_obj_inv * err_vector;

            info!("Targeter -- Iteration #{}", it);
            for obj in &objmsg {
                info!("{}", obj);
            }
            info!("Mapped {:?} error = {:.3}", self.variables, delta.norm());

            // And finally apply it to the xi
            for (i, var) in self.variables.iter().enumerate() {
                match var {
                    Vary::PositionX => {
                        xi.orbit.x += delta[0];
                        total_correction[i] += delta[0];
                    }
                    Vary::PositionY => {
                        xi.orbit.y += delta[1];
                        total_correction[i] += delta[1];
                    }
                    Vary::PositionZ => {
                        xi.orbit.z += delta[2];
                        total_correction[i] += delta[2];
                    }
                    Vary::VelocityX => {
                        xi.orbit.vx += delta[3];
                        total_correction[i] += delta[3];
                    }
                    Vary::VelocityY => {
                        xi.orbit.vy += delta[4];
                        total_correction[i] += delta[4];
                    }
                    Vary::VelocityZ => {
                        xi.orbit.vz += delta[5];
                        total_correction[i] += delta[5];
                    }
                }
            }
            // if self.variables == Vary::Position {
            //     xi.orbit.x += delta[0];
            //     xi.orbit.y += delta[1];
            //     xi.orbit.z += delta[2];

            //     total_correction[0] += delta[0];
            //     total_correction[1] += delta[1];
            //     total_correction[2] += delta[2];
            // } else {
            //     xi.orbit.vx += delta[3];
            //     xi.orbit.vy += delta[4];
            //     xi.orbit.vz += delta[5];

            //     total_correction[0] += delta[3];
            //     total_correction[1] += delta[4];
            //     total_correction[2] += delta[5];
            // }
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

            if param_err.abs() > 10.0 * obj.tolerance {
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
