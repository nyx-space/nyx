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
    pub value_precision: f64,
}

/// Defines the kind of correction to apply in the targeter
/// Note that the targeter will error if the variables do not affect the correction desired.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Corrector {
    Position,
    Velocity,
    All,
}

/// Defines a targeter. The targeter will use a differential corrector if the problem is fully defined.
/// Otherwise, it will use PANOC (TODO -- might be in another setup unsure yet).
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
    /// The list of variables
    pub variables: Vec<StateParameter>,
    /// The kind of correction to apply to achieve the objectives
    pub corrector: Corrector,
    /// Maximum number of iterations (defaults to 25)
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
                obj.parameter, obj.desired_value, obj.value_precision
            ));
        }

        write!(
            f,
            "Targeter:\n\tVary: {:?}\n\tObjectives: {}\n\tCorrect: {:?}",
            self.variables, objmsg, self.corrector
        )
    }
}

impl<'a, D: Dynamics<StateType = Spacecraft>, E: ErrorCtrl> Targeter<'a, D, E>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::PropVecSize>
        + Allocator<f64, <D::StateType as State>::Size>,
{
    pub fn try_achieve_from(
        &self,
        initial_state: Spacecraft,
        correction_epoch: Epoch,
        achivement_epoch: Epoch,
    ) -> Result<Vector3<f64>, NyxError> {
        if self.variables.is_empty() {
            return Err(NyxError::UnderdeterminedProblem);
        } else if self.variables.len() > 6 {
            return Err(NyxError::OverdeterminedProblem);
        }

        // Let's immediately make sure that the correction can be applied.
        let orbit_dual = OrbitDual::from(initial_state.orbit);

        if self.corrector == Corrector::Position {
            if self.variables.len() > 3 {
                return Err(NyxError::OverdeterminedProblem);
            }
            // Check that at least one of the parameters changes the position
            let mut will_affect = false;
            for variable in &self.variables {
                let partials = orbit_dual.partial_for(variable)?;
                if partials.wtr_x().abs() > 0.0
                    || partials.wtr_y().abs() > 0.0
                    || partials.wtr_z().abs() > 0.0
                {
                    will_affect = true;
                } else {
                    warn!("Variable {:?} has no effect on position", variable);
                }
            }

            if !will_affect {
                return Err(NyxError::CorrectionIneffective(format!(
                    "Correcting {:?} has no effect on position",
                    self.variables
                )));
            }
        } else if self.corrector == Corrector::Velocity {
            if self.variables.len() > 3 {
                return Err(NyxError::OverdeterminedProblem);
            }
            // Check that at least one of the parameters changes the position
            let mut will_affect = false;
            for variable in &self.variables {
                let partials = orbit_dual.partial_for(variable)?;
                if partials.wtr_vx().abs() > 0.0
                    || partials.wtr_vy().abs() > 0.0
                    || partials.wtr_vz().abs() > 0.0
                {
                    will_affect = true;
                } else {
                    warn!("Variable {:?} has no effect on velocity", variable);
                }
            }

            if !will_affect {
                return Err(NyxError::CorrectionIneffective(format!(
                    "Correcting {:?} has no effect on velocity",
                    self.variables
                )));
            }
        }

        // Now we know that the problem is correctly defined, so let's propagate as is to the epoch
        // where the correction should be applied.
        let mut xi = self
            .prop
            .with(initial_state)
            .until_epoch(correction_epoch)?;

        // Store the total correction in Vector3
        // TODO: Branch for a "Both" correction
        let mut total_correction = Vector3::zeros();

        // STM index to split on
        let split_on = if self.corrector == Corrector::Velocity {
            3
        } else {
            0
        };

        let mut prev_err_norm = std::f64::INFINITY;

        for _ in 0..=self.iterations {
            // Now, enable the trajectory STM for this state so we can apply the correction
            xi.enable_traj_stm();

            let xf = self.prop.with(xi).until_epoch(achivement_epoch)?;

            let phi_k_to_0 = xf.stm();
            let phi_inv = match phi_k_to_0
                .fixed_slice::<U3, U3>(0, split_on)
                .to_owned()
                .try_inverse()
            {
                Some(inv) => inv,
                None => return Err(NyxError::SingularStateTransitionMatrix),
            };

            // Build the partials
            let xf_dual = OrbitDual::from(xf.orbit);

            // Build the error vector
            let mut param_errors = Vec::new();
            let mut jac_rows = Vec::new();
            let mut converged = true;

            for obj in &self.objectives {
                let partial = xf_dual.partial_for(&obj.parameter)?;
                let param_err = -partial.real() + obj.desired_value;

                if param_err.abs() > obj.value_precision {
                    converged = false;
                }
                param_errors.push(param_err);

                if self.corrector == Corrector::Position {
                    jac_rows.push(vec![partial.wtr_x(), partial.wtr_y(), partial.wtr_z()]);
                } else {
                    jac_rows.push(vec![partial.wtr_vx(), partial.wtr_vy(), partial.wtr_vz()]);
                }
            }
            if converged {
                return Ok(total_correction);
            }

            // We haven't converged yet, so let's build the error vector
            let err_vector = DVector::from(param_errors);

            if err_vector.norm() >= prev_err_norm {
                return Err(NyxError::CorrectionIneffective(
                    format!("Targeting is ineffective at reducing the error:\nprev err norm: {:.3} km\tcur err norm: {:.3} km", prev_err_norm, err_vector.norm())
                ));
            }

            prev_err_norm = err_vector.norm();

            // And the Jacobian
            let mut jacobian = DMatrix::from_element(self.objectives.len(), 3, 0.0);
            for (i, row) in jac_rows.iter().enumerate() {
                jacobian[(i, 0)] = row[0];
                jacobian[(i, 1)] = row[1];
                jacobian[(i, 2)] = row[2];
            }

            // Compute the correction at xf
            // XXX: Do I need to invert it??
            let delta = jacobian * err_vector;

            // Apply this correction back to the correction epoch using the STM
            let delta_next = phi_inv * delta;

            // And finally apply it to the xi
            if self.corrector == Corrector::Position {
                xi.orbit.x += delta_next[0];
                xi.orbit.y += delta_next[1];
                xi.orbit.z += delta_next[2];
            } else {
                xi.orbit.vx += delta_next[0];
                xi.orbit.vy += delta_next[1];
                xi.orbit.vz += delta_next[2];
            }

            total_correction += delta_next;
        }

        Err(NyxError::MaxIterReached(format!(
            "Failed after {} iterations: {}",
            self.iterations, self
        )))
    }
}
