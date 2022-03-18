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

use crate::errors::TargetingError;
use crate::md::objective::Objective;
use crate::md::ui::*;
use crate::md::StateParameter;
pub use crate::md::{Variable, Vary};
use crate::propagators::error_ctrl::ErrorCtrl;
use std::convert::TryInto;
use std::fmt;

use super::solution::TargeterSolution;

/// An optimizer structure with V control variables and O objectives.
#[derive(Clone)]
pub struct Optimizer<'a, E: ErrorCtrl, const V: usize, const O: usize> {
    /// The propagator setup (kind, stages, etc.)
    pub prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
    /// The list of objectives of this targeter
    pub objectives: [Objective; O],
    /// An optional frame (and Cosm) to compute the objectives in.
    /// Needed if the propagation frame is separate from objectives frame (e.g. for B Plane targeting).
    pub objective_frame: Option<(Frame, Arc<Cosm>)>,
    /// The kind of correction to apply to achieve the objectives
    pub variables: [Variable; V],
    /// The frame in which the correction should be applied, must be either a local frame or inertial
    pub correction_frame: Option<Frame>,
    /// Maximum number of iterations
    pub iterations: usize,
}

impl<'a, E: ErrorCtrl, const V: usize, const O: usize> fmt::Display for Optimizer<'a, E, V, O> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut objmsg = String::from("");
        for obj in &self.objectives {
            objmsg.push_str(&format!("{}; ", obj));
        }

        let mut varmsg = String::from("");
        for var in &self.variables {
            varmsg.push_str(&format!("{}; ", var));
        }

        write!(
            f,
            "Targeter:\n\tObjectives: {}\n\tCorrect: {}",
            objmsg, varmsg
        )
    }
}

impl<'a, E: ErrorCtrl, const O: usize> Optimizer<'a, E, 3, O> {
    /// Create a new Targeter which will apply an impulsive delta-v correction.
    pub fn delta_v(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        objectives: [Objective; O],
    ) -> Self {
        Self {
            prop,
            objectives,
            variables: [
                Vary::VelocityX.try_into().unwrap(),
                Vary::VelocityY.try_into().unwrap(),
                Vary::VelocityZ.try_into().unwrap(),
            ],
            iterations: 100,
            objective_frame: None,
            correction_frame: None,
        }
    }

    /// Create a new Targeter which will MOVE the position of the spacecraft at the correction epoch
    pub fn delta_r(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        objectives: [Objective; O],
    ) -> Self {
        Self {
            prop,
            objectives,
            variables: [
                Vary::PositionX.try_into().unwrap(),
                Vary::PositionY.try_into().unwrap(),
                Vary::PositionZ.try_into().unwrap(),
            ],
            iterations: 100,
            objective_frame: None,
            correction_frame: None,
        }
    }

    /// Create a new Targeter which will apply an impulsive delta-v correction on all components of the VNC frame. By default, max step is 0.5 km/s.
    pub fn vnc(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        objectives: [Objective; O],
    ) -> Self {
        Self {
            prop,
            objectives,
            variables: [
                Vary::VelocityX.try_into().unwrap(),
                Vary::VelocityY.try_into().unwrap(),
                Vary::VelocityZ.try_into().unwrap(),
            ],
            iterations: 100,
            objective_frame: None,
            correction_frame: Some(Frame::VNC),
        }
    }
}

impl<'a, E: ErrorCtrl, const O: usize> Optimizer<'a, E, 4, O> {
    /// Create a new Targeter which will apply a continuous thrust for the whole duration of the segment
    pub fn thrust_dir(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        objectives: [Objective; O],
    ) -> Self {
        Self {
            prop,
            objectives,
            variables: [
                Variable::from(Vary::ThrustX),
                Variable::from(Vary::ThrustY),
                Variable::from(Vary::ThrustZ),
                Variable::from(Vary::ThrustLevel),
            ],
            iterations: 20,
            objective_frame: None,
            correction_frame: None,
        }
    }
}

impl<'a, E: ErrorCtrl, const O: usize> Optimizer<'a, E, 7, O> {
    /// Create a new Targeter which will apply a continuous thrust for the whole duration of the segment
    pub fn thrust_dir_rate(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        objectives: [Objective; O],
    ) -> Self {
        Self {
            prop,
            objectives,
            variables: [
                Variable::from(Vary::ThrustX),
                Variable::from(Vary::ThrustY),
                Variable::from(Vary::ThrustZ),
                Variable::from(Vary::ThrustLevel),
                Variable::from(Vary::ThrustRateX),
                Variable::from(Vary::ThrustRateY),
                Variable::from(Vary::ThrustRateZ),
            ],
            iterations: 50,
            objective_frame: None,
            correction_frame: None,
        }
    }
}

impl<'a, E: ErrorCtrl, const O: usize> Optimizer<'a, E, 10, O> {
    /// Create a new Targeter which will apply a continuous thrust for the whole duration of the segment
    pub fn thrust_profile(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        objectives: [Objective; O],
    ) -> Self {
        Self {
            prop,
            objectives,
            variables: [
                Variable::from(Vary::ThrustX),
                Variable::from(Vary::ThrustY),
                Variable::from(Vary::ThrustZ),
                Variable::from(Vary::ThrustLevel),
                Variable::from(Vary::ThrustRateX),
                Variable::from(Vary::ThrustRateY),
                Variable::from(Vary::ThrustRateZ),
                Variable::from(Vary::ThrustAccelX),
                Variable::from(Vary::ThrustAccelY),
                Variable::from(Vary::ThrustAccelZ),
            ],
            iterations: 50,
            objective_frame: None,
            correction_frame: None,
        }
    }
}

impl<'a, E: ErrorCtrl, const V: usize, const O: usize> Optimizer<'a, E, V, O> {
    /// Create a new Targeter which will apply an impulsive delta-v correction.
    pub fn new(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        variables: [Variable; V],
        objectives: [Objective; O],
    ) -> Self {
        Self {
            prop,
            objectives,
            variables,
            iterations: 100,
            objective_frame: None,
            correction_frame: None,
        }
    }

    /// Create a new Targeter which will apply an impulsive delta-v correction.
    pub fn in_frame(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        variables: [Variable; V],
        objectives: [Objective; O],
        objective_frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Self {
        Self {
            prop,
            objectives,
            variables,
            iterations: 100,
            objective_frame: Some((objective_frame, cosm)),
            correction_frame: None,
        }
    }

    /// Create a new Targeter which will apply an impulsive delta-v correction on the specified components of the VNC frame.
    pub fn vnc_with_components(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        variables: [Variable; V],
        objectives: [Objective; O],
    ) -> Self {
        Self {
            prop,
            objectives,
            variables,
            iterations: 100,
            objective_frame: None,
            correction_frame: Some(Frame::VNC),
        }
    }

    /// Runs the targeter using finite differencing (for now).
    #[allow(clippy::identity_op)]
    pub fn try_achieve_from(
        &self,
        initial_state: Spacecraft,
        correction_epoch: Epoch,
        achievement_epoch: Epoch,
    ) -> Result<TargeterSolution<V, O>, NyxError> {
        self.try_achieve_fd(initial_state, correction_epoch, achievement_epoch)
    }

    /// Apply a correction and propagate to achievement epoch. Also checks that the objectives are indeed matched
    pub fn apply(&self, solution: &TargeterSolution<V, O>) -> Result<Spacecraft, NyxError> {
        let (xf, _) = self.apply_with_traj(solution)?;
        Ok(xf)
    }

    /// Apply a correction and propagate to achievement epoch, return the final state and trajectory.
    /// Also checks that the objectives are indeed matched.
    pub fn apply_with_traj(
        &self,
        solution: &TargeterSolution<V, O>,
    ) -> Result<(Spacecraft, Traj<Spacecraft>), NyxError> {
        let (xf, traj) = match solution.to_mnvr() {
            Ok(mnvr) => {
                println!("{}", mnvr);
                let mut prop = self.prop.clone();
                prop.dynamics = prop.dynamics.with_guidance_law(Arc::new(mnvr));
                prop.with(solution.corrected_state)
                    .until_epoch_with_traj(solution.achieved_state.epoch())?
            }
            Err(_) => {
                // This isn't a finite burn maneuver, let's just apply the correction
                // Propagate until achievement epoch
                self.prop
                    .with(solution.corrected_state)
                    .until_epoch_with_traj(solution.achieved_state.epoch())?
            }
        };

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
                    "{:?} = {:.3} BUT should be {:.3} (± {:.1e}) (error = {:.3})",
                    obj.parameter,
                    obj.desired_value + param_errors[i],
                    obj.desired_value,
                    obj.tolerance,
                    param_errors[i]
                ));
            }
            Err(NyxError::Targeter(TargetingError::Verification(objmsg)))
        }
    }
}
