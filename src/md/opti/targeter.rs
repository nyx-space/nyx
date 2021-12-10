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

use crate::errors::TargetingError;
use crate::linalg::DVector;
use crate::md::objective::Objective;
use crate::md::ui::*;
use crate::md::StateParameter;
pub use crate::md::{Variable, Vary};
use crate::propagators::error_ctrl::ErrorCtrl;
use std::convert::TryInto;
use std::fmt;
use std::time::Duration;

/// Defines a targeter solution
#[derive(Clone, Debug)]
pub struct TargeterSolution {
    /// The corrected spacecraft state at the correction epoch
    pub corrected_state: Spacecraft,
    /// The state at which the objectives are achieved
    pub achieved_state: Spacecraft,
    /// The correction vector applied
    pub correction: DVector<f64>,
    /// The kind of correction (position or velocity)
    pub variables: Vec<Variable>,
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

        let mut corrmsg = format!(
            "Correction @ {}:",
            self.corrected_state.epoch().as_gregorian_utc_str()
        );
        let mut is_only_position = true;
        let mut is_only_velocity = true;
        for (i, var) in self.variables.iter().enumerate() {
            let unit = match var.component {
                Vary::PositionX | Vary::PositionY | Vary::PositionZ => {
                    is_only_velocity = false;
                    "m"
                }
                Vary::VelocityX | Vary::VelocityY | Vary::VelocityZ => {
                    is_only_position = false;
                    "m/s"
                }
                _ => {
                    is_only_position = false;
                    ""
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

        writeln!(
            f,
            "Targeter solution correcting {:?} (converged in {:.3} seconds, {} iterations):\n\t{}\n\tAchieved @ {}:{}\n\tCorrected state:\n\t\t{}\n\t\t{:x}\n\tAchieved state:\n\t\t{}\n\t\t{:x}",
            self.variables.iter().map(|v| format!("{:?}", v.component)).collect::<Vec<String>>(),
            self.computation_dur.as_secs_f64(), self.iterations, corrmsg, self.achieved_state.epoch().as_gregorian_utc_str(), objmsg, self.corrected_state, self.corrected_state, self.achieved_state, self.achieved_state
        )
    }
}

/// The target is a differential corrector.
#[derive(Clone)]
pub struct Targeter<'a, E: ErrorCtrl> {
    /// The propagator setup (kind, stages, etc.)
    // pub prop: &'a Propagator<'a, D, E>,
    pub prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
    /// The list of objectives of this targeter
    pub objectives: Vec<Objective>,
    /// An optional frame (and Cosm) to compute the objectives in.
    /// Needed if the propagation frame is separate from objectives frame (e.g. for B Plane targeting).
    pub objective_frame: Option<(Frame, Arc<Cosm>)>,
    /// The kind of correction to apply to achieve the objectives
    pub variables: Vec<Variable>,
    /// The frame in which the correction should be applied, must be either a local frame or inertial
    pub correction_frame: Option<Frame>,
    /// Maximum number of iterations
    pub iterations: usize,
}

impl<'a, E: ErrorCtrl> fmt::Display for Targeter<'a, E> {
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

impl<'a, E: ErrorCtrl> Targeter<'a, E> {
    /// Create a new Targeter which will apply an impulsive delta-v correction.
    pub fn new(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        variables: Vec<Variable>,
        objectives: Vec<Objective>,
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
            correction_frame: None,
        }
    }

    /// Create a new Targeter which will apply an impulsive delta-v correction on all components of the VNC frame. By default, max step is 0.5 km/s.
    pub fn vnc(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        objectives: Vec<Objective>,
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
            objective_frame: None,
            correction_frame: Some(Frame::VNC),
        }
    }

    /// Create a new Targeter which will apply an impulsive delta-v correction on the specified components of the VNC frame.
    pub fn vnc_with_components(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        variables: Vec<Variable>,
        objectives: Vec<Objective>,
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

    /// Create a new Targeter which will apply an impulsive delta-v correction.
    pub fn delta_v(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        objectives: Vec<Objective>,
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
            objective_frame: None,
            correction_frame: None,
        }
    }

    /// Create a new Targeter which will MOVE the position of the spacecraft at the correction epoch
    pub fn delta_r(
        prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
        objectives: Vec<Objective>,
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
            objective_frame: None,
            correction_frame: None,
        }
    }

    /// Runs the targeter using finite differencing (for now).
    #[allow(clippy::identity_op)]
    pub fn try_achieve_from(
        &self,
        initial_state: Spacecraft,
        correction_epoch: Epoch,
        achievement_epoch: Epoch,
    ) -> Result<TargeterSolution, NyxError> {
        self.try_achieve_fd(initial_state, correction_epoch, achievement_epoch)
    }

    /// Apply a correction and propagate to achievement epoch. Also checks that the objectives are indeed matched
    pub fn apply(&self, solution: &TargeterSolution) -> Result<Spacecraft, NyxError> {
        let (xf, _) = self.apply_with_traj(solution)?;
        Ok(xf)
    }

    /// Apply a correction and propagate to achievement epoch, return the final state and trajectory.
    /// Also checks that the objectives are indeed matched.
    pub fn apply_with_traj(
        &self,
        solution: &TargeterSolution,
    ) -> Result<(Spacecraft, Traj<Spacecraft>), NyxError> {
        // Propagate until achievement epoch
        let (xf, traj) = self
            .prop
            .with(solution.corrected_state)
            .until_epoch_with_traj(solution.achieved_state.epoch())?;

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
            Err(NyxError::Targeter(TargetingError::Verification(objmsg)))
        }
    }
}
