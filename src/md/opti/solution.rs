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

use crate::linalg::SVector;
use crate::md::objective::Objective;
use crate::md::ui::*;
pub use crate::md::{Variable, Vary};
use std::fmt;
use std::time::Duration;

/// Defines a targeter solution
#[derive(Clone, Debug)]
pub struct TargeterSolution<const V: usize, const O: usize> {
    /// The corrected spacecraft state at the correction epoch
    pub corrected_state: Spacecraft,
    /// The state at which the objectives are achieved
    pub achieved_state: Spacecraft,
    /// The correction vector applied
    pub correction: SVector<f64, V>,
    /// The kind of correction (position or velocity)
    pub variables: [Variable; V],
    /// The errors achieved
    pub achieved_errors: SVector<f64, O>,
    /// The objectives set in the targeter
    pub achieved_objectives: [Objective; O],
    /// The number of iterations required
    pub iterations: usize,
    /// Computation duration
    pub computation_dur: Duration,
}

impl<const V: usize, const O: usize> fmt::Display for TargeterSolution<V, O> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut objmsg = String::from("");
        for (i, obj) in self.achieved_objectives.iter().enumerate() {
            objmsg.push_str(&format!(
                "\n\t\t{:?} = {:.3} (wanted {:.3} ± {:.1e})",
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
