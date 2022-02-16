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

use hifitime::TimeUnits;

use crate::dynamics::guidance::Mnvr;
use crate::linalg::SVector;
use crate::md::objective::Objective;
use crate::md::ui::*;
pub use crate::md::{Variable, Vary};
use crate::polyfit::CommonPolynomial;
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

impl<const V: usize, const O: usize> TargeterSolution<V, O> {
    /// Returns whether this solution is a finite burn solution or not
    pub fn is_finite_burn(&self) -> bool {
        for var in &self.variables {
            if var.component.is_finite_burn() {
                return true;
            }
        }
        false
    }

    /// Returns a maneuver if targeter solution was a finite burn maneuver
    pub fn to_mnvr(&self) -> Result<Mnvr, NyxError> {
        if !self.is_finite_burn() {
            Err(NyxError::CustomError(
                "Not a finite burn solution".to_string(),
            ))
        } else {
            let correction_epoch = self.corrected_state.epoch();
            let achievement_epoch = self.achieved_state.epoch();
            let mut mnvr = Mnvr {
                start: correction_epoch,
                end: achievement_epoch,
                thrust_lvl: 1.0,
                alpha_inplane_radians: CommonPolynomial::Quadratic(0.0, 0.0, 0.0),
                delta_outofplane_radians: CommonPolynomial::Quadratic(0.0, 0.0, 0.0),
                frame: Frame::RCN,
            };

            for (i, var) in self.variables.iter().enumerate() {
                let corr = self.correction[i];

                // Modify the maneuver, but do not change the epochs of the maneuver unless the change is greater than one millisecond
                match var.component {
                    Vary::Duration => {
                        if corr.abs() > 1e-3 {
                            // Check that we are within the bounds
                            let init_duration_s =
                                (correction_epoch - achievement_epoch).in_seconds();
                            let acceptable_corr = var.apply_bounds(init_duration_s).seconds();
                            mnvr.end = mnvr.start + acceptable_corr;
                        }
                    }
                    Vary::EndEpoch => {
                        if corr.abs() > 1e-3 {
                            // Check that we are within the bounds
                            let total_end_corr =
                                (mnvr.end + corr.seconds() - achievement_epoch).in_seconds();
                            let acceptable_corr = var.apply_bounds(total_end_corr).seconds();
                            mnvr.end += acceptable_corr;
                        }
                    }
                    Vary::StartEpoch => {
                        if corr.abs() > 1e-3 {
                            // Check that we are within the bounds
                            let total_start_corr =
                                (mnvr.start + corr.seconds() - correction_epoch).in_seconds();
                            let acceptable_corr = var.apply_bounds(total_start_corr).seconds();
                            mnvr.end += acceptable_corr;

                            mnvr.start += corr.seconds()
                        }
                    }
                    Vary::MnvrAlpha | Vary::MnvrAlphaDot | Vary::MnvrAlphaDDot => {
                        mnvr.alpha_inplane_radians = mnvr
                            .alpha_inplane_radians
                            .add_val_in_order(corr, var.component.vec_index())
                            .unwrap();
                    }
                    Vary::MnvrDelta | Vary::MnvrDeltaDot | Vary::MnvrDeltaDDot => {
                        mnvr.delta_outofplane_radians = mnvr
                            .delta_outofplane_radians
                            .add_val_in_order(corr, var.component.vec_index())
                            .unwrap();
                    }
                    Vary::ThrustX | Vary::ThrustY | Vary::ThrustZ => {
                        let mut vector = mnvr.direction();
                        vector[var.component.vec_index()] += corr;
                        var.ensure_bounds(&mut vector[var.component.vec_index()]);
                        mnvr.set_direction(vector)?;
                    }
                    Vary::ThrustRateX | Vary::ThrustRateY | Vary::ThrustRateZ => {
                        let mut vector = mnvr.rate();
                        let idx = (var.component.vec_index() - 1) % 3;
                        vector[idx] += corr;
                        var.ensure_bounds(&mut vector[idx]);
                        mnvr.set_rate(vector)?;
                    }
                    Vary::ThrustAccelX | Vary::ThrustAccelY | Vary::ThrustAccelZ => {
                        let mut vector = mnvr.accel();
                        let idx = (var.component.vec_index() - 1) % 3;
                        vector[idx] += corr;
                        var.ensure_bounds(&mut vector[idx]);
                        mnvr.set_accel(vector)?;
                    }
                    Vary::ThrustLevel => {
                        mnvr.thrust_lvl += corr;
                        var.ensure_bounds(&mut mnvr.thrust_lvl);
                    }
                    _ => unreachable!(),
                }
            }

            Ok(mnvr)
        }
    }
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
                    is_only_velocity = false;
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
        } else if self.is_finite_burn() {
            let mnvr = self.to_mnvr().unwrap();
            corrmsg.push_str(&format!("\n\t\t{}\n", mnvr));
        }

        writeln!(
            f,
            "Targeter solution correcting {:?} (converged in {:.3} seconds, {} iterations):\n\t{}\n\tAchieved @ {}:{}\n\tCorrected state:\n\t\t{}\n\t\t{:x}\n\tAchieved state:\n\t\t{}\n\t\t{:x}",
            self.variables.iter().map(|v| format!("{:?}", v.component)).collect::<Vec<String>>(),
            self.computation_dur.as_secs_f64(), self.iterations, corrmsg, self.achieved_state.epoch().as_gregorian_utc_str(), objmsg, self.corrected_state, self.corrected_state, self.achieved_state, self.achieved_state
        )
    }
}
