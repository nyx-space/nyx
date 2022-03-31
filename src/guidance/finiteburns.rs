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
extern crate rayon;

use super::{GuidanceLaw, Mnvr};
use crate::cosmic::{Frame, GuidanceMode, Spacecraft};
use crate::linalg::Vector3;
use crate::State;
use std::fmt;
use std::sync::Arc;
/// A controller for a set of pre-determined maneuvers.
#[derive(Clone, Debug)]
pub struct FiniteBurns {
    /// Maneuvers should be provided in chronological order, first maneuver first in the list
    pub mnvrs: Vec<Mnvr>,
}

impl FiniteBurns {
    /// Builds a schedule from the vector of maneuvers, must be provided in chronological order.
    pub fn from_mnvrs(mnvrs: Vec<Mnvr>) -> Arc<Self> {
        Arc::new(Self { mnvrs })
    }
}

impl fmt::Display for FiniteBurns {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FiniteBurns with {} maneuvers", self.mnvrs.len())
    }
}

impl GuidanceLaw<GuidanceMode> for FiniteBurns {
    fn direction(&self, osc: &Spacecraft) -> Vector3<f64> {
        // NOTE: We do not increment the mnvr number here. The power function is called first,
        // so we let that function handle starting and stopping of the maneuver.
        match osc.mode() {
            GuidanceMode::Custom(mnvr_no) => {
                let next_mnvr = self.mnvrs[mnvr_no as usize];
                if next_mnvr.start <= osc.epoch() {
                    if matches!(next_mnvr.frame, Frame::Inertial) {
                        next_mnvr.vector(osc.epoch())
                    } else {
                        osc.orbit.dcm_from_traj_frame(next_mnvr.frame).unwrap()
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
        match osc.mode() {
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

    fn next(&self, sc: &mut Spacecraft) {
        // Here, we're using the Custom field of the mode to store the current maneuver number we're executing
        let next_mode = match sc.mode() {
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
        };
        sc.mut_mode(next_mode);
    }
}
