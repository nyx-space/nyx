/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use hifitime::Epoch;

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

    fn maneuver_at(&self, epoch: Epoch) -> Option<&Mnvr> {
        // Find the maneuver with the closest start epoch that is less than or equal to the current epoch
        let index = self.mnvrs.binary_search_by_key(&epoch, |mnvr| mnvr.start);
        match index {
            Err(0) => None, // No maneuvers start before the current epoch
            Ok(index) => Some(&self.mnvrs[index]),
            Err(index) => Some(&self.mnvrs[index - 1]), // Return the maneuver with the closest start epoch
        }
    }
}

impl fmt::Display for FiniteBurns {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "FiniteBurns with {} maneuvers", self.mnvrs.len())
    }
}

impl GuidanceLaw for FiniteBurns {
    fn direction(&self, osc: &Spacecraft) -> Vector3<f64> {
        // NOTE: We do not increment the mnvr number here. The power function is called first,
        // so we let that function handle starting and stopping of the maneuver.
        match osc.mode() {
            GuidanceMode::Thrust => {
                if let Some(next_mnvr) = self.maneuver_at(osc.epoch()) {
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
                } else {
                    Vector3::zeros()
                }
            }
            _ => Vector3::zeros(),
        }
    }

    fn throttle(&self, osc: &Spacecraft) -> f64 {
        match osc.mode() {
            GuidanceMode::Thrust => {
                if let Some(next_mnvr) = self.maneuver_at(osc.epoch()) {
                    if next_mnvr.start <= osc.epoch() {
                        next_mnvr.thrust_prct
                    } else {
                        0.0
                    }
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
        // Only switch to thrusting mode is the thrusting is not enabled and not inhibited.
        if sc.mode() == GuidanceMode::Coast {
            // If we haven't started the maneuvers yet, let's get ready to do so by switching to the mode
            sc.mut_mode(GuidanceMode::Thrust)
        };
    }
}
