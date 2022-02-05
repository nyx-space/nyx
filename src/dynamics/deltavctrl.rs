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

use crate::cosmic::{Frame, Orbit};
use crate::linalg::Vector3;
use crate::State;

pub use super::guidance::Mnvr;

/// The `DeltaVctrl` trait handles control laws, optimizations, and other such methods for
/// controlling the change in velocity of a point mass during a mission arc (`MissionArc`).
pub trait DeltaVctrl
where
    Self: Clone + Sized,
{
    /// Returns the control vector corresponding to the change in velocity direction in the inertial frame.
    fn ctrl_vector(&self, state: &Orbit) -> Vector3<f64>;

    /// Prepares the controller for the next maneuver (called from set_state of the dynamics).
    fn next(&mut self, state: &Orbit);
}

#[derive(Clone, Debug)]
pub struct ImpulsiveBurns {
    /// Maneuvers should be provided in chronological order, first maneuver first in the list
    pub mnvrs: Vec<Mnvr>,
    pub mnvr_no: usize,
}

impl ImpulsiveBurns {
    /// Builds a schedule from the vector of maneuvers, must be provided in chronological order.
    pub fn from_mnvrs(mnvrs: Vec<Mnvr>) -> Self {
        Self { mnvrs, mnvr_no: 0 }
    }
}

impl DeltaVctrl for ImpulsiveBurns {
    fn ctrl_vector(&self, state: &Orbit) -> Vector3<f64> {
        if self.mnvr_no >= self.mnvrs.len() {
            Vector3::zeros()
        } else {
            let next_mnvr = self.mnvrs[self.mnvr_no];
            if next_mnvr.start <= state.dt && next_mnvr.end >= state.dt {
                state.dcm_from_traj_frame(Frame::VNC).unwrap() * next_mnvr.vector(state.epoch())
            } else {
                Vector3::zeros()
            }
        }
    }

    fn next(&mut self, state: &Orbit) {
        if self.mnvr_no < self.mnvrs.len() {
            let cur_mnvr = self.mnvrs[self.mnvr_no];
            if state.dt >= cur_mnvr.end {
                self.mnvr_no += 1;
            }
        }
    }
}
