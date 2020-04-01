use super::na::Vector3;
use celestia::{Frame, State};

pub use super::thrustctrl::Mnvr;

/// The `DeltaVctrl` trait handles control laws, optimizations, and other such methods for
/// controlling the change in velocity of a point mass during a mission arc (`MissionArc`).
pub trait DeltaVctrl
where
    Self: Clone + Sized,
{
    /// Returns the control vector corresponding to the change in velocity direction in the inertial frame.
    fn ctrl_vector(&self, state: &State) -> Vector3<f64>;

    /// Prepares the controller for the next maneuver (called from set_state of the dynamics).
    fn next(&mut self, state: &State);
}

#[derive(Clone, Debug)]
pub struct InstantBurns {
    /// Maneuvers should be provided in chronological order, first maneuver first in the list
    pub mnvrs: Vec<Mnvr>,
    pub mnvr_no: usize,
}

impl InstantBurns {
    /// Builds a schedule from the vector of maneuvers, must be provided in chronological order.
    pub fn from_mnvrs(mnvrs: Vec<Mnvr>) -> Self {
        Self { mnvrs, mnvr_no: 0 }
    }
}

impl DeltaVctrl for InstantBurns {
    fn ctrl_vector(&self, state: &State) -> Vector3<f64> {
        if self.mnvr_no >= self.mnvrs.len() {
            Vector3::zeros()
        } else {
            let next_mnvr = self.mnvrs[self.mnvr_no];
            if next_mnvr.start <= state.dt && next_mnvr.end >= state.dt {
                state.dcm_to_inertial(Frame::VNC) * next_mnvr.vector
            } else {
                Vector3::zeros()
            }
        }
    }

    fn next(&mut self, state: &State) {
        if self.mnvr_no < self.mnvrs.len() {
            let cur_mnvr = self.mnvrs[self.mnvr_no];
            if state.dt >= cur_mnvr.end {
                self.mnvr_no += 1;
            }
        }
    }
}
