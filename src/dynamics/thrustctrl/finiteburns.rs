use super::ThrustControl;
use crate::celestia::{Frame, GuidanceMode, SpacecraftState};
use crate::dimensions::Vector3;
use crate::state::TimeTagged;
use crate::time::{Epoch, TimeUnit};
use std::sync::Arc;

/// Mnvr defined a single maneuver. Direction MUST be in the VNC frame (Velocity / Normal / Cross).
/// It may be used with a maneuver scheduler.
#[derive(Copy, Clone, Debug)]
pub struct Mnvr {
    /// Start epoch of the maneuver
    pub start: Epoch,
    /// End epoch of the maneuver
    pub end: Epoch,
    /// Thrust level, if 1.0 use all thruster available at full power
    pub thrust_lvl: f64,
    /// Direction of the thrust in the VNC frame
    pub vector: Vector3<f64>,
}

impl Mnvr {
    /// Creates an instantaneous maneuver whose vector is the deltaV.
    pub fn instantaneous(dt: Epoch, vector: Vector3<f64>) -> Self {
        Self {
            start: dt,
            end: dt + TimeUnit::Microsecond,
            thrust_lvl: 1.0,
            vector,
        }
    }
}

/// A controller for a set of pre-determined maneuvers.
#[derive(Clone, Debug)]
pub struct FiniteBurns {
    /// Maneuvers should be provided in chronological order, first maneuver first in the list
    pub mnvrs: Vec<Mnvr>,
    /// The frame in which the maneuvers are defined.
    pub frame: Frame,
}

impl FiniteBurns {
    /// Builds a schedule from the vector of maneuvers, must be provided in chronological order.
    pub fn from_mnvrs(mnvrs: Vec<Mnvr>, frame: Frame) -> Arc<Self> {
        assert!(
            matches!(frame, Frame::Inertial | Frame::VNC),
            "Maneuvers must be either in the inertial frame or in a body frame"
        );
        Arc::new(Self { mnvrs, frame })
    }
}

impl ThrustControl for FiniteBurns {
    fn direction(&self, osc: &SpacecraftState) -> Vector3<f64> {
        // NOTE: We do not increment the mnvr number here. The power function is called first,
        // so we let that function handle starting and stopping of the maneuver.
        match osc.mode {
            GuidanceMode::Custom(mnvr_no) => {
                let next_mnvr = self.mnvrs[mnvr_no as usize];
                if next_mnvr.start <= osc.epoch() {
                    if matches!(self.frame, Frame::Inertial) {
                        next_mnvr.vector
                    } else {
                        osc.orbit.dcm_to_inertial(self.frame) * next_mnvr.vector
                    }
                } else {
                    Vector3::zeros()
                }
            }
            _ => Vector3::zeros(),
        }
    }

    fn throttle(&self, osc: &SpacecraftState) -> f64 {
        match osc.mode {
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

    fn next(&self, sc: &SpacecraftState) -> GuidanceMode {
        // Here, we're using the Custom field of the mode to store the current maneuver number we're executing
        match sc.mode {
            GuidanceMode::Custom(mnvr_no) => {
                if (mnvr_no as usize) < self.mnvrs.len() {
                    let cur_mnvr = self.mnvrs[mnvr_no as usize];
                    if sc.epoch() >= cur_mnvr.end {
                        GuidanceMode::Custom(mnvr_no + 1)
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
        }
    }
}
