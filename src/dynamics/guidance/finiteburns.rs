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

use super::GuidanceLaw;
use crate::cosmic::{Frame, GuidanceMode, Spacecraft};
use crate::linalg::Vector3;
use crate::md::ui::{Propagator, SpacecraftDynamics};
use crate::propagators::ErrorCtrl;
use crate::time::{Epoch, TimeUnit};
use crate::{NyxError, State};
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
    /// Creates an impulsive maneuver whose vector is the deltaV.
    pub fn from_impulsive(dt: Epoch, vector: Vector3<f64>) -> Self {
        Self {
            start: dt,
            end: dt + TimeUnit::Microsecond,
            thrust_lvl: 1.0,
            vector,
        }
    }

    /// Creates a manneuver from the provided time-invariant delta-v, in km/s
    pub fn from_time_invariant(
        start: Epoch,
        end: Epoch,
        thrust_lvl: f64,
        dvector: Vector3<f64>,
    ) -> Self {
        Self {
            start,
            end,
            thrust_lvl,
            vector: dvector,
        }
    }

    /// Converts the input delta-v vector in km/s at the provided Epoch to a finite burn
    /// Uses Copernicus algorithm as described in "AAS 12-236: Recent Improvements to the Copernicus Trajectory Design and Optimization System" by Williams et al.
    /// The vector is expected to be in the same frame as the spaceraft's orbit.
    /// Convergence criteria:
    ///     1. If the change in the duration of the burn is less than 1 seconds; or
    ///     2. If the magnitude of the thrust vector matches the impulsive maneuver to less than 1 mm/s
    pub fn impulsive_to_finite<'a, E: ErrorCtrl>(
        epoch: Epoch,
        dv: Vector3<f64>,
        spacecraft: Spacecraft,
        prop: &'a Propagator<'a, SpacecraftDynamics, E>,
    ) -> Result<Self, NyxError> {
        if spacecraft.thruster.is_none() {
            // Can't do any conversion to finite burns without a thruster
            return Err(NyxError::CtrlExistsButNoThrusterAvail);
        }

        // Clone the dynamics
        let mut prop = prop.clone();
        // Propagate to the dv epoch
        prop.dynamics = prop.dynamics.without_ctrl();
        let sc_at_dv_epoch = prop.with(spacecraft).until_epoch(epoch)?;
        // Calculate the u, dot u (=0) and ddot u from this state
        let u = dv / dv.norm();
        let r = sc_at_dv_epoch.orbit.radius();
        let rmag = sc_at_dv_epoch.orbit.rmag();
        let u_ddot = (3.0 * sc_at_dv_epoch.orbit.frame.gm() / rmag.powi(5))
            * (r.dot(&u) * r - (r.dot(&u).powi(2) * u));
        // Compute the control rates

        // Compute a few thruster parameters
        let thruster = spacecraft.thruster.as_ref().unwrap();
        let c = thruster.exhaust_velocity();

        let mut delta_tfb =
            ((c * spacecraft.mass_kg()) / thruster.thrust) * (1.0 - (-dv.norm() / c).exp());

        // Start iteration
        let max_iter = 25;
        let mut converged = false;

        for _ in 0..max_iter {
            let fb_guess = FiniteBurns {
                mnvrs: vec![Self {
                    start: epoch - 0.5 * delta_tfb * TimeUnit::Second,
                    end: epoch + 0.5 * delta_tfb * TimeUnit::Second,
                    thrust_lvl: 1.0,
                    vector: u,
                }],
                frame: spacecraft.orbit.frame,
            };
            // Set

            prop.dynamics = prop.dynamics.with_ctrl(Arc::new(fb_guess));
        }

        if !converged {
            return Err(NyxError::MaxIterReached(format!(
                "Finite burn failed to converge after {} iterations.",
                max_iter
            )));
        }

        todo!()
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

impl GuidanceLaw for FiniteBurns {
    fn direction(&self, osc: &Spacecraft) -> Vector3<f64> {
        // NOTE: We do not increment the mnvr number here. The power function is called first,
        // so we let that function handle starting and stopping of the maneuver.
        match osc.mode {
            GuidanceMode::Custom(mnvr_no) => {
                let next_mnvr = self.mnvrs[mnvr_no as usize];
                if next_mnvr.start <= osc.epoch() {
                    if matches!(self.frame, Frame::Inertial) {
                        next_mnvr.vector
                    } else {
                        osc.orbit.dcm_from_traj_frame(self.frame).unwrap() * next_mnvr.vector
                    }
                } else {
                    Vector3::zeros()
                }
            }
            _ => Vector3::zeros(),
        }
    }

    fn throttle(&self, osc: &Spacecraft) -> f64 {
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

    fn next(&self, sc: &Spacecraft) -> GuidanceMode {
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
