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

use super::{plane_angles_from_unit_vector, unit_vector_from_plane_angles, GuidanceLaw};
use crate::cosmic::{Frame, GuidanceMode, Spacecraft};
use crate::linalg::Vector3;
use crate::polyfit::CommonPolynomial;
use crate::time::{Epoch, TimeUnit};
use crate::State;
use std::fmt;
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
    /// The interpolation polynomial for the in-plane angle
    pub alpha_inplane: CommonPolynomial,
    /// The interpolation polynomial for the out-of-plane angle
    pub beta_outofplane: CommonPolynomial,
}

impl fmt::Display for Mnvr {
    /// Prints the polynomial with the least significant coefficients first
    #[allow(clippy::identity_op)]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.end - self.start >= 1 * TimeUnit::Millisecond {
            write!(f, "Finite burn maneuver @ {} for {}\n\tin-plane angle a(t)={}\n\tin-plane angle a(t)={}", self.start, self.end-self.start, self.alpha_inplane, self.beta_outofplane)
        } else {
            write!(
                f,
                "Impulsive maneuver @ {}\n\tin-plane angle a(t)={}\n\tout-plane angle b(t)={}",
                self.start, self.alpha_inplane, self.beta_outofplane
            )
        }
    }
}

impl Mnvr {
    /// Creates an impulsive maneuver whose vector is the deltaV.
    /// TODO: This should use William's algorithm
    pub fn from_impulsive(dt: Epoch, vector: Vector3<f64>) -> Self {
        Self::from_time_invariant(dt, dt + TimeUnit::Millisecond, 1.0, vector)
    }

    /// Creates a manneuver from the provided time-invariant delta-v, in km/s
    pub fn from_time_invariant(
        start: Epoch,
        end: Epoch,
        thrust_lvl: f64,
        vector: Vector3<f64>,
    ) -> Self {
        // Convert to angles
        let (alpha, beta) = plane_angles_from_unit_vector(vector);
        Self {
            start,
            end,
            thrust_lvl,
            alpha_inplane: CommonPolynomial::Constant(alpha),
            beta_outofplane: CommonPolynomial::Constant(beta),
        }
    }

    /// Converts the input delta-v vector in km/s at the provided Epoch to a finite burn
    /// Uses Copernicus algorithm as described in "AAS 12-236: Recent Improvements to the Copernicus Trajectory Design and Optimization System" by Williams et al.
    /// The vector is expected to be in the same frame as the spaceraft's orbit.
    /// Convergence criteria:
    ///     1. If the change in the duration of the burn is less than 1 seconds; or
    ///     2. If the magnitude of the thrust vector matches the impulsive maneuver to less than 1 mm/s
    #[cfg(feature = "broken-donotuse")]
    pub fn impulsive_to_finite<'a, E: ErrorCtrl>(
        epoch: Epoch,
        dv: Vector3<f64>,
        spacecraft: Spacecraft,
        prop: &'a Propagator<'a, SpacecraftDynamics, E>,
    ) -> Result<Self, NyxError> {
        // Imports specific to this broken feature
        use crate::md::ui::{Propagator, SpacecraftDynamics};
        use crate::propagators::ErrorCtrl;
        use crate::NyxError;

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
        // Compute the control rates at the time of the impulsive maneuver (tdv)
        let (alpha_tdv, beta_tdv) = plane_angles_from_unit_vector(u);
        let (alpha_ddot_tdv, beta_ddot_tdv) = plane_angles_from_unit_vector(u_ddot);
        // Build the maneuver polynomial angles from these
        let alpha_inplane = CommonPolynomial::Quadratic(alpha_ddot_tdv, 0.0, alpha_tdv);
        let beta_outofplane = CommonPolynomial::Quadratic(beta_ddot_tdv, 0.0, beta_tdv);

        // Compute a few thruster parameters
        let thruster = spacecraft.thruster.as_ref().unwrap();
        let c_km_s = thruster.exhaust_velocity() * 1e-3;

        let mut delta_tfb = 1000.0
            * ((c_km_s * spacecraft.mass_kg()) / thruster.thrust)
            * (1.0 - (-dv.norm() / c_km_s).exp());

        // Build the estimated maneuver
        let mut mnvr_guess = Self {
            start: epoch - 0.5 * delta_tfb * TimeUnit::Second,
            end: epoch + 0.5 * delta_tfb * TimeUnit::Second,
            thrust_lvl: 1.0,
            alpha_inplane,
            beta_outofplane,
        };

        // Start iteration
        let max_iter = 25;
        let mut converged = false;

        for _ in 0..max_iter {
            println!("{}", mnvr_guess);
            // Start by propagating the dynamics back to the estimated start of the maneuver without thrust.
            prop.dynamics = prop.dynamics.without_ctrl();
            let sc_start = prop
                .with(sc_at_dv_epoch.with_guidance_mode(GuidanceMode::Coast))
                .for_duration(-(1.0 + delta_tfb) * TimeUnit::Second)?;
            // Propagate for the duration of the burn without maneuver enabled to compute the baseline
            let sc_post_no_mnvr = prop.with(sc_start).until_epoch(mnvr_guess.end)?;

            // Build the finite burn
            let fb_guess = FiniteBurns {
                mnvrs: vec![mnvr_guess],
                frame: Frame::Inertial,
            };

            // Propagate for the duration of the burn
            prop.set_max_step(mnvr_guess.end - mnvr_guess.start);
            prop.dynamics = prop.dynamics.with_ctrl(Arc::new(fb_guess));
            let sc_post_mnvr = prop
                .with(sc_start.with_guidance_mode(GuidanceMode::Custom(0)))
                .until_epoch(mnvr_guess.end)?;
            // Compute the applied delta-v
            let accumulated_dv = sc_post_mnvr.orbit.velocity() - sc_post_no_mnvr.orbit.velocity();
            println!(
                "Wanted {:.3} m/s {}\nGot {:.3} m/s {}\nError = {:.3} m/s",
                dv.norm() * 1e3,
                dv,
                accumulated_dv.norm() * 1e3,
                accumulated_dv,
                (accumulated_dv - dv).norm() * 1e3
            );
            if (accumulated_dv - dv).norm() * 1e3 < 10.0 {
                // Less than 10cm/s error, we're good
                converged = true;
                break;
            } else {
                if accumulated_dv.norm() > dv.norm() {
                    // Burn was too long let's decrease the burn duration
                    delta_tfb = 0.25 * delta_tfb;
                } else {
                    // Burn was too short, let's icnrease burn duration
                    delta_tfb = 1.25 * delta_tfb;
                }
                // Update the maneuver
                mnvr_guess = Self {
                    start: epoch - delta_tfb * TimeUnit::Second,
                    end: epoch + delta_tfb * TimeUnit::Second,
                    thrust_lvl: 1.0,
                    alpha_inplane,
                    beta_outofplane,
                };
            }
        }

        if !converged {
            return Err(NyxError::MaxIterReached(format!(
                "Finite burn failed to converge after {} iterations.",
                max_iter
            )));
        }
        Ok(mnvr_guess)
    }

    /// Return the thrust vector computed at the provided epoch
    pub fn vector(&self, epoch: Epoch) -> Vector3<f64> {
        let t = (epoch - self.start).in_seconds();
        let alpha = self.alpha_inplane.eval(t);
        let beta = self.beta_outofplane.eval(t);
        unit_vector_from_plane_angles(alpha, beta)
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
                        next_mnvr.vector(osc.epoch())
                    } else {
                        osc.orbit.dcm_from_traj_frame(self.frame).unwrap()
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

#[cfg(feature = "broken-donotuse")]
#[test]
fn name_impulsive_to_finite() {
    use crate::cosmic::Cosm;
    use crate::dynamics::guidance::Thruster;
    use crate::dynamics::OrbitalDynamics;
    use crate::Orbit;

    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    /* Define the parking orbit */
    let epoch = Epoch::from_gregorian_utc_at_noon(2022, 03, 04);
    let start = Orbit::keplerian_altitude(300.0, 0.01, 30.0, 90.0, 90.0, 60.0, epoch, eme2k);

    /* Build the spacecraft -- really only the mass is needed here */
    let sc = Spacecraft {
        orbit: start,
        dry_mass_kg: 1000.0,
        fuel_mass_kg: 2000.0,
        thruster: Some(Thruster {
            thrust: 3000.0,
            isp: 300.0,
        }),
        mode: GuidanceMode::Thrust,

        ..Default::default()
    };

    let prop = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));

    let epoch = Epoch::from_gregorian_utc_hms(2022, 03, 04, 12, 00, 00);
    let dv = Vector3::new(-0.002, -0.011, 0.001);

    let m = Mnvr::impulsive_to_finite(epoch, dv, sc, &prop).unwrap();
    println!("{}", m);
}
