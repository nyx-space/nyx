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

use super::{ra_dec_from_unit_vector, GuidanceErrors, GuidanceLaw};
use crate::cosmic::{Frame, GuidanceMode, Spacecraft};
use crate::dynamics::guidance::unit_vector_from_ra_dec;
use crate::linalg::Vector3;
use crate::polyfit::CommonPolynomial;
use crate::time::{Epoch, Unit};
use crate::State;
use hifitime::{Duration, TimeUnits};
use std::fmt;

/// Mnvr defined a single maneuver. Direction MUST be in the VNC frame (Velocity / Normal / Cross).
/// It may be used with a maneuver scheduler.
#[derive(Copy, Clone, Debug)]
pub struct Mnvr {
    /// Start epoch of the maneuver
    pub start: Epoch,
    /// End epoch of the maneuver
    pub end: Epoch,
    /// TODO: Add a thruster group set to specify which set of thrusters to use for this maneuver, should be a key to a thruster (maybe change thruster to a hashmap actually now that I don't care about embedded stuff).
    /// Thrust level, if 1.0 use all thruster available at full power
    /// TODO: Convert this to a common polynomial as well to optimize throttle, throttle rate (and accel?)
    pub thrust_prct: f64,
    /// The interpolation polynomial for the in-plane angle
    pub alpha_inplane_radians: CommonPolynomial,
    /// The interpolation polynomial for the out-of-plane angle
    pub delta_outofplane_radians: CommonPolynomial,
    /// The frame in which the maneuvers are defined.
    pub frame: Frame,
}

impl fmt::Display for Mnvr {
    /// Prints the polynomial with the least significant coefficients first
    #[allow(clippy::identity_op)]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.end != self.start {
            let start_vec = self.vector(self.start);
            let end_vec = self.vector(self.end);
            write!(
                f,
                "Finite burn maneuver @ {:.2}% on {} for {} (ending on {})",
                100.0 * self.thrust_prct,
                self.start,
                self.end - self.start,
                self.end,
            )?;
            write!(
                f,
                "\n\tin-plane angle α: {}\n\tout-of-plane angle β: {}",
                self.alpha_inplane_radians, self.delta_outofplane_radians
            )?;
            write!(
                f,
                "\n\tinitial dir: [{:.6}, {:.6}, {:.6}]\n\tfinal dir  : [{:.6}, {:.6}, {:.6}]",
                start_vec[0], start_vec[1], start_vec[2], end_vec[0], end_vec[1], end_vec[2]
            )
        } else {
            write!(
                f,
                "Impulsive maneuver @ {}\n\tin-plane angle α: {}\n\tout-of-plane angle β: {}",
                self.start, self.alpha_inplane_radians, self.delta_outofplane_radians
            )
        }
    }
}

impl Mnvr {
    /// Creates an impulsive maneuver whose vector is the deltaV.
    /// TODO: This should use William's algorithm
    pub fn from_impulsive(dt: Epoch, vector: Vector3<f64>, frame: Frame) -> Self {
        Self::from_time_invariant(dt, dt + Unit::Millisecond, 1.0, vector, frame)
    }

    /// Creates a maneuver from the provided time-invariant delta-v, in km/s
    pub fn from_time_invariant(
        start: Epoch,
        end: Epoch,
        thrust_lvl: f64,
        vector: Vector3<f64>,
        frame: Frame,
    ) -> Self {
        // Convert to angles
        let (alpha, delta) = ra_dec_from_unit_vector(vector);
        Self {
            start,
            end,
            thrust_prct: thrust_lvl,
            alpha_inplane_radians: CommonPolynomial::Constant(alpha),
            delta_outofplane_radians: CommonPolynomial::Constant(delta),
            frame,
        }
    }

    /// Return the thrust vector computed at the provided epoch
    pub fn vector(&self, epoch: Epoch) -> Vector3<f64> {
        let t = (epoch - self.start).to_seconds();
        let alpha = self.alpha_inplane_radians.eval(t);
        let delta = self.delta_outofplane_radians.eval(t);
        unit_vector_from_ra_dec(alpha, delta)
    }

    /// Return the duration of this maneuver
    pub fn duration(&self) -> Duration {
        self.end - self.start
    }

    /// Return whether this is an antichronological maneuver
    pub fn antichronological(&self) -> bool {
        self.duration().abs() > 1.microseconds() && self.duration() < 1.microseconds()
    }

    /// Returns the direction of the burn at the start of the burn, useful for setting new angles
    pub fn direction(&self) -> Vector3<f64> {
        let alpha = self.alpha_inplane_radians.coeff_in_order(0).unwrap();
        let delta = self.delta_outofplane_radians.coeff_in_order(0).unwrap();
        unit_vector_from_ra_dec(alpha, delta)
    }

    /// Set the time-invariant direction for this finite burn while keeping the other components as they are
    pub fn set_direction(&mut self, vector: Vector3<f64>) -> Result<(), GuidanceErrors> {
        self.set_direction_and_rates(vector, self.rate(), self.accel())
    }

    /// Returns the rate of direction of the burn at the start of the burn, useful for setting new angles
    pub fn rate(&self) -> Vector3<f64> {
        match self.alpha_inplane_radians.coeff_in_order(1) {
            Ok(alpha) => {
                let delta = self.delta_outofplane_radians.coeff_in_order(1).unwrap();
                unit_vector_from_ra_dec(alpha, delta)
            }
            Err(_) => Vector3::zeros(),
        }
    }

    /// Set the rate of direction for this finite burn while keeping the other components as they are
    pub fn set_rate(&mut self, rate: Vector3<f64>) -> Result<(), GuidanceErrors> {
        self.set_direction_and_rates(self.direction(), rate, self.accel())
    }

    /// Returns the acceleration of the burn at the start of the burn, useful for setting new angles
    pub fn accel(&self) -> Vector3<f64> {
        match self.alpha_inplane_radians.coeff_in_order(2) {
            Ok(alpha) => {
                let delta = self.delta_outofplane_radians.coeff_in_order(2).unwrap();
                unit_vector_from_ra_dec(alpha, delta)
            }
            Err(_) => Vector3::zeros(),
        }
    }

    /// Set the acceleration of the direction of this finite burn while keeping the other components as they are
    pub fn set_accel(&mut self, accel: Vector3<f64>) -> Result<(), GuidanceErrors> {
        self.set_direction_and_rates(self.direction(), self.rate(), accel)
    }

    /// Set the initial direction, direction rate, and direction acceleration for this finite burn
    pub fn set_direction_and_rates(
        &mut self,
        dir: Vector3<f64>,
        rate: Vector3<f64>,
        accel: Vector3<f64>,
    ) -> Result<(), GuidanceErrors> {
        let (alpha, delta) = ra_dec_from_unit_vector(dir);
        if alpha.is_nan() || delta.is_nan() {
            return Err(GuidanceErrors::InvalidDirection {
                x: dir[0],
                y: dir[1],
                z: dir[2],
                in_plane_deg: alpha.to_degrees(),
                out_of_plane_deg: delta.to_degrees(),
            });
        }
        if rate.norm() < 2e-16 && accel.norm() < 2e-16 {
            self.alpha_inplane_radians = CommonPolynomial::Constant(alpha);
            self.delta_outofplane_radians = CommonPolynomial::Constant(delta);
        } else {
            let (alpha_dt, delta_dt) = ra_dec_from_unit_vector(rate);
            if alpha_dt.is_nan() || delta_dt.is_nan() {
                return Err(GuidanceErrors::InvalidRate {
                    x: rate[0],
                    y: rate[1],
                    z: rate[2],
                    in_plane_deg_s: alpha_dt.to_degrees(),
                    out_of_plane_deg_s: delta_dt.to_degrees(),
                });
            }
            if accel.norm() < 2e-16 {
                self.alpha_inplane_radians = CommonPolynomial::Linear(alpha_dt, alpha);
                self.delta_outofplane_radians = CommonPolynomial::Linear(delta_dt, delta);
            } else {
                let (alpha_ddt, delta_ddt) = ra_dec_from_unit_vector(accel);
                if alpha_ddt.is_nan() || delta_ddt.is_nan() {
                    return Err(GuidanceErrors::InvalidAcceleration {
                        x: accel[0],
                        y: accel[1],
                        z: accel[2],
                        in_plane_deg_s2: alpha_ddt.to_degrees(),
                        out_of_plane_deg_s2: delta_ddt.to_degrees(),
                    });
                }
                self.alpha_inplane_radians =
                    CommonPolynomial::Quadratic(alpha_ddt, alpha_dt, alpha);
                self.delta_outofplane_radians =
                    CommonPolynomial::Quadratic(delta_ddt, delta_dt, delta);
            }
        }
        Ok(())
    }
}

impl GuidanceLaw for Mnvr {
    fn direction(&self, osc: &Spacecraft) -> Vector3<f64> {
        match osc.mode() {
            GuidanceMode::Thrust => {
                if matches!(self.frame, Frame::Inertial) {
                    self.vector(osc.epoch())
                } else {
                    osc.orbit.dcm_from_traj_frame(self.frame).unwrap() * self.vector(osc.epoch())
                }
            }
            _ => Vector3::zeros(),
        }
    }

    fn throttle(&self, osc: &Spacecraft) -> f64 {
        // match self.next(osc) {
        match osc.mode() {
            GuidanceMode::Thrust => self.thrust_prct,
            _ => {
                // We aren't in maneuver mode, so return 0% throttle
                0.0
            }
        }
        // self.thrust_lvl
    }

    fn next(&self, sc: &mut Spacecraft) {
        let next_mode = if sc.epoch() >= self.start && sc.epoch() < self.end {
            GuidanceMode::Thrust
        } else {
            GuidanceMode::Coast
        };
        sc.mut_mode(next_mode);
    }
}
