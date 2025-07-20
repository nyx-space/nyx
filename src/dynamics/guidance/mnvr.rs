/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::{
    ra_dec_from_unit_vector, GuidanceError, GuidanceLaw, GuidancePhysicsSnafu, LocalFrame,
};
use crate::cosmic::{GuidanceMode, Spacecraft};
use crate::dynamics::guidance::unit_vector_from_ra_dec;
use crate::linalg::Vector3;
use crate::polyfit::CommonPolynomial;
use crate::time::{Epoch, Unit};
use crate::State;
use anise::prelude::Almanac;
use hifitime::{Duration, TimeUnits};
use serde::{Deserialize, Serialize};
use snafu::ResultExt;
use std::fmt;
use std::sync::Arc;

/// Mnvr defined a single maneuver. Direction MUST be in the VNC frame (Velocity / Normal / Cross).
/// It may be used with a maneuver scheduler.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Maneuver {
    /// Start epoch of the maneuver
    pub start: Epoch,
    /// End epoch of the maneuver
    pub end: Epoch,
    /// TODO: Add a thruster group set to specify which set of thrusters to use for this maneuver, should be a key to a thruster (maybe change thruster to a hashmap actually now that I don't care about embedded stuff).
    /// Thrust level, if 1.0 use all thruster available at full power
    /// TODO: Convert this to a common polynomial as well to optimize throttle, throttle rate (and accel?)
    pub thrust_prct: f64,
    /// The representation of this maneuver.
    pub representation: MnvrRepr,
    /// The frame in which the maneuvers are defined.
    pub frame: LocalFrame,
}

impl fmt::Display for Maneuver {
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
            write!(f, "\n{}", self.representation)?;
            write!(
                f,
                "\n\tinitial dir: [{:.6}, {:.6}, {:.6}]\n\tfinal dir  : [{:.6}, {:.6}, {:.6}]",
                start_vec[0], start_vec[1], start_vec[2], end_vec[0], end_vec[1], end_vec[2]
            )
        } else {
            write!(
                f,
                "Impulsive maneuver @ {}\n{}",
                self.start, self.representation
            )
        }
    }
}

/// Defines the available maneuver representations.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum MnvrRepr {
    /// Represents the maneuver as a fixed vector in the local frame.
    Vector(Vector3<f64>),
    /// Represents the maneuver as a polynominal of azimuth (right ascension / in-plane) and elevation (declination / out of plane)
    Angles {
        azimuth: CommonPolynomial,
        elevation: CommonPolynomial,
    },
}

impl fmt::Display for MnvrRepr {
    /// Prints the polynomial with the least significant coefficients first
    #[allow(clippy::identity_op)]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            MnvrRepr::Vector(vector) => write!(f, "{vector}"),
            MnvrRepr::Angles { azimuth, elevation } => write!(
                f,
                "\tazimuth (in-plane) α: {azimuth}\n\televation (out-of-plane) β: {elevation}"
            ),
        }
    }
}

impl Maneuver {
    /// Creates an impulsive maneuver whose vector is the deltaV.
    /// TODO: This should use William's algorithm
    pub fn from_impulsive(dt: Epoch, vector: Vector3<f64>, frame: LocalFrame) -> Self {
        Self::from_time_invariant(dt, dt + Unit::Millisecond, 1.0, vector, frame)
    }

    /// Creates a maneuver from the provided time-invariant delta-v, in km/s
    pub fn from_time_invariant(
        start: Epoch,
        end: Epoch,
        thrust_lvl: f64,
        vector: Vector3<f64>,
        frame: LocalFrame,
    ) -> Self {
        Self {
            start,
            end,
            thrust_prct: thrust_lvl,
            representation: MnvrRepr::Vector(vector),
            frame,
        }
    }

    /// Return the thrust vector computed at the provided epoch
    pub fn vector(&self, epoch: Epoch) -> Vector3<f64> {
        match self.representation {
            MnvrRepr::Vector(vector) => vector,
            MnvrRepr::Angles { azimuth, elevation } => {
                let t = (epoch - self.start).to_seconds();
                let alpha = azimuth.eval(t);
                let delta = elevation.eval(t);
                unit_vector_from_ra_dec(alpha, delta)
            }
        }
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
        match self.representation {
            MnvrRepr::Vector(vector) => vector / vector.norm(),
            MnvrRepr::Angles { azimuth, elevation } => {
                let alpha = azimuth.coeff_in_order(0).unwrap();
                let delta = elevation.coeff_in_order(0).unwrap();
                unit_vector_from_ra_dec(alpha, delta)
            }
        }
    }

    /// Set the time-invariant direction for this finite burn while keeping the other components as they are
    pub fn set_direction(&mut self, vector: Vector3<f64>) -> Result<(), GuidanceError> {
        self.set_direction_and_rates(vector, self.rate(), self.accel())
    }

    /// Returns the rate of direction of the burn at the start of the burn, useful for setting new angles
    pub fn rate(&self) -> Vector3<f64> {
        match self.representation {
            MnvrRepr::Vector(_) => Vector3::zeros(),
            MnvrRepr::Angles { azimuth, elevation } => match azimuth.coeff_in_order(1) {
                Ok(alpha) => {
                    let delta = elevation.coeff_in_order(1).unwrap();
                    unit_vector_from_ra_dec(alpha, delta)
                }
                Err(_) => Vector3::zeros(),
            },
        }
    }

    /// Set the rate of direction for this finite burn while keeping the other components as they are
    pub fn set_rate(&mut self, rate: Vector3<f64>) -> Result<(), GuidanceError> {
        self.set_direction_and_rates(self.direction(), rate, self.accel())
    }

    /// Returns the acceleration of the burn at the start of the burn, useful for setting new angles
    pub fn accel(&self) -> Vector3<f64> {
        match self.representation {
            MnvrRepr::Vector(_) => Vector3::zeros(),
            MnvrRepr::Angles { azimuth, elevation } => match azimuth.coeff_in_order(2) {
                Ok(alpha) => {
                    let delta = elevation.coeff_in_order(2).unwrap();
                    unit_vector_from_ra_dec(alpha, delta)
                }
                Err(_) => Vector3::zeros(),
            },
        }
    }

    /// Set the acceleration of the direction of this finite burn while keeping the other components as they are
    pub fn set_accel(&mut self, accel: Vector3<f64>) -> Result<(), GuidanceError> {
        self.set_direction_and_rates(self.direction(), self.rate(), accel)
    }

    /// Set the initial direction, direction rate, and direction acceleration for this finite burn
    pub fn set_direction_and_rates(
        &mut self,
        dir: Vector3<f64>,
        rate: Vector3<f64>,
        accel: Vector3<f64>,
    ) -> Result<(), GuidanceError> {
        if rate.norm() < f64::EPSILON && accel.norm() < f64::EPSILON {
            // Set as a vector
            self.representation = MnvrRepr::Vector(dir)
        } else {
            let (alpha, delta) = ra_dec_from_unit_vector(dir);
            if alpha.is_nan() || delta.is_nan() {
                return Err(GuidanceError::InvalidDirection {
                    x: dir[0],
                    y: dir[1],
                    z: dir[2],
                    in_plane_deg: alpha.to_degrees(),
                    out_of_plane_deg: delta.to_degrees(),
                });
            }
            if rate.norm() < f64::EPSILON && accel.norm() < f64::EPSILON {
                self.representation = MnvrRepr::Angles {
                    azimuth: CommonPolynomial::Constant(alpha),
                    elevation: CommonPolynomial::Constant(delta),
                };
            } else {
                let (alpha_dt, delta_dt) = ra_dec_from_unit_vector(rate);
                if alpha_dt.is_nan() || delta_dt.is_nan() {
                    return Err(GuidanceError::InvalidRate {
                        x: rate[0],
                        y: rate[1],
                        z: rate[2],
                        in_plane_deg_s: alpha_dt.to_degrees(),
                        out_of_plane_deg_s: delta_dt.to_degrees(),
                    });
                }
                if accel.norm() < f64::EPSILON {
                    self.representation = MnvrRepr::Angles {
                        azimuth: CommonPolynomial::Linear(alpha_dt, alpha),
                        elevation: CommonPolynomial::Linear(delta_dt, delta),
                    };
                } else {
                    let (alpha_ddt, delta_ddt) = ra_dec_from_unit_vector(accel);
                    if alpha_ddt.is_nan() || delta_ddt.is_nan() {
                        return Err(GuidanceError::InvalidAcceleration {
                            x: accel[0],
                            y: accel[1],
                            z: accel[2],
                            in_plane_deg_s2: alpha_ddt.to_degrees(),
                            out_of_plane_deg_s2: delta_ddt.to_degrees(),
                        });
                    }

                    self.representation = MnvrRepr::Angles {
                        azimuth: CommonPolynomial::Quadratic(alpha_ddt, alpha_dt, alpha),
                        elevation: CommonPolynomial::Quadratic(delta_ddt, delta_dt, delta),
                    };
                }
            }
        }
        Ok(())
    }
}

impl GuidanceLaw for Maneuver {
    fn direction(&self, osc: &Spacecraft) -> Result<Vector3<f64>, GuidanceError> {
        match osc.mode() {
            GuidanceMode::Thrust => match self.frame {
                LocalFrame::Inertial => Ok(self.vector(osc.epoch())),
                _ => Ok(self.frame.dcm_to_inertial(osc.orbit).context({
                    GuidancePhysicsSnafu {
                        action: "computing RCN frame",
                    }
                })? * self.vector(osc.epoch())),
            },
            _ => Ok(Vector3::zeros()),
        }
    }

    fn throttle(&self, osc: &Spacecraft) -> Result<f64, GuidanceError> {
        // match self.next(osc) {
        match osc.mode() {
            GuidanceMode::Thrust => Ok(self.thrust_prct),
            _ => {
                // We aren't in maneuver mode, so return 0% throttle
                Ok(0.0)
            }
        }
    }

    fn next(&self, sc: &mut Spacecraft, _almanac: Arc<Almanac>) {
        let next_mode = if sc.epoch() >= self.start && sc.epoch() < self.end {
            GuidanceMode::Thrust
        } else {
            GuidanceMode::Coast
        };
        sc.mut_mode(next_mode);
    }
}

#[cfg(test)]
mod ut_mnvr {
    use hifitime::Epoch;
    use nalgebra::Vector3;

    use crate::dynamics::guidance::LocalFrame;

    use super::Maneuver;

    #[test]
    fn serde_mnvr() {
        let epoch = Epoch::from_gregorian_utc_at_midnight(2012, 2, 29);
        let mnvr = Maneuver::from_impulsive(epoch, Vector3::new(1.0, 1.0, 0.0), LocalFrame::RCN);

        let mnvr_yml = serde_yml::to_string(&mnvr).unwrap();
        println!("{mnvr_yml}");

        let mnvr2 = serde_yml::from_str(&mnvr_yml).unwrap();
        assert_eq!(mnvr, mnvr2);
    }
}
