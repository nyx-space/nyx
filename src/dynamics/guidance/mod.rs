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

use crate::cosmic::{Frame, GuidanceMode, Orbit, Spacecraft, STD_GRAVITY};
use crate::errors::NyxError;
use crate::linalg::Vector3;

mod finiteburns;
pub use finiteburns::{FiniteBurns, Mnvr};

mod ruggiero;
pub use ruggiero::Ruggiero;

/// Defines a thruster with a maximum isp and a maximum thrust.
#[derive(Copy, Clone, Debug)]
pub struct Thruster {
    /// The thrust is to be provided in Newtons
    pub thrust: f64,
    /// The Isp is to be provided in seconds
    pub isp: f64,
}

impl Thruster {
    /// Returns the exhaust velocity v_e in meters per second
    pub fn exhaust_velocity(&self) -> f64 {
        self.isp * STD_GRAVITY
    }
}

/// The `GuidanceLaw` trait handles guidance laws, optimizations, and other such methods for
/// controlling the overall thrust direction when tied to a `Spacecraft`. For delta V control,
/// tie the DeltaVctrl to a MissionArc.
pub trait GuidanceLaw: Send + Sync {
    /// Returns a unit vector corresponding to the thrust direction in the inertial frame.
    fn direction(&self, state: &Spacecraft) -> Vector3<f64>;

    /// Returns a number between [0;1] corresponding to the engine throttle level.
    /// For example, 0 means coasting, i.e. no thrusting, and 1 means maximum thrusting.
    fn throttle(&self, state: &Spacecraft) -> f64;

    /// Prepares the controller for the next maneuver by returning the next guidance mode.
    fn next(&self, state: &Spacecraft) -> GuidanceMode;

    /// Returns whether this thrust control has been achieved, if it has an objective
    fn achieved(&self, _state: &Spacecraft) -> Result<bool, NyxError> {
        Err(NyxError::NoObjectiveDefined)
    }
}

/// Goals used for sub-optimal controls
#[derive(Copy, Clone, Debug)]
pub enum Achieve {
    Sma { target: f64, tol: f64 },
    Ecc { target: f64, tol: f64 },
    Inc { target: f64, tol: f64 },
    Raan { target: f64, tol: f64 },
    Aop { target: f64, tol: f64 },
}

impl Achieve {
    pub fn achieved(&self, state: &Orbit) -> bool {
        match *self {
            Achieve::Sma { target, tol } => (state.sma() - target).abs() < tol,
            Achieve::Ecc { target, tol } => (state.ecc() - target).abs() < tol,
            Achieve::Inc { target, tol } => (state.inc() - target).abs() < tol,
            Achieve::Raan { target, tol } => (state.raan() - target).abs() < tol,
            Achieve::Aop { target, tol } => (state.aop() - target).abs() < tol,
        }
    }
}

/// Converts the alpha (in-plane) and beta (out-of-plane) angles in the RCN frame to the unit vector in the RCN frame
fn unit_vector_from_plane_angles(alpha: f64, beta: f64) -> Vector3<f64> {
    Vector3::new(
        alpha.sin() * beta.cos(),
        alpha.cos() * beta.cos(),
        beta.sin(),
    )
}

/// Converts the provided unit vector into in-plane and out-of-plane angles in the RCN frame, returned in radians
fn plane_angles_from_unit_vector(vhat: Vector3<f64>) -> (f64, f64) {
    (vhat[0].atan2(vhat[1]), vhat[2].asin())
}

#[test]
fn name() {
    use std::f64::consts::FRAC_PI_6;
    let alpha = 5.0 * FRAC_PI_6;
    let beta = 2.0 * FRAC_PI_6;
    let (alpha2, beta2) = plane_angles_from_unit_vector(unit_vector_from_plane_angles(alpha, beta));
    assert!((alpha2 - alpha).abs() < 2e-16);
    assert!((beta2 - beta).abs() < 2e-16);
}
