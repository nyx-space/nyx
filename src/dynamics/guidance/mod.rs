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

use crate::cosmic::{GuidanceMode, Orbit, Spacecraft, STD_GRAVITY};
use crate::errors::{NyxError, StateError};
use crate::linalg::Vector3;
use anise::astro::PhysicsResult;
use anise::errors::PhysicsError;
use anise::math::rotation::DCM;
use anise::prelude::Almanac;
use serde::{Deserialize, Serialize};

mod finiteburns;
pub use finiteburns::FiniteBurns;

mod mnvr;
pub use mnvr::{Maneuver, MnvrRepr};

mod ruggiero;
pub use ruggiero::{Objective, Ruggiero, StateParameter};
use snafu::Snafu;

use std::fmt;
use std::sync::Arc;

/// Defines a thruster with a maximum isp and a maximum thrust.
#[allow(non_snake_case)]
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct Thruster {
    /// The thrust is to be provided in Newtons
    pub thrust_N: f64,
    /// The Isp is to be provided in seconds
    pub isp_s: f64,
}

impl Thruster {
    /// Returns the exhaust velocity v_e in meters per second
    pub fn exhaust_velocity_m_s(&self) -> f64 {
        self.isp_s * STD_GRAVITY
    }
}

/// The `GuidanceLaw` trait handles guidance laws, optimizations, and other such methods for
/// controlling the overall thrust direction when tied to a `Spacecraft`. For delta V control,
/// tie the DeltaVctrl to a MissionArc.
pub trait GuidanceLaw: fmt::Display + Send + Sync {
    /// Returns a unit vector corresponding to the thrust direction in the inertial frame.
    fn direction(&self, osc_state: &Spacecraft) -> Result<Vector3<f64>, GuidanceError>;

    /// Returns a number between [0;1] corresponding to the engine throttle level.
    /// For example, 0 means coasting, i.e. no thrusting, and 1 means maximum thrusting.
    fn throttle(&self, osc_state: &Spacecraft) -> Result<f64, GuidanceError>;

    /// Updates the state of the BaseSpacecraft for the next maneuver, e.g. prepares the controller for the next maneuver
    fn next(&self, next_state: &mut Spacecraft, almanac: Arc<Almanac>);

    /// Returns whether this thrust control has been achieved, if it has an objective
    fn achieved(&self, _osc_state: &Spacecraft) -> Result<bool, GuidanceError> {
        Err(GuidanceError::NoGuidanceObjectiveDefined)
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
pub fn plane_angles_from_unit_vector(vhat: Vector3<f64>) -> (f64, f64) {
    (vhat[1].atan2(vhat[0]), vhat[2].asin())
}

/// Converts the alpha (in-plane) and beta (out-of-plane) angles in the RCN frame to the unit vector in the RCN frame
pub(crate) fn unit_vector_from_ra_dec(alpha: f64, delta: f64) -> Vector3<f64> {
    Vector3::new(
        delta.cos() * alpha.cos(),
        delta.cos() * alpha.sin(),
        delta.sin(),
    )
}

/// Converts the provided unit vector into in-plane and out-of-plane angles in the RCN frame, returned in radians
pub(crate) fn ra_dec_from_unit_vector(vhat: Vector3<f64>) -> (f64, f64) {
    let alpha = vhat[1].atan2(vhat[0]);
    let delta = vhat[2].asin();
    (alpha, delta)
}

#[derive(Debug, PartialEq, Snafu)]
pub enum GuidanceError {
    #[snafu(display("No thruster attached to spacecraft"))]
    NoThrustersDefined,
    #[snafu(display("Throttle is not between 0.0 and 1.0: {ratio}"))]
    ThrottleRatio { ratio: f64 },
    #[snafu(display("Invalid finite burn control direction u = [{x}, {y}, {z}] => i-plane = {in_plane_deg} deg, Delta = {out_of_plane_deg} deg",))]
    InvalidDirection {
        x: f64,
        y: f64,
        z: f64,
        in_plane_deg: f64,
        out_of_plane_deg: f64,
    },
    #[snafu(display("Invalid finite burn control rate u = [{x}, {y}, {z}] => in-plane = {in_plane_deg_s} deg/s, out of plane = {out_of_plane_deg_s} deg/s",))]
    InvalidRate {
        x: f64,
        y: f64,
        z: f64,
        in_plane_deg_s: f64,
        out_of_plane_deg_s: f64,
    },
    #[snafu(display("Invalid finite burn control acceleration u = [{x}, {y}, {z}] => in-plane = {in_plane_deg_s2} deg/s^2, out of plane = {out_of_plane_deg_s2} deg/s^2",))]
    InvalidAcceleration {
        x: f64,
        y: f64,
        z: f64,
        in_plane_deg_s2: f64,
        out_of_plane_deg_s2: f64,
    },
    #[snafu(display("when {action} encountered {source}"))]
    GuidancePhysicsError {
        action: &'static str,
        source: PhysicsError,
    },
    #[snafu(display(
        "An objective based analysis or control was attempted, but no objective was defined"
    ))]
    NoGuidanceObjectiveDefined,
    #[snafu(display("{param} is not a control variable in this guidance law"))]
    InvalidControl { param: StateParameter },
    #[snafu(display("guidance encountered {source}"))]
    GuidState { source: StateError },
}

/// Local frame options, used notably for guidance laws.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum LocalFrame {
    Inertial,
    RIC,
    VNC,
    RCN,
}

impl LocalFrame {
    pub fn dcm_to_inertial(&self, state: Orbit) -> PhysicsResult<DCM> {
        match self {
            LocalFrame::Inertial => Ok(DCM::identity(
                state.frame.orientation_id,
                state.frame.orientation_id,
            )),
            LocalFrame::RIC => state.dcm_from_ric_to_inertial(),
            LocalFrame::VNC => state.dcm_from_vnc_to_inertial(),
            LocalFrame::RCN => state.dcm_from_rcn_to_inertial(),
        }
    }
}

#[test]
fn ra_dec_from_vec() {
    use std::f64::consts::{FRAC_PI_2, PI, TAU};
    let mut delta = -FRAC_PI_2;
    let mut alpha = 0.0;
    loop {
        loop {
            let unit_v = unit_vector_from_ra_dec(alpha, delta);
            let (alpha2, delta2) = ra_dec_from_unit_vector(unit_v);
            assert!((alpha - alpha2).abs() < f64::EPSILON);
            assert!((delta - delta2).abs() < f64::EPSILON);
            alpha += TAU * 0.1; // Increment right ascension by one tenth of a circle
            if alpha > PI {
                alpha = 0.0;
                break;
            }
        }
        delta += TAU * 0.1; // Increment declination by one tenth of a circle
        if delta > FRAC_PI_2 {
            break;
        }
    }
}
