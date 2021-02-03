use crate::celestia::{Frame, GuidanceMode, Orbit, SpacecraftState};
use crate::dimensions::Vector3;
use crate::errors::NyxError;

mod finiteburns;
pub use finiteburns::{FiniteBurns, Mnvr};

mod ruggiero;
pub use ruggiero::Ruggiero;

/// Defines a thruster with a maximum isp and a maximum thrust.
#[derive(Copy, Clone, Debug)]
pub struct Thruster {
    pub thrust: f64,
    pub isp: f64,
}

/// The `ThrustControl` trait handles control laws, optimizations, and other such methods for
/// controlling the overall thrust direction when tied to a `Spacecraft`. For delta V control,
/// tie the DeltaVctrl to a MissionArc.
pub trait ThrustControl: Send + Sync {
    /// Returns a unit vector corresponding to the thrust direction in the inertial frame.
    fn direction(&self, state: &SpacecraftState) -> Vector3<f64>;

    /// Returns a number between [0;1] corresponding to the engine throttle level.
    /// For example, 0 means coasting, i.e. no thrusting, and 1 means maximum thrusting.
    fn throttle(&self, state: &SpacecraftState) -> f64;

    /// Prepares the controller for the next maneuver by returning the next guidance mode.
    fn next(&self, state: &SpacecraftState) -> GuidanceMode;

    /// Returns whether this thrust control has been achieved, if it has an objective
    fn achieved(&self, _state: &SpacecraftState) -> Result<bool, NyxError> {
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

fn unit_vector_from_angles(alpha: f64, beta: f64) -> Vector3<f64> {
    Vector3::new(
        alpha.sin() * beta.cos(),
        alpha.cos() * beta.cos(),
        beta.sin(),
    )
}
