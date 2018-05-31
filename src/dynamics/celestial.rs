use super::Dynamics;
use celestia::CelestialBody;

extern crate nalgebra as na;
use self::na::{U3, U6, Vector6, VectorN};
use std::f64;

/// `TwoBody` exposes the equations of motion for a simple two body propagation.
#[derive(Copy, Clone)]
pub struct TwoBody {
    time: f64,
    pos_vel: Vector6<f64>,
    mu: f64,
}

impl TwoBody {
    /// Initialize TwoBody dynamics given a provided gravitional parameter (as `mu`)
    pub fn from_state_vec_with_gm(state: &Vector6<f64>, mu: f64) -> TwoBody {
        TwoBody {
            time: 0.0,
            pos_vel: *state,
            mu: mu,
        }
    }

    /// Initialize TwoBody dynamics around a provided `CelestialBody` from the provided state vector (cf. nyx::celestia).
    pub fn from_state_vec<B: CelestialBody>(state: &Vector6<f64>) -> TwoBody {
        TwoBody {
            time: 0.0,
            pos_vel: *state,
            mu: B::gm(),
        }
    }
}

impl Dynamics for TwoBody {
    type StateSize = U6;
    fn time(&self) -> f64 {
        self.time
    }

    fn state(&self) -> VectorN<f64, Self::StateSize> {
        self.pos_vel
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        self.time = new_t;
        self.pos_vel = *new_state;
    }

    fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let radius = state.fixed_rows::<U3>(0).into_owned();
        let velocity = state.fixed_rows::<U3>(3).into_owned();
        let body_acceleration = (-self.mu / radius.norm().powi(3)) * radius;
        Vector6::from_iterator(velocity.iter().chain(body_acceleration.iter()).cloned())
    }
}
