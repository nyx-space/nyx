use super::Dynamics;
use celestia::{Body, Planet};

extern crate nalgebra as na;
use std::f64;
use self::na::{U1, U3, U6, Vector6, VectorN};

/// TwoBody exposes the equations of motion for a simple two body propagation.
#[derive(Copy, Clone)]
pub struct TwoBody {
    time: f64, // XXX: Ugh, if I have this here, it means that each Dynamics implementor will have its own time too. Sounds like extra storage, but then again it's only an f64.
    pos_vel: Vector6<f64>,
    mu: f64,
}

impl TwoBody {
    /// Initialize TwoBody dynamics given a provided gravitional parameter (as `mu`)
    pub fn with_gm(state: &Vector6<f64>, mu: f64) -> TwoBody {
        TwoBody {
            time: 0.0,
            pos_vel: *state,
            mu: mu,
        }
    }

    /// Initialize TwoBody dynamics given a provided a reference to a `Planet` (cf. nyx::celestia).
    pub fn around(state: &Vector6<f64>, obj: &Planet) -> TwoBody {
        TwoBody {
            time: 0.0,
            pos_vel: *state,
            mu: obj.gm(),
        }
    }
}

impl Dynamics for TwoBody {
    type StateSize = U6;
    fn time(&self) -> f64 {
        self.time
    }

    fn state(&self) -> &VectorN<f64, Self::StateSize> {
        &self.pos_vel
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        self.time = new_t;
        self.pos_vel = *new_state;
    }

    fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let radius = state.fixed_slice::<U3, U1>(0, 0);
        let velocity = state.fixed_slice::<U3, U1>(3, 0);
        let body_acceleration = (-self.mu / radius.norm().powi(3)) * radius;
        Vector6::from_iterator(velocity.iter().chain(body_acceleration.iter()).cloned())
    }
}
