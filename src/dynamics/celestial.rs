use super::Dynamics;

extern crate nalgebra as na;
use std::f64;
use self::na::{U1, U3, U6, Vector6, VectorN};

pub struct TwoBody {
    pos_vel: Vector6<f64>,
    mu: f64,
}

impl TwoBody {
    pub fn new(state: &Vector6<f64>, mu: f64) -> TwoBody {
        TwoBody {
            pos_vel: *state,
            mu: mu,
        }
    }
}

impl Dynamics for TwoBody {
    type StateSize = U6;
    fn time(&self) -> f64 {
        0.0
    }

    fn state(&self) -> &VectorN<f64, Self::StateSize> {
        &self.pos_vel
    }

    fn set_state(&mut self, _new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        self.pos_vel = *new_state;
    }

    fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let radius = state.fixed_slice::<U3, U1>(0, 0);
        let velocity = state.fixed_slice::<U3, U1>(3, 0);
        let body_acceleration = (-self.mu / radius.norm().powi(3)) * radius;
        Vector6::from_iterator(velocity.iter().chain(body_acceleration.iter()).cloned())
    }
}
