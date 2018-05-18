use super::Dynamics;
use celestia::CelestialBody;

extern crate nalgebra as na;
use std::f64;
use self::na::{U1, U3, U6, Vector6, VectorN};

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

    /// This is a state error estimator.
    ///
    /// It returns the largest error between the position and the velocity errors using the step error from the integration (`prop_err`)
    /// given the difference in the candidate state and the previous state (`state_delta`).
    /// This error estimator is from the physical model estimator of GMAT
    /// https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/PhysicalModel.cpp#L987
    pub fn max_pos_vel_error_estimator(
        prop_err: &VectorN<f64, <TwoBody as Dynamics>::StateSize>,
        state_delta: &VectorN<f64, <TwoBody as Dynamics>::StateSize>,
    ) -> f64 {
        let delta_radius = state_delta.fixed_slice::<U3, U1>(0, 0).norm();
        let prop_err_radius = prop_err.fixed_slice::<U3, U1>(0, 0).norm();
        let err_radius = if delta_radius > 0.0 {
            prop_err_radius / delta_radius
        } else {
            prop_err_radius
        };

        let delta_velocity = state_delta.fixed_slice::<U3, U1>(3, 0).norm();
        let prop_err_velocity = state_delta.fixed_slice::<U3, U1>(3, 0).norm();
        let err_velocity = if delta_velocity > 0.0 {
            prop_err_velocity / delta_velocity
        } else {
            prop_err_velocity
        };

        // Always return the largest of both errors
        if err_radius > err_velocity {
            err_radius
        } else {
            err_velocity
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
        let radius = state.fixed_slice::<U3, U1>(0, 0);
        let velocity = state.fixed_slice::<U3, U1>(3, 0);
        let body_acceleration = (-self.mu / radius.norm().powi(3)) * radius;
        Vector6::from_iterator(velocity.iter().chain(body_acceleration.iter()).cloned())
    }
}
