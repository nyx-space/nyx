extern crate nalgebra as na;
use self::na::{Matrix6, MatrixMN, U1, U3, U36, U42, U6, Vector6, VectorN};
use super::Dynamics;
use celestia::{CelestialBody, CoordinateFrame, State};
use od::Linearization;
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

/// `TwoBody` exposes the equations of motion for a simple two body propagation. It inherently supports
/// the State Transition Matrix for orbital determination.
#[derive(Copy, Clone)]
pub struct TwoBodyWithStm {
    pub compute_stm: bool,
    pub pos_vel: Vector6<f64>,
    pub stm: Matrix6<f64>,
    time: f64,
    mu: f64,
}

impl TwoBodyWithStm {
    /// Initialize TwoBody dynamics given a provided gravitional parameter (as `mu`)
    pub fn from_state_vec_with_gm(state: &Vector6<f64>, mu: f64) -> TwoBodyWithStm {
        TwoBodyWithStm {
            compute_stm: false,
            time: 0.0,
            pos_vel: *state,
            mu: mu,
            stm: Matrix6::identity(),
        }
    }

    /// Initialize TwoBody dynamics around a provided `CelestialBody` from the provided state vector (cf. nyx::celestia).
    pub fn from_state_vec<B: CelestialBody>(state: &Vector6<f64>) -> TwoBodyWithStm {
        TwoBodyWithStm {
            compute_stm: false,
            time: 0.0,
            pos_vel: *state,
            mu: B::gm(),
            stm: Matrix6::identity(),
        }
    }

    /// Initialize TwoBody dynamics around a provided `CelestialBody` from the provided position and velocity state (cf. nyx::celestia).
    pub fn from_state<B: CelestialBody, F: CoordinateFrame>(state: State<F>) -> TwoBodyWithStm {
        TwoBodyWithStm {
            compute_stm: false,
            time: 0.0,
            pos_vel: state.to_cartesian_vec(),
            mu: B::gm(),
            stm: Matrix6::identity(),
        }
    }
}

impl Dynamics for TwoBodyWithStm {
    type StateSize = U42;
    fn time(&self) -> f64 {
        self.time
    }

    fn state(&self) -> VectorN<f64, Self::StateSize> {
        let stm_as_vec = self.stm.fixed_resize::<U1, U36>(0.0);
        VectorN::<f64, Self::StateSize>::from_iterator(self.pos_vel.iter().chain(stm_as_vec.iter()).cloned())
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        self.time = new_t;
        self.pos_vel = new_state.fixed_rows::<U6>(0).into_owned();
    }

    fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let radius = state.fixed_rows::<U3>(0).into_owned();
        let velocity = state.fixed_rows::<U3>(3).into_owned();
        let body_acceleration = (-self.mu / radius.norm().powi(3)) * radius;
        VectorN::<f64, Self::StateSize>::from_iterator(velocity.iter().chain(body_acceleration.iter()).cloned())
    }
}

impl Linearization for TwoBodyWithStm {
    type StateSize = U6;

    fn gradient(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> MatrixMN<f64, Self::StateSize, Self::StateSize> {
        let mut grad = MatrixMN::<f64, U6, U6>::zeros();
        // Top right is Identity 3x3
        grad[(3, 0)] = 1.0;
        grad[(4, 1)] = 1.0;
        grad[(2, 5)] = 1.0;
        // Bottom left is where the magic happens.
        let x = state[(0, 0)];
        let y = state[(1, 0)];
        let z = state[(2, 0)];
        let x2 = x.powi(2);
        let y2 = y.powi(2);
        let z2 = z.powi(2);
        let vx = state[(3, 0)];
        let vy = state[(4, 0)];
        let vz = state[(5, 0)];
        let r232 = (x2 + y2 + z2).powf(3.0 / 2.0);
        let r252 = (x2 + y2 + z2).powf(5.0 / 2.0);

        // Add the body perturbations
        let dAxDx = 3.0 * self.mu * x2 / r252 - self.mu / r232;
        let dAxDy = 3.0 * self.mu * x * y / r252;
        let dAxDz = 3.0 * self.mu * x * z / r252;
        let dAyDx = 3.0 * self.mu * x * y / r252;
        let dAyDy = 3.0 * self.mu * y2 / r252 - self.mu / r232;
        let dAyDz = 3.0 * self.mu * y * z / r252;
        let dAzDx = 3.0 * self.mu * x * z / r252;
        let dAzDy = 3.0 * self.mu * y * z / r252;
        let dAzDz = 3.0 * self.mu * z2 / r252 - self.mu / r232;

        // Let the gradient
        grad[(0, 3)] = dAxDx;
        grad[(0, 4)] = dAyDx;
        grad[(0, 5)] = dAzDx;
        grad[(1, 3)] = dAxDy;
        grad[(1, 4)] = dAyDy;
        grad[(1, 5)] = dAzDy;
        grad[(2, 3)] = dAxDz;
        grad[(2, 4)] = dAyDz;
        grad[(2, 5)] = dAzDz;

        grad
    }
}
