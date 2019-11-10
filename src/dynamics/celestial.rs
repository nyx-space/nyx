extern crate hyperdual;

use self::hyperdual::linalg::norm;
use self::hyperdual::{hyperspace_from_vector, Float, Hyperdual};
use super::hifitime::Epoch;
use super::na::{DimName, Matrix6, Vector3, Vector6, VectorN, U3, U36, U42, U6, U7};
use super::Dynamics;
use celestia::{Cosm, Geoid, State};
use od::AutoDiffDynamics;
use std::f64;

/// `CelestialDynamics` provides the equations of motion for any celestial dynamic, without state transition matrix computation.
pub struct CelestialDynamics<'a> {
    pub state: State<Geoid>,
    pub bodies: Vec<i32>,
    // Loss in precision is avoided by using a relative time parameter initialized to zero
    relative_time: f64,
    // Allows us to rebuilt the true epoch
    init_tai_secs: f64,
    pub cosm: Option<&'a Cosm>,
}

impl<'a> CelestialDynamics<'a> {
    /// Initialize third body dynamics given the EXB IDs and a Cosm
    pub fn new(state: State<Geoid>, bodies: Vec<i32>, cosm: &'a Cosm) -> Self {
        for exb_id in &bodies {
            cosm.try_geoid_from_id(*exb_id)
                .expect("unknown EXB ID in list of third bodies");
        }
        Self {
            state,
            bodies,
            relative_time: 0.0,
            init_tai_secs: state.dt.as_tai_seconds(),
            cosm: Some(cosm),
        }
    }

    /// Initializes a CelestialDynamics which does not simulate the gravity pull of other celestial objects but the primary one.
    pub fn two_body(state: State<Geoid>) -> Self {
        Self {
            state,
            bodies: Vec::new(),
            relative_time: 0.0,
            init_tai_secs: state.dt.as_tai_seconds(),
            cosm: None,
        }
    }

    pub fn state_ctor(&self, rel_time: f64, state_vec: &Vector6<f64>) -> State<Geoid> {
        State::<Geoid>::from_cartesian(
            state_vec[0],
            state_vec[1],
            state_vec[2],
            state_vec[3],
            state_vec[4],
            state_vec[5],
            Epoch::from_tai_seconds(self.init_tai_secs + rel_time),
            self.state.frame,
        )
    }
}

impl<'a> Dynamics for CelestialDynamics<'a> {
    type StateSize = U6;
    type StateType = State<Geoid>;
    /// Returns the relative time to the propagator. Use prop.dynamics.state.dt for absolute time
    fn time(&self) -> f64 {
        self.relative_time
    }

    fn state_vector(&self) -> VectorN<f64, Self::StateSize> {
        self.state.to_cartesian_vec()
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        self.relative_time = new_t;
        self.state.dt = Epoch::from_tai_seconds(self.init_tai_secs + new_t);
        self.state.x = new_state[0];
        self.state.y = new_state[1];
        self.state.z = new_state[2];
        self.state.vx = new_state[3];
        self.state.vy = new_state[4];
        self.state.vz = new_state[5];
    }

    fn state(&self) -> State<Geoid> {
        self.state
    }

    fn eom(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let radius = state.fixed_rows::<U3>(0).into_owned();
        let velocity = state.fixed_rows::<U3>(3).into_owned();
        let body_acceleration = (-self.state.frame.gm / radius.norm().powi(3)) * radius;
        let mut d_x =
            Vector6::from_iterator(velocity.iter().chain(body_acceleration.iter()).cloned());

        // Get all of the position vectors between the center body and the third bodies
        let jde = Epoch::from_tai_seconds(self.init_tai_secs + t);
        for exb_id in &self.bodies {
            let third_body = self.cosm.unwrap().geoid_from_id(*exb_id);
            // State of j-th body as seen from primary body
            let st_ij = self
                .cosm
                .unwrap()
                .celestial_state(*exb_id, jde, self.state.frame.id);

            let r_ij = st_ij.radius();
            let r_ij3 = st_ij.rmag().powi(3);
            let r_j = radius - r_ij; // sc as seen from 3rd body
            let r_j3 = r_j.norm().powi(3);
            let third_body_acc = -third_body.gm * (r_j / r_j3 + r_ij / r_ij3);

            d_x[3] += third_body_acc[0];
            d_x[4] += third_body_acc[1];
            d_x[5] += third_body_acc[2];
        }

        d_x
    }
}

/// `CelestialDynamicsStm` provides the equations of motion for any celestial dynamic, **with** state transition matrix computation.
pub struct CelestialDynamicsStm<'a> {
    pub state: State<Geoid>,
    pub bodies: Vec<i32>,
    pub stm: Matrix6<f64>,
    // Loss in precision is avoided by using a relative time parameter initialized to zero
    relative_time: f64,
    // Allows us to rebuilt the true epoch
    init_tai_secs: f64,
    pub cosm: Option<&'a Cosm>,
}

impl<'a> CelestialDynamicsStm<'a> {
    /// Initialize third body dynamics given the EXB IDs and a Cosm
    pub fn new(state: State<Geoid>, bodies: Vec<i32>, cosm: &'a Cosm) -> Self {
        // Check that these bodies are present in the EXB.
        for exb_id in &bodies {
            cosm.try_geoid_from_id(*exb_id)
                .expect("unknown EXB ID in list of third bodies");
        }
        Self {
            state,
            bodies,
            stm: Matrix6::identity(),
            relative_time: 0.0,
            init_tai_secs: state.dt.as_tai_seconds(),
            cosm: Some(cosm),
        }
    }

    /// Initializes a CelestialDynamicsStm which does not simulate the gravity pull of other celestial objects but the primary one.
    pub fn two_body(state: State<Geoid>) -> Self {
        Self {
            state,
            bodies: Vec::new(),
            stm: Matrix6::identity(),
            relative_time: 0.0,
            init_tai_secs: state.dt.as_tai_seconds(),
            cosm: None,
        }
    }

    /// Provides a copy to the state.
    pub fn as_state(&self) -> State<Geoid> {
        self.state
    }

    /// Used only to set the orbital state, useful for Extended Kalman Filters.
    pub fn set_orbital_state(&mut self, new_t: f64, new_state: &Vector6<f64>) {
        self.relative_time = new_t;
        self.state.dt = Epoch::from_tai_seconds(self.init_tai_secs + new_t);
        self.state.x = new_state[0];
        self.state.y = new_state[1];
        self.state.z = new_state[2];
        self.state.vx = new_state[3];
        self.state.vy = new_state[4];
        self.state.vz = new_state[5];
    }
}

impl<'a> AutoDiffDynamics for CelestialDynamicsStm<'a> {
    type HyperStateSize = U7;
    type STMSize = U6;

    fn dual_eom(
        &self,
        t: f64,
        state: &VectorN<Hyperdual<f64, U7>, U6>,
    ) -> (Vector6<f64>, Matrix6<f64>) {
        // Extract data from hyperspace
        let radius = state.fixed_rows::<U3>(0).into_owned();
        let velocity = state.fixed_rows::<U3>(3).into_owned();

        // Code up math as usual
        let rmag = norm(&radius);
        let body_acceleration =
            radius * (Hyperdual::<f64, U7>::from_real(-self.state.frame.gm) / rmag.powi(3));

        // Extract result into Vector6 and Matrix6
        let mut fx = Vector6::zeros();
        let mut grad = Matrix6::zeros();
        for i in 0..U6::dim() {
            fx[i] = if i < 3 {
                velocity[i].real()
            } else {
                body_acceleration[i - 3].real()
            };
            for j in 1..U7::dim() {
                grad[(i, j - 1)] = if i < 3 {
                    velocity[i][j]
                } else {
                    body_acceleration[i - 3][j]
                };
            }
        }

        // Get all of the position vectors between the center body and the third bodies
        let jde = Epoch::from_tai_seconds(self.init_tai_secs + t);
        for exb_id in &self.bodies {
            let third_body = self.cosm.unwrap().geoid_from_id(*exb_id);
            let gm_d = Hyperdual::<f64, U7>::from_real(-third_body.gm);

            // State of j-th body as seen from primary body
            let st_ij = self
                .cosm
                .unwrap()
                .celestial_state(*exb_id, jde, self.state.frame.id);

            let r_ij: Vector3<Hyperdual<f64, U7>> = hyperspace_from_vector(&st_ij.radius());
            let r_ij3 = norm(&r_ij).powi(3) / gm_d;
            // The difference leads to the dual parts nulling themselves out, so let's fix that.
            let mut r_j = radius - r_ij; // sc as seen from 3rd body
            r_j[0][1] = 1.0;
            r_j[1][2] = 1.0;
            r_j[2][3] = 1.0;

            let r_j3 = norm(&r_j).powi(3) / gm_d;
            let third_body_acc_d = r_j / r_j3 + r_ij / r_ij3;

            for i in 0..U3::dim() {
                fx[i + 3] += third_body_acc_d[i][0];
                for j in 1..U7::dim() {
                    grad[(i + 3, j - 1)] += third_body_acc_d[i][j];
                }
            }
        }

        (fx, grad)
    }
}

impl<'a> Dynamics for CelestialDynamicsStm<'a> {
    type StateSize = U42;
    type StateType = (State<Geoid>, Matrix6<f64>);
    /// Returns the relative time to the propagator. Use prop.dynamics.state.dt for absolute time
    fn time(&self) -> f64 {
        self.relative_time
    }

    fn state_vector(&self) -> VectorN<f64, Self::StateSize> {
        let mut stm_as_vec = VectorN::<f64, U36>::zeros();
        let mut stm_idx = 0;
        for i in 0..6 {
            for j in 0..6 {
                stm_as_vec[(stm_idx, 0)] = self.stm[(i, j)];
                stm_idx += 1;
            }
        }
        VectorN::<f64, Self::StateSize>::from_iterator(
            self.state
                .to_cartesian_vec()
                .iter()
                .chain(stm_as_vec.iter())
                .cloned(),
        )
    }

    /// Returns the celestial state and the state transition matrix
    fn state(&self) -> Self::StateType {
        (self.state, self.stm)
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        self.relative_time = new_t;
        self.state.dt = Epoch::from_tai_seconds(self.init_tai_secs + new_t);
        self.state.x = new_state[0];
        self.state.y = new_state[1];
        self.state.z = new_state[2];
        self.state.vx = new_state[3];
        self.state.vy = new_state[4];
        self.state.vz = new_state[5];
        // And update the STM
        let mut stm_k_to_0 = Matrix6::zeros();
        let mut stm_idx = 6;
        for i in 0..6 {
            for j in 0..6 {
                stm_k_to_0[(i, j)] = new_state[(stm_idx, 0)];
                stm_idx += 1;
            }
        }

        let mut stm_prev = self.stm;
        if !stm_prev.try_inverse_mut() {
            println!("{}", self.stm);
            panic!("STM not invertible");
        }
        self.stm = stm_k_to_0 * stm_prev;
    }

    fn eom(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let pos_vel = state.fixed_rows::<U6>(0).into_owned();
        let (state, grad) = self.compute(t, &pos_vel);
        let stm_dt = self.stm * grad;
        // Rebuild the STM as a vector.
        let mut stm_as_vec = VectorN::<f64, U36>::zeros();
        let mut stm_idx = 0;
        for i in 0..6 {
            for j in 0..6 {
                stm_as_vec[(stm_idx, 0)] = stm_dt[(i, j)];
                stm_idx += 1;
            }
        }
        VectorN::<f64, Self::StateSize>::from_iterator(
            state.iter().chain(stm_as_vec.iter()).cloned(),
        )
    }
}
