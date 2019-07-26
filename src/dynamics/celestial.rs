extern crate hifitime;
extern crate hyperdual;
extern crate nalgebra as na;

use self::hifitime::instant::{Era, Instant};
use self::hifitime::julian::ModifiedJulian;
use self::hifitime::TimeSystem;
use self::hifitime::SECONDS_PER_DAY;
use self::hyperdual::linalg::norm;
use self::hyperdual::{Float, Hyperdual};
use self::na::{DimName, Matrix6, MatrixMN, Vector6, VectorN, U3, U36, U42, U6, U7};
use super::Dynamics;
use celestia::{Cosm, Geoid, State};
use od::{AutoDiffDynamics, Linearization};
use std::f64;

/// `CelestialDynamics` provides the equations of motion for any celestial dynamic, without state transition matrix computation.
pub struct CelestialDynamics<'a> {
    pub state: State<Geoid>,
    pub bodies: Vec<i32>,
    init_mjd: ModifiedJulian,
    time: f64,
    cosm: Option<&'a Cosm>,
}

impl<'a> CelestialDynamics<'a> {
    /// Initialize third body dynamics given the EXB IDs and a Cosm
    pub fn new(state: State<Geoid>, bodies: Vec<i32>, cosm: &'a Cosm) -> Self {
        Self {
            state,
            bodies,
            init_mjd: state.dt_as_modified_julian(),
            time: 0.0,
            cosm: Some(cosm),
        }
    }

    /// Initializes a CelestialDynamics which does not simulate the gravity pull of other celestial objects but the primary one.
    pub fn two_body(state: State<Geoid>) -> Self {
        Self {
            state,
            bodies: Vec::new(),
            init_mjd: state.dt_as_modified_julian(),
            time: 0.0,
            cosm: None,
        }
    }

    /// Provides a copy to the state.
    pub fn as_state(&self) -> State<Geoid> {
        self.state
    }
}

impl<'a> Dynamics for CelestialDynamics<'a> {
    type StateSize = U6;
    /// Returns modified Julian seconds since 01 Jan 2000 00:00
    fn time(&self) -> f64 {
        self.time
    }

    fn state(&self) -> VectorN<f64, Self::StateSize> {
        self.state.to_cartesian_vec()
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        self.time = new_t;
        self.state.x = new_state[0];
        self.state.y = new_state[1];
        self.state.z = new_state[2];
        self.state.vx = new_state[3];
        self.state.vy = new_state[4];
        self.state.vz = new_state[5];
        let mjd = ModifiedJulian {
            days: self.init_mjd.days + self.time / SECONDS_PER_DAY,
        };
        debug!("set state at MJD: {}", mjd);
        self.state.dt = mjd.into_instant();
    }

    fn eom(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let radius = state.fixed_rows::<U3>(0).into_owned();
        let velocity = state.fixed_rows::<U3>(3).into_owned();
        let body_acceleration = (-self.state.frame.gm / radius.norm().powi(3)) * radius;
        let mut dX = Vector6::from_iterator(velocity.iter().chain(body_acceleration.iter()).cloned());

        // Get all of the position vectors between the center body and the third bodies
        let jde = ModifiedJulian {
            days: self.init_mjd.days + t / SECONDS_PER_DAY,
        }
        .julian_days();
        for exb_id in &self.bodies {
            let third_body = self
                .cosm
                .unwrap()
                .geoid_from_id(*exb_id)
                .expect("unknown EXB ID in list of third bodies");
            let st_ij = self.cosm.unwrap().celestial_state(self.state.frame.id, jde, *exb_id).unwrap(); // frame center to 3rd body
            let r_ij = st_ij.radius();
            let r_ij3 = st_ij.rmag().powi(3);
            let r_j = radius - r_ij; // sc to 3rd body
            let r_j3 = r_j.norm().powi(3);
            let third_body_acc = -third_body.gm * (r_j / r_j3 + r_ij / r_ij3);
            dX[3] += third_body_acc[0];
            dX[4] += third_body_acc[1];
            dX[5] += third_body_acc[2];
        }

        dX
    }
}

/// `TwoBodyWithStm` exposes the equations of motion for a simple two body propagation. It inherently supports
/// the State Transition Matrix for orbital determination.
pub struct TwoBodyWithStm<'a> {
    pub stm: Matrix6<f64>,
    pub two_body_dyn: CelestialDynamics<'a>,
    pub geoid: Geoid,
    pub initial_instant: Instant,
}

impl<'a> TwoBodyWithStm<'a> {
    /// Initialize TwoBody dynamics around a provided `CelestialBody` from the provided position and velocity state (cf. nyx::celestia).
    pub fn from_state(state: State<Geoid>) -> TwoBodyWithStm<'a> {
        TwoBodyWithStm {
            stm: Matrix6::identity(),
            two_body_dyn: CelestialDynamics::two_body(state),
            geoid: state.frame,
            initial_instant: state.dt,
        }
    }

    pub fn to_state(&self) -> State<Geoid> {
        // Compute the new time
        let instant_in_secs = self.initial_instant.secs() as f64 + (f64::from(self.initial_instant.nanos())) * 1e-9;
        State::<Geoid>::from_cartesian_vec(
            &self.two_body_dyn.state(),
            ModifiedJulian::from_instant(Instant::from_precise_seconds(instant_in_secs, Era::Present)),
            self.geoid,
        )
    }
}

impl<'a> Dynamics for TwoBodyWithStm<'a> {
    type StateSize = U42;
    fn time(&self) -> f64 {
        self.two_body_dyn.time()
    }

    fn state(&self) -> VectorN<f64, Self::StateSize> {
        let mut stm_as_vec = VectorN::<f64, U36>::zeros();
        let mut stm_idx = 0;
        for i in 0..6 {
            for j in 0..6 {
                stm_as_vec[(stm_idx, 0)] = self.stm[(i, j)];
                stm_idx += 1;
            }
        }
        VectorN::<f64, Self::StateSize>::from_iterator(self.two_body_dyn.state().iter().chain(stm_as_vec.iter()).cloned())
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        let pos_vel = new_state.fixed_rows::<U6>(0).into_owned();
        self.two_body_dyn.set_state(new_t, &pos_vel);
        let mut stm_k_to_0 = Matrix6::zeros();
        let mut stm_idx = 6;
        for i in 0..6 {
            for j in 0..6 {
                stm_k_to_0[(i, j)] = new_state[stm_idx];
                stm_idx += 1;
            }
        }

        let mut stm_prev = self.stm;
        if !stm_prev.try_inverse_mut() {
            panic!("STM not invertible");
        }
        self.stm = stm_k_to_0 * stm_prev;
    }

    fn eom(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let pos_vel = state.fixed_rows::<U6>(0).into_owned();
        let two_body_dt = self.two_body_dyn.eom(t, &pos_vel);
        let stm_dt = self.stm * self.gradient(t, &pos_vel);
        // Rebuild the STM as a vector.
        let mut stm_as_vec = VectorN::<f64, U36>::zeros();
        let mut stm_idx = 0;
        for i in 0..6 {
            for j in 0..6 {
                stm_as_vec[(stm_idx, 0)] = stm_dt[(i, j)];
                stm_idx += 1;
            }
        }
        VectorN::<f64, Self::StateSize>::from_iterator(two_body_dt.iter().chain(stm_as_vec.iter()).cloned())
    }
}

impl<'a> Linearization for TwoBodyWithStm<'a> {
    type StateSize = U6;

    fn gradient(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> MatrixMN<f64, Self::StateSize, Self::StateSize> {
        let mut grad = MatrixMN::<f64, U6, U6>::zeros();
        // Top right is Identity 3x3
        grad[(0, 3)] = 1.0;
        grad[(1, 4)] = 1.0;
        grad[(2, 5)] = 1.0;
        // Bottom left is where the magic happens.
        let x = state[(0, 0)];
        let y = state[(1, 0)];
        let z = state[(2, 0)];
        let x2 = x.powi(2);
        let y2 = y.powi(2);
        let z2 = z.powi(2);
        let r232 = (x2 + y2 + z2).powf(3.0 / 2.0);
        let r252 = (x2 + y2 + z2).powf(5.0 / 2.0);

        // Add the body perturbations
        let dax_dx = 3.0 * self.geoid.gm * x2 / r252 - self.geoid.gm / r232;
        let dax_dy = 3.0 * self.geoid.gm * x * y / r252;
        let dax_dz = 3.0 * self.geoid.gm * x * z / r252;
        let day_dx = 3.0 * self.geoid.gm * x * y / r252;
        let day_dy = 3.0 * self.geoid.gm * y2 / r252 - self.geoid.gm / r232;
        let day_dz = 3.0 * self.geoid.gm * y * z / r252;
        let daz_dx = 3.0 * self.geoid.gm * x * z / r252;
        let daz_dy = 3.0 * self.geoid.gm * y * z / r252;
        let daz_dz = 3.0 * self.geoid.gm * z2 / r252 - self.geoid.gm / r232;

        // Let the gradient
        grad[(3, 0)] = dax_dx;
        grad[(4, 0)] = day_dx;
        grad[(5, 0)] = daz_dx;
        grad[(3, 1)] = dax_dy;
        grad[(4, 1)] = day_dy;
        grad[(5, 1)] = daz_dy;
        grad[(3, 2)] = dax_dz;
        grad[(4, 2)] = day_dz;
        grad[(5, 2)] = daz_dz;

        grad
    }
}

#[derive(Copy, Clone)]
pub struct TwoBodyWithDualStm {
    pub mu: f64,
    pub stm: Matrix6<f64>,
    pub pos_vel: Vector6<f64>,
    time: f64,
    geoid: Geoid,
    initial_instant: Instant,
}

impl TwoBodyWithDualStm {
    /// Initialize TwoBody dynamics around a provided `CelestialBody` from the provided position and velocity state (cf. nyx::celestia).
    pub fn from_state(state: &State<Geoid>) -> TwoBodyWithDualStm {
        TwoBodyWithDualStm {
            mu: state.frame.gm,
            stm: Matrix6::identity(),
            pos_vel: state.to_cartesian_vec(),
            time: 0.0,
            geoid: state.frame,
            initial_instant: state.dt,
        }
    }

    pub fn to_state(&self) -> State<Geoid> {
        // Compute the new time
        let instant_in_secs = self.initial_instant.secs() as f64 + (f64::from(self.initial_instant.nanos())) * 1e-9;
        State::<Geoid>::from_cartesian_vec(
            &self.pos_vel,
            ModifiedJulian::from_instant(Instant::from_precise_seconds(instant_in_secs, Era::Present)),
            self.geoid,
        )
    }
}

impl AutoDiffDynamics for TwoBodyWithDualStm {
    type HyperStateSize = U7;
    type STMSize = U6;

    fn dual_eom(&self, _t: f64, state: &VectorN<Hyperdual<f64, U7>, U6>) -> (Vector6<f64>, Matrix6<f64>) {
        // Extract data from hyperspace
        let radius = state.fixed_rows::<U3>(0).into_owned();
        let velocity = state.fixed_rows::<U3>(3).into_owned();

        // Code up math as usual
        let rmag = norm(&radius);
        let body_acceleration = radius * (Hyperdual::<f64, U7>::from_real(-self.mu) / rmag.powi(3));

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
                grad[(i, j - 1)] = if i < 3 { velocity[i][j] } else { body_acceleration[i - 3][j] };
            }
        }

        (fx, grad)
    }
}

impl Dynamics for TwoBodyWithDualStm {
    type StateSize = U42;
    fn time(&self) -> f64 {
        self.time
    }

    fn state(&self) -> VectorN<f64, Self::StateSize> {
        let mut stm_as_vec = VectorN::<f64, U36>::zeros();
        let mut stm_idx = 0;
        for i in 0..6 {
            for j in 0..6 {
                stm_as_vec[(stm_idx, 0)] = self.stm[(i, j)];
                stm_idx += 1;
            }
        }
        VectorN::<f64, Self::StateSize>::from_iterator(self.pos_vel.iter().chain(stm_as_vec.iter()).cloned())
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        self.time = new_t;
        self.pos_vel = new_state.fixed_rows::<U6>(0).into_owned();

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
        VectorN::<f64, Self::StateSize>::from_iterator(state.iter().chain(stm_as_vec.iter()).cloned())
    }
}
