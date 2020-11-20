use super::hyperdual::linalg::norm;
use super::hyperdual::{extract_jacobian_and_result, hyperspace_from_vector, Float, Hyperdual};
use super::{AccelModel, AutoDiff, Dynamics, NyxError};
use crate::celestia::{Cosm, Frame, LTCorr, Orbit};
use crate::dimensions::{
    allocator::Allocator, DefaultAllocator, DimName, Matrix3, Matrix6, Vector3, Vector6, VectorN,
    U3, U36, U4, U42, U6, U7,
};
use crate::time::Epoch;
use crate::{State, TimeTagged};
use od::Estimable;
use std::f64;

pub use super::sph_harmonics::{Harmonics, HarmonicsDiff};

// TODO: I don't think I need this anymore, I think I can simply initialize the Dynamics as a <T1, T2>
// like a Result since it doesn't store the state anymore
pub trait OrbitalDynamicsT: Dynamics<StateType = Orbit> {}
// D: Estimable<MsrIn, LinStateSize = Msr::StateSize>,
// pub trait OrbitalDynamicsT: Dynamics {
//     fn stm(&self) -> Option<Matrix6<f64>>;

//     fn orbital_state_ctor(&self, rel_time: f64, state_vec: &VectorN<f64, Self::StateSize>) -> Orbit
//     where
//         DefaultAllocator: Allocator<f64, Self::StateSize>;
// }

/// `OrbitalDynamics` provides the equations of motion for any celestial dynamic, without state transition matrix computation.
pub struct OrbitalDynamics<'a> {
    /// Current state of these dynamics
    // pub integr_frame: Frame,
    // /// Loss in precision is avoided by using a relative time parameter initialized to zero
    // relative_time: f64,
    // /// Allows us to rebuilt the true epoch
    // init_tai_secs: f64,
    pub accel_models: Vec<Box<dyn AccelModel + Sync + 'a>>,
}

impl<'a> OrbitalDynamicsT for OrbitalDynamics<'a> {
    //     fn stm(&self) -> Option<Matrix6<f64>> {
    //         None
    //     }

    //     fn orbital_state_ctor(
    //         &self,
    //         rel_time: f64,
    //         state_vec: &VectorN<f64, Self::StateSize>,
    //     ) -> Orbit {
    //         self.state_ctor(rel_time, state_vec)
    //     }
}

impl<'a> OrbitalDynamics<'a> {
    /// Initialize point mass dynamics given the EXB IDs and a Cosm
    pub fn point_masses(integr_frame: Frame, bodies: &[i32], cosm: &'a Cosm) -> Self {
        // Create the point masses
        let pts = PointMasses::new(integr_frame, bodies, cosm);
        Self {
            // integr_frame,
            // relative_time: 0.0,
            // init_tai_secs: state.dt.as_tai_seconds(),
            accel_models: vec![Box::new(pts)],
        }
    }

    /// Initializes a OrbitalDynamics which does not simulate the gravity pull of other celestial objects but the primary one.
    pub fn two_body(state: Orbit) -> Self {
        Self {
            // state,
            // relative_time: 0.0,
            // init_tai_secs: state.dt.as_tai_seconds(),
            accel_models: Vec::new(),
        }
    }

    /// Initialize orbital dynamics with a list of acceleration models
    pub fn new(state: Orbit, accel_models: Vec<Box<dyn AccelModel + Sync + 'a>>) -> Self {
        Self {
            // state,
            // relative_time: 0.0,
            // init_tai_secs: state.dt.as_tai_seconds(),
            accel_models,
        }
    }

    pub fn add_model(&mut self, accel_model: Box<dyn AccelModel + Sync + 'a>) {
        self.accel_models.push(accel_model);
    }
}

impl<'a> Dynamics for OrbitalDynamics<'a> {
    type StateSize = U6;

    fn eom(
        &self,
        delta_t_s: f64,
        state: &VectorN<f64, Self::StateSize>,
        ctx: &Orbit,
    ) -> VectorN<f64, Self::StateSize> {
        let osc = ctx.ctor_from(delta_t_s, state);
        // TODO: Speed check this with the PointMasses only, including the integration frame point mass
        let body_acceleration = (-osc.frame.gm() / osc.rmag().powi(3)) * osc.radius();
        let mut d_x = Vector6::from_iterator(
            osc.velocity()
                .iter()
                .chain(body_acceleration.iter())
                .cloned(),
        );

        // Apply the acceleration models
        for model in &self.accel_models {
            let model_acc = model.eom(&osc);
            for i in 0..3 {
                d_x[i + 3] += model_acc[i];
            }
        }

        d_x
    }
}

/// `OrbitalDynamicsStm` provides the equations of motion for any celestial dynamic, **with** state transition matrix computation.
pub struct OrbitalDynamicsStm<'a> {
    // pub state: Orbit,
    // // Loss in precision is avoided by using a relative time parameter initialized to zero
    // relative_time: f64,
    // // Allows us to rebuilt the true epoch
    // init_tai_secs: f64,
    pub accel_models: Vec<Box<dyn AutoDiff<U7, U3, CtxType = Orbit> + 'a>>,
}

impl<'a> OrbitalDynamicsT for OrbitalDynamicsStm<'a> {
    //     fn orbital_state(&self) -> Orbit {
    //         self.state
    //     }

    //     fn stm(&self) -> Option<Matrix6<f64>> {
    //         self.state.stm
    //     }

    //     fn orbital_state_ctor(
    //         &self,
    //         rel_time: f64,
    //         state_vec: &VectorN<f64, Self::StateSize>,
    //     ) -> Orbit {
    //         self.state_ctor(rel_time, state_vec).0
    //     }
}

impl<'a> OrbitalDynamicsStm<'a> {
    /// Initialize third body dynamics given the EXB IDs and a Cosm
    pub fn point_masses(state: Orbit, bodies: &[i32], cosm: &'a Cosm) -> Self {
        let pts = PointMasses::new(state.frame, bodies, cosm);
        let mut state = state;
        state.stm_identity();
        Self {
            // state,
            // relative_time: 0.0,
            // init_tai_secs: state.dt.as_tai_seconds(),
            accel_models: vec![Box::new(pts)],
        }
    }

    /// Initializes a OrbitalDynamicsStm which does not simulate the gravity pull of other celestial objects but the primary one.
    pub fn two_body(state: Orbit) -> Self {
        let mut state = state;
        state.stm_identity();
        Self {
            // state,
            // relative_time: 0.0,
            // init_tai_secs: state.dt.as_tai_seconds(),
            accel_models: vec![],
        }
    }

    /// Add a model to these celestial dynamics (must be differentiable by automatic differentiation)
    pub fn add_model(&mut self, accel_model: Box<dyn AutoDiff<U7, U3, CtxType = Orbit> + 'a>) {
        self.accel_models.push(accel_model);
    }

    // /// Provides a copy to the state.
    // pub fn as_state(&self) -> Orbit {
    //     self.state
    // }

    // TODO: Figure this out
    // / Used only to set the orbital state, useful for Extended Kalman Filters.
    // pub fn set_orbital_state(&mut self, new_t: f64, new_state: &Vector6<f64>) {
    //     self.relative_time = new_t;
    //     self.state.dt = Epoch::from_tai_seconds(self.init_tai_secs + new_t);
    //     self.state.x = new_state[0];
    //     self.state.y = new_state[1];
    //     self.state.z = new_state[2];
    //     self.state.vx = new_state[3];
    //     self.state.vy = new_state[4];
    //     self.state.vz = new_state[5];
    // }

    // / Rebuild the state and STM from the provided vector
    // pub fn state_ctor(&self, t: f64, in_state: &VectorN<f64, U42>) -> (Orbit, Matrix6<f64>) {
    //     // Copy the current state, and then modify it
    //     let mut state = self.state();
    //     state.dt = Epoch::from_tai_seconds(self.init_tai_secs + t);
    //     state.x = in_state[0];
    //     state.y = in_state[1];
    //     state.z = in_state[2];
    //     state.vx = in_state[3];
    //     state.vy = in_state[4];
    //     state.vz = in_state[5];

    //     let mut stm = state.stm();
    //     for i in 0..6 {
    //         for j in 0..6 {
    //             stm[(i, j)] = in_state[i * 6 + j + 6];
    //         }
    //     }

    //     (state, stm)
    // }
}

impl<'a> AutoDiff<U7, U6> for OrbitalDynamicsStm<'a> {
    // type STMRows = U6;
    // type HyperStateSize = U7;
    type CtxType = Orbit;

    fn dual_eom(
        &self,
        // epoch: Epoch,
        // integr_frame: Frame,
        state: &VectorN<Hyperdual<f64, U7>, U6>,
        ctx: &Orbit,
    ) -> (Vector6<f64>, Matrix6<f64>) {
        // Extract data from hyperspace
        let radius = state.fixed_rows::<U3>(0).into_owned();
        let velocity = state.fixed_rows::<U3>(3).into_owned();

        // Code up math as usual
        let rmag = norm(&radius);
        let body_acceleration =
            radius * (Hyperdual::<f64, U7>::from_real(-ctx.frame.gm()) / rmag.powi(3));

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

        // Apply the acceleration models
        for model in &self.accel_models {
            let (model_acc, model_grad) = model.dual_eom(&radius, ctx);
            for i in 0..U3::dim() {
                fx[i + 3] += model_acc[i];
                for j in 1..U4::dim() {
                    grad[(i + 3, j - 1)] += model_grad[(i, j - 1)];
                }
            }
        }

        (fx, grad)
    }
}

impl<'a> Dynamics for OrbitalDynamicsStm<'a> {
    type StateSize = U42;

    // fn state_vector(&self) -> VectorN<f64, Self::StateSize> {
    //     let mut stm_as_vec = VectorN::<f64, U36>::zeros();
    //     let mut stm_idx = 0;
    //     let stm = self.state.stm();
    //     for i in 0..6 {
    //         for j in 0..6 {
    //             stm_as_vec[(stm_idx, 0)] = stm[(i, j)];
    //             stm_idx += 1;
    //         }
    //     }
    //     VectorN::<f64, Self::StateSize>::from_iterator(
    //         self.state
    //             .to_cartesian_vec()
    //             .iter()
    //             .chain(stm_as_vec.iter())
    //             .cloned(),
    //     )
    // }

    // fn set_state(
    //     &mut self,
    //     new_t: f64,
    //     new_state: &VectorN<f64, Self::StateSize>,
    // ) -> Result<(), NyxError> {
    //     self.relative_time = new_t;
    //     self.state.dt = Epoch::from_tai_seconds(self.init_tai_secs + new_t);
    //     self.state.x = new_state[0];
    //     self.state.y = new_state[1];
    //     self.state.z = new_state[2];
    //     self.state.vx = new_state[3];
    //     self.state.vy = new_state[4];
    //     self.state.vz = new_state[5];
    //     // And update the STM
    //     let mut stm_k_to_0 = Matrix6::zeros();
    //     let mut stm_idx = 6;
    //     for i in 0..6 {
    //         for j in 0..6 {
    //             stm_k_to_0[(i, j)] = new_state[(stm_idx, 0)];
    //             stm_idx += 1;
    //         }
    //     }

    //     let mut stm_prev = self.state.stm();
    //     if !stm_prev.try_inverse_mut() {
    //         error!("STM not invertible: {}", stm_prev);
    //         return Err(NyxError::SingularStateTransitionMatrix);
    //     }
    //     self.state.stm = Some(stm_k_to_0 * stm_prev);
    //     // self.state.stm = Some(stm_k_to_0);
    //     Ok(())
    // }

    fn eom(
        &self,
        delta_t_s: f64,
        state: &VectorN<f64, Self::StateSize>,
        ctx: &Orbit,
    ) -> VectorN<f64, Self::StateSize> {
        let pos_vel = state.fixed_rows::<U6>(0).into_owned();
        let epoch = Epoch::from_tai_seconds(self.init_tai_secs + t);
        let (state, grad) = self.eom_grad(epoch, self.state.frame, &pos_vel);
        let stm_dt = self.state.stm() * grad;
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

impl<'a> Estimable<Orbit> for OrbitalDynamicsStm<'a> {
    type LinStateSize = U6;

    fn to_measurement(&self, prop_state: &Orbit) -> Orbit {
        *prop_state
    }

    fn extract_stm(&self, prop_state: &Orbit) -> Matrix6<f64> {
        prop_state.stm()
    }

    fn extract_estimated_state(&self, prop_state: &Orbit) -> VectorN<f64, Self::LinStateSize> {
        prop_state.to_cartesian_vec()
    }

    /// Returns the estimated state
    fn set_estimated_state(&mut self, new_state: VectorN<f64, Self::LinStateSize>) {
        self.state.x = new_state[0];
        self.state.y = new_state[1];
        self.state.z = new_state[2];
        self.state.vx = new_state[3];
        self.state.vy = new_state[4];
        self.state.vz = new_state[5];
    }
}

/// PointMasses model
pub struct PointMasses<'a> {
    /// The propagation frame
    pub frame: Frame,
    pub bodies: Vec<i32>,
    /// Optional point to a Cosm, needed if extra point masses are needed
    pub cosm: &'a Cosm,
    /// Light-time correction computation if extra point masses are needed
    pub correction: LTCorr,
}

impl<'a> PointMasses<'a> {
    pub fn new(propagation_frame: Frame, bodies: &[i32], cosm: &'a Cosm) -> Self {
        Self::with_correction(propagation_frame, bodies.to_vec(), cosm, LTCorr::None)
    }

    pub fn with_correction(
        propagation_frame: Frame,
        bodies: Vec<i32>,
        cosm: &'a Cosm,
        correction: LTCorr,
    ) -> Self {
        // Check that these celestial bodies exist
        for exb_id in &bodies {
            cosm.try_frame_by_exb_id(*exb_id)
                .expect("unknown EXB ID in list of third bodies");
        }

        Self {
            frame: propagation_frame,
            bodies,
            cosm,
            correction,
        }
    }
}

impl<'a> AccelModel for PointMasses<'a> {
    fn eom(&self, osc: &Orbit) -> Vector3<f64> {
        let mut d_x = Vector3::zeros();
        // Get all of the position vectors between the center body and the third bodies
        for exb_id in &self.bodies {
            let third_body = self.cosm.frame_by_exb_id(*exb_id);
            // Orbit of j-th body as seen from primary body
            let st_ij = self
                .cosm
                .celestial_state(*exb_id, osc.dt, self.frame, self.correction);

            let r_ij = st_ij.radius();
            let r_ij3 = st_ij.rmag().powi(3);
            let r_j = osc.radius() - r_ij; // sc as seen from 3rd body
            let r_j3 = r_j.norm().powi(3);
            d_x += -third_body.gm() * (r_j / r_j3 + r_ij / r_ij3);
        }
        d_x
    }
}

impl<'a> AutoDiff<U7, U3> for PointMasses<'a> {
    // type HyperStateSize = U7;
    // type STMRows = U3;
    type CtxType = Orbit;

    fn dual_eom(
        &self,
        state: &VectorN<Hyperdual<f64, U7>, U3>,
        ctx: &Orbit,
    ) -> (Vector3<f64>, Matrix3<f64>) {
        // Extract data from hyperspace
        let radius = state.fixed_rows::<U3>(0).into_owned();
        // Extract result into Vector6 and Matrix6
        let mut fx = Vector3::zeros();
        let mut grad = Matrix3::zeros();

        // Get all of the position vectors between the center body and the third bodies
        for exb_id in &self.bodies {
            let third_body = self.cosm.frame_by_exb_id(*exb_id);
            let gm_d = Hyperdual::<f64, U7>::from_real(-third_body.gm());

            // Orbit of j-th body as seen from primary body
            let st_ij =
                self.cosm
                    .celestial_state(*exb_id, ctx.epoch(), self.frame, self.correction);

            let r_ij: Vector3<Hyperdual<f64, U7>> = hyperspace_from_vector(&st_ij.radius());
            let r_ij3 = norm(&r_ij).powi(3) / gm_d;
            // The difference leads to the dual parts nulling themselves out, so let's fix that.
            let mut r_j = radius - r_ij; // sc as seen from 3rd body
            r_j[0][1] = 1.0;
            r_j[1][2] = 1.0;
            r_j[2][3] = 1.0;

            let r_j3 = norm(&r_j).powi(3) / gm_d;
            let third_body_acc_d = r_j / r_j3 + r_ij / r_ij3;

            let (fxp, gradp) = extract_jacobian_and_result::<_, U3, U3, _>(&third_body_acc_d);
            fx += fxp;
            grad += gradp;
        }

        (fx, grad)
    }
}
