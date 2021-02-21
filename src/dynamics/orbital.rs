/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::hyperdual::linalg::norm;
use super::hyperdual::{extract_jacobian_and_result, hyperspace_from_vector, Float, Hyperdual};
use super::{AccelModel, Dynamics, NyxError};
use crate::celestia::{Bodies, Cosm, Frame, LTCorr, Orbit};
use crate::dimensions::{
    DimName, Matrix3, Matrix6, Vector3, Vector6, VectorN, U3, U36, U4, U42, U6, U7,
};
use crate::{State, TimeTagged};
// use od::Estimable;
use std::f64;
use std::sync::Arc;

pub use super::sph_harmonics::Harmonics;

/// `OrbitalDynamics` provides the equations of motion for any celestial dynamic, without state transition matrix computation.
#[derive(Clone)]
pub struct OrbitalDynamics<'a> {
    pub accel_models: Vec<Arc<dyn AccelModel + Sync + 'a>>,
}

impl<'a> OrbitalDynamics<'a> {
    /// Initialize point mass dynamics given the EXB IDs and a Cosm
    pub fn point_masses(integr_frame: Frame, bodies: &[Bodies], cosm: Arc<Cosm>) -> Arc<Self> {
        // Create the point masses
        Self::new(vec![PointMasses::new(integr_frame, bodies, cosm)])
    }

    /// Initializes a OrbitalDynamics which does not simulate the gravity pull of other celestial objects but the primary one.
    pub fn two_body() -> Arc<Self> {
        Self::new(vec![])
    }

    /// Initialize orbital dynamics with a list of acceleration models
    pub fn new(accel_models: Vec<Arc<dyn AccelModel + Sync + 'a>>) -> Arc<Self> {
        Arc::new(Self::new_raw(accel_models))
    }

    /// Initialize orbital dynamics with a list of acceleration models, _without_ encapsulating it in an Arc
    /// Use this only if you need to mutate the dynamics as you'll need to wrap it in an Arc before propagation.
    pub fn new_raw(accel_models: Vec<Arc<dyn AccelModel + Sync + 'a>>) -> Self {
        Self { accel_models }
    }

    /// Initialize new orbital mechanics with the provided model.
    /// **Note:** Orbital dynamics _always_ include two body dynamics, these cannot be turned off.
    pub fn with_model(accel_model: Arc<dyn AccelModel + Sync + 'a>) -> Arc<Self> {
        Self::new(vec![accel_model])
    }

    /// Add a model to the currently defined orbital dynamics
    pub fn add_model(&mut self, accel_model: Arc<dyn AccelModel + Sync + 'a>) {
        self.accel_models.push(accel_model);
    }
}

impl<'a> Dynamics for OrbitalDynamics<'a> {
    type HyperdualSize = U7;
    type StateType = Orbit;

    fn eom(
        &self,
        delta_t_s: f64,
        state: &VectorN<f64, U42>,
        ctx: &Orbit,
    ) -> Result<VectorN<f64, U42>, NyxError> {
        let (new_state, new_stm) = if ctx.stm.is_some() {
            // Then call the dual_eom with the correct state size
            let pos_vel = state.fixed_rows::<U6>(0).into_owned();
            let (state, grad) = self.eom_grad(delta_t_s, &pos_vel, ctx)?;
            let stm_dt = ctx.stm() * grad;
            // Rebuild the STM as a vector.
            let mut stm_as_vec = VectorN::<f64, U36>::zeros();
            let mut stm_idx = 0;
            for i in 0..U6::dim() {
                for j in 0..U6::dim() {
                    stm_as_vec[(stm_idx, 0)] = stm_dt[(i, j)];
                    stm_idx += 1;
                }
            }
            (state, stm_as_vec)
        } else {
            // Still return something of size 42, but the STM will be zeros.

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
                let model_acc = model.eom(&osc)?;
                for i in 0..3 {
                    d_x[i + 3] += model_acc[i];
                }
            }

            (d_x, VectorN::<f64, U36>::zeros())
        };
        Ok(VectorN::<f64, U42>::from_iterator(
            new_state.iter().chain(new_stm.iter()).cloned(),
        ))
    }

    fn dual_eom(
        &self,
        _delta_t_s: f64,
        state: &VectorN<Hyperdual<f64, U7>, U6>,
        ctx: &Orbit,
    ) -> Result<(Vector6<f64>, Matrix6<f64>), NyxError> {
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
            let (model_acc, model_grad) = model.dual_eom(&radius, ctx)?;
            for i in 0..U3::dim() {
                fx[i + 3] += model_acc[i];
                for j in 1..U4::dim() {
                    grad[(i + 3, j - 1)] += model_grad[(i, j - 1)];
                }
            }
        }

        Ok((fx, grad))
    }
}

/// Stores the reference to the third body, including the gm to avoid having to fetch it every time
pub struct ThirdBodyRef {
    ephem: Vec<usize>,
    gm: f64,
}

/// PointMasses model
pub struct PointMasses {
    /// The propagation frame
    pub frame: Frame,
    pub bodies: Vec<ThirdBodyRef>,
    /// Optional point to a Cosm, needed if extra point masses are needed
    pub cosm: Arc<Cosm>,
    /// Light-time correction computation if extra point masses are needed
    pub correction: LTCorr,
}

impl PointMasses {
    /// Initializes the multibody point mass dynamics with the provided list of bodies
    pub fn new(propagation_frame: Frame, body_names: &[Bodies], cosm: Arc<Cosm>) -> Arc<Self> {
        Arc::new(Self::with_correction(
            propagation_frame,
            body_names,
            cosm,
            LTCorr::None,
        ))
    }

    /// Initializes the multibody point mass dynamics with the provided list of bodies, and accounting for some light time correction
    pub fn with_correction(
        propagation_frame: Frame,
        bodies: &[Bodies],
        cosm: Arc<Cosm>,
        correction: LTCorr,
    ) -> Self {
        let mut refs = Vec::new();
        // Check that these celestial bodies exist and build their references
        for body in bodies {
            refs.push(ThirdBodyRef {
                ephem: body.ephem_path().to_vec(),
                gm: cosm.frame_from_ephem_path(&body.ephem_path()).gm(),
            });
        }

        Self {
            frame: propagation_frame,
            bodies: refs,
            cosm,
            correction,
        }
    }

    /// Allows using bodies by name, defined in the non-default XB
    pub fn specific(
        propagation_frame: Frame,
        bodies: &[String],
        cosm: Arc<Cosm>,
        correction: LTCorr,
    ) -> Self {
        let mut refs = Vec::new();
        // Check that these celestial bodies exist and build their references
        for body in bodies {
            let path = cosm.xb.ephemeris_find_path(body.to_string()).unwrap();
            refs.push(ThirdBodyRef {
                ephem: path.to_vec(),
                gm: cosm.frame_from_ephem_path(&path).gm(),
            });
        }

        Self {
            frame: propagation_frame,
            bodies: refs,
            cosm,
            correction,
        }
    }
}

impl AccelModel for PointMasses {
    fn eom(&self, osc: &Orbit) -> Result<Vector3<f64>, NyxError> {
        let mut d_x = Vector3::zeros();
        // Get all of the position vectors between the center body and the third bodies
        for third_body in &self.bodies {
            // Orbit of j-th body as seen from primary body
            let st_ij =
                self.cosm
                    .celestial_state(&third_body.ephem, osc.dt, self.frame, self.correction);

            let r_ij = st_ij.radius();
            let r_ij3 = st_ij.rmag().powi(3);
            let r_j = osc.radius() - r_ij; // sc as seen from 3rd body
            let r_j3 = r_j.norm().powi(3);
            d_x += -third_body.gm * (r_j / r_j3 + r_ij / r_ij3);
        }
        Ok(d_x)
    }

    fn dual_eom(
        &self,
        state: &VectorN<Hyperdual<f64, U7>, U3>,
        ctx: &Orbit,
    ) -> Result<(Vector3<f64>, Matrix3<f64>), NyxError> {
        // Extract data from hyperspace
        let radius = state.fixed_rows::<U3>(0).into_owned();
        // Extract result into Vector6 and Matrix6
        let mut fx = Vector3::zeros();
        let mut grad = Matrix3::zeros();

        // Get all of the position vectors between the center body and the third bodies
        for third_body in &self.bodies {
            let gm_d = Hyperdual::<f64, U7>::from_real(-third_body.gm);

            // Orbit of j-th body as seen from primary body
            let st_ij = self.cosm.try_celestial_state(
                &third_body.ephem,
                ctx.epoch(),
                self.frame,
                self.correction,
            )?;

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

        Ok((fx, grad))
    }
}
