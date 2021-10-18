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
use crate::cosmic::{Bodies, Cosm, Frame, LightTimeCalc, Orbit};
use crate::linalg::{Const, Matrix3, Matrix6, OVector, Vector3, Vector6};
use crate::State;
use std::f64;
use std::fmt;
use std::sync::Arc;

pub use super::sph_harmonics::Harmonics;

/// `OrbitalDynamics` provides the equations of motion for any celestial dynamic, without state transition matrix computation.
#[derive(Clone)]
pub struct OrbitalDynamics<'a> {
    pub accel_models: Vec<Arc<dyn AccelModel + Sync + 'a>>,
}

impl<'a> OrbitalDynamics<'a> {
    /// Initialize point mass dynamics given the EXB IDs and a Cosm
    pub fn point_masses(bodies: &[Bodies], cosm: Arc<Cosm>) -> Arc<Self> {
        // Create the point masses
        Self::new(vec![PointMasses::new(bodies, cosm)])
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
    pub fn from_model(accel_model: Arc<dyn AccelModel + Sync + 'a>) -> Arc<Self> {
        Self::new(vec![accel_model])
    }

    /// Add a model to the currently defined orbital dynamics
    pub fn add_model(&mut self, accel_model: Arc<dyn AccelModel + Sync + 'a>) {
        self.accel_models.push(accel_model);
    }

    /// Clone these dynamics and add a model to the currently defined orbital dynamics
    pub fn with_model(self, accel_model: Arc<dyn AccelModel + Sync + 'a>) -> Arc<Self> {
        let mut me = self.clone();
        me.add_model(accel_model);
        Arc::new(me)
    }
}

impl<'a> fmt::Display for OrbitalDynamics<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let as_string: String = self
            .accel_models
            .iter()
            .map(|x| format!("{}; ", x))
            .collect();
        write!(f, "Orbital dynamics: {}", as_string)
    }
}

impl<'a> Dynamics for OrbitalDynamics<'a> {
    type HyperdualSize = Const<7>;
    type StateType = Orbit;

    fn eom(
        &self,
        delta_t_s: f64,
        state: &OVector<f64, Const<42>>,
        ctx: &Orbit,
    ) -> Result<OVector<f64, Const<42>>, NyxError> {
        let osc = ctx.set_with_delta_seconds(delta_t_s, state);
        let (new_state, new_stm) = if ctx.stm.is_some() {
            let (state, grad) = self.dual_eom(delta_t_s, &osc)?;

            let stm_dt = ctx.stm()? * grad;
            // Rebuild the STM as a vector.
            let stm_as_vec = OVector::<f64, Const<36>>::from_column_slice(stm_dt.as_slice());
            (state, stm_as_vec)
        } else {
            // Still return something of size 42, but the STM will be zeros.
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

            (d_x, OVector::<f64, Const<36>>::zeros())
        };
        Ok(OVector::<f64, Const<42>>::from_iterator(
            new_state.iter().chain(new_stm.iter()).cloned(),
        ))
    }

    fn dual_eom(
        &self,
        _delta_t_s: f64,
        osc: &Orbit,
    ) -> Result<(Vector6<f64>, Matrix6<f64>), NyxError> {
        // Extract data from hyperspace
        // Build full state vector with partials in the right position (hence building with all six components)
        let state: Vector6<Hyperdual<f64, Const<7>>> =
            hyperspace_from_vector(&osc.to_cartesian_vec());

        let radius = state.fixed_rows::<3>(0).into_owned();
        let velocity = state.fixed_rows::<3>(3).into_owned();

        // Code up math as usual
        let rmag = norm(&radius);
        let body_acceleration =
            radius * (Hyperdual::<f64, Const<7>>::from_real(-osc.frame.gm()) / rmag.powi(3));

        // Extract result into Vector6 and Matrix6
        let mut dx = Vector6::zeros();
        let mut grad = Matrix6::zeros();
        for i in 0..6 {
            dx[i] = if i < 3 {
                velocity[i].real()
            } else {
                body_acceleration[i - 3].real()
            };
            for j in 1..7 {
                grad[(i, j - 1)] = if i < 3 {
                    velocity[i][j]
                } else {
                    body_acceleration[i - 3][j]
                };
            }
        }

        // Apply the acceleration models
        for model in &self.accel_models {
            // let (model_acc, model_grad) = model.dual_eom(&radius, osc)?;
            let (model_acc, model_grad) = model.dual_eom(osc)?;
            for i in 0..3 {
                dx[i + 3] += model_acc[i];
                for j in 1..4 {
                    grad[(i + 3, j - 1)] += model_grad[(i, j - 1)];
                }
            }
        }

        // This function returns the time derivative of each function. The propagator will add this to the state vector (which has the previous STM).
        // This is why we don't multiply the gradient (A matrix) with the previous STM
        Ok((dx, grad))
    }
}

/// PointMasses model
pub struct PointMasses {
    pub bodies: Vec<Frame>,
    /// Optional point to a Cosm, needed if extra point masses are needed
    pub cosm: Arc<Cosm>,
    /// Light-time correction computation if extra point masses are needed
    pub correction: LightTimeCalc,
}

impl PointMasses {
    /// Initializes the multibody point mass dynamics with the provided list of bodies
    pub fn new(bodies: &[Bodies], cosm: Arc<Cosm>) -> Arc<Self> {
        Arc::new(Self::with_correction(bodies, cosm, LightTimeCalc::None))
    }

    /// Initializes the multibody point mass dynamics with the provided list of bodies, and accounting for some light time correction
    pub fn with_correction(bodies: &[Bodies], cosm: Arc<Cosm>, correction: LightTimeCalc) -> Self {
        let mut refs = Vec::with_capacity(bodies.len());
        // Check that these celestial bodies exist and build their references
        for body in bodies {
            refs.push(cosm.frame_from_ephem_path(body.ephem_path()));
        }

        Self {
            bodies: refs,
            cosm,
            correction,
        }
    }

    /// Allows using bodies by name, defined in the non-default XB
    pub fn specific(body_names: &[String], cosm: Arc<Cosm>, correction: LightTimeCalc) -> Self {
        let mut refs = Vec::with_capacity(body_names.len());
        // Check that these celestial bodies exist and build their references
        for body in body_names {
            refs.push(cosm.frame(body));
        }

        Self {
            bodies: refs,
            cosm,
            correction,
        }
    }
}

impl fmt::Display for PointMasses {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let body_list: String = self.bodies.iter().map(|x| format!("{}; ", x)).collect();
        write!(f, "Point masses of {}", body_list)
    }
}

impl AccelModel for PointMasses {
    fn eom(&self, osc: &Orbit) -> Result<Vector3<f64>, NyxError> {
        let mut d_x = Vector3::zeros();
        // Get all of the position vectors between the center body and the third bodies
        for third_body in &self.bodies {
            if third_body == &osc.frame {
                // Ignore the contribution of the integration frame, that's handled by OrbitalDynamics
                continue;
            }
            // Orbit of j-th body as seen from primary body
            let st_ij = self.cosm.celestial_state(
                &third_body.ephem_path(),
                osc.dt,
                osc.frame,
                self.correction,
            );

            let r_ij = st_ij.radius();
            let r_ij3 = st_ij.rmag().powi(3);
            let r_j = osc.radius() - r_ij; // sc as seen from 3rd body
            let r_j3 = r_j.norm().powi(3);
            d_x += -third_body.gm() * (r_j / r_j3 + r_ij / r_ij3);
        }
        Ok(d_x)
    }

    fn dual_eom(&self, osc: &Orbit) -> Result<(Vector3<f64>, Matrix3<f64>), NyxError> {
        // Build the hyperdual space of the radius vector
        let radius: Vector3<Hyperdual<f64, Const<7>>> = hyperspace_from_vector(&osc.radius());
        // Extract result into Vector6 and Matrix6
        let mut fx = Vector3::zeros();
        let mut grad = Matrix3::zeros();

        // Get all of the position vectors between the center body and the third bodies
        for third_body in &self.bodies {
            if third_body == &osc.frame {
                // Ignore the contribution of the integration frame, that's handled by OrbitalDynamics
                continue;
            }
            let gm_d = Hyperdual::<f64, Const<7>>::from_real(-third_body.gm());

            // Orbit of j-th body as seen from primary body
            let st_ij = self.cosm.celestial_state(
                &third_body.ephem_path(),
                osc.dt,
                osc.frame,
                self.correction,
            );

            let r_ij: Vector3<Hyperdual<f64, Const<7>>> = hyperspace_from_vector(&st_ij.radius());
            let r_ij3 = norm(&r_ij).powi(3);

            // The difference leads to the dual parts nulling themselves out, so let's fix that.
            let mut r_j = radius - r_ij; // sc as seen from 3rd body
            r_j[0][1] = 1.0;
            r_j[1][2] = 1.0;
            r_j[2][3] = 1.0;

            let r_j3 = norm(&r_j).powi(3);
            let mut third_body_acc_d = r_j / r_j3 + r_ij / r_ij3;
            third_body_acc_d[0] *= gm_d;
            third_body_acc_d[1] *= gm_d;
            third_body_acc_d[2] *= gm_d;

            let (fxp, gradp) =
                extract_jacobian_and_result::<_, Const<3>, Const<3>, _>(&third_body_acc_d);
            fx += fxp;
            grad += gradp;
        }

        Ok((fx, grad))
    }
}
