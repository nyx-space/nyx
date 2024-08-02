/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::{
    AccelModel, DynamicsAlmanacSnafu, DynamicsAstroSnafu, DynamicsError, DynamicsPlanetarySnafu,
};
use crate::cosmic::{AstroPhysicsSnafu, Frame, Orbit};
use crate::linalg::{Const, Matrix3, Matrix6, OVector, Vector3, Vector6};

use anise::almanac::Almanac;
use anise::astro::Aberration;
use hyperdual::linalg::norm;
use hyperdual::{extract_jacobian_and_result, hyperspace_from_vector, Float, OHyperdual};
use snafu::ResultExt;
use std::f64;
use std::fmt;
use std::sync::Arc;

pub use super::sph_harmonics::Harmonics;

/// `OrbitalDynamics` provides the equations of motion for any celestial dynamic, without state transition matrix computation.
#[derive(Clone)]

pub struct OrbitalDynamics {
    pub accel_models: Vec<Arc<dyn AccelModel + Sync>>,
}

impl OrbitalDynamics {
    /// Initializes the point masses gravities with the provided list of bodies
    pub fn point_masses(celestial_objects: Vec<i32>) -> Self {
        // Create the point masses
        Self::new(vec![PointMasses::new(celestial_objects)])
    }

    /// Initializes a OrbitalDynamics which does not simulate the gravity pull of other celestial objects but the primary one.
    pub fn two_body() -> Self {
        Self::new(vec![])
    }

    /// Initialize orbital dynamics with a list of acceleration models
    pub fn new(accel_models: Vec<Arc<dyn AccelModel + Sync>>) -> Self {
        Self { accel_models }
    }

    /// Initialize new orbital mechanics with the provided model.
    /// **Note:** Orbital dynamics _always_ include two body dynamics, these cannot be turned off.
    pub fn from_model(accel_model: Arc<dyn AccelModel + Sync>) -> Self {
        Self::new(vec![accel_model])
    }
}

impl fmt::Display for OrbitalDynamics {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let models: Vec<String> = self.accel_models.iter().map(|x| format!("{x}")).collect();
        write!(f, "Orbital dynamics: {}", models.join("; "))
    }
}

impl OrbitalDynamics {
    pub(crate) fn eom(
        &self,
        osc: &Orbit,
        almanac: Arc<Almanac>,
    ) -> Result<OVector<f64, Const<42>>, DynamicsError> {
        // Still return something of size 42, but the STM will be zeros.
        let body_acceleration = (-osc
            .frame
            .mu_km3_s2()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?
            / osc.rmag_km().powi(3))
            * osc.radius_km;

        let mut d_x = Vector6::from_iterator(
            osc.velocity_km_s
                .iter()
                .chain(body_acceleration.iter())
                .cloned(),
        );

        // Apply the acceleration models
        for model in &self.accel_models {
            let model_acc = model.eom(osc, almanac.clone())?;
            for i in 0..3 {
                d_x[i + 3] += model_acc[i];
            }
        }

        Ok(OVector::<f64, Const<42>>::from_iterator(
            d_x.iter()
                .chain(OVector::<f64, Const<36>>::zeros().iter())
                .cloned(),
        ))
    }

    pub fn dual_eom(
        &self,
        _delta_t_s: f64,
        osc: &Orbit,
        almanac: Arc<Almanac>,
    ) -> Result<(Vector6<f64>, Matrix6<f64>), DynamicsError> {
        // Extract data from hyperspace
        // Build full state vector with partials in the right position (hence building with all six components)
        let state: Vector6<OHyperdual<f64, Const<7>>> =
            hyperspace_from_vector(&osc.to_cartesian_pos_vel());

        let radius = state.fixed_rows::<3>(0).into_owned();
        let velocity = state.fixed_rows::<3>(3).into_owned();

        // Code up math as usual
        let rmag = norm(&radius);
        let body_acceleration = radius
            * (OHyperdual::<f64, Const<7>>::from_real(
                -osc.frame
                    .mu_km3_s2()
                    .context(AstroPhysicsSnafu)
                    .context(DynamicsAstroSnafu)?,
            ) / rmag.powi(3));

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
            let (model_acc, model_grad) = model.dual_eom(osc, almanac.clone())?;
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
    pub celestial_objects: Vec<i32>,
    /// Light-time correction computation if extra point masses are needed
    pub correction: Option<Aberration>,
}

impl PointMasses {
    /// Initializes the point masses gravities with the provided list of bodies
    pub fn new(celestial_objects: Vec<i32>) -> Arc<Self> {
        Arc::new(Self {
            celestial_objects,
            correction: None,
        })
    }

    /// Initializes the point masses gravities with the provided list of bodies, and accounting for some light time correction
    pub fn with_correction(celestial_objects: Vec<i32>, correction: Aberration) -> Self {
        Self {
            celestial_objects,
            correction: Some(correction),
        }
    }
}

impl fmt::Display for PointMasses {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let masses: Vec<String> = self
            .celestial_objects
            .iter()
            .map(|third_body| format!("{}", Frame::from_ephem_j2000(*third_body)))
            .collect();
        write!(f, "Point masses of {}", masses.join(", "))
    }
}

impl AccelModel for PointMasses {
    fn eom(&self, osc: &Orbit, almanac: Arc<Almanac>) -> Result<Vector3<f64>, DynamicsError> {
        let mut d_x = Vector3::zeros();
        // Get all of the position vectors between the center body and the third bodies
        for third_body in self.celestial_objects.iter().copied() {
            if osc.frame.ephem_origin_id_match(third_body) {
                // Ignore the contribution of the integration frame, that's handled by OrbitalDynamics
                continue;
            }

            let third_body_frame = almanac
                .frame_from_uid(osc.frame.with_ephem(third_body))
                .context(DynamicsPlanetarySnafu {
                    action: "planetary data from third body not loaded",
                })?;

            // Orbit of j-th body as seen from primary body
            let st_ij = almanac
                .transform(third_body_frame, osc.frame, osc.epoch, self.correction)
                .context(DynamicsAlmanacSnafu {
                    action: "computing third body gravitational pull",
                })?;

            let r_ij = st_ij.radius_km;
            let r_ij3 = st_ij.rmag_km().powi(3);
            let r_j = osc.radius_km - r_ij; // sc as seen from 3rd body
            let r_j3 = r_j.norm().powi(3);
            d_x += -third_body_frame
                .mu_km3_s2()
                .context(AstroPhysicsSnafu)
                .context(DynamicsAstroSnafu)?
                * (r_j / r_j3 + r_ij / r_ij3);
        }
        Ok(d_x)
    }

    fn dual_eom(
        &self,
        osc: &Orbit,
        almanac: Arc<Almanac>,
    ) -> Result<(Vector3<f64>, Matrix3<f64>), DynamicsError> {
        // Build the hyperdual space of the radius vector
        let radius: Vector3<OHyperdual<f64, Const<7>>> = hyperspace_from_vector(&osc.radius_km);
        // Extract result into Vector6 and Matrix6
        let mut fx = Vector3::zeros();
        let mut grad = Matrix3::zeros();

        // Get all of the position vectors between the center body and the third bodies
        for third_body in &self.celestial_objects {
            let third_body_frame = almanac
                .frame_from_uid(Frame::from_ephem_j2000(*third_body))
                .context(DynamicsPlanetarySnafu {
                    action: "planetary data from third body not loaded",
                })?;

            if osc.frame.ephem_origin_match(third_body_frame) {
                // Ignore the contribution of the integration frame, that's handled by OrbitalDynamics
                continue;
            }

            let gm_d = OHyperdual::<f64, Const<7>>::from_real(
                -third_body_frame
                    .mu_km3_s2()
                    .context(AstroPhysicsSnafu)
                    .context(DynamicsAstroSnafu)?,
            );

            // Orbit of j-th body as seen from primary body
            let st_ij = almanac
                .transform(third_body_frame, osc.frame, osc.epoch, self.correction)
                .context(DynamicsAlmanacSnafu {
                    action: "computing third body gravitational pull",
                })?;

            let r_ij: Vector3<OHyperdual<f64, Const<7>>> = hyperspace_from_vector(&st_ij.radius_km);
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

            let (fxp, gradp) = extract_jacobian_and_result::<_, 3, 3, 7>(&third_body_acc_d);
            fx += fxp;
            grad += gradp;
        }

        Ok((fx, grad))
    }
}
