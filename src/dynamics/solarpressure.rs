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

use super::{DynamicsAlmanacSnafu, DynamicsError, DynamicsPlanetarySnafu, ForceModel};
use crate::cosmic::eclipse::EclipseLocator;
use crate::cosmic::{Frame, Spacecraft, AU, SPEED_OF_LIGHT_M_S};
use crate::linalg::{Const, Matrix4x3, Vector3};
use anise::almanac::Almanac;
use anise::constants::frames::SUN_J2000;
use hyperdual::{hyperspace_from_vector, linalg::norm, Float, OHyperdual};
use log::warn;
use snafu::ResultExt;
use std::fmt;
use std::sync::Arc;

// Default solar flux in W/m^2
#[allow(non_upper_case_globals)]
pub const SOLAR_FLUX_W_m2: f64 = 1367.0;

/// Computation of solar radiation pressure is based on STK: http://help.agi.com/stk/index.htm#gator/eq-solar.htm .
#[derive(Clone)]
pub struct SolarPressure {
    /// solar flux at 1 AU, in W/m^2
    pub phi: f64,
    pub e_loc: EclipseLocator,
    /// Set to true to estimate the coefficient of reflectivity
    pub estimate: bool,
}

impl SolarPressure {
    /// Will set the solar flux at 1 AU to: Phi = 1367.0
    pub fn default_raw(
        shadow_bodies: Vec<Frame>,
        almanac: Arc<Almanac>,
    ) -> Result<Self, DynamicsError> {
        let e_loc = EclipseLocator {
            light_source: almanac.frame_from_uid(SUN_J2000).context({
                DynamicsPlanetarySnafu {
                    action: "planetary data from third body not loaded",
                }
            })?,
            shadow_bodies: shadow_bodies
                .iter()
                .filter_map(|object| match almanac.frame_from_uid(object) {
                    Ok(loaded_obj) => Some(loaded_obj),
                    Err(e) => {
                        warn!("when initializing SRP model for {object}, {e}");
                        None
                    }
                })
                .collect(),
        };
        Ok(Self {
            phi: SOLAR_FLUX_W_m2,
            e_loc,
            estimate: true,
        })
    }

    /// Accounts for the shadowing of only one body and will set the solar flux at 1 AU to: Phi = 1367.0
    pub fn default(shadow_body: Frame, almanac: Arc<Almanac>) -> Result<Arc<Self>, DynamicsError> {
        Ok(Arc::new(Self::default_raw(vec![shadow_body], almanac)?))
    }

    /// Accounts for the shadowing of only one body and will set the solar flux at 1 AU to: Phi = 1367.0
    pub fn default_no_estimation(
        shadow_bodies: Vec<Frame>,
        almanac: Arc<Almanac>,
    ) -> Result<Arc<Self>, DynamicsError> {
        let mut srp = Self::default_raw(shadow_bodies, almanac)?;
        srp.estimate = false;
        Ok(Arc::new(srp))
    }

    /// Must provide the flux in W/m^2
    pub fn with_flux(
        flux_w_m2: f64,
        shadow_bodies: Vec<Frame>,
        almanac: Arc<Almanac>,
    ) -> Result<Arc<Self>, DynamicsError> {
        let mut me = Self::default_raw(shadow_bodies, almanac)?;
        me.phi = flux_w_m2;
        Ok(Arc::new(me))
    }

    /// Solar radiation pressure force model accounting for the provided shadow bodies.
    pub fn new(
        shadow_bodies: Vec<Frame>,
        almanac: Arc<Almanac>,
    ) -> Result<Arc<Self>, DynamicsError> {
        Ok(Arc::new(Self::default_raw(shadow_bodies, almanac)?))
    }
}

impl ForceModel for SolarPressure {
    fn estimation_index(&self) -> Option<usize> {
        if self.estimate {
            Some(6)
        } else {
            None
        }
    }

    fn eom(&self, ctx: &Spacecraft, almanac: Arc<Almanac>) -> Result<Vector3<f64>, DynamicsError> {
        let osc = ctx.orbit;
        // Compute the position of the Sun as seen from the spacecraft
        let r_sun = almanac
            .transform_to(ctx.orbit, self.e_loc.light_source, None)
            .context(DynamicsAlmanacSnafu {
                action: "transforming state to vector seen from Sun",
            })?
            .radius_km;

        let r_sun_unit = r_sun / r_sun.norm();

        // ANISE returns the occultation percentage (or factor), which is the opposite as the illumination factor.
        let occult = self
            .e_loc
            .compute(osc, almanac)
            .context(DynamicsAlmanacSnafu {
                action: "solar radiation pressure computation",
            })?
            .factor();

        // Compute the illumination factor.
        let k: f64 = (occult - 1.0).abs();

        let r_sun_au = r_sun.norm() / AU;
        // in N/(m^2)
        let flux_pressure = (k * self.phi / SPEED_OF_LIGHT_M_S) * (1.0 / r_sun_au).powi(2);

        // Note the 1e-3 is to convert the SRP from m/s^2 to km/s^2
        Ok(1e-3 * ctx.srp.cr * ctx.srp.area_m2 * flux_pressure * r_sun_unit)
    }

    fn dual_eom(
        &self,
        ctx: &Spacecraft,
        almanac: Arc<Almanac>,
    ) -> Result<(Vector3<f64>, Matrix4x3<f64>), DynamicsError> {
        let osc = ctx.orbit;

        // Compute the position of the Sun as seen from the spacecraft
        let r_sun = almanac
            .transform_to(ctx.orbit, self.e_loc.light_source, None)
            .context(DynamicsAlmanacSnafu {
                action: "transforming state to vector seen from Sun",
            })?
            .radius_km;

        let r_sun_d: Vector3<OHyperdual<f64, Const<9>>> = hyperspace_from_vector(&r_sun);
        let r_sun_unit = r_sun_d / norm(&r_sun_d);

        // ANISE returns the occultation percentage (or factor), which is the opposite as the illumination factor.
        let occult = self
            .e_loc
            .compute(osc, almanac.clone())
            .context(DynamicsAlmanacSnafu {
                action: "solar radiation pressure computation",
            })?
            .factor();

        // Compute the illumination factor.
        let k: f64 = (occult - 1.0).abs();

        let r_sun_au = norm(&r_sun_d) / AU;
        let inv_r_sun_au = OHyperdual::<f64, Const<9>>::from_real(1.0) / (r_sun_au);
        let inv_r_sun_au_p2 = inv_r_sun_au.powi(2);
        // in N/(m^2)
        let flux_pressure =
            OHyperdual::<f64, Const<9>>::from_real(k * self.phi / SPEED_OF_LIGHT_M_S)
                * inv_r_sun_au_p2;

        // Note the 1e-3 is to convert the SRP from m/s^2 to km/s^2
        let dual_force_scalar =
            OHyperdual::<f64, Const<9>>::from_real(1e-3 * ctx.srp.cr * ctx.srp.area_m2);
        let mut dual_force: Vector3<OHyperdual<f64, Const<9>>> = Vector3::zeros();
        dual_force[0] = dual_force_scalar * flux_pressure * r_sun_unit[0];
        dual_force[1] = dual_force_scalar * flux_pressure * r_sun_unit[1];
        dual_force[2] = dual_force_scalar * flux_pressure * r_sun_unit[2];

        // Extract result into Vector6 and Matrix6
        let mut dx = Vector3::zeros();
        let mut grad = Matrix4x3::zeros();
        for i in 0..3 {
            dx[i] += dual_force[i].real();
            // NOTE: Although the hyperdual state is of size 7, we're only setting the values up to 3 (Matrix3)
            for j in 0..3 {
                grad[(i, j)] += dual_force[i][j + 1];
            }
        }

        // Compute the partial wrt to Cr.
        let wrt_cr = self.eom(ctx, almanac)? / ctx.srp.cr;
        for j in 0..3 {
            grad[(3, j)] = wrt_cr[j];
        }

        Ok((dx, grad))
    }
}

impl fmt::Display for SolarPressure {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "SRP with φ = {} W/m^2 and eclipse {}",
            self.phi, self.e_loc
        )
    }
}
