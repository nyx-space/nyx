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

use super::hyperdual::{hyperspace_from_vector, linalg::norm, Hyperdual};
use super::ForceModel;
use crate::celestia::eclipse::{EclipseLocator, EclipseState};
use crate::celestia::{Cosm, Frame, LTCorr, Spacecraft, AU, SPEED_OF_LIGHT};
use crate::dimensions::{DimName, Matrix3, Vector3, U3, U7};
use crate::errors::NyxError;
use std::sync::Arc;

/// Computation of solar radiation pressure is based on STK: http://help.agi.com/stk/index.htm#gator/eq-solar.htm .
#[derive(Clone)]
pub struct SolarPressure {
    /// solar flux at 1 AU, in W/m^2
    pub phi: f64,
    pub e_loc: EclipseLocator,
}

impl<'a> SolarPressure {
    /// Will set the solar flux at 1 AU to: Phi = 1367.0
    pub fn default_raw(shadow_bodies: Vec<Frame>, cosm: Arc<Cosm>) -> Self {
        let e_loc = EclipseLocator {
            light_source: cosm.frame("Sun J2000"),
            shadow_bodies,
            cosm,
            correction: LTCorr::None,
        };
        Self { phi: 1367.0, e_loc }
    }

    /// Accounts for the shadowing of only one body and will set the solar flux at 1 AU to: Phi = 1367.0
    pub fn default(shadow_body: Frame, cosm: Arc<Cosm>) -> Arc<Self> {
        Arc::new(Self::default_raw(vec![shadow_body], cosm))
    }

    /// Must provide the flux in W/m^2
    pub fn with_flux(flux_w_m2: f64, shadow_bodies: Vec<Frame>, cosm: Arc<Cosm>) -> Arc<Self> {
        let mut me = Self::default_raw(shadow_bodies, cosm);
        me.phi = flux_w_m2;
        Arc::new(me)
    }
}

impl ForceModel for SolarPressure {
    fn eom(&self, ctx: &Spacecraft) -> Result<Vector3<f64>, NyxError> {
        let osc = &ctx.orbit;
        // Compute the position of the Sun as seen from the spacecraft
        let r_sun = self
            .e_loc
            .cosm
            .frame_chg(osc, self.e_loc.light_source)
            .radius();
        let r_sun_unit = r_sun / r_sun.norm();

        // Compute the shaddowing factor.
        let k = match self.e_loc.compute(osc) {
            EclipseState::Umbra => 0.0,
            EclipseState::Visibilis => 1.0,
            EclipseState::Penumbra(val) => val,
        };

        let r_sun_au = r_sun.norm() / AU;
        // in N/(m^2)
        let flux_pressure = (k * self.phi / SPEED_OF_LIGHT) * (1.0 / r_sun_au).powi(2);

        // Note the 1e-3 is to convert the SRP from m/s^2 to km/s^2
        Ok(-1e-3 * ctx.cr * ctx.srp_area_m2 * flux_pressure * r_sun_unit)
    }

    fn dual_eom(
        &self,
        _radius: &Vector3<Hyperdual<f64, U7>>,
        ctx: &Spacecraft,
    ) -> Result<(Vector3<f64>, Matrix3<f64>), NyxError> {
        let osc = ctx.orbit;

        // Compute the position of the Sun as seen from the spacecraft
        let r_sun = self
            .e_loc
            .cosm
            .frame_chg(&osc, self.e_loc.light_source)
            .radius();

        let r_sun_d: Vector3<Hyperdual<f64, U7>> = hyperspace_from_vector(&r_sun);
        let r_sun_unit = r_sun_d / norm(&r_sun_d);

        // Compute the shaddowing factor.
        let k = match self.e_loc.compute(&osc) {
            EclipseState::Umbra => 0.0,
            EclipseState::Visibilis => 1.0,
            EclipseState::Penumbra(val) => val,
        };

        let inv_r_sun_au = Hyperdual::<f64, U7>::from_real(1.0) / (norm(&r_sun_d) / AU);
        let inv_r_sun_au_p2 = inv_r_sun_au * inv_r_sun_au;
        // in N/(m^2)
        let flux_pressure =
            Hyperdual::<f64, U7>::from_real(k * self.phi / SPEED_OF_LIGHT) * inv_r_sun_au_p2;

        // Note the 1e-3 is to convert the SRP from m/s^2 to km/s^2
        let dual_force_scalar =
            Hyperdual::<f64, U7>::from_real(-1e-3 * ctx.cr * ctx.srp_area_m2) * flux_pressure;
        let mut dual_force: Vector3<Hyperdual<f64, U7>> = Vector3::zeros();
        dual_force[0] = dual_force_scalar * r_sun_unit[0];
        dual_force[1] = dual_force_scalar * r_sun_unit[1];
        dual_force[2] = dual_force_scalar * r_sun_unit[2];

        // Extract result into Vector6 and Matrix6
        let mut fx = Vector3::zeros();
        let mut grad = Matrix3::zeros();
        for i in 0..U3::dim() {
            fx[i] += dual_force[i][0];
            // NOTE: Although the hyperdual state is of size 7, we're only setting the values up to 3 (Matrix3)
            for j in 0..U3::dim() {
                grad[(i, j)] += dual_force[i][j + 1];
            }
        }

        Ok((fx, grad))
    }
}
