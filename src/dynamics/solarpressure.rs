use super::hyperdual::{hyperspace_from_vector, linalg::norm, vector_from_hyperspace, Hyperdual};
use super::AutoDiff;
use super::ForceModel;
use crate::dimensions::{DimName, Matrix3, Vector3, U3, U7};
use crate::time::Epoch;
use celestia::eclipse::{EclipseLocator, EclipseState};
use celestia::{Cosm, Frame, LTCorr, Orbit, AU, SPEED_OF_LIGHT};

/// Computation of solar radiation pressure is based on STK: http://help.agi.com/stk/index.htm#gator/eq-solar.htm .
#[derive(Clone)]
pub struct SolarPressure<'a> {
    /// in m^2
    pub sc_area: f64,
    /// coefficient of reflectivity, must be between 0.0 (translucent) and 2.0 (all radiation absorbed and twice the force is transmitted back).
    pub cr: f64,
    /// solar flux at 1 AU, in W/m^2
    pub phi: f64,
    pub e_loc: EclipseLocator<'a>,
}

impl<'a> SolarPressure<'a> {
    /// Will use Cr = 1.8, Phi = 1367.0
    pub fn default(sc_area: f64, shadow_bodies: Vec<Frame>, cosm: &'a Cosm) -> Self {
        let e_loc = EclipseLocator {
            light_source: cosm.frame("Sun J2000"),
            shadow_bodies,
            cosm: &cosm,
            correction: LTCorr::None,
        };
        Self {
            sc_area,
            cr: 1.8,
            phi: 1367.0,
            e_loc,
        }
    }
}

impl<'a> ForceModel for SolarPressure<'a> {
    fn eom(&self, osc: &Orbit) -> Vector3<f64> {
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
        -1e-3 * self.cr * self.sc_area * flux_pressure * r_sun_unit
    }
}

impl<'a> AutoDiff for SolarPressure<'a> {
    type HyperStateSize = U7;
    type STMSize = U3;

    fn dual_eom(
        &self,
        dt: Epoch,
        integr_frame: Frame,
        radius: &Vector3<Hyperdual<f64, U7>>,
    ) -> (Vector3<f64>, Matrix3<f64>) {
        // Extract data from hyperspace
        let cart_r = vector_from_hyperspace(&radius.fixed_rows::<U3>(0).into_owned());
        // Recreate the osculating state
        let osc = Orbit::cartesian(
            cart_r[0],
            cart_r[1],
            cart_r[2],
            0.0,
            0.0,
            0.0,
            dt,
            integr_frame,
        );

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
            Hyperdual::<f64, U7>::from_real(-1e-3 * self.cr * self.sc_area) * flux_pressure;
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

        (fx, grad)
    }
}
