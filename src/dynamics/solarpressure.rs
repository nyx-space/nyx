use super::na::Vector3;
use super::ForceModel;
use celestia::eclipse::{EclipseLocator, EclipseState};
use celestia::{bodies, Cosm, Geoid, LTCorr, State, AU, SPEED_OF_LIGHT};

/// Computation of solar radiation pressure is based on STK: http://help.agi.com/stk/index.htm#gator/eq-solar.htm .
#[derive(Clone)]
pub struct SolarPressure<'a> {
    /// in kg, set in the Spacecraft's eom.
    pub sc_mass: f64,
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
    pub fn default(sc_area: f64, shadow_bodies: Vec<Geoid>, cosm: &'a Cosm) -> Self {
        let sun = cosm.geoid_from_id(bodies::SUN);
        let e_loc = EclipseLocator {
            light_source: sun,
            shadow_bodies,
            cosm: &cosm,
            correction: LTCorr::None,
        };
        Self {
            sc_mass: 0.0,
            sc_area,
            cr: 1.8,
            phi: 1367.0,
            e_loc,
        }
    }
}

impl<'a> ForceModel<Geoid> for SolarPressure<'a> {
    fn eom(&self, osc: &State<Geoid>) -> Vector3<f64> {
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
