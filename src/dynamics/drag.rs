use super::ForceModel;
use crate::dimensions::Vector3;
use celestia::{Cosm, State};

/// `ConstantDrag` implements a constant drag model as defined in Vallado, 4th ed., page 551, with an important caveat.
///
/// **WARNING:** This basic model assumes that the velocity of the spacecraft is identical to the velocity of the upper atmosphere,
/// This is a **bad** assumption and **should not** be used for high fidelity simulations.
/// This will be resolved after https://gitlab.com/chrisrabotin/nyx/issues/93 is implemented.
#[derive(Clone)]
pub struct ConstantDrag<'a> {
    /// in m^2
    pub sc_area: f64,
    /// coefficient of drag; (spheres are between 2.0 and 2.1, use 2.2 in Earth's atmosphere).
    pub cd: f64,
    /// atmospheric density in kg/m^3
    pub rho: f64,
    /// Geoid causing the drag
    pub drag_frame_id: i32,
    /// a Cosm reference is needed to convert to the state around the correct planet
    pub cosm: &'a Cosm,
}

impl<'a> ForceModel for ConstantDrag<'a> {
    fn eom(&self, osc: &State) -> Vector3<f64> {
        let osc = self.cosm.frame_chg_by_id(&osc, self.drag_frame_id);
        let velocity = osc.velocity();
        -0.5 * self.rho * self.cd * self.sc_area * velocity.norm() * velocity
    }
}

/// `ExpEarthDrag` implements an exponential decay drag model.
///
/// **WARNING:** This model assumes that the velocity of the spacecraft is identical to the velocity of the upper atmosphere,
/// This is a **bad** assumption and **should not** be used for high fidelity simulations.
/// /// This will be resolved after https://gitlab.com/chrisrabotin/nyx/issues/93 is implemented.
#[derive(Clone)]
pub struct ExpEarthDrag<'a> {
    /// in m^2
    pub sc_area: f64,
    /// coefficient of drag; (spheres are between 2.0 and 2.1, use 2.2 in Earth's atmosphere).
    pub cd: f64,
    /// a Cosm reference is needed to convert to the state around the correct planet
    pub cosm: &'a Cosm,
}

impl<'a> ForceModel for ExpEarthDrag<'a> {
    fn eom(&self, osc: &State) -> Vector3<f64> {
        let eme2k = self.cosm.frame("EME2000");
        // Compute the density
        let rho0 = 3.614e-13; // # kg/m^3
        let r0 = 700_000.0 + eme2k.equatorial_radius();
        let h = 88_667.0; // m
        let rho = rho0 * (-(osc.rmag() - r0) / h).exp(); // # Exponential decay model for density

        let osc = self.cosm.frame_chg(&osc, eme2k);

        // Incorrectly transform to some ECEF frame
        let earth_rot = 7.292_115_855_3e-5;

        let velocity = osc.velocity() - Vector3::new(earth_rot * osc.y, -earth_rot * osc.x, 0.0);
        -0.5 * rho * self.cd * self.sc_area * velocity.norm() * velocity
    }
}
