/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::ForceModel;
use crate::cosmic::{BaseSpacecraft, Cosm, Frame, SpacecraftExt};
use crate::errors::NyxError;
use crate::linalg::{Matrix3, Vector3};
use std::fmt;
use std::sync::Arc;

/// Density in kg/m^3 and altitudes in meters, not kilometers!
#[derive(Clone, Copy, Debug)]
pub enum AtmDensity {
    Constant(f64),
    Exponential { rho0: f64, r0: f64, ref_alt_m: f64 },
    StdAtm { max_alt_m: f64 },
}

/// `ConstantDrag` implements a constant drag model as defined in Vallado, 4th ed., page 551, with an important caveat.
///
/// **WARNING:** This basic model assumes that the velocity of the spacecraft is identical to the velocity of the upper atmosphere,
/// This is a **bad** assumption and **should not** be used for high fidelity simulations.
/// This will be resolved after https://gitlab.com/chrisrabotin/nyx/issues/93 is implemented.
#[derive(Clone)]
pub struct ConstantDrag {
    /// atmospheric density in kg/m^3
    pub rho: f64,
    /// Geoid causing the drag
    pub drag_frame: Frame,
    /// a Cosm reference is needed to convert to the state around the correct planet
    pub cosm: Arc<Cosm>,
}

impl fmt::Display for ConstantDrag {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\tConstant Drag rho = {} kg/m^3 in frame {}",
            self.rho, self.drag_frame
        )
    }
}

impl<X: SpacecraftExt> ForceModel<X> for ConstantDrag {
    fn eom(&self, ctx: &BaseSpacecraft<X>) -> Result<Vector3<f64>, NyxError> {
        let osc = self.cosm.frame_chg(&ctx.orbit, self.drag_frame);
        let velocity = osc.velocity();
        Ok(-0.5 * self.rho * ctx.cd * ctx.drag_area_m2 * velocity.norm() * velocity)
    }

    fn dual_eom(
        &self,
        _osc_ctx: &BaseSpacecraft<X>,
    ) -> Result<(Vector3<f64>, Matrix3<f64>), NyxError> {
        Err(NyxError::PartialsUndefined)
    }
}

/// `Drag` implements all three drag models.
#[derive(Clone)]
pub struct Drag {
    /// Density computation method
    pub density: AtmDensity,
    /// Frame to compute the drag in
    pub drag_frame: Frame,
    /// a Cosm reference is needed to convert to the state around the correct planet
    pub cosm: Arc<Cosm>,
}

impl Drag {
    /// Common exponential drag model for the Earth
    pub fn earth_exp(cosm: Arc<Cosm>) -> Arc<Self> {
        Arc::new(Self {
            density: AtmDensity::Exponential {
                rho0: 3.614e-13,
                r0: 700_000.0,
                ref_alt_m: 88_667.0,
            },
            drag_frame: cosm.frame("IAU Earth"),
            cosm,
        })
    }

    /// Drag model which uses the standard atmosphere 1976 model for atmospheric density
    pub fn std_atm1976(cosm: Arc<Cosm>) -> Arc<Self> {
        Arc::new(Self {
            density: AtmDensity::StdAtm {
                max_alt_m: 1_000_000.0,
            },
            drag_frame: cosm.frame("IAU Earth"),
            cosm,
        })
    }
}

impl fmt::Display for Drag {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\tDrag density {:?} in frame {}",
            self.density, self.drag_frame
        )
    }
}

impl<X: SpacecraftExt> ForceModel<X> for Drag {
    fn eom(&self, ctx: &BaseSpacecraft<X>) -> Result<Vector3<f64>, NyxError> {
        let integration_frame = ctx.orbit.frame;
        let osc = self.cosm.frame_chg(&ctx.orbit, self.drag_frame);
        match self.density {
            AtmDensity::Constant(rho) => {
                let velocity = osc.velocity();
                Ok(-0.5 * rho * ctx.cd * ctx.drag_area_m2 * velocity.norm() * velocity)
            }

            AtmDensity::Exponential {
                rho0,
                r0,
                ref_alt_m,
            } => {
                let rho = rho0
                    * (-(osc.rmag() - (r0 + self.drag_frame.equatorial_radius())) / ref_alt_m)
                        .exp();

                let velocity_integr_frame = self.cosm.frame_chg(&osc, integration_frame).velocity();

                let velocity = velocity_integr_frame - osc.velocity();
                Ok(-0.5 * rho * ctx.cd * ctx.drag_area_m2 * velocity.norm() * velocity)
            }

            AtmDensity::StdAtm { max_alt_m } => {
                let altitude_km = osc.rmag() - self.drag_frame.equatorial_radius();
                let rho = if altitude_km > max_alt_m / 1_000.0 {
                    // Use a constant density
                    10.0_f64.powf((-7e-5) * altitude_km - 14.464)
                } else {
                    // Code from AVS/Schaub's Basilisk
                    // Calculating the density based on a scaled 6th order polynomial fit to the log of density
                    let scale = (altitude_km - 526.8000) / 292.8563;
                    let logdensity =
                        0.34047 * scale.powi(6) - 0.5889 * scale.powi(5) - 0.5269 * scale.powi(4)
                            + 1.0036 * scale.powi(3)
                            + 0.60713 * scale.powi(2)
                            - 2.3024 * scale
                            - 12.575;

                    /* Calculating density by raising 10 to the log of density */
                    10.0_f64.powf(logdensity)
                };

                let velocity_integr_frame = self.cosm.frame_chg(&osc, integration_frame).velocity();

                let velocity = velocity_integr_frame - osc.velocity();
                Ok(-0.5 * rho * ctx.cd * ctx.drag_area_m2 * velocity.norm() * velocity)
            }
        }
    }

    fn dual_eom(
        &self,
        _osc_ctx: &BaseSpacecraft<X>,
    ) -> Result<(Vector3<f64>, Matrix3<f64>), NyxError> {
        Err(NyxError::PartialsUndefined)
    }
}
