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
    DynamicsAlmanacSnafu, DynamicsAstroSnafu, DynamicsError, DynamicsPlanetarySnafu, ForceModel,
};
use crate::cosmic::{AstroPhysicsSnafu, Frame, Spacecraft};
use crate::dynamics::nrlmsise00::msise00_density;
use crate::io::space_weather::SpaceWeatherData;
use crate::linalg::{Const, Matrix4x3, Vector3};
use anise::almanac::Almanac;
use anise::constants::frames::IAU_EARTH_FRAME;
use anise::errors::OrientationSnafu;
use hifitime::Unit;
use hyperdual::{hyperspace_from_vector, linalg::norm, Float, OHyperdual};
use serde::{Deserialize, Serialize};
use serde_dhall::StaticType;
use snafu::ResultExt;
use std::fmt;
use std::sync::Arc;

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyType;

pub mod nrlmsise00;

/// Density in kg/m^3 and altitudes in kilometers
#[derive(Clone, Debug, Serialize, Deserialize, StaticType)]
#[cfg_attr(feature = "python", pyclass(from_py_object, get_all, set_all))]
pub enum AtmDensity {
    /// Homogeneous, static atmospheric mass density ($\text{kg/m}^3$).
    ///
    /// Ignores altitude, spatial position, and temporal variations. Useful for analytical
    /// baseline tests, sanity-checking drag accelerations, or short propagation steps.
    Constant(f64),

    /// Barometric scale-height density model.
    ///
    /// Evaluates atmospheric density using a single-layer exponential decay:
    /// $$\rho(h) = \rho_0 \exp\left(-\frac{h - h_0}{H}\right)$$
    /// where $h$ is geodetic altitude ($\text{m}$), $h_0$ (`ref_alt_m`) is the reference altitude,
    /// $\rho_0$ (`rho0`) is reference density ($\text{kg/m}^3$), and $H$ (`scale_height_m`) is
    /// the density scale height.
    ///
    /// **Limitations:** Ignores solar/geomagnetic activity and diurnal variations. Accuracy
    /// degrades rapidly outside a narrow altitude band around $h_0$.
    Exponential {
        /// Reference atmospheric density $\rho_0$ at altitude $h_0$ [kg/m³].
        rho0: f64,
        /// Reference geodetic altitude $h_0$ [km].
        ref_alt_km: f64,
        /// Atmospheric scale height $H = \frac{R T}{M g}$ [km].
        scale_height_km: f64,
    },

    /// U.S. Standard Atmosphere 1976 (USSA76) empirical density model.
    ///
    /// Evaluates piecewise atmospheric temperature and pressure profiles up to $1,000\text{ km}$ ($10^6\text{ m}$).
    /// Assumes hydrostatic equilibrium and perfect gas behavior across defined atmospheric layers.
    ///
    /// **Limitations:** Static global average model. Does not capture solar EUV heating cycles,
    /// geomagnetic storm surges, or diurnal day/night atmospheric expansion.
    StdAtm {
        /// Maximum operational altitude [km]. Above this threshold, density returns 0.0 kg/m³.
        max_alt_km: f64,
    },

    NRLMSISE00 {
        weather: SpaceWeatherData,
    },
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl AtmDensity {
    /// Constructs a standard exponential drag model for Earth orbiters.
    ///
    /// Configured with nominal LEO reference parameters at $h_0 = 700\text{ km}$:
    /// * $\rho_0 = 3.614 \times 10^{-13}\text{ kg/m}^3$
    /// * $H = 88.667\text{ km}$ ($88,667\text{ m}$)
    #[classmethod]
    fn earth_exponential(_cls: &Bound<'_, PyType>) -> Self {
        AtmDensity::Exponential {
            rho0: 3.614e-13,
            ref_alt_km: 700.000,
            scale_height_km: 88.667,
        }
    }
}

/// `Drag` implements all three drag models.
#[derive(Clone, Debug, Serialize, Deserialize, StaticType)]
#[cfg_attr(feature = "python", pyclass(from_py_object, get_all, set_all))]
pub struct Drag {
    /// Density computation method
    pub density: AtmDensity,
    /// Frame to compute the drag in
    pub frame: Frame,
    /// Set to true to estimate the coefficient of drag
    pub estimate: bool,
}

impl Drag {
    /// Constructs a standard exponential drag model for Earth orbiters.
    ///
    /// Configured with nominal LEO reference parameters at $h_0 = 700\text{ km}$:
    /// * $\rho_0 = 3.614 \times 10^{-13}\text{ kg/m}^3$
    /// * $H = 88.667\text{ km}$ ($88,667\text{ m}$)
    pub fn earth_exp(almanac: &Almanac) -> Result<Arc<Self>, DynamicsError> {
        Ok(Arc::new(Self {
            density: AtmDensity::Exponential {
                rho0: 3.614e-13,
                ref_alt_km: 700.000,
                scale_height_km: 88.667,
            },
            frame: almanac
                .frame_info(IAU_EARTH_FRAME)
                .context(DynamicsPlanetarySnafu {
                    action: "planetary data from third body not loaded",
                })?,
            estimate: false,
        }))
    }

    /// Constructs a U.S. Standard Atmosphere 1976 drag model for Earth orbiters.
    ///
    /// Valid for altitudes up to $1,000\text{ km}$ ($1,000,000\text{ m}$). Suitable for general
    /// trajectory analysis where space weather data ($F_{10.7}$, $A_p$) is unavailable.
    pub fn std_atm1976(almanac: &Almanac) -> Result<Arc<Self>, DynamicsError> {
        Ok(Arc::new(Self {
            density: AtmDensity::StdAtm {
                max_alt_km: 1_000.0,
            },
            frame: almanac
                .frame_info(IAU_EARTH_FRAME)
                .context(DynamicsPlanetarySnafu {
                    action: "planetary data from third body not loaded",
                })?,
            estimate: false,
        }))
    }
}

impl fmt::Display for Drag {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "\tDrag density {:?} in frame {}",
            self.density, self.frame
        )
    }
}

impl ForceModel for Drag {
    fn estimation_index(&self) -> Option<usize> {
        if self.estimate {
            Some(7)
        } else {
            None
        }
    }

    fn eom(&self, ctx: &Spacecraft, almanac: &Almanac) -> Result<Vector3<f64>, DynamicsError> {
        let integration_frame = ctx.orbit.frame;

        let drag_frame = almanac
            .frame_info(self.frame)
            .context(DynamicsPlanetarySnafu {
                action: "fetching drag frame information",
            })?;

        let osc_drag_frame =
            almanac
                .transform_to(ctx.orbit, self.frame, None)
                .context(DynamicsAlmanacSnafu {
                    action: "transforming into drag frame",
                })?;

        let rho_kg_m3 = match &self.density {
            AtmDensity::Constant(rho) => *rho,

            AtmDensity::Exponential {
                rho0,
                scale_height_km,
                ref_alt_km,
            } => {
                let altitude_km = osc_drag_frame
                    .altitude_km()
                    .context(AstroPhysicsSnafu)
                    .context(DynamicsAstroSnafu)?;
                rho0 * (-(altitude_km - ref_alt_km) / scale_height_km).exp()
            }

            AtmDensity::StdAtm { max_alt_km } => {
                let altitude_km = osc_drag_frame
                    .altitude_km()
                    .context(AstroPhysicsSnafu)
                    .context(DynamicsAstroSnafu)?;

                if altitude_km > *max_alt_km {
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

                    // Calculating density by raising 10 to the log of density
                    10.0_f64.powf(logdensity)
                }
            }

            AtmDensity::NRLMSISE00 { weather } => {
                let (lat_deg, long_deg, alt_km) = osc_drag_frame
                    .latlongalt()
                    .context(AstroPhysicsSnafu)
                    .context(DynamicsAstroSnafu)?;

                let lst_h = almanac.local_solar_time(osc_drag_frame, None).context(
                    DynamicsAlmanacSnafu {
                        action: "computing local solar time",
                    },
                )?;

                let epoch = ctx.orbit.epoch;

                let sw = weather.msise_weather(epoch);

                msise00_density(
                    sw,
                    lst_h.to_unit(Unit::Hour),
                    lat_deg,
                    long_deg,
                    alt_km,
                    ctx.orbit.epoch,
                )?
                .total_mass_density_kg_m3
            }
        };

        let v_km_s = osc_drag_frame.velocity_km_s;

        // Note this is in kg*km/s^2 (or kN) because the vehicle mass has not yet been divided.
        let accel_drag_frame_kg_km_s2 = -0.5
            * 1e3
            * rho_kg_m3
            * ctx.drag.coeff_drag
            * ctx.drag.area_m2
            * v_km_s.norm()
            * v_km_s;

        let accel_integr_frame = almanac
            .rotate(drag_frame, integration_frame, ctx.orbit.epoch)
            .context(OrientationSnafu {
                action: "rotating drafg force into integration frame",
            })
            .context(DynamicsAlmanacSnafu {
                action: "rotating drag force into integration frame",
            })?
            * accel_drag_frame_kg_km_s2;

        // Finally, apply the drag model.
        Ok(accel_integr_frame)
    }

    fn gradient(
        &self,
        ctx: &Spacecraft,
        almanac: &Almanac,
    ) -> Result<(Vector3<f64>, Matrix4x3<f64>), DynamicsError> {
        todo!()
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl Drag {
    #[pyo3(signature = (density, frame, estimate=true))]
    #[new]
    fn py_new(density: AtmDensity, frame: Frame, estimate: bool) -> Self {
        Self {
            density,
            frame,
            estimate,
        }
    }

    fn __str__(&self) -> String {
        format!("{self}")
    }

    fn __repr__(&self) -> String {
        format!("{self} @ {self:p}")
    }
}
