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
pub use crate::dynamics::nrlmsise00::{GeomagneticMode, Nrlmsise00Flags, OutputUnits};
use crate::io::space_weather::SpaceWeatherData;
use crate::linalg::{Matrix4x3, Vector3};
use anise::almanac::Almanac;
use anise::constants::frames::IAU_EARTH_FRAME;
use anise::errors::OrientationSnafu;
use hifitime::Unit;
use serde::{Deserialize, Serialize};
use serde_dhall::StaticType;
use snafu::ResultExt;
use std::fmt;
use std::sync::Arc;
use trig_const::{ln, sqrt};

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyType;

pub mod nrlmsise00;

const SIN_30_DEG: f64 = 0.5;
const COS_30_DEG: f64 = sqrt(3.0_f64) / 2.0_f64;

/// Standard Harris-Priester altitude profile node (100 km to 1000 km)
/// Densities in kg/m^3
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
struct HpNode {
    alt_km: f64,
    min_density_kg_m3: f64,
    max_density_kg_m3: f64,
}

impl HpNode {
    const fn ln_min_density(&self) -> f64 {
        ln(self.min_density_kg_m3)
    }
    const fn ln_max_density(&self) -> f64 {
        ln(self.max_density_kg_m3)
    }
}

/// Baseline Harris-Priester reference table for mean solar activity (~F10.7 = 150)
const HP_TABLE: &[HpNode] = &[
    HpNode {
        alt_km: 100.0,
        min_density_kg_m3: 4.974e-07,
        max_density_kg_m3: 4.974e-07,
    },
    HpNode {
        alt_km: 120.0,
        min_density_kg_m3: 2.490e-08,
        max_density_kg_m3: 2.490e-08,
    },
    HpNode {
        alt_km: 140.0,
        min_density_kg_m3: 3.840e-09,
        max_density_kg_m3: 3.840e-09,
    },
    HpNode {
        alt_km: 160.0,
        min_density_kg_m3: 1.170e-09,
        max_density_kg_m3: 1.170e-09,
    },
    HpNode {
        alt_km: 180.0,
        min_density_kg_m3: 4.820e-10,
        max_density_kg_m3: 5.220e-10,
    },
    HpNode {
        alt_km: 200.0,
        min_density_kg_m3: 2.260e-10,
        max_density_kg_m3: 2.620e-10,
    },
    HpNode {
        alt_km: 240.0,
        min_density_kg_m3: 6.880e-11,
        max_density_kg_m3: 9.380e-11,
    },
    HpNode {
        alt_km: 280.0,
        min_density_kg_m3: 2.570e-11,
        max_density_kg_m3: 4.180e-11,
    },
    HpNode {
        alt_km: 320.0,
        min_density_kg_m3: 1.090e-11,
        max_density_kg_m3: 2.060e-11,
    },
    HpNode {
        alt_km: 360.0,
        min_density_kg_m3: 4.980e-12,
        max_density_kg_m3: 1.070e-11,
    },
    HpNode {
        alt_km: 400.0,
        min_density_kg_m3: 2.380e-12,
        max_density_kg_m3: 5.820e-12,
    },
    HpNode {
        alt_km: 440.0,
        min_density_kg_m3: 1.180e-12,
        max_density_kg_m3: 3.250e-12,
    },
    HpNode {
        alt_km: 480.0,
        min_density_kg_m3: 6.020e-13,
        max_density_kg_m3: 1.860e-12,
    },
    HpNode {
        alt_km: 520.0,
        min_density_kg_m3: 3.150e-13,
        max_density_kg_m3: 1.080e-12,
    },
    HpNode {
        alt_km: 560.0,
        min_density_kg_m3: 1.680e-13,
        max_density_kg_m3: 6.400e-13,
    },
    HpNode {
        alt_km: 600.0,
        min_density_kg_m3: 9.100e-14,
        max_density_kg_m3: 3.830e-13,
    },
    HpNode {
        alt_km: 680.0,
        min_density_kg_m3: 2.820e-14,
        max_density_kg_m3: 1.440e-13,
    },
    HpNode {
        alt_km: 760.0,
        min_density_kg_m3: 9.200e-15,
        max_density_kg_m3: 5.760e-14,
    },
    HpNode {
        alt_km: 840.0,
        min_density_kg_m3: 3.100e-15,
        max_density_kg_m3: 2.400e-14,
    },
    HpNode {
        alt_km: 920.0,
        min_density_kg_m3: 1.100e-15,
        max_density_kg_m3: 1.050e-14,
    },
    HpNode {
        alt_km: 1000.0,
        min_density_kg_m3: 4.000e-16,
        max_density_kg_m3: 4.800e-15,
    },
];

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
        rho0_kg_m3: f64,
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

    /// NRLMSISE-00 empirical atmosphere model.
    ///
    /// Computes neutral atmospheric density and composition from 0 to ~1000 km altitude
    /// as a function of location, time, solar activity (F10.7), and geomagnetic
    /// activity (Ap).
    NRLMSISE00 {
        weather: SpaceWeatherData,
        #[serde(default)]
        flags: Option<Nrlmsise00Flags>,
    },

    /// Harris-Priester atmospheric density model.
    ///
    /// Computes density accounting for diurnal atmospheric expansion using tabular
    /// min/max density profiles interpolated exponentially across altitude bands,
    /// modified by the diurnal bulge angle offset.
    HarrisPriester {
        /// Diurnal bulge parameter $n$ (typically 2 for low inclination, 6 for polar).
        n_parameter: usize,
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
            rho0_kg_m3: 3.614e-13,
            ref_alt_km: 700.000,
            scale_height_km: 88.667,
        }
    }

    /// Constructs an NRLMSISE-00 atmospheric density model with optional flags.
    #[classmethod]
    #[pyo3(name = "NRLMSISE00")]
    #[pyo3(signature = (weather, flags = None))]
    fn py_nrlmsise00(
        _cls: &Bound<'_, PyType>,
        weather: SpaceWeatherData,
        flags: Option<Nrlmsise00Flags>,
    ) -> Self {
        AtmDensity::NRLMSISE00 { weather, flags }
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
                rho0_kg_m3: 3.614e-13,
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

    /// Calculate the density as a private function, since it's duplicated in the EOM and Gradient
    fn rho_kg_m3(&self, ctx: &Spacecraft, almanac: &Almanac) -> Result<f64, DynamicsError> {
        let osc_drag_frame =
            almanac
                .transform_to(ctx.orbit, self.frame, None)
                .context(DynamicsAlmanacSnafu {
                    action: "transforming into drag frame",
                })?;

        let rho_kg_m3 = match &self.density {
            AtmDensity::Constant(rho) => *rho,

            AtmDensity::Exponential {
                rho0_kg_m3,
                scale_height_km,
                ref_alt_km,
            } => {
                let altitude_km = osc_drag_frame
                    .altitude_km()
                    .context(AstroPhysicsSnafu)
                    .context(DynamicsAstroSnafu)?;
                rho0_kg_m3 * (-(altitude_km - ref_alt_km) / scale_height_km).exp()
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

            AtmDensity::NRLMSISE00 { weather, flags } => {
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
                    flags.unwrap_or(Nrlmsise00Flags {
                        units: OutputUnits::Si,
                        geomagnetic: GeomagneticMode::StandardDailyAp,
                        ..Default::default()
                    }),
                )?
                .total_mass_density_kg_m3
            }

            AtmDensity::HarrisPriester { n_parameter } => {
                let altitude_km = osc_drag_frame
                    .altitude_km()
                    .context(AstroPhysicsSnafu)
                    .context(DynamicsAstroSnafu)?;

                if altitude_km < HP_TABLE[0].alt_km
                    || altitude_km > HP_TABLE[HP_TABLE.len() - 1].alt_km
                {
                    0.0
                } else {
                    // Find altitude layer
                    let idx = HP_TABLE
                        .windows(2)
                        .position(|w| altitude_km >= w[0].alt_km && altitude_km <= w[1].alt_km)
                        .unwrap_or(0);

                    let n0 = &HP_TABLE[idx];
                    let n1 = &HP_TABLE[idx + 1];

                    // Scale height interpolation for min and max density
                    let h_min =
                        (n0.alt_km - n1.alt_km) / (n1.ln_min_density() - n0.ln_min_density());
                    let h_max =
                        (n0.alt_km - n1.alt_km) / (n1.ln_max_density() - n0.ln_max_density());

                    let rho_min = n0.min_density_kg_m3 * (-(altitude_km - n0.alt_km) / h_min).exp();
                    let rho_max = n0.max_density_kg_m3 * (-(altitude_km - n0.alt_km) / h_max).exp();

                    // Compute Sun unit vector in drag frame
                    let u_sun = almanac
                        .sun_unit_vector(ctx.orbit.epoch, self.frame, None)
                        .context(DynamicsAlmanacSnafu {
                            action: "fetching sun position for Harris-Priester model",
                        })?;

                    // Diurnal bulge apex: lagging Sun by ~30 deg in Right Ascension
                    let u_bulge = Vector3::new(
                        u_sun.x * COS_30_DEG - u_sun.y * SIN_30_DEG,
                        u_sun.x * SIN_30_DEG + u_sun.y * COS_30_DEG,
                        u_sun.z,
                    );

                    // Cosine of angle between spacecraft position vector and diurnal bulge apex
                    let u_pos = osc_drag_frame.r_hat();
                    let cos_psi = u_pos.dot(&u_bulge).clamp(-1.0, 1.0);

                    // Diurnal variation modifier: cos^(n)(psi / 2)
                    let cos_half_psi = ((1.0 + cos_psi) / 2.0).sqrt();
                    let mod_factor = cos_half_psi.powi(*n_parameter as i32);

                    rho_min + (rho_max - rho_min) * mod_factor
                }
            }
        };

        Ok(rho_kg_m3)
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
        if self.estimate { Some(7) } else { None }
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

        let rho_kg_m3 = self.rho_kg_m3(ctx, almanac)?;

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

    /// This model uses central differencing for gradient computation instead of hyperdual numbers.
    /// This is required given the complexity of the NRLMSISE00 model.
    fn gradient(
        &self,
        ctx: &Spacecraft,
        almanac: &Almanac,
    ) -> Result<(Vector3<f64>, Matrix4x3<f64>), DynamicsError> {
        let dx = self.eom(ctx, almanac)?;

        let mut grad = Matrix4x3::zeros();

        // Central differencing: 6 EOM evaluations, O(h^2) error
        for j in 0..3 {
            // Optimal step size for central diff: h ~ eps^(1/3) * |r|
            let h = 6.0e-6 * ctx.orbit.radius_km[j].abs().max(1.0);

            let mut ctx_plus = *ctx;
            ctx_plus.orbit.radius_km[j] += h;
            let f_plus = self.eom(&ctx_plus, almanac)?;

            let mut ctx_minus = *ctx;
            ctx_minus.orbit.radius_km[j] -= h;
            let f_minus = self.eom(&ctx_minus, almanac)?;

            let df_dr = (f_plus - f_minus) / (2.0 * h);
            for i in 0..3 {
                grad[(i, j)] = df_dr[i];
            }
        }

        // Partial wrt Cd: drag acceleration is exactly linear in coeff_drag,
        // so this is computed analytically rather than by finite differencing.
        // (This is d(accel)/d(Cd), not d(Cd)/d(Cd) -- the latter is the separate,
        // legitimately-zero term for Cd's own dynamics under a constant model.)
        let wrt_cd = dx / ctx.drag.coeff_drag;
        for j in 0..3 {
            grad[(3, j)] = wrt_cd[j];
        }

        Ok((dx, grad))
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
