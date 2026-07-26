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

//! NRLMSISE-00 empirical atmosphere model.
//!
//! Clean-room implementation based on the following references:
//! - Picone, J.M. et al. (2002), "NRLMSISE-00 empirical model of the atmosphere:
//!   Statistical comparisons and scientific issues", J. Geophys. Res., 107(A12), 1468,
//!   doi:10.1029/2002JA009430
//! - Hedin, A.E. (1991), "Extension of the MSIS thermosphere model into the middle
//!   and lower atmosphere", J. Geophys. Res., 96(A2), 1159-1172.
//! - Hedin, A.E. (1987), "MSIS-86 thermospheric model",
//!   J. Geophys. Res., 92(A5), 4649-4662.
//!
//! Coefficient values are from the official NRL distribution, treated as published data.
//! NRLMSISE-00 is believed to be in the public domain as a U.S. Government work
//! (17 U.S.C. § 105), though no explicit license was provided by NRL.
//!
//! Note: MSIS is a registered trademark. This module uses the name "NRLMSISE-00"
//! for nominative fair use (identifying compatibility with the NRL model).
//!
//! Validated against `pymsis` (official NRL Fortran wrapper, `version=0`).

use crate::cosmic::{AstroError, AstroPhysicsSnafu, Spacecraft};
use crate::dynamics::{DynamicsAlmanacSnafu, DynamicsAstroSnafu, DynamicsError, ForceModel};
pub use crate::io::space_weather::Msise00DailyWeather;
use crate::linalg::{Matrix4x3, Vector3};
use anise::almanac::Almanac;
use anise::frames::Frame;
use core::fmt;
use hifitime::{Epoch, Unit};
use snafu::ResultExt;
use std::collections::BTreeMap;

pub mod coefficients;
mod model;

/// Full output of the NRLMSISE-00 model.
///
/// Includes temperatures and all species number densities.
#[derive(Debug, Clone)]
pub struct Nrlmsise00Output {
    /// Exospheric temperature [K].
    pub temp_exo_k: f64,
    /// Temperature at altitude [K].
    pub temp_alt_k: f64,
    /// He number density [cm⁻³].
    pub density_he_per_cm3: f64,
    /// O number density [cm⁻³].
    pub density_o_per_cm3: f64,
    /// N₂ number density [cm⁻³].
    pub density_n2_per_cm3: f64,
    /// O₂ number density [cm⁻³].
    pub density_o2_per_cm3: f64,
    /// Ar number density [cm⁻³].
    pub density_ar_per_cm3: f64,
    /// H number density [cm⁻³].
    pub density_h_per_cm3: f64,
    /// N number density [cm⁻³].
    pub density_n_per_cm3: f64,
    /// Anomalous oxygen number density [cm⁻³].
    pub density_anomalous_o_per_cm3: f64,
    /// Total mass density [kg/m³].
    pub total_mass_density_kg_m3: f64,
}

/// Input parameters for a single NRLMSISE-00 evaluation.
#[derive(Debug, Clone)]
pub struct Nrlmsise00Input {
    /// Day of year [1-366].
    pub day_of_year: u32,
    /// Universal time [seconds since midnight].
    pub ut_seconds: f64,
    /// Geodetic altitude [km].
    pub altitude_km: f64,
    /// Geodetic latitude [degrees, -90 to 90].
    pub latitude_deg: f64,
    /// Geodetic longitude [degrees, 0 to 360 or -180 to 180].
    pub longitude_deg: f64,
    /// Local apparent solar time [hours, 0-24].
    pub local_solar_time_hours: f64,
    /// Previous day's F10.7 [SFU].
    pub f107_daily: f64,
    /// 81-day centered average F10.7 [SFU].
    pub f107_avg: f64,
    /// Daily Ap index.
    pub ap_daily: f64,
    /// 7-element Ap array for magnetic activity variations.
    pub ap_array: [f64; 7],
}

/// NRLMSISE-00 empirical atmosphere model.
///
/// Computes neutral atmospheric density and composition from 0 to ~1000 km altitude
/// as a function of location, time, solar activity (F10.7), and geomagnetic
/// activity (Ap).
///
/// Generic over the space weather provider `P`. Use [`ConstantWeather`] for
/// fixed conditions, or [`CssiSpaceWeather`] for time-varying data.

#[derive(Clone, Debug)]
pub struct Nrlmsise00 {
    /// Frame causing the drag
    pub frame: Frame,
    pub weather: BTreeMap<Epoch, Msise00DailyWeather>,
}

impl Nrlmsise00 {
    /// Create a new NRLMSISE-00 model with the given space weather provider.
    pub fn new(weather: BTreeMap<Epoch, Msise00DailyWeather>, frame: Frame) -> Self {
        Self { weather, frame }
    }

    /// Compute full NRLMSISE-00 output for the given input parameters.
    ///
    /// Returns temperatures and all species number densities.
    pub fn calculate(&self, input: &Nrlmsise00Input) -> Nrlmsise00Output {
        let (d, temp_exo, temp_alt) = model::compute(input);
        // d[0..8]: He, O, N2, O2, Ar, total_mass(g/cm³), H, N, anomO
        // Total mass density: d[5] is in g/cm³, convert to kg/m³ (* 1000)
        Nrlmsise00Output {
            temp_exo_k: temp_exo,
            temp_alt_k: temp_alt,
            density_he_per_cm3: d[0],
            density_o_per_cm3: d[1],
            density_n2_per_cm3: d[2],
            density_o2_per_cm3: d[3],
            density_ar_per_cm3: d[4],
            density_h_per_cm3: d[6],
            density_n_per_cm3: d[7],
            density_anomalous_o_per_cm3: d[8],
            // Convert g/cm^3 to kg/m^3: g->kg <=> 1e-3; cm^3 -> m^3 <-> 1e6 => 1e3
            total_mass_density_kg_m3: d[5] * 1e3,
        }
    }

    /// Compute full atmospheric composition from geodetic coordinates and epoch.
    ///
    /// Returns the complete NRLMSISE-00 output including:
    /// - Total mass density \[kg/m³\]
    /// - Number densities \[cm⁻³\] for 9 species: He, O, N₂, O₂, Ar, H, N, anomalous O
    /// - Exospheric and local temperatures \[K\]
    ///
    /// This is the high-level API that takes pre-computed geodetic coordinates.
    /// For direct low-level access with explicit NRLMSISE-00 input parameters,
    /// use [`Nrlmsise00::calculate()`].
    pub fn density_with_composition(
        &self,
        lst_h: f64,
        latitude_deg: f64,
        longitude_deg: f64,
        altitude_km: f64,
        epoch: Epoch,
    ) -> Result<Nrlmsise00Output, DynamicsError> {
        // THis is an O(log N) look up
        let sw = self
            .weather
            .range(..=epoch)
            .next_back()
            .map(|(_, v)| v)
            .or_else(|| self.weather.values().next())
            .copied()
            .unwrap_or_default();

        let at_midnight = epoch.with_hms(0, 0, 0);
        let ut_seconds = (epoch - at_midnight).to_seconds();

        let input = Nrlmsise00Input {
            day_of_year: at_midnight.day_of_year() as u32,
            ut_seconds,
            altitude_km,
            latitude_deg,
            longitude_deg,
            local_solar_time_hours: lst_h,
            f107_daily: sw.f107_daily_sfu,
            f107_avg: sw.f107_avg_sfu,
            ap_daily: sw.ap_daily,
            ap_array: sw.ap_3hour_history,
        };

        Ok(self.calculate(&input))
    }
}

impl fmt::Display for Nrlmsise00 {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "NRLMSISE-00 drag model")
    }
}

impl ForceModel for Nrlmsise00 {
    fn estimation_index(&self) -> Option<usize> {
        None
    }

    fn eom(&self, ctx: &Spacecraft, almanac: &Almanac) -> Result<Vector3<f64>, DynamicsError> {
        let integration_frame = ctx.orbit.frame;
        let earth_frame = almanac.frame_info(self.frame).unwrap();

        let osc_drag_frame =
            almanac
                .transform_to(ctx.orbit, earth_frame, None)
                .context(DynamicsAlmanacSnafu {
                    action: "transforming into drag frame",
                })?;

        let lat_deg = osc_drag_frame.latitude_deg().unwrap();
        let lon_deg = osc_drag_frame.longitude_deg();

        let alt_km = osc_drag_frame
            .altitude_km()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;

        let lst_h =
            almanac
                .local_solar_time(osc_drag_frame, None)
                .context(DynamicsAlmanacSnafu {
                    action: "computing local solar time",
                })?;

        let out = self.density_with_composition(
            lst_h.to_unit(Unit::Hour),
            lat_deg,
            lon_deg,
            alt_km,
            ctx.orbit.epoch,
        )?;

        let rho_kg_m3 = out.total_mass_density_kg_m3;

        let velocity_integr_frame = almanac
            .transform_to(osc_drag_frame, integration_frame, None)
            .context(DynamicsAlmanacSnafu {
                action: "rotating into the integration frame",
            })?
            .velocity_km_s;

        let v_km_s = velocity_integr_frame - osc_drag_frame.velocity_km_s;

        // Finally, apply the drag model.
        Ok(
            -0.5 * 1e3
                * rho_kg_m3
                * ctx.drag.coeff_drag
                * ctx.drag.area_m2
                * v_km_s.norm()
                * v_km_s,
        )
    }

    fn gradient(
        &self,
        _osc_ctx: &Spacecraft,
        _almanac: &Almanac,
    ) -> Result<(Vector3<f64>, Matrix4x3<f64>), DynamicsError> {
        Err(DynamicsError::DynamicsAstro {
            source: AstroError::PartialsUndefined,
        })
    }
}
