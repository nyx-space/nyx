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

use crate::dynamics::DynamicsError;
pub use crate::io::space_weather::Msise00DailyWeather;
use hifitime::Epoch;

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

/// Compute full NRLMSISE-00 output for the given input parameters.
///
/// Returns temperatures and all species number densities.
fn calculate(input: &Nrlmsise00Input) -> Nrlmsise00Output {
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
pub fn msise00_density(
    sw: Msise00DailyWeather,
    lst_h: f64,
    latitude_deg: f64,
    longitude_deg: f64,
    altitude_km: f64,
    epoch: Epoch,
) -> Result<Nrlmsise00Output, DynamicsError> {
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

    Ok(calculate(&input))
}
