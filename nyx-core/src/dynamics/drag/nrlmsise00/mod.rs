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
use serde::{Deserialize, Serialize};
use serde_dhall::StaticType;

#[cfg(feature = "python")]
use pyo3::prelude::*;

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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, StaticType)]
#[cfg_attr(feature = "python", pyclass(from_py_object, eq, eq_int))]
pub enum GeomagneticMode {
    Off,
    StandardDailyAp,
    ExtendedHistory57h, // Sets sw[9] = -1.0
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, StaticType)]
#[cfg_attr(feature = "python", pyclass(from_py_object, get_all, set_all))]
pub struct Nrlmsise00Flags {
    pub geomagnetic: GeomagneticMode,
    pub f107_solar_flux: bool,
    pub time_independent: bool,
    pub annual_harmonics: bool,
    pub semiannual_harmonics: bool,
    pub diurnal_tides: bool,
    pub semidiurnal_tides: bool,
    pub terdiurnal_tides: bool,
    pub ut_and_longitude: bool,
    pub exospheric_temp_variations: bool,
    pub lower_boundary_temp_variations: bool,
    pub gradient_variations: bool,
    pub departures_from_diffusive_equilibrium: bool,
    pub lower_thermosphere_temp_variations: bool,
    pub upper_stratosphere_temp_variations: bool,
    pub boundary_density_variations: bool,
    pub lower_mesosphere_temp_variations: bool,
    pub turbopause_scale_height_variations: bool,
}

impl Default for Nrlmsise00Flags {
    fn default() -> Self {
        Self {
            geomagnetic: GeomagneticMode::StandardDailyAp,
            f107_solar_flux: true,
            time_independent: true,
            annual_harmonics: true,
            semiannual_harmonics: true,
            diurnal_tides: true,
            semidiurnal_tides: true,
            terdiurnal_tides: true,
            ut_and_longitude: true,
            exospheric_temp_variations: true,
            lower_boundary_temp_variations: true,
            gradient_variations: true,
            departures_from_diffusive_equilibrium: true,
            lower_thermosphere_temp_variations: true,
            upper_stratosphere_temp_variations: true,
            boundary_density_variations: true,
            lower_mesosphere_temp_variations: true,
            turbopause_scale_height_variations: true,
        }
    }
}

impl Nrlmsise00Flags {
    /// Compiles the high-level flags into the raw 24-element float array consumed by the kernel.
    pub(crate) fn to_switches(self) -> [f64; 24] {
        let mut sw = [1.0f64; 24];

        // NOTE Unit selection is ALWAYS set to 1.0 because the calculation code does not even check it.

        // Geomagnetic storm mode
        sw[9] = match self.geomagnetic {
            GeomagneticMode::Off => 0.0,
            GeomagneticMode::StandardDailyAp => 1.0,
            GeomagneticMode::ExtendedHistory57h => -1.0,
        };

        // Specific feature toggles
        if !self.f107_solar_flux {
            sw[1] = 0.0;
        }
        if !self.time_independent {
            sw[2] = 0.0;
        }
        if !self.annual_harmonics {
            sw[3] = 0.0;
            sw[5] = 0.0;
        }
        if !self.semiannual_harmonics {
            sw[4] = 0.0;
            sw[6] = 0.0;
        }
        if !self.diurnal_tides {
            sw[7] = 0.0;
        }
        if !self.semidiurnal_tides {
            sw[8] = 0.0;
        }
        if !self.terdiurnal_tides {
            sw[14] = 0.0;
        }
        if !self.ut_and_longitude {
            sw[10] = 0.0;
            sw[11] = 0.0;
            sw[12] = 0.0;
            sw[13] = 0.0;
        }
        if !self.exospheric_temp_variations {
            sw[16] = 0.0;
        }
        if !self.lower_boundary_temp_variations {
            sw[17] = 0.0;
        }
        if !self.gradient_variations {
            sw[19] = 0.0;
        }
        if !self.departures_from_diffusive_equilibrium {
            sw[15] = 0.0;
        }
        if !self.lower_thermosphere_temp_variations {
            sw[18] = 0.0;
        }
        if !self.upper_stratosphere_temp_variations {
            sw[20] = 0.0;
        }
        if !self.boundary_density_variations {
            sw[21] = 0.0;
        }
        if !self.lower_mesosphere_temp_variations {
            sw[22] = 0.0;
        }
        if !self.turbopause_scale_height_variations {
            sw[23] = 0.0;
        }

        sw
    }
}

#[cfg(feature = "python")]
#[pymethods]
impl Nrlmsise00Flags {
    #[new]
    #[pyo3(signature = (
        geomagnetic = None,
        f107_solar_flux = true,
        time_independent = true,
        annual_harmonics = true,
        semiannual_harmonics = true,
        diurnal_tides = true,
        semidiurnal_tides = true,
        terdiurnal_tides = true,
        ut_and_longitude = true,
        exospheric_temp_variations = true,
        lower_boundary_temp_variations = true,
        gradient_variations = true,
        departures_from_diffusive_equilibrium = true,
        lower_thermosphere_temp_variations = true,
        upper_stratosphere_temp_variations = true,
        boundary_density_variations = true,
        lower_mesosphere_temp_variations = true,
        turbopause_scale_height_variations = true,
    ))]
    #[allow(clippy::too_many_arguments)]
    fn py_new(
        geomagnetic: Option<GeomagneticMode>,
        f107_solar_flux: bool,
        time_independent: bool,
        annual_harmonics: bool,
        semiannual_harmonics: bool,
        diurnal_tides: bool,
        semidiurnal_tides: bool,
        terdiurnal_tides: bool,
        ut_and_longitude: bool,
        exospheric_temp_variations: bool,
        lower_boundary_temp_variations: bool,
        gradient_variations: bool,
        departures_from_diffusive_equilibrium: bool,
        lower_thermosphere_temp_variations: bool,
        upper_stratosphere_temp_variations: bool,
        boundary_density_variations: bool,
        lower_mesosphere_temp_variations: bool,
        turbopause_scale_height_variations: bool,
    ) -> Self {
        Self {
            geomagnetic: geomagnetic.unwrap_or(GeomagneticMode::StandardDailyAp),
            f107_solar_flux,
            time_independent,
            annual_harmonics,
            semiannual_harmonics,
            diurnal_tides,
            semidiurnal_tides,
            terdiurnal_tides,
            ut_and_longitude,
            exospheric_temp_variations,
            lower_boundary_temp_variations,
            gradient_variations,
            departures_from_diffusive_equilibrium,
            lower_thermosphere_temp_variations,
            upper_stratosphere_temp_variations,
            boundary_density_variations,
            lower_mesosphere_temp_variations,
            turbopause_scale_height_variations,
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self)
    }

    fn __str__(&self) -> String {
        format!("{:?} @ {self:p}", self)
    }
}

/// Compute full NRLMSISE-00 output for the given input parameters.
///
/// Returns temperatures and all species number densities.
fn calculate(input: &Nrlmsise00Input, flags: Nrlmsise00Flags) -> Nrlmsise00Output {
    let sw = flags.to_switches();
    let (d, temp_exo, temp_alt) = model::compute(input, &sw);
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
    flags: Nrlmsise00Flags,
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

    Ok(calculate(&input, flags))
}

#[cfg(test)]
mod tests {
    use super::*;
    // use crate::dynamics::nrlmsise00::{compute, Nrlmsise00Input};
    use approx::assert_relative_eq;
    use serde::Deserialize;
    use std::fs;
    use std::path::PathBuf;

    #[test]
    fn test_nrlmsise00_flags() {
        let default_flags = Nrlmsise00Flags::default();
        let default_switches = default_flags.to_switches();

        // Spot-check standard switches mapping
        assert_eq!(default_switches[0], 1.0);
        assert_eq!(default_switches[9], 1.0);
        assert_eq!(default_switches[1], 1.0);
        assert_eq!(default_switches[2], 1.0);

        // Customize some flags
        let mut custom_flags = Nrlmsise00Flags {
            geomagnetic: GeomagneticMode::StandardDailyAp,
            f107_solar_flux: false,
            time_independent: false,
            annual_harmonics: false,
            semiannual_harmonics: false,
            diurnal_tides: false,
            semidiurnal_tides: false,
            terdiurnal_tides: false,
            ut_and_longitude: false,
            exospheric_temp_variations: false,
            lower_boundary_temp_variations: false,
            gradient_variations: false,
            departures_from_diffusive_equilibrium: false,
            lower_thermosphere_temp_variations: false,
            upper_stratosphere_temp_variations: false,
            boundary_density_variations: false,
            lower_mesosphere_temp_variations: false,
            turbopause_scale_height_variations: false,
        };

        let custom_switches = custom_flags.to_switches();
        assert_eq!(custom_switches[0], 1.0); // Always 1.0
        assert_eq!(custom_switches[9], 1.0); // GeomagneticMode::StandardDailyAp -> 1.0
        assert_eq!(custom_switches[1], 0.0);
        assert_eq!(custom_switches[2], 0.0);
        assert_eq!(custom_switches[3], 0.0);
        assert_eq!(custom_switches[5], 0.0);
        assert_eq!(custom_switches[4], 0.0);
        assert_eq!(custom_switches[6], 0.0);
        assert_eq!(custom_switches[7], 0.0);
        assert_eq!(custom_switches[8], 0.0);
        assert_eq!(custom_switches[14], 0.0);
        assert_eq!(custom_switches[10], 0.0);
        assert_eq!(custom_switches[11], 0.0);
        assert_eq!(custom_switches[12], 0.0);
        assert_eq!(custom_switches[13], 0.0);
        assert_eq!(custom_switches[16], 0.0);
        assert_eq!(custom_switches[17], 0.0);
        assert_eq!(custom_switches[19], 0.0);
        assert_eq!(custom_switches[15], 0.0);
        assert_eq!(custom_switches[18], 0.0);
        assert_eq!(custom_switches[20], 0.0);
        assert_eq!(custom_switches[21], 0.0);
        assert_eq!(custom_switches[22], 0.0);
        assert_eq!(custom_switches[23], 0.0);

        // Test GeomagneticMode::Off
        custom_flags.geomagnetic = GeomagneticMode::Off;
        let switches_off = custom_flags.to_switches();
        assert_eq!(switches_off[9], 0.0);
    }

    #[derive(Deserialize)]
    struct MsisTestCase {
        altitude_km: f64,
        latitude_deg: f64,
        longitude_deg: f64,
        day_of_year: u32,
        ut_seconds: f64,
        f107_daily: f64,
        f107_avg: f64,
        ap_array: [f64; 7],
        is_storm: bool,
        expected_total_density_kg_m3: f64,
        expected_temperature_k: f64,
    }

    #[test]
    fn pymsis_validation() {
        let test_data: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "../data/03_tests/nrlmsise00_validation.json",
        ]
        .iter()
        .collect();

        let data = fs::read_to_string(test_data).expect("Failed to read validation JSON");
        let test_cases: Vec<MsisTestCase> =
            serde_json::from_str(&data).expect("Failed to deserialize test cases");

        for (i, tc) in test_cases.iter().enumerate() {
            let input = Nrlmsise00Input {
                day_of_year: tc.day_of_year,
                ut_seconds: tc.ut_seconds,
                altitude_km: tc.altitude_km,
                latitude_deg: tc.latitude_deg,
                longitude_deg: tc.longitude_deg,
                local_solar_time_hours: (tc.ut_seconds / 3600.0 + tc.longitude_deg / 15.0)
                    .rem_euclid(24.0),
                f107_daily: tc.f107_daily,
                f107_avg: tc.f107_avg,
                ap_daily: tc.ap_array[0],
                ap_array: tc.ap_array,
            };

            // Standard switches: 57-hour history enabled
            let sw = Nrlmsise00Flags {
                geomagnetic: if tc.is_storm {
                    GeomagneticMode::ExtendedHistory57h
                } else {
                    GeomagneticMode::StandardDailyAp
                },
                ..Default::default()
            };

            let output = calculate(&input, sw);

            let total_density_kg_m3 = output.total_mass_density_kg_m3;
            let t_alt = output.temp_alt_k;

            // Verify temperature at altitude with a 0.5% relative tolerance
            println!(
                "[storm={}] Temperature #{i}: alt={} km, lat={} deg. Rust: {t_alt}, Pymsis: {}",
                tc.is_storm, tc.altitude_km, tc.latitude_deg, tc.expected_temperature_k
            );
            assert_relative_eq!(
                t_alt,
                tc.expected_temperature_k,
                max_relative = 0.005,
                epsilon = 1e-5
            );

            // Verify total mass density with a 1.0% relative tolerance
            // (Density integration magnifies the floating point differences in the exponential term)
            println!(
                "[storm={}] Density #{i}: alt={} km. Rust: {total_density_kg_m3:.6e}, Pymsis: {:.6e}",
                tc.is_storm, tc.altitude_km, tc.expected_total_density_kg_m3
            );
            assert_relative_eq!(
                total_density_kg_m3,
                tc.expected_total_density_kg_m3,
                max_relative = 0.01,
                epsilon = 1e-18,
            );
        }
    }
}
