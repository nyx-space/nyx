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

pub mod coefficients;
pub mod geo;
mod model;

use hifitime::Epoch;
use snafu::ResultExt;
use std::collections::BTreeMap;

/// Full output of the NRLMSISE-00 model.
///
/// Includes temperatures and all species number densities.
#[derive(Debug, Clone)]
pub struct Nrlmsise00Output {
    /// Exospheric temperature [K].
    pub temp_exo: f64,
    /// Temperature at altitude [K].
    pub temp_alt: f64,
    /// He number density [cm⁻³].
    pub density_he: f64,
    /// O number density [cm⁻³].
    pub density_o: f64,
    /// N₂ number density [cm⁻³].
    pub density_n2: f64,
    /// O₂ number density [cm⁻³].
    pub density_o2: f64,
    /// Ar number density [cm⁻³].
    pub density_ar: f64,
    /// H number density [cm⁻³].
    pub density_h: f64,
    /// N number density [cm⁻³].
    pub density_n: f64,
    /// Anomalous oxygen number density [cm⁻³].
    pub density_anomalous_o: f64,
    /// Total mass density [kg/m³].
    pub total_mass_density: f64,
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

/// Daily space weather data.
#[derive(Debug, Clone, Copy, Default)]
pub struct DailySpaceWeather {
    /// Daily F10.7 [SFU].
    pub f107_daily: f64,
    /// 81-day centered average F10.7 [SFU].
    pub f107_avg: f64,
    /// Daily Ap index.
    pub ap_daily: f64,
    /// 7-element Ap array for magnetic activity variations.
    pub ap_3hour_history: [f64; 7],
}

#[derive(Clone, Debug)]
pub struct Nrlmsise00 {
    pub weather: BTreeMap<Epoch, DailySpaceWeather>,
}

impl Nrlmsise00 {
    /// Create a new NRLMSISE-00 model with the given space weather provider.
    pub fn new(weather: BTreeMap<Epoch, DailySpaceWeather>) -> Self {
        Self { weather }
    }

    fn get_weather(&self, epoch: Epoch) -> DailySpaceWeather {
        if self.weather.is_empty() {
            return DailySpaceWeather::default();
        }
        let mut prev = self.weather.keys().next().unwrap();
        for k in self.weather.keys() {
            if *k > epoch {
                break;
            }
            prev = k;
        }
        self.weather[prev]
    }

    /// Compute full NRLMSISE-00 output for the given input parameters.
    ///
    /// Returns temperatures and all species number densities.
    pub fn calculate(&self, input: &Nrlmsise00Input) -> Nrlmsise00Output {
        let (d, temp_exo, temp_alt) = model::compute(input);
        // d[0..8]: He, O, N2, O2, Ar, total_mass(g/cm³), H, N, anomO
        // Total mass density: d[5] is in g/cm³, convert to kg/m³ (* 1000)
        Nrlmsise00Output {
            temp_exo,
            temp_alt,
            density_he: d[0],
            density_o: d[1],
            density_n2: d[2],
            density_o2: d[3],
            density_ar: d[4],
            density_h: d[6],
            density_n: d[7],
            density_anomalous_o: d[8],
            total_mass_density: d[5] * 1000.0,
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
        latitude_deg: f64,
        longitude_deg: f64,
        altitude_km: f64,
        epoch: Epoch,
    ) -> Nrlmsise00Output {
        let (doy, ut_seconds) = geo::epoch_to_day_of_year_and_ut(&epoch);
        let lst = geo::local_solar_time(ut_seconds, longitude_deg, &epoch);
        let sw = self.get_weather(epoch);

        let input = Nrlmsise00Input {
            day_of_year: doy,
            ut_seconds,
            altitude_km,
            latitude_deg,
            longitude_deg,
            local_solar_time_hours: lst,
            f107_daily: sw.f107_daily,
            f107_avg: sw.f107_avg,
            ap_daily: sw.ap_daily,
            ap_array: sw.ap_3hour_history,
        };

        self.calculate(&input)
    }
}

use crate::cosmic::AstroError;
use crate::cosmic::AstroPhysicsSnafu;
use crate::cosmic::Spacecraft;
use crate::dynamics::{DynamicsAlmanacSnafu, DynamicsAstroSnafu, DynamicsError, ForceModel};
use crate::linalg::{Matrix4x3, Vector3};
use anise::almanac::Almanac;
use anise::constants::frames::IAU_EARTH_FRAME;
use core::fmt;

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
        let earth_frame = almanac.frame_info(IAU_EARTH_FRAME).unwrap();

        let osc_drag_frame =
            almanac
                .transform_to(ctx.orbit, earth_frame, None)
                .context(DynamicsAlmanacSnafu {
                    action: "transforming into drag frame",
                })?;

        let lat_deg = osc_drag_frame.latitude_deg().unwrap();
        let lon_deg = osc_drag_frame.longitude_deg();
        let r0 = earth_frame
            .mean_equatorial_radius_km()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;
        let alt_km = osc_drag_frame.rmag_km() - r0;

        let out = self.density_with_composition(lat_deg, lon_deg, alt_km, ctx.orbit.epoch);
        let rho = out.total_mass_density;

        let velocity_integr_frame = almanac
            .transform_to(osc_drag_frame, integration_frame, None)
            .context(DynamicsAlmanacSnafu {
                action: "rotating into the integration frame",
            })?
            .velocity_km_s;

        let velocity = velocity_integr_frame - osc_drag_frame.velocity_km_s;

        Ok(-0.5 * 1e3 * rho * ctx.drag.coeff_drag * ctx.drag.area_m2 * velocity.norm() * velocity)
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
