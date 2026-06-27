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
    SolarPressure,
};
use crate::cosmic::{AU, AstroPhysicsSnafu, Frame, SPEED_OF_LIGHT_M_S, Spacecraft};
use crate::linalg::Matrix4x3;
use anise::almanac::Almanac;
use anise::constants::frames::SUN_J2000;
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use serde_dhall::StaticType;
use snafu::ResultExt;
use std::fmt;
use std::sync::Arc;

/// `Albedo` implements the albedo radiation pressure force model.
#[derive(Clone, Debug, Serialize, Deserialize, StaticType)]
pub struct Albedo {
    /// The frame of the planet reflecting the light.
    pub planet_frame: Frame,
    /// Average albedo of the planet (0.0 to 1.0).
    pub average_albedo: f64,
}

impl Albedo {
    /// Earth albedo model with average value of 0.306
    pub fn earth(almanac: Arc<Almanac>) -> Result<Arc<Self>, DynamicsError> {
        use anise::constants::frames::IAU_EARTH_FRAME;
        Ok(Arc::new(Self {
            planet_frame: almanac.frame_info(IAU_EARTH_FRAME).context(DynamicsPlanetarySnafu {
                action: "fetching Earth frame",
            })?,
            average_albedo: 0.306,
        }))
    }

    /// Moon albedo model with average value of 0.11
    pub fn moon(almanac: Arc<Almanac>) -> Result<Arc<Self>, DynamicsError> {
        use anise::constants::frames::IAU_MOON_FRAME;
        Ok(Arc::new(Self {
            planet_frame: almanac.frame_info(IAU_MOON_FRAME).context(DynamicsPlanetarySnafu {
                action: "fetching Moon frame",
            })?,
            average_albedo: 0.11,
        }))
    }

    /// Mars albedo model with average value of 0.25
    pub fn mars(almanac: Arc<Almanac>) -> Result<Arc<Self>, DynamicsError> {
        use anise::constants::frames::IAU_MARS_FRAME;
        Ok(Arc::new(Self {
            planet_frame: almanac.frame_info(IAU_MARS_FRAME).context(DynamicsPlanetarySnafu {
                action: "fetching Mars frame",
            })?,
            average_albedo: 0.25,
        }))
    }

    /// Venus albedo model with average value of 0.75
    pub fn venus(almanac: Arc<Almanac>) -> Result<Arc<Self>, DynamicsError> {
        use anise::constants::frames::IAU_VENUS_FRAME;
        Ok(Arc::new(Self {
            planet_frame: almanac.frame_info(IAU_VENUS_FRAME).context(DynamicsPlanetarySnafu {
                action: "fetching Venus frame",
            })?,
            average_albedo: 0.75,
        }))
    }
}

impl ForceModel for Albedo {
    fn estimation_index(&self) -> Option<usize> {
        // Albedo Cr is at index 9
        Some(9)
    }

    fn eom(&self, ctx: &Spacecraft, almanac: Arc<Almanac>) -> Result<Vector3<f64>, DynamicsError> {
        // 1. Get position of Sun and Spacecraft relative to the planet
        let sc_state_in_planet = almanac
            .transform_to(ctx.orbit, self.planet_frame, None)
            .context(DynamicsAlmanacSnafu {
                action: "transforming SC state to planet frame",
            })?;
        let r_sc = sc_state_in_planet.radius_km;
        let r_sc_mag = r_sc.norm();
        let r_sc_unit = r_sc / r_sc_mag;

        let sun_state_in_planet = almanac
            .transform(SUN_J2000, self.planet_frame, ctx.orbit.epoch, None)
            .context(DynamicsAlmanacSnafu {
                action: "fetching Sun position in planet frame",
            })?;
        let r_sun = sun_state_in_planet.radius_km;
        let r_sun_mag = r_sun.norm();
        let r_sun_unit = r_sun / r_sun_mag;

        // 2. Compute Sun elevation at nadir (angle Sun-Planet-Probe)
        // cos_theta is the cosine of the angle between Planet-Sun and Planet-Probe vectors
        let cos_theta = r_sun_unit.dot(&r_sc_unit);

        // Requirement: Sun-Planet-Probe angle should be < 90 degrees (noon is 0), else nadir is in the dark
        if cos_theta <= 0.0 {
            return Ok(Vector3::zeros());
        }

        // 3. Fetch albedo value
        let albedo = self.average_albedo;

        // 4. Compute reflected flux
        // Incident solar flux at the planet
        let r_sun_au = r_sun_mag / AU;
        let phi_inc = crate::dynamics::solarpressure::SOLAR_FLUX_W_m2 / r_sun_au.powi(2);

        // Planet radius
        let planet_radius = self
            .planet_frame
            .mean_equatorial_radius_km()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;

        // Reflected pressure (Lambertian model)
        // P_refl = (Phi_inc / c) * albedo * cos_theta * (R/r)^2
        let flux_pressure = (phi_inc / SPEED_OF_LIGHT_M_S)
            * albedo
            * cos_theta
            * (planet_radius / r_sc_mag).powi(2);

        // 5. Compute attenuation (placeholder for now)
        let attenuation = 1.0;

        // 6. Compute force using refactored SRP helper
        // Direction is upward (away from planet)
        Ok(SolarPressure::compute_force(
            flux_pressure * attenuation,
            ctx.albedo.area_m2,
            ctx.albedo.coeff_reflectivity,
            r_sc_unit,
        ))
    }

    fn gradient(
        &self,
        _osc_ctx: &Spacecraft,
        _almanac: Arc<Almanac>,
    ) -> Result<(Vector3<f64>, Matrix4x3<f64>), DynamicsError> {
        // Partials for albedo are complex and not yet implemented
        Err(DynamicsError::DynamicsAstro {
            source: crate::cosmic::AstroError::PartialsUndefined,
        })
    }
}

impl fmt::Display for Albedo {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Albedo for {} with average = {}",
            self.planet_frame, self.average_albedo
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cosmic::{Epoch, Orbit};
    use anise::constants::frames::{EARTH_J2000, IAU_EARTH_FRAME};
    use std::path::PathBuf;
    use std::str::FromStr;

    #[test]
    fn test_albedo_earth() {
        let data_folder: PathBuf = [env!("CARGO_MANIFEST_DIR"), "../data/01_planetary"]
            .iter()
            .collect();
        let mut almanac = Almanac::default();
        // Load kernels
        almanac = almanac
            .load(data_folder.join("de440s.bsp").to_str().unwrap())
            .unwrap();
        almanac = almanac
            .load(data_folder.join("pck08.pca").to_str().unwrap())
            .unwrap();
        let almanac = Arc::new(almanac);

        // Noon over Earth (Sun-Earth-Probe angle ~ 0)
        let epoch = Epoch::from_str("2024-03-20T12:00:00 UTC").unwrap();
        let frame = Frame::from(EARTH_J2000);

        // Transform Sun to Earth frame at this epoch to find its direction
        let sun_state = almanac.transform(SUN_J2000, IAU_EARTH_FRAME, epoch, None).unwrap();
        let sun_unit = sun_state.radius_km / sun_state.radius_km.norm();

        // Place spacecraft in same direction as Sun (noon)
        let sc_pos = sun_unit * 7000.0;
        let sc_orbit = Orbit::cartesian(sc_pos.x, sc_pos.y, sc_pos.z, 0.0, 7.5, 0.0, epoch, IAU_EARTH_FRAME);
        let sc = Spacecraft::from_srp_defaults(sc_orbit, 1000.0, 10.0).with_albedo(10.0, 1.0);

        let albedo_model = Albedo::earth(almanac.clone()).unwrap();
        let acc = albedo_model.eom(&sc, almanac.clone()).unwrap();

        println!("Albedo acceleration at noon: {:?}", acc);
        assert!(acc.norm() > 0.0);
        // Upward force: direction should be similar to sc_unit
        let sc_unit = sc_pos / sc_pos.norm();
        let acc_unit = acc / acc.norm();
        assert!(acc_unit.dot(&sc_unit) > 0.999);

        // Test midnight (nadir in dark)
        let sc_midnight_pos = -sun_unit * 7000.0;
        let sc_midnight_orbit = Orbit::cartesian(sc_midnight_pos.x, sc_midnight_pos.y, sc_midnight_pos.z, 0.0, 7.5, 0.0, epoch, IAU_EARTH_FRAME);
        let sc_midnight = sc.with_orbit(sc_midnight_orbit);
        let acc_midnight = albedo_model.eom(&sc_midnight, almanac.clone()).unwrap();
        assert_eq!(acc_midnight.norm(), 0.0);
    }
}
