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

use anise::astro::PhysicsResult;
use anise::constants::frames::EARTH_J2000;
pub use anise::prelude::Orbit;

pub use anise::structure::spacecraft::{DragData, Mass, SRPData};
use nalgebra::Vector3;
use serde::{Deserialize, Serialize};
use snafu::ResultExt;
use typed_builder::TypedBuilder;

use super::{AstroPhysicsSnafu, BPlane, State};
use crate::dynamics::guidance::Thruster;
use crate::dynamics::DynamicsError;
use crate::errors::{StateAstroSnafu, StateError};
use crate::io::ConfigRepr;
use crate::linalg::{Const, DimName, OMatrix, OVector};
use crate::md::StateParameter;
use crate::time::Epoch;
use crate::utils::{cartesian_to_spherical, spherical_to_cartesian};

use std::default::Default;
use std::fmt;
use std::ops::Add;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]

pub enum GuidanceMode {
    /// Guidance is turned off and Guidance Law may switch mode to Thrust for next call
    Coast,
    /// Guidance is turned on and Guidance Law may switch mode to Coast for next call
    Thrust,
    /// Guidance is turned off and Guidance Law may not change its mode (will need to be done externally to the guidance law).
    Inhibit,
}

impl Default for GuidanceMode {
    fn default() -> Self {
        Self::Coast
    }
}

impl From<f64> for GuidanceMode {
    fn from(value: f64) -> Self {
        if value >= 1.0 {
            Self::Thrust
        } else if value < 0.0 {
            Self::Inhibit
        } else {
            Self::Coast
        }
    }
}

impl From<GuidanceMode> for f64 {
    fn from(mode: GuidanceMode) -> f64 {
        match mode {
            GuidanceMode::Coast => 0.0,
            GuidanceMode::Thrust => 1.0,
            GuidanceMode::Inhibit => -1.0,
        }
    }
}

/// A spacecraft state, composed of its orbit, its masses (dry, prop, extra, all in kg), its SRP configuration, its drag configuration, its thruster configuration, and its guidance mode.
///
/// Optionally, the spacecraft state can also store the state transition matrix from the start of the propagation until the current time (i.e. trajectory STM, not step-size STM).
#[derive(Clone, Copy, Debug, Serialize, Deserialize, TypedBuilder)]
pub struct Spacecraft {
    /// Initial orbit of the vehicle
    pub orbit: Orbit,
    /// Dry, propellant, and extra masses
    #[builder(default)]
    pub mass: Mass,
    /// Solar Radiation Pressure configuration for this spacecraft
    #[builder(default)]
    #[serde(default)]
    pub srp: SRPData,
    #[builder(default)]
    #[serde(default)]
    pub drag: DragData,
    #[builder(default, setter(strip_option))]
    pub thruster: Option<Thruster>,
    /// Any extra information or extension that is needed for specific guidance laws
    #[builder(default)]
    #[serde(default)]
    pub mode: GuidanceMode,
    /// Optionally stores the state transition matrix from the start of the propagation until the current time (i.e. trajectory STM, not step-size STM)
    /// STM is contains position and velocity, Cr, Cd, prop mass
    #[builder(default, setter(strip_option))]
    #[serde(skip)]
    pub stm: Option<OMatrix<f64, Const<9>, Const<9>>>,
}

impl Default for Spacecraft {
    fn default() -> Self {
        Self {
            orbit: Orbit::zero(EARTH_J2000),
            mass: Mass::default(),
            srp: SRPData::default(),
            drag: DragData::default(),
            thruster: None,
            mode: GuidanceMode::default(),
            stm: None,
        }
    }
}

impl From<Orbit> for Spacecraft {
    fn from(orbit: Orbit) -> Self {
        Self::builder().orbit(orbit).build()
    }
}

impl Spacecraft {
    /// Initialize a spacecraft state from all of its parameters
    pub fn new(
        orbit: Orbit,
        dry_mass_kg: f64,
        prop_mass_kg: f64,
        srp_area_m2: f64,
        drag_area_m2: f64,
        coeff_reflectivity: f64,
        coeff_drag: f64,
    ) -> Self {
        Self {
            orbit,
            mass: Mass::from_dry_and_prop_masses(dry_mass_kg, prop_mass_kg),
            srp: SRPData {
                area_m2: srp_area_m2,
                coeff_reflectivity,
            },
            drag: DragData {
                area_m2: drag_area_m2,
                coeff_drag,
            },
            ..Default::default()
        }
    }

    /// Initialize a spacecraft state from only a thruster and mass. Use this when designing guidance laws while ignoring drag and SRP.
    pub fn from_thruster(
        orbit: Orbit,
        dry_mass_kg: f64,
        prop_mass_kg: f64,
        thruster: Thruster,
        mode: GuidanceMode,
    ) -> Self {
        Self {
            orbit,
            mass: Mass::from_dry_and_prop_masses(dry_mass_kg, prop_mass_kg),
            thruster: Some(thruster),
            mode,
            ..Default::default()
        }
    }

    /// Initialize a spacecraft state from the SRP default 1.8 for coefficient of reflectivity (prop mass and drag parameters nullified!)
    pub fn from_srp_defaults(orbit: Orbit, dry_mass_kg: f64, srp_area_m2: f64) -> Self {
        Self {
            orbit,
            mass: Mass::from_dry_mass(dry_mass_kg),
            srp: SRPData::from_area(srp_area_m2),
            ..Default::default()
        }
    }

    /// Initialize a spacecraft state from the SRP default 1.8 for coefficient of drag (prop mass and SRP parameters nullified!)
    pub fn from_drag_defaults(orbit: Orbit, dry_mass_kg: f64, drag_area_m2: f64) -> Self {
        Self {
            orbit,
            mass: Mass::from_dry_mass(dry_mass_kg),
            drag: DragData::from_area(drag_area_m2),
            ..Default::default()
        }
    }

    pub fn with_dv_km_s(mut self, dv_km_s: Vector3<f64>) -> Self {
        self.orbit.apply_dv_km_s(dv_km_s);
        self
    }

    /// Returns a copy of the state with a new dry mass
    pub fn with_dry_mass(mut self, dry_mass_kg: f64) -> Self {
        self.mass.dry_mass_kg = dry_mass_kg;
        self
    }

    /// Returns a copy of the state with a new prop mass
    pub fn with_prop_mass(mut self, prop_mass_kg: f64) -> Self {
        self.mass.prop_mass_kg = prop_mass_kg;
        self
    }

    /// Returns a copy of the state with a new SRP area and CR
    pub fn with_srp(mut self, srp_area_m2: f64, coeff_reflectivity: f64) -> Self {
        self.srp = SRPData {
            area_m2: srp_area_m2,
            coeff_reflectivity,
        };

        self
    }

    /// Returns a copy of the state with a new SRP area
    pub fn with_srp_area(mut self, srp_area_m2: f64) -> Self {
        self.srp.area_m2 = srp_area_m2;
        self
    }

    /// Returns a copy of the state with a new coefficient of reflectivity
    pub fn with_cr(mut self, coeff_reflectivity: f64) -> Self {
        self.srp.coeff_reflectivity = coeff_reflectivity;
        self
    }

    /// Returns a copy of the state with a new drag area and CD
    pub fn with_drag(mut self, drag_area_m2: f64, coeff_drag: f64) -> Self {
        self.drag = DragData {
            area_m2: drag_area_m2,
            coeff_drag,
        };
        self
    }

    /// Returns a copy of the state with a new SRP area
    pub fn with_drag_area(mut self, drag_area_m2: f64) -> Self {
        self.drag.area_m2 = drag_area_m2;
        self
    }

    /// Returns a copy of the state with a new coefficient of drag
    pub fn with_cd(mut self, coeff_drag: f64) -> Self {
        self.drag.coeff_drag = coeff_drag;
        self
    }

    /// Returns a copy of the state with a new orbit
    pub fn with_orbit(mut self, orbit: Orbit) -> Self {
        self.orbit = orbit;
        self
    }

    /// Returns the root sum square error between this spacecraft and the other, in kilometers for the position, kilometers per second in velocity, and kilograms in prop
    pub fn rss(&self, other: &Self) -> PhysicsResult<(f64, f64, f64)> {
        let rss_p_km = self.orbit.rss_radius_km(&other.orbit)?;
        let rss_v_km_s = self.orbit.rss_velocity_km_s(&other.orbit)?;
        let rss_prop_kg = (self.mass.prop_mass_kg - other.mass.prop_mass_kg)
            .powi(2)
            .sqrt();

        Ok((rss_p_km, rss_v_km_s, rss_prop_kg))
    }

    /// Sets the STM of this state of identity, which also enables computation of the STM for spacecraft navigation
    pub fn enable_stm(&mut self) {
        self.stm = Some(OMatrix::<f64, Const<9>, Const<9>>::identity());
    }

    /// Returns the total mass in kilograms
    pub fn mass_kg(&self) -> f64 {
        self.mass.total_mass_kg()
    }

    /// Returns a copy of the state with the provided guidance mode
    pub fn with_guidance_mode(mut self, mode: GuidanceMode) -> Self {
        self.mode = mode;
        self
    }

    pub fn mode(&self) -> GuidanceMode {
        self.mode
    }

    pub fn mut_mode(&mut self, mode: GuidanceMode) {
        self.mode = mode;
    }
}

impl PartialEq for Spacecraft {
    fn eq(&self, other: &Self) -> bool {
        let mass_tol = 1e-6; // milligram
        self.orbit.eq_within(&other.orbit, 1e-9, 1e-12)
            && (self.mass - other.mass).abs().total_mass_kg() < mass_tol
            && self.srp == other.srp
            && self.drag == other.drag
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::Display for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*}", mass_prec, self.mass.total_mass_kg()),
            format!("{:.*}", orbit_prec, self.orbit),
            self.mode,
        )
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::LowerExp for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*e}", mass_prec, self.mass.total_mass_kg()),
            format!("{:.*e}", orbit_prec, self.orbit),
            self.mode,
        )
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::LowerHex for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*}", mass_prec, self.mass.total_mass_kg()),
            format!("{:.*x}", orbit_prec, self.orbit),
            self.mode,
        )
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::UpperHex for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*e}", mass_prec, self.mass.total_mass_kg()),
            format!("{:.*X}", orbit_prec, self.orbit),
            self.mode,
        )
    }
}

impl State for Spacecraft {
    type Size = Const<9>;
    type VecLength = Const<90>;

    /// Copies the current state but sets the STM to identity
    fn with_stm(mut self) -> Self {
        self.enable_stm();
        self
    }

    fn reset_stm(&mut self) {
        self.stm = Some(OMatrix::<f64, Const<9>, Const<9>>::identity());
    }

    fn zeros() -> Self {
        Self::default()
    }

    /// The vector is organized as such:
    /// [X, Y, Z, Vx, Vy, Vz, Cr, Cd, Fuel mass, STM(9x9)]
    fn to_vector(&self) -> OVector<f64, Const<90>> {
        let mut vector = OVector::<f64, Const<90>>::zeros();
        // Set the orbit state info
        for (i, val) in self.orbit.radius_km.iter().enumerate() {
            // Place the orbit state first, then skip three (Cr, Cd, Fuel), then copy orbit STM
            vector[i] = *val;
        }
        for (i, val) in self.orbit.velocity_km_s.iter().enumerate() {
            // Place the orbit state first, then skip three (Cr, Cd, Fuel), then copy orbit STM
            vector[i + 3] = *val;
        }
        // Set the spacecraft parameters
        vector[6] = self.srp.coeff_reflectivity;
        vector[7] = self.drag.coeff_drag;
        vector[8] = self.mass.prop_mass_kg;
        // Add the STM to the vector
        if let Some(stm) = self.stm {
            for (idx, stm_val) in stm.as_slice().iter().enumerate() {
                vector[idx + Self::Size::dim()] = *stm_val;
            }
        }
        vector
    }

    /// Vector is expected to be organized as such:
    /// [X, Y, Z, Vx, Vy, Vz, Cr, Cd, Fuel mass, STM(9x9)]
    fn set(&mut self, epoch: Epoch, vector: &OVector<f64, Const<90>>) {
        let sc_state =
            OVector::<f64, Self::Size>::from_column_slice(&vector.as_slice()[..Self::Size::dim()]);

        if self.stm.is_some() {
            let sc_full_stm = OMatrix::<f64, Self::Size, Self::Size>::from_column_slice(
                &vector.as_slice()[Self::Size::dim()..],
            );

            self.stm = Some(sc_full_stm);
        }

        let radius_km = sc_state.fixed_rows::<3>(0).into_owned();
        let vel_km_s = sc_state.fixed_rows::<3>(3).into_owned();
        self.orbit.epoch = epoch;
        self.orbit.radius_km = radius_km;
        self.orbit.velocity_km_s = vel_km_s;
        self.srp.coeff_reflectivity = sc_state[6].clamp(0.0, 2.0);
        self.drag.coeff_drag = sc_state[7];
        self.mass.prop_mass_kg = sc_state[8];
    }

    /// diag(STM) = [X,Y,Z,Vx,Vy,Vz,Cr,Cd,Fuel]
    /// WARNING: Currently the STM assumes that the prop mass is constant at ALL TIMES!
    fn stm(&self) -> Result<OMatrix<f64, Self::Size, Self::Size>, DynamicsError> {
        match self.stm {
            Some(stm) => Ok(stm),
            None => Err(DynamicsError::StateTransitionMatrixUnset),
        }
    }

    fn epoch(&self) -> Epoch {
        self.orbit.epoch
    }

    fn set_epoch(&mut self, epoch: Epoch) {
        self.orbit.epoch = epoch
    }

    fn add(self, other: OVector<f64, Self::Size>) -> Self {
        self + other
    }

    fn value(&self, param: StateParameter) -> Result<f64, StateError> {
        match param {
            StateParameter::Cd => Ok(self.drag.coeff_drag),
            StateParameter::Cr => Ok(self.srp.coeff_reflectivity),
            StateParameter::DryMass => Ok(self.mass.dry_mass_kg),
            StateParameter::PropMass => Ok(self.mass.prop_mass_kg),
            StateParameter::TotalMass => Ok(self.mass.total_mass_kg()),
            StateParameter::Isp => match self.thruster {
                Some(thruster) => Ok(thruster.isp_s),
                None => Err(StateError::NoThrusterAvail),
            },
            StateParameter::Thrust => match self.thruster {
                Some(thruster) => Ok(thruster.thrust_N),
                None => Err(StateError::NoThrusterAvail),
            },
            StateParameter::GuidanceMode => Ok(self.mode.into()),
            StateParameter::ApoapsisRadius => self
                .orbit
                .apoapsis_km()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::AoL => self
                .orbit
                .aol_deg()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::AoP => self
                .orbit
                .aop_deg()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::BdotR => Ok(BPlane::new(self.orbit)
                .context(StateAstroSnafu { param })?
                .b_r_km
                .real()),
            StateParameter::BdotT => Ok(BPlane::new(self.orbit)
                .context(StateAstroSnafu { param })?
                .b_t_km
                .real()),
            StateParameter::BLTOF => Ok(BPlane::new(self.orbit)
                .context(StateAstroSnafu { param })?
                .ltof_s
                .real()),
            StateParameter::C3 => self
                .orbit
                .c3_km2_s2()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::Declination => Ok(self.orbit.declination_deg()),
            StateParameter::EccentricAnomaly => self
                .orbit
                .ea_deg()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::Eccentricity => self
                .orbit
                .ecc()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::Energy => self
                .orbit
                .energy_km2_s2()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::FlightPathAngle => self
                .orbit
                .fpa_deg()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::Height => self
                .orbit
                .height_km()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::Latitude => self
                .orbit
                .latitude_deg()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::Longitude => Ok(self.orbit.longitude_deg()),
            StateParameter::Hmag => self
                .orbit
                .hmag()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::HX => self
                .orbit
                .hx()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::HY => self
                .orbit
                .hy()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::HZ => self
                .orbit
                .hz()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::HyperbolicAnomaly => self
                .orbit
                .hyperbolic_anomaly_deg()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::Inclination => self
                .orbit
                .inc_deg()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::MeanAnomaly => self
                .orbit
                .ma_deg()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::PeriapsisRadius => self
                .orbit
                .periapsis_km()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::Period => Ok(self
                .orbit
                .period()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param })?
                .to_seconds()),
            StateParameter::RightAscension => Ok(self.orbit.right_ascension_deg()),
            StateParameter::RAAN => self
                .orbit
                .raan_deg()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::Rmag => Ok(self.orbit.rmag_km()),
            StateParameter::SemiMinorAxis => self
                .orbit
                .semi_minor_axis_km()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::SemiParameter => self
                .orbit
                .semi_parameter_km()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::SMA => self
                .orbit
                .sma_km()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::TrueAnomaly => self
                .orbit
                .ta_deg()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::TrueLongitude => self
                .orbit
                .tlong_deg()
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param }),
            StateParameter::VelocityDeclination => Ok(self.orbit.velocity_declination_deg()),
            StateParameter::Vmag => Ok(self.orbit.vmag_km_s()),
            StateParameter::X => Ok(self.orbit.radius_km.x),
            StateParameter::Y => Ok(self.orbit.radius_km.y),
            StateParameter::Z => Ok(self.orbit.radius_km.z),
            StateParameter::VX => Ok(self.orbit.velocity_km_s.x),
            StateParameter::VY => Ok(self.orbit.velocity_km_s.y),
            StateParameter::VZ => Ok(self.orbit.velocity_km_s.z),
            _ => Err(StateError::Unavailable { param }),
        }
    }

    fn set_value(&mut self, param: StateParameter, val: f64) -> Result<(), StateError> {
        match param {
            StateParameter::Cd => self.drag.coeff_drag = val,
            StateParameter::Cr => self.srp.coeff_reflectivity = val,
            StateParameter::PropMass => self.mass.prop_mass_kg = val,
            StateParameter::DryMass => self.mass.dry_mass_kg = val,
            StateParameter::Isp => match self.thruster {
                Some(ref mut thruster) => thruster.isp_s = val,
                None => return Err(StateError::NoThrusterAvail),
            },
            StateParameter::Thrust => match self.thruster {
                Some(ref mut thruster) => thruster.thrust_N = val,
                None => return Err(StateError::NoThrusterAvail),
            },
            StateParameter::AoP => self
                .orbit
                .set_aop_deg(val)
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param })?,
            StateParameter::Eccentricity => self
                .orbit
                .set_ecc(val)
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param })?,
            StateParameter::Inclination => self
                .orbit
                .set_inc_deg(val)
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param })?,
            StateParameter::RAAN => self
                .orbit
                .set_raan_deg(val)
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param })?,
            StateParameter::SMA => self
                .orbit
                .set_sma_km(val)
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param })?,
            StateParameter::TrueAnomaly => self
                .orbit
                .set_ta_deg(val)
                .context(AstroPhysicsSnafu)
                .context(StateAstroSnafu { param })?,
            StateParameter::X => self.orbit.radius_km.x = val,
            StateParameter::Y => self.orbit.radius_km.y = val,
            StateParameter::Z => self.orbit.radius_km.z = val,
            StateParameter::Rmag => {
                // Convert the position to spherical coordinates
                let (_, θ, φ) = cartesian_to_spherical(&self.orbit.radius_km);
                // Convert back to cartesian after setting the new range value
                self.orbit.radius_km = spherical_to_cartesian(val, θ, φ);
            }
            StateParameter::VX => self.orbit.velocity_km_s.x = val,
            StateParameter::VY => self.orbit.velocity_km_s.y = val,
            StateParameter::VZ => self.orbit.velocity_km_s.z = val,
            StateParameter::Vmag => {
                // Convert the velocity to spherical coordinates
                let (_, θ, φ) = cartesian_to_spherical(&self.orbit.velocity_km_s);
                // Convert back to cartesian after setting the new range value
                self.orbit.velocity_km_s = spherical_to_cartesian(val, θ, φ);
            }
            _ => return Err(StateError::ReadOnly { param }),
        }
        Ok(())
    }

    fn unset_stm(&mut self) {
        self.stm = None;
    }

    fn orbit(&self) -> Orbit {
        self.orbit
    }

    fn set_orbit(&mut self, orbit: Orbit) {
        self.orbit = orbit;
    }
}

impl Add<OVector<f64, Const<6>>> for Spacecraft {
    type Output = Self;

    /// Adds the provided state deviation to this orbit
    fn add(mut self, other: OVector<f64, Const<6>>) -> Self {
        let radius_km = other.fixed_rows::<3>(0).into_owned();
        let vel_km_s = other.fixed_rows::<3>(3).into_owned();

        self.orbit.radius_km += radius_km;
        self.orbit.velocity_km_s += vel_km_s;

        self
    }
}

impl Add<OVector<f64, Const<9>>> for Spacecraft {
    type Output = Self;

    /// Adds the provided state deviation to this orbit
    fn add(mut self, other: OVector<f64, Const<9>>) -> Self {
        let radius_km = other.fixed_rows::<3>(0).into_owned();
        let vel_km_s = other.fixed_rows::<3>(3).into_owned();

        self.orbit.radius_km += radius_km;
        self.orbit.velocity_km_s += vel_km_s;
        self.srp.coeff_reflectivity = (self.srp.coeff_reflectivity + other[6]).clamp(0.0, 2.0);
        self.drag.coeff_drag += other[7];
        self.mass.prop_mass_kg += other[8];

        self
    }
}

impl ConfigRepr for Spacecraft {}

#[test]
fn test_serde() {
    use serde_yml;
    use std::str::FromStr;

    use anise::constants::frames::EARTH_J2000;

    let orbit = Orbit::new(
        -9042.862234,
        18536.333069,
        6999.957069,
        -3.288789,
        -2.226285,
        1.646738,
        Epoch::from_str("2018-09-15T00:15:53.098 UTC").unwrap(),
        EARTH_J2000,
    );

    let sc = Spacecraft::new(orbit, 500.0, 159.0, 2.0, 2.0, 1.8, 2.2);

    let serialized_sc = serde_yml::to_string(&sc).unwrap();
    println!("{serialized_sc}");

    let deser_sc: Spacecraft = serde_yml::from_str(&serialized_sc).unwrap();

    assert_eq!(sc, deser_sc);

    // Check that we can omit the thruster info entirely.
    let s = r#"
orbit:
    radius_km:
        - -9042.862234
        - 18536.333069
        - 6999.957069
    velocity_km_s:
        - -3.288789
        - -2.226285
        - 1.646738
    epoch: 2018-09-15T00:15:53.098000000 UTC
    frame:
      ephemeris_id: 399
      orientation_id: 1
      mu_km3_s2: null
      shape: null
mass:
    dry_mass_kg: 500.0
    prop_mass_kg: 159.0
    extra_mass_kg: 0.0
srp:
    area_m2: 2.0
    coeff_reflectivity: 1.8
drag:
    area_m2: 2.0
    coeff_drag: 2.2
    "#;

    let deser_sc: Spacecraft = serde_yml::from_str(s).unwrap();
    assert_eq!(sc, deser_sc);

    // Check that we can specify a thruster info entirely.
    let s = r#"
orbit:
    radius_km:
    - -9042.862234
    - 18536.333069
    - 6999.957069
    velocity_km_s:
    - -3.288789
    - -2.226285
    - 1.646738
    epoch: 2018-09-15T00:15:53.098000000 UTC
    frame:
        ephemeris_id: 399
        orientation_id: 1
        mu_km3_s2: null
        shape: null
mass:
    dry_mass_kg: 500.0
    prop_mass_kg: 159.0
    extra_mass_kg: 0.0
srp:
    area_m2: 2.0
    coeff_reflectivity: 1.8
drag:
    area_m2: 2.0
    coeff_drag: 2.2
thruster:
    thrust_N: 1e-5
    isp_s: 300.5
    "#;

    let mut sc_thruster = sc;
    sc_thruster.thruster = Some(Thruster {
        isp_s: 300.5,
        thrust_N: 1e-5,
    });
    let deser_sc: Spacecraft = serde_yml::from_str(s).unwrap();
    assert_eq!(sc_thruster, deser_sc);

    // Tests the minimum definition which will set all of the defaults too
    let s = r#"
orbit:
    radius_km:
    - -9042.862234
    - 18536.333069
    - 6999.957069
    velocity_km_s:
    - -3.288789
    - -2.226285
    - 1.646738
    epoch: 2018-09-15T00:15:53.098000000 UTC
    frame:
        ephemeris_id: 399
        orientation_id: 1
        mu_km3_s2: null
        shape: null
mass:
    dry_mass_kg: 500.0
    prop_mass_kg: 159.0
    extra_mass_kg: 0.0
"#;

    let deser_sc: Spacecraft = serde_yml::from_str(s).unwrap();

    let sc = Spacecraft::new(orbit, 500.0, 159.0, 0.0, 0.0, 1.8, 2.2);
    assert_eq!(sc, deser_sc);
}
