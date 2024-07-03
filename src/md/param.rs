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

use super::NyxError;
use arrow::datatypes::{DataType, Field};
use core::fmt;
use enum_iterator::Sequence;
#[cfg(feature = "python")]
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};
use std::{collections::HashMap, str::FromStr};

/// Common state parameters
#[allow(non_camel_case_types, clippy::upper_case_acronyms)]
#[derive(Copy, Clone, Debug, PartialEq, Sequence, Serialize, Deserialize)]
#[cfg_attr(feature = "python", pyclass)]
pub enum StateParameter {
    /// Argument of Latitude (deg)
    AoL,
    /// Argument of Periapse (deg)
    AoP,
    /// Apoapsis, shortcut for TA == 180.0
    Apoapsis,
    /// Radius of apoapsis (km)
    ApoapsisRadius,
    /// B-Plane B⋅R
    BdotR,
    /// B-Plane B⋅T
    BdotT,
    /// B-Plane LTOF
    BLTOF,
    /// C_3 in (km/s)^2
    C3,
    /// Coefficient of drag
    Cd,
    /// Coefficient of reflectivity
    Cr,
    /// Declination (deg) (also called elevation if in a body fixed frame)
    Declination,
    /// Dry mass (kg)
    DryMass,
    /// The epoch of the state
    Epoch,
    /// Eccentric anomaly (deg)
    EccentricAnomaly,
    /// Eccentricity (no unit)
    Eccentricity,
    /// Specific energy
    Energy,
    /// Flight path angle (deg)
    FlightPathAngle,
    /// fuel mass in kilograms
    FuelMass,
    /// Geodetic height (km)
    Height,
    /// Geodetic latitude (deg)
    Latitude,
    /// Geodetic longitude (deg)
    Longitude,
    /// Return the guidance mode of the spacecraft
    GuidanceMode,
    /// Orbital momentum
    Hmag,
    /// X component of the orbital momentum vector
    HX,
    /// Y component of the orbital momentum vector
    HY,
    /// Z component of the orbital momentum vector
    HZ,
    /// Hyperbolic anomaly (deg), only valid for hyperbolic orbits
    HyperbolicAnomaly,
    /// Inclination (deg)
    Inclination,
    /// Specific impulse (isp) in seconds
    Isp,
    /// Mean anomaly (deg)
    MeanAnomaly,
    /// Periapsis, shortcut for TA == 0.0
    Periapsis,
    /// Radius of periapse (km)
    PeriapsisRadius,
    /// Orbital period (s)
    Period,
    /// Right ascension (deg)
    RightAscension,
    /// Right ascension of the ascending node (deg)
    RAAN,
    /// Norm of the radius vector
    Rmag,
    /// Semi parameter (km)
    SemiParameter,
    /// Semi major axis (km)
    SMA,
    /// Semi minor axis (km)
    SemiMinorAxis,
    /// Thrust (Newtons)
    Thrust,
    /// True anomaly
    TrueAnomaly,
    /// True longitude
    TrueLongitude,
    /// Velocity declination (deg)
    VelocityDeclination,
    /// Norm of the velocity vector (km/s)
    Vmag,
    /// X component of the radius (km)
    X,
    /// Y component of the radius (km)
    Y,
    /// Z component of the radius (km)
    Z,
    /// X component of the velocity (km/s)
    VX,
    /// Y component of the velocity (km/s)
    VY,
    /// Z component of the velocity (km/s)
    VZ,
}

#[cfg_attr(feature = "python", pymethods)]
impl StateParameter {
    /// Returns the default event finding precision in the unit of that parameter
    pub fn default_event_precision(&self) -> f64 {
        match self {
            Self::Eccentricity => 1e-5,
            // Non anomaly angles
            Self::AoL
            | Self::AoP
            | Self::Declination
            | Self::Latitude
            | Self::Longitude
            | Self::FlightPathAngle
            | Self::Inclination
            | Self::RightAscension
            | Self::RAAN
            | Self::TrueLongitude
            | Self::VelocityDeclination => 1e-1,

            // Anomaly angles
            Self::Apoapsis
            | Self::Periapsis
            | Self::MeanAnomaly
            | Self::EccentricAnomaly
            | Self::HyperbolicAnomaly
            | Self::TrueAnomaly => 1e-3,

            // Distances
            Self::ApoapsisRadius
            | Self::BdotR
            | Self::BdotT
            | Self::Height
            | Self::Hmag
            | Self::HX
            | Self::HY
            | Self::HZ
            | Self::PeriapsisRadius
            | Self::Rmag
            | Self::SemiParameter
            | Self::SMA
            | Self::SemiMinorAxis
            | Self::X
            | Self::Y
            | Self::Z => 1e-3,

            // Velocities
            Self::C3 | Self::VX | Self::VY | Self::VZ | Self::Vmag => 1e-3,

            // Special
            Self::Energy => 1e-3,
            Self::DryMass | Self::FuelMass => 1e-3,
            Self::Period => 1e-1,
            _ => unimplemented!("{self} cannot be used for event finding"),
        }
    }

    /// Returns whether this parameter is of the B-Plane kind
    pub const fn is_b_plane(&self) -> bool {
        matches!(&self, Self::BdotR | Self::BdotT | Self::BLTOF)
    }

    /// Returns whether this is an orbital parameter
    pub const fn is_orbital(&self) -> bool {
        !self.is_for_spacecraft() && !matches!(self, Self::Apoapsis | Self::Periapsis | Self::Epoch)
    }

    /// Returns whether this parameter is only applicable to a spacecraft state
    pub const fn is_for_spacecraft(&self) -> bool {
        matches!(
            &self,
            Self::DryMass
                | Self::FuelMass
                | Self::Cr
                | Self::Cd
                | Self::Isp
                | Self::GuidanceMode
                | Self::Thrust
        )
    }

    pub const fn unit(&self) -> &'static str {
        match self {
            // Angles
            Self::AoL
            | Self::AoP
            | Self::Declination
            | Self::Latitude
            | Self::Longitude
            | Self::FlightPathAngle
            | Self::Inclination
            | Self::RightAscension
            | Self::RAAN
            | Self::TrueLongitude
            | Self::VelocityDeclination
            | Self::Apoapsis
            | Self::Periapsis
            | Self::MeanAnomaly
            | Self::EccentricAnomaly
            | Self::HyperbolicAnomaly
            | Self::TrueAnomaly => "deg",

            // Distances
            Self::ApoapsisRadius
            | Self::BdotR
            | Self::BdotT
            | Self::Height
            | Self::Hmag
            | Self::HX
            | Self::HY
            | Self::HZ
            | Self::PeriapsisRadius
            | Self::Rmag
            | Self::SemiParameter
            | Self::SMA
            | Self::SemiMinorAxis
            | Self::X
            | Self::Y
            | Self::Z => "km",

            // Velocities
            Self::VX | Self::VY | Self::VZ | Self::Vmag => "km/s",

            Self::C3 | Self::Energy => "km^2/s^2",

            Self::DryMass | Self::FuelMass => "kg",
            Self::Isp => "isp",
            Self::Thrust => "N",
            _ => "",
        }
    }

    /// Prints this orbit in Keplerian form
    #[cfg(feature = "python")]
    fn __str__(&self) -> String {
        format!("{self}")
    }

    #[cfg(feature = "python")]
    #[new]
    fn py_new(name: String) -> Result<Self, NyxError> {
        Self::from_str(&name)
    }
}

impl StateParameter {
    /// Returns the parquet field of this parameter
    pub(crate) fn to_field(self, more_meta: Option<Vec<(String, String)>>) -> Field {
        let mut meta = HashMap::new();
        meta.insert("unit".to_string(), self.unit().to_string());
        if let Some(more_data) = more_meta {
            for (k, v) in more_data {
                meta.insert(k, v);
            }
        }

        Field::new(
            format!("{self}"),
            if self == Self::GuidanceMode {
                DataType::Utf8
            } else {
                DataType::Float64
            },
            false,
        )
        .with_metadata(meta)
    }
}

impl FromStr for StateParameter {
    type Err = NyxError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let keyword = s.split_whitespace().next().ok_or(NyxError::LoadingError {
            msg: format!("Unknown state parameter: {s}"),
        })?;

        match keyword.to_lowercase().as_str() {
            "apoapsis" => Ok(Self::Apoapsis),
            "periapsis" => Ok(Self::Periapsis),
            "aol" => Ok(Self::AoL),
            "aop" => Ok(Self::AoP),
            "bltof" => Ok(Self::BLTOF),
            "bdotr" => Ok(Self::BdotR),
            "bdott" => Ok(Self::BdotT),
            "c3" => Ok(Self::C3),
            "cd" => Ok(Self::Cd),
            "cr" => Ok(Self::Cr),
            "declin" => Ok(Self::Declination),
            "dry_mass" => Ok(Self::DryMass),
            "apoapsis_radius" => Ok(Self::ApoapsisRadius),
            "ea" => Ok(Self::EccentricAnomaly),
            "ecc" => Ok(Self::Eccentricity),
            "energy" => Ok(Self::Energy),
            "fpa" => Ok(Self::FlightPathAngle),
            "fuel_mass" => Ok(Self::FuelMass),
            "guidance_mode" | "mode" => Ok(Self::GuidanceMode),
            "geodetic_height" => Ok(Self::Height),
            "geodetic_latitude" => Ok(Self::Latitude),
            "geodetic_longitude" => Ok(Self::Longitude),
            "ha" => Ok(Self::HyperbolicAnomaly),
            "hmag" => Ok(Self::Hmag),
            "hx" => Ok(Self::HX),
            "hy" => Ok(Self::HY),
            "hz" => Ok(Self::HZ),
            "inc" => Ok(Self::Inclination),
            "isp" => Ok(Self::Isp),
            "ma" => Ok(Self::MeanAnomaly),
            "periapsis_radius" => Ok(Self::PeriapsisRadius),
            "period" => Ok(Self::Period),
            "right_asc" => Ok(Self::RightAscension),
            "raan" => Ok(Self::RAAN),
            "rmag" => Ok(Self::Rmag),
            "semi_parameter" => Ok(Self::SemiParameter),
            "semi_minor" => Ok(Self::SemiMinorAxis),
            "sma" => Ok(Self::SMA),
            "ta" => Ok(Self::TrueAnomaly),
            "tlong" => Ok(Self::TrueLongitude),
            "thrust" => Ok(Self::Thrust),
            "vdeclin" => Ok(Self::VelocityDeclination),
            "vmag" => Ok(Self::Vmag),
            "x" => Ok(Self::X),
            "y" => Ok(Self::Y),
            "z" => Ok(Self::Z),
            "vx" => Ok(Self::VX),
            "vy" => Ok(Self::VY),
            "vz" => Ok(Self::VZ),
            _ => Err(NyxError::LoadingError {
                msg: format!("Unknown state parameter: {s}"),
            }),
        }
    }
}

impl fmt::Display for StateParameter {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let repr = match *self {
            Self::Apoapsis => "apoapsis",
            Self::Periapsis => "periapsis",
            Self::AoL => "aol",
            Self::AoP => "aop",
            Self::BLTOF => "BLToF",
            Self::BdotR => "BdotR",
            Self::BdotT => "BdotT",
            Self::C3 => "c3",
            Self::Cd => "cd",
            Self::Cr => "cr",
            Self::Declination => "declin",
            Self::DryMass => "dry_mass",
            Self::Epoch => "epoch",
            Self::ApoapsisRadius => "apoapsis_radius",
            Self::EccentricAnomaly => "ea",
            Self::Eccentricity => "ecc",
            Self::Energy => "energy",
            Self::FlightPathAngle => "fpa",
            Self::FuelMass => "fuel_mass",
            Self::GuidanceMode => "guidance_mode",
            Self::Height => "geodetic_height",
            Self::Latitude => "geodetic_latitude",
            Self::Longitude => "geodetic_longitude",
            Self::HyperbolicAnomaly => "ha",
            Self::Hmag => "hmag",
            Self::HX => "hx",
            Self::HY => "hy",
            Self::HZ => "hz",
            Self::Inclination => "inc",
            Self::Isp => "isp",
            Self::MeanAnomaly => "ma",
            Self::PeriapsisRadius => "periapsis_radius",
            Self::Period => "period",
            Self::RightAscension => "right_asc",
            Self::RAAN => "raan",
            Self::Rmag => "rmag",
            Self::SemiParameter => "semi_parameter",
            Self::SemiMinorAxis => "semi_minor",
            Self::SMA => "sma",
            Self::Thrust => "thrust",
            Self::TrueAnomaly => "ta",
            Self::TrueLongitude => "tlong",
            Self::VelocityDeclination => "vdeclin",
            Self::Vmag => "vmag",
            Self::X => "x",
            Self::Y => "y",
            Self::Z => "z",
            Self::VX => "vx",
            Self::VY => "vy",
            Self::VZ => "vz",
            // _ => &default,
        };
        let unit = if self.unit().is_empty() {
            String::new()
        } else {
            format!(" ({})", self.unit())
        };
        write!(f, "{repr}{unit}")
    }
}

#[cfg(test)]
mod ut_state_param {
    use super::{FromStr, StateParameter};
    #[test]
    fn test_str_to_from() {
        for s in [
            StateParameter::Apoapsis,
            StateParameter::Periapsis,
            StateParameter::AoL,
            StateParameter::AoP,
            StateParameter::BdotR,
            StateParameter::BdotT,
            StateParameter::BLTOF,
            StateParameter::C3,
            StateParameter::Cd,
            StateParameter::Cr,
            StateParameter::Declination,
            StateParameter::DryMass,
            StateParameter::ApoapsisRadius,
            StateParameter::EccentricAnomaly,
            StateParameter::Eccentricity,
            StateParameter::Energy,
            StateParameter::FlightPathAngle,
            StateParameter::FuelMass,
            StateParameter::GuidanceMode,
            StateParameter::Height,
            StateParameter::Latitude,
            StateParameter::Longitude,
            StateParameter::HyperbolicAnomaly,
            StateParameter::Hmag,
            StateParameter::HX,
            StateParameter::HY,
            StateParameter::HZ,
            StateParameter::Inclination,
            StateParameter::Isp,
            StateParameter::MeanAnomaly,
            StateParameter::PeriapsisRadius,
            StateParameter::Period,
            StateParameter::RightAscension,
            StateParameter::RAAN,
            StateParameter::Rmag,
            StateParameter::SemiParameter,
            StateParameter::SemiMinorAxis,
            StateParameter::SMA,
            StateParameter::Thrust,
            StateParameter::TrueAnomaly,
            StateParameter::TrueLongitude,
            StateParameter::VelocityDeclination,
            StateParameter::Vmag,
            StateParameter::X,
            StateParameter::Y,
            StateParameter::Z,
            StateParameter::VX,
            StateParameter::VY,
            StateParameter::VZ,
        ] {
            let as_str = format!("{s}");
            println!("{as_str}");
            let loaded = StateParameter::from_str(&as_str).unwrap();

            assert_eq!(loaded, s);
        }
    }
}
