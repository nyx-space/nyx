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

use super::NyxError;
use core::fmt;
use std::str::FromStr;

/// Common state parameters
#[allow(non_camel_case_types, clippy::upper_case_acronyms)]
#[derive(Copy, Clone, Debug, PartialEq)]
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
    /// Declination (deg)
    Declination,
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
    GeodeticHeight,
    /// Geodetic latitude (deg)
    GeodeticLatitude,
    /// Geodetic longitude (deg)
    GeodeticLongitude,
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
    /// Computes the slant angle by returning the arccos of the dot product between state's radius vector and the provided vector coordinates. The vector will be normalized if needed.
    SlantAngle { x: f64, y: f64, z: f64 },
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
    /// Allows creating new custom events by matching on the value of the field
    Custom { mapping: usize },
}

impl StateParameter {
    /// Returns the default event finding precision in the unit of that parameter
    pub fn default_event_precision(self) -> f64 {
        match self {
            Self::Eccentricity => 1e-5,
            // Non anomaly angles
            Self::AoL
            | Self::AoP
            | Self::Declination
            | Self::GeodeticLatitude
            | Self::GeodeticLongitude
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
            | Self::SlantAngle { .. }
            | Self::TrueAnomaly => 1e-3,

            // Distances
            Self::ApoapsisRadius
            | Self::BdotR
            | Self::BdotT
            | Self::GeodeticHeight
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
            Self::FuelMass => 1e-3,
            Self::Period => 1e-1,
            _ => unimplemented!(),
        }
    }

    /// Returns whether this parameter is of the B-Plane kind
    pub fn is_b_plane(self) -> bool {
        self == Self::BdotR || self == Self::BdotT || self == Self::BLTOF
    }

    pub fn unit(self) -> &'static str {
        match self {
            // Angles
            Self::AoL
            | Self::AoP
            | Self::Declination
            | Self::GeodeticLatitude
            | Self::GeodeticLongitude
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
            | Self::SlantAngle { .. }
            | Self::TrueAnomaly => "deg",

            // Distances
            Self::ApoapsisRadius
            | Self::BdotR
            | Self::BdotT
            | Self::GeodeticHeight
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
            Self::C3 | Self::VX | Self::VY | Self::VZ | Self::Vmag => "km/s",

            Self::FuelMass => "kg",
            Self::Isp => "isp",
            Self::Thrust => "N",
            _ => "",
        }
    }
}

impl FromStr for StateParameter {
    type Err = NyxError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let keyword = s.trim().to_lowercase();

        match keyword.as_str() {
            "apoapsis" => Ok(Self::Apoapsis),
            "periapsis" => Ok(Self::Periapsis),
            "aol" => Ok(Self::AoL),
            "aop" => Ok(Self::AoP),
            "c3" => Ok(Self::C3),
            "cd" => Ok(Self::Cd),
            "cr" => Ok(Self::Cr),
            "declin" => Ok(Self::Declination),
            "apoapsis_radius" => Ok(Self::ApoapsisRadius),
            "ea" => Ok(Self::EccentricAnomaly),
            "ecc" => Ok(Self::Eccentricity),
            "energy" => Ok(Self::Energy),
            "fpa" => Ok(Self::FlightPathAngle),
            "fuel_mass" => Ok(Self::FuelMass),
            "geodetic_height" => Ok(Self::GeodeticHeight),
            "geodetic_latitude" => Ok(Self::GeodeticLatitude),
            "geodetic_longitude" => Ok(Self::GeodeticLongitude),
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
            _ => Err(NyxError::LoadingError(format!(
                "Unknown state parameter: {}",
                s
            ))),
        }
    }
}

impl fmt::Display for StateParameter {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let default = format!("{:?}", self);
        let repr = match *self {
            Self::Apoapsis => "apoapsis",
            Self::Periapsis => "periapsis",
            Self::AoL => "aol",
            Self::AoP => "aop",
            Self::C3 => "c3",
            Self::Cd => "cd",
            Self::Cr => "cr",
            Self::Declination => "declin",
            Self::ApoapsisRadius => "apoapsis_radius",
            Self::EccentricAnomaly => "ea",
            Self::Eccentricity => "ecc",
            Self::Energy => "energy",
            Self::FlightPathAngle => "fpa",
            Self::FuelMass => "fuel_mass",
            Self::GeodeticHeight => "geodetic_height",
            Self::GeodeticLatitude => "geodetic_latitude",
            Self::GeodeticLongitude => "geodetic_longitude",
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
            _ => &default,
        };
        write!(f, "{}", repr)
    }
}
