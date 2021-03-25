/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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

extern crate bacon_sci;
extern crate crossbeam;
extern crate csv;
extern crate rayon;

use crate::errors::NyxError;
use crate::io::formatter::StateFormatter;
use crate::{Orbit, Spacecraft};
use std::fs::File;
use std::str::FromStr;

pub mod ui;

pub mod trajectory;

mod events;
pub use events::{Event, EventEvaluator};

pub mod targeter;

pub type ScTraj = trajectory::Traj<Spacecraft>;
pub type Ephemeris = trajectory::Traj<Orbit>;

/// A Mission Design handler
pub trait MdHdlr<StateType: Copy>: Send + Sync {
    fn handle(&mut self, state: &StateType);
}

pub struct OrbitStateOutput {
    csv_out: csv::Writer<File>,
    fmtr: StateFormatter,
}

impl OrbitStateOutput {
    pub fn new(fmtr: StateFormatter) -> Result<Self, NyxError> {
        match csv::Writer::from_path(fmtr.filename.clone()) {
            Ok(mut wtr) => {
                wtr.serialize(&fmtr.headers)
                    .expect("could not write headers");
                info!("Saving output to {}", fmtr.filename);

                Ok(Self { csv_out: wtr, fmtr })
            }
            Err(e) => Err(NyxError::ExportError(e.to_string())),
        }
    }
}

impl MdHdlr<Spacecraft> for OrbitStateOutput {
    fn handle(&mut self, state: &Spacecraft) {
        self.csv_out
            .serialize(self.fmtr.fmt(&state.orbit))
            .expect("could not format state");
    }
}

impl MdHdlr<Orbit> for OrbitStateOutput {
    fn handle(&mut self, state: &Orbit) {
        self.csv_out
            .serialize(self.fmtr.fmt(&state))
            .expect("could not format state");
    }
}

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
    /// C_3 in (km/s)^2
    C3,
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
            | Self::TrueAnomaly => 1e-3,

            // Distances
            Self::ApoapsisRadius
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
