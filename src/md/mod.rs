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
use crate::{Orbit, SpacecraftState};
use std::fs::File;
use std::str::FromStr;

pub mod ui;

pub mod trajectory;

mod events;
pub use events::{Event, EventEvaluator};

pub type ScTraj = trajectory::Traj<SpacecraftState>;
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
    pub fn new(fmtr: StateFormatter) -> Self {
        let mut wtr = csv::Writer::from_path(fmtr.filename.clone()).expect("could not create file");
        wtr.serialize(&fmtr.headers)
            .expect("could not write headers");
        info!("Saving output to {}", fmtr.filename);

        Self { csv_out: wtr, fmtr }
    }
}

impl MdHdlr<SpacecraftState> for OrbitStateOutput {
    fn handle(&mut self, state: &SpacecraftState) {
        self.csv_out
            .serialize(self.fmtr.fmt(&state.orbit))
            .expect("could not format state");
    }
}

/// Common state parameters
#[allow(non_camel_case_types)]
#[derive(Clone, Debug, PartialEq)]
pub enum StateParameter {
    /// Argument of Latitude (deg)
    AoL(f64),
    /// Argument of Periapse (deg)
    AoP(f64),
    /// Apoapsis, shortcut for TA == 180.0
    Apoapsis,
    /// Radius of apoapsis (km)
    ApoapsisRadius(f64),
    /// Declination (deg)
    Declination(f64),
    /// Eccentric anomaly (deg)
    EccentricAnomaly(f64),
    /// Eccentricity (no unit)
    Eccentricity(f64),
    /// Specific energy
    Energy(f64),
    /// fuel mass in kilograms
    FuelMass(f64),
    /// Geodetic height (km)
    GeodeticHeight(f64),
    /// Geodetic latitude (deg)
    GeodeticLatitude(f64),
    /// Geodetic longitude (deg)
    GeodeticLongitude(f64),
    /// Orbital momentum
    Hmag(f64),
    /// X component of the orbital momentum vector
    HX(f64),
    /// Y component of the orbital momentum vector
    HY(f64),
    /// Z component of the orbital momentum vector
    HZ(f64),
    /// Inclination (deg)
    Inclination(f64),
    /// Mean anomaly (deg)
    MeanAnomaly(f64),
    /// Periapsis, shortcut for TA == 0.0
    Periapsis,
    /// Radius of periapse (km)
    PeriapsisRadius(f64),
    /// Orbital period (s)
    Period(f64),
    /// Right ascension (deg)
    RightAscension(f64),
    /// Right ascension of the ascending node (deg)
    RAAN(f64),
    /// Norm of the radius vector
    Rmag(f64),
    /// Semi parameter (km)
    SemiParameter(f64),
    /// Semi major axis (km)
    SMA(f64),
    /// Semi minor axis (km)
    SemiMinorAxis(f64),
    /// True anomaly
    TrueAnomaly(f64),
    /// True longitude
    TrueLongitude(f64),
    /// Norm of the velocity vector (km/s)
    Vmag(f64),
    /// X component of the radius (km)
    X(f64),
    /// Y component of the radius (km)
    Y(f64),
    /// Z component of the radius (km)
    Z(f64),
    /// X component of the velocity (km/s)
    VX(f64),
    /// Y component of the velocity (km/s)
    VY(f64),
    /// Z component of the velocity (km/s)
    VZ(f64),
    /// Allows creating new custom events by matching on the value of the field
    Custom(usize, f64),
}

impl FromStr for StateParameter {
    type Err = NyxError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let rplt = s.to_lowercase().replace("=", "");
        let parts: Vec<&str> = rplt.split(' ').collect();
        if parts.len() == 1 {
            match parts[0].trim() {
                "apoapsis" => Ok(Self::Apoapsis),
                "periapsis" => Ok(Self::Periapsis),
                _ => Err(NyxError::LoadingError(format!(
                    "Unknown event state parameter without value: {}",
                    s
                ))),
            }
        } else if parts.len() == 2 {
            match parts[1].trim().parse::<f64>() {
                Ok(value) => match parts[0].trim() {
                    "aol" => Ok(Self::AoL(value)),
                    "aop" => Ok(Self::AoP(value)),
                    "declin" => Ok(Self::Declination(value)),
                    "apoapsisradius" => Ok(Self::ApoapsisRadius(value)),
                    "eccentricanomaly" => Ok(Self::EccentricAnomaly(value)),
                    "eccentricity" => Ok(Self::Eccentricity(value)),
                    "energy" => Ok(Self::Energy(value)),
                    "fuelmass" => Ok(Self::FuelMass(value)),
                    "geodeticheight" => Ok(Self::GeodeticHeight(value)),
                    "geodeticlatitude" => Ok(Self::GeodeticLatitude(value)),
                    "geodeticlongitude" => Ok(Self::GeodeticLongitude(value)),
                    "hmag" => Ok(Self::Hmag(value)),
                    "hx" => Ok(Self::HX(value)),
                    "hy" => Ok(Self::HY(value)),
                    "hz" => Ok(Self::HZ(value)),
                    "inclination" => Ok(Self::Inclination(value)),
                    "meananomaly" => Ok(Self::MeanAnomaly(value)),
                    "periapsisradius" => Ok(Self::PeriapsisRadius(value)),
                    "period" => Ok(Self::Period(value)),
                    "right_asc" => Ok(Self::RightAscension(value)),
                    "raan" => Ok(Self::RAAN(value)),
                    "rmag" => Ok(Self::Rmag(value)),
                    "semiparameter" => Ok(Self::SemiParameter(value)),
                    "semiminor" => Ok(Self::SemiMinorAxis(value)),
                    "sma" => Ok(Self::SMA(value)),
                    "trueanomaly" => Ok(Self::TrueAnomaly(value)),
                    "truelongitude" => Ok(Self::TrueLongitude(value)),
                    "vmag" => Ok(Self::Vmag(value)),
                    "x" => Ok(Self::X(value)),
                    "y" => Ok(Self::Y(value)),
                    "z" => Ok(Self::Z(value)),
                    "vx" => Ok(Self::VX(value)),
                    "vy" => Ok(Self::VY(value)),
                    "vz" => Ok(Self::VZ(value)),
                    _ => Err(NyxError::LoadingError(format!(
                        "Unknown event state parameter: {}",
                        s
                    ))),
                },
                Err(_) => Err(NyxError::LoadingError(format!(
                    "Could not parse value in event: {}",
                    s
                ))),
            }
        } else {
            Err(NyxError::LoadingError(format!(
                "Could not parse event: {}",
                s
            )))
        }
    }
}
