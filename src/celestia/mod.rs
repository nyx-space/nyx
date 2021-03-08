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

extern crate hyperdual;
extern crate nalgebra as na;
extern crate prost;

use std::convert::TryFrom;
use std::fs::File;
use std::io::Read;
use std::time::Instant;

pub use self::xb::Xb;
use self::xb::{Ephemeris, Epoch as XbEpoch};
pub use crate::errors::NyxError;
use crate::time::{Epoch, TimeUnit, SECONDS_PER_DAY};

impl XbEpoch {
    /// Returns the epoch as a raw f64, allows for speed ups if you know what is the stored time system
    pub fn as_raw(&self) -> f64 {
        f64::from(self.days) + self.seconds / SECONDS_PER_DAY
    }

    pub fn to_epoch(&self) -> Epoch {
        let epoch_delta = self.days * TimeUnit::Day + self.seconds * TimeUnit::Second;
        match self.ts {
            0 => {
                unimplemented!("TAI")
            }
            1 => match self.repr {
                5 => Epoch::from_jde_et(epoch_delta.in_unit_f64(TimeUnit::Day)),
                _ => unimplemented!("ET"),
            },
            2 => match self.repr {
                0 => Epoch::from_tt_seconds(epoch_delta.in_seconds()),
                _ => unimplemented!("TT"),
            },
            3 => {
                unimplemented!("UTC")
            }
            4 => match self.repr {
                2 => Epoch::from_tdb_seconds(epoch_delta.in_seconds()),
                5 => Epoch::from_jde_tdb(epoch_delta.in_unit_f64(TimeUnit::Day)),
                _ => unimplemented!("TDB"),
            },
            _ => unimplemented!(),
        }
    }
}

impl Xb {
    /// Loads the provided input_filename as an XB
    pub fn from_file(input_filename: &str) -> Result<Self, NyxError> {
        let mut input_xb_buf = Vec::new();

        match File::open(input_filename) {
            Err(e) => return Err(NyxError::LoadingError(format!("{}", e))),
            Ok(mut f) => {
                if f.read_to_end(&mut input_xb_buf).is_err() {
                    return Err(NyxError::LoadingError("Could not read buffer".to_string()));
                }
            }
        };

        Self::from_buffer(&input_xb_buf)
    }

    /// Loads the provided input buffer as an XB
    pub fn from_buffer(input_xb_buf: &[u8]) -> Result<Self, NyxError> {
        use self::prost::Message;
        if input_xb_buf.is_empty() {
            return Err(NyxError::LoadingError("XB buffer is empty".to_string()));
        }

        let decode_start = Instant::now();

        match Self::decode(&*input_xb_buf) {
            Ok(xb) => {
                info!("Loaded XB in {} ms.", decode_start.elapsed().as_millis());
                Ok(xb)
            }
            Err(e) => Err(NyxError::LoadingError(format!(
                "Could not decode XB: {}",
                e
            ))),
        }
    }

    /// Finds the ephemeris provided the path as usize, e.g. [3,1] would return the Moon with any DE xb.
    pub fn ephemeris_from_path<'a>(&'a self, path: &[usize]) -> Result<&'a Ephemeris, NyxError> {
        match &self.ephemeris_root {
            None => Err(NyxError::ObjectNotFound("not ephemeris root".to_string())),
            Some(root) => {
                if path.is_empty() {
                    return Ok(&root);
                }
                for pos in path {
                    if root.children.get(*pos).is_none() {
                        let hpath: String =
                            path.iter().map(|p| format!("{}", p)).collect::<String>();
                        return Err(NyxError::ObjectNotFound(hpath));
                    }
                }

                // This is absolutely terrible, and there must be a better way to do it, but it's late.
                match path.len() {
                    1 => Ok(&self.ephemeris_root.as_ref().unwrap().children[path[0]]),
                    2 => Ok(
                        &self.ephemeris_root.as_ref().unwrap().children[path[0]].children[path[1]],
                    ),
                    3 => Ok(
                        &self.ephemeris_root.as_ref().unwrap().children[path[0]].children[path[1]]
                            .children[path[2]],
                    ),
                    _ => unimplemented!(),
                }
            }
        }
    }

    /// Seek an ephemeris from its celestial name (e.g. Earth Moon Barycenter)
    fn ephemeris_seek_by_name(
        name: &str,
        cur_path: &mut Vec<usize>,
        e: &Ephemeris,
    ) -> Result<Vec<usize>, NyxError> {
        if e.name == name {
            Ok(cur_path.to_vec())
        } else if e.children.is_empty() {
            Err(NyxError::ObjectNotFound(name.to_string()))
        } else {
            for (cno, child) in e.children.iter().enumerate() {
                let mut this_path = cur_path.clone();
                this_path.push(cno);
                let child_attempt = Self::ephemeris_seek_by_name(name, &mut this_path, child);
                if let Ok(found_path) = child_attempt {
                    return Ok(found_path);
                }
            }
            // Could not find name in iteration, fail
            Err(NyxError::ObjectNotFound(name.to_string()))
        }
    }

    /// Returns the machine path of the requested ephemeris
    pub fn ephemeris_find_path(&self, name: String) -> Result<Vec<usize>, NyxError> {
        match &self.ephemeris_root {
            None => Err(NyxError::ObjectNotFound("No root!".to_string())),
            Some(root) => {
                if root.name == name {
                    // Return an empty vector (but OK because we're asking for the root)
                    Ok(Vec::new())
                } else {
                    let mut path = Vec::new();
                    Self::ephemeris_seek_by_name(&name, &mut path, &root)
                }
            }
        }
    }

    fn ephemeris_names(mut names: &mut Vec<String>, e: &Ephemeris) {
        names.push(e.name.clone());
        for child in &e.children {
            Self::ephemeris_names(&mut names, child);
        }
    }

    pub fn ephemeris_get_names(&self) -> Vec<String> {
        let mut names = Vec::new();
        if let Some(root) = &self.ephemeris_root {
            Self::ephemeris_names(&mut names, &root);
        }
        names
    }
}

/// Known orientation IDs defined for ease of access. All Cosm objects may be accessed via Cosm directly.
pub mod orientations {
    /// J2000 orientation frame
    pub const J2000: i32 = 1;
}

/// Defines the default celestial bodies in the provided de438 XB.
#[derive(Copy, Clone, Debug)]
pub enum Bodies {
    SSB,
    Sun,
    MercuryBarycenter,
    Mercury,
    VenusBarycenter,
    Venus,
    EarthBarycenter,
    Earth,
    Luna,
    MarsBarycenter,
    JupiterBarycenter,
    SaturnBarycenter,
    UranusBarycenter,
    NeptuneBarycenter,
}

impl Bodies {
    /// Returns the path in the standard de438 XB
    pub fn ephem_path(&self) -> &'static [usize] {
        match *self {
            Self::SSB => &[],
            Self::Sun => &[0],
            Self::MercuryBarycenter => &[1],
            Self::Mercury => &[1],
            Self::VenusBarycenter => &[2],
            Self::Venus => &[2],
            Self::EarthBarycenter => &[3],
            Self::Earth => &[3, 0],
            Self::Luna => &[3, 1],
            Self::MarsBarycenter => &[4],
            Self::JupiterBarycenter => &[5],
            Self::SaturnBarycenter => &[6],
            Self::UranusBarycenter => &[7],
            Self::NeptuneBarycenter => &[8],
        }
    }

    /// Returns the human name
    pub fn name(&self) -> String {
        match *self {
            Self::SSB => "Solar System Barycenter".to_string(),
            Self::Sun => "Sun".to_string(),
            Self::MercuryBarycenter => "Mercury".to_string(),
            Self::Mercury => "Mercury".to_string(),
            Self::VenusBarycenter => "Venus".to_string(),
            Self::Venus => "Venus".to_string(),
            Self::EarthBarycenter => "Earth Moon Barycenter".to_string(),
            Self::Earth => "Earth".to_string(),
            Self::Luna => "Moon".to_string(),
            Self::MarsBarycenter => "Mars".to_string(),
            Self::JupiterBarycenter => "Jupiter Barycenter".to_string(),
            Self::SaturnBarycenter => "Saturn Barycenter".to_string(),
            Self::UranusBarycenter => "Uranus Barycenter".to_string(),
            Self::NeptuneBarycenter => "Neptune Barycenter".to_string(),
        }
    }
}

impl TryFrom<String> for Bodies {
    type Error = NyxError;

    fn try_from(name: String) -> Result<Self, Self::Error> {
        match name.to_lowercase().as_str() {
            "solar system barycenter" | "ssb" => Ok(Self::SSB),
            "sun" => Ok(Self::Sun),
            "mercury" => Ok(Self::Mercury),
            "venus" => Ok(Self::Venus),
            "earth moon barycenter" => Ok(Self::EarthBarycenter),
            "earth" => Ok(Self::Earth),
            "moon" | "luna" => Ok(Self::Luna),
            "mars" | "mars barycenter" => Ok(Self::MarsBarycenter),
            "jupiter" | "jupiter barycenter" => Ok(Self::JupiterBarycenter),
            "saturn" | "saturn barycenter" => Ok(Self::SaturnBarycenter),
            "uranus" | "uranus barycenter" => Ok(Self::UranusBarycenter),
            "neptune" | "neptune barycenter" => Ok(Self::NeptuneBarycenter),
            _ => Err(NyxError::ObjectNotFound(name)),
        }
    }
}

impl TryFrom<Vec<usize>> for Bodies {
    type Error = NyxError;

    fn try_from(ephem_path: Vec<usize>) -> Result<Self, Self::Error> {
        match ephem_path.len() {
            0 => Ok(Self::SSB),
            1 => match ephem_path[0] {
                0 => Ok(Self::Sun),
                1 => Ok(Self::Mercury),
                2 => Ok(Self::Venus),
                3 => Ok(Self::EarthBarycenter),
                4 => Ok(Self::MarsBarycenter),
                5 => Ok(Self::JupiterBarycenter),
                6 => Ok(Self::SaturnBarycenter),
                7 => Ok(Self::UranusBarycenter),
                8 => Ok(Self::NeptuneBarycenter),
                _ => Err(NyxError::ObjectNotFound(format!("{:?}", ephem_path))),
            },
            2 if ephem_path[0] == 3 => match ephem_path[1] {
                // This only support the Earth system
                0 => Ok(Self::Earth),
                1 => Ok(Self::Luna),
                _ => Err(NyxError::ObjectNotFound(format!("{:?}", ephem_path))),
            },
            _ => Err(NyxError::ObjectNotFound(format!("{:?}", ephem_path))),
        }
    }
}

// Re-Export orbit
mod orbit;
pub use self::orbit::*;

// Re-Export OrbitDual
mod orbitdual;
pub use self::orbitdual::*;

// Re-Export B Plane
mod bplane;
pub use self::bplane::*;

// Re-Export spacecraft
mod spacecraft;
pub use self::spacecraft::*;

// Re-Export frames
mod frames;
pub use self::frames::*;

mod rotations;
pub use self::rotations::*;

mod cosm;
mod xb;
pub use self::cosm::*;

/// The eclipse module allows finding eclipses and (conversely) visibility between a state and another one (e.g. a planet or the Sun).
pub mod eclipse;

/// Speed of light in meters per second
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;
/// Speed of light in kilometers per second
pub const SPEED_OF_LIGHT_KMS: f64 = 299_792.458;

/// Astronomical unit, in kilometers, according to the [IAU](https://www.iau.org/public/themes/measuring/).
pub const AU: f64 = 149_597_870.700;
