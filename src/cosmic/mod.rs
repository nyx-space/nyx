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

extern crate hyperdual;
extern crate nalgebra as na;
extern crate prost;

pub use self::xb::Xb;
use self::xb::{Ephemeris, Epoch as XbEpoch};
pub use crate::cosmic::{Frame, GuidanceMode, Orbit, Spacecraft};
pub use crate::errors::NyxError;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector};
use crate::md::StateParameter;
use crate::time::{Duration, Epoch, Unit, SECONDS_PER_DAY};
use std::fmt;
use std::fs::File;
use std::io::Read;
use std::time::Instant;

/// A trait allowing for something to have an epoch
pub trait TimeTagged {
    /// Retrieve the Epoch
    fn epoch(&self) -> Epoch;
    /// Set the Epoch
    fn set_epoch(&mut self, epoch: Epoch);

    /// Shift this epoch by a duration (can be negative)
    fn shift_by(&mut self, duration: Duration) {
        self.set_epoch(self.epoch() + duration);
    }
}

/// A trait for generate propagation and estimation state.
/// The first parameter is the size of the state, the second is the size of the propagated state including STM and extra items.
pub trait State: Copy + PartialEq + fmt::Display + fmt::LowerExp + Send + Sync
where
    Self: Sized,
    DefaultAllocator: Allocator<f64, Self::Size>
        + Allocator<f64, Self::Size, Self::Size>
        + Allocator<f64, Self::VecLength>,
{
    /// Size of the state and its STM
    type Size: DimName;
    type VecLength: DimName;

    /// Initialize an empty state
    /// By default, this is not implemented. This function must be implemented when filtering on this state.
    fn zeros() -> Self {
        unimplemented!()
    }

    /// Return this state as a vector for the propagation/estimation
    fn as_vector(&self) -> Result<OVector<f64, Self::VecLength>, NyxError>;

    /// Return this state as a vector for the propagation/estimation
    /// By default, this is not implemented. This function must be implemented when filtering on this state.
    fn stm(&self) -> Result<OMatrix<f64, Self::Size, Self::Size>, NyxError> {
        unimplemented!()
    }

    /// Return this state as a vector for the propagation/estimation
    /// By default, this is not implemented. This function must be implemented when filtering on this state.
    fn reset_stm(&mut self) {
        unimplemented!()
    }

    /// Set this state
    fn set(&mut self, epoch: Epoch, vector: &OVector<f64, Self::VecLength>)
        -> Result<(), NyxError>;

    /// Reconstruct a new State from the provided delta time in seconds compared to the current state
    /// and with the provided vector.
    fn set_with_delta_seconds(self, delta_t_s: f64, vector: &OVector<f64, Self::VecLength>) -> Self
    where
        DefaultAllocator: Allocator<f64, Self::VecLength>,
    {
        let mut me = self;
        me.set(me.epoch() + delta_t_s, vector).unwrap();
        me
    }

    /// Retrieve the Epoch
    fn epoch(&self) -> Epoch;
    /// Set the Epoch
    fn set_epoch(&mut self, epoch: Epoch);

    /// By default, this is not implemented. This function must be implemented when filtering on this state.
    fn add(self, _other: OVector<f64, Self::Size>) -> Self {
        unimplemented!()
    }

    /// Return the value of the parameter, returns an error by default
    fn value(&self, _param: &StateParameter) -> Result<f64, NyxError> {
        Err(NyxError::StateParameterUnavailable)
    }

    /// Allows setting the value of the given parameter.
    /// NOTE: Most paramaters where the `value` is available CANNOT be also set for that parameter (it's a much harder problem!)
    fn set_value(&mut self, _param: &StateParameter, _val: f64) -> Result<(), NyxError> {
        Err(NyxError::StateParameterUnavailable)
    }
}

impl XbEpoch {
    /// Returns the epoch as a raw f64, allows for speed ups if you know what is the stored time system
    pub fn as_raw(&self) -> f64 {
        f64::from(self.days) + self.seconds / SECONDS_PER_DAY
    }

    pub fn to_epoch(&self) -> Epoch {
        let epoch_delta = i64::from(self.days) * Unit::Day + self.seconds * Unit::Second;
        match self.ts {
            0 => {
                unimplemented!("TAI")
            }
            1 => match self.repr {
                5 => Epoch::from_jde_et(epoch_delta.in_unit(Unit::Day)),
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
                5 => Epoch::from_jde_tdb(epoch_delta.in_unit(Unit::Day)),
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
                    return Ok(root);
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
                    Self::ephemeris_seek_by_name(&name, &mut path, root)
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
            Self::ephemeris_names(&mut names, root);
        }
        names
    }
}

/// Known orientation IDs defined for ease of access. All Cosm objects may be accessed via Cosm directly.
pub mod orientations {
    /// J2000 orientation frame
    pub const J2000: i32 = 1;
}

// Re-Export bodies
mod bodies;
pub use self::bodies::*;

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
pub const SPEED_OF_LIGHT_KMS: f64 = SPEED_OF_LIGHT / 1000.0;

/// Astronomical unit, in kilometers, according to the [IAU](https://www.iau.org/public/themes/measuring/).
pub const AU: f64 = 149_597_870.700;

/// From NIST special publication 330, 2008 edition, in meters per second squared
pub const STD_GRAVITY: f64 = 9.80665;
