/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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

pub use self::xb::Xb;
use self::xb::{Ephemeris, Epoch as XbEpoch};
pub use crate::cosmic::{Frame, GuidanceMode, Orbit, Spacecraft};
use crate::dynamics::DynamicsError;
pub use crate::errors::NyxError;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector};
use crate::md::StateParameter;
use crate::time::{Duration, Epoch, Unit};
use hifitime::SECONDS_PER_DAY;
use snafu::Snafu;
use std::fmt::{self, Write};
use std::fs::File;
use std::io::Read;
#[cfg(not(target_arch = "wasm32"))]
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
pub trait State: Default + Copy + PartialEq + fmt::Display + fmt::LowerExp + Send + Sync
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
    fn as_vector(&self) -> OVector<f64, Self::VecLength>;

    /// Return this state as a vector for the propagation/estimation
    /// By default, this is not implemented. This function must be implemented when filtering on this state.
    fn stm(&self) -> Result<OMatrix<f64, Self::Size, Self::Size>, DynamicsError> {
        unimplemented!()
    }

    /// Return this state as a vector for the propagation/estimation
    /// By default, this is not implemented. This function must be implemented when filtering on this state.
    fn reset_stm(&mut self) {
        unimplemented!()
    }

    /// Unsets the STM for this state
    fn unset_stm(&mut self);

    /// Set this state
    fn set(&mut self, epoch: Epoch, vector: &OVector<f64, Self::VecLength>);

    /// Reconstruct a new State from the provided delta time in seconds compared to the current state
    /// and with the provided vector.
    fn set_with_delta_seconds(self, delta_t_s: f64, vector: &OVector<f64, Self::VecLength>) -> Self
    where
        DefaultAllocator: Allocator<f64, Self::VecLength>,
    {
        let mut me = self;
        me.set(me.epoch() + delta_t_s, vector);
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
    fn value(&self, param: StateParameter) -> Result<f64, NyxError> {
        Err(NyxError::StateParameterUnavailable {
            param,
            msg: "unimplemented in State trait".to_string(),
        })
    }

    /// Allows setting the value of the given parameter.
    /// NOTE: Most parameters where the `value` is available CANNOT be also set for that parameter (it's a much harder problem!)
    fn set_value(&mut self, param: StateParameter, _val: f64) -> Result<(), NyxError> {
        Err(NyxError::StateParameterUnavailable {
            param,
            msg: "unimplemented in State trait".to_string(),
        })
    }
}

#[derive(Copy, Clone, Debug, PartialEq, Snafu)]
pub enum AstroError {
    #[snafu(display("B Plane jacobian invariant must be either VX, VY or VZ"))]
    BPlaneInvariant,
    #[snafu(display("operation requires a local frame"))]
    NotLocalFrame,
    #[snafu(display("partial derivatives not defined for this parameter"))]
    PartialsUndefined,
    #[snafu(display("Orbit is not hyperbolic so there is no hyperbolic anomaly."))]
    NotHyperbolic,
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
                5 => Epoch::from_jde_et(epoch_delta.to_unit(Unit::Day)),
                _ => unimplemented!("ET"),
            },
            2 => match self.repr {
                0 => Epoch::from_tt_seconds(epoch_delta.to_seconds()),
                _ => unimplemented!("TT"),
            },
            3 => {
                unimplemented!("UTC")
            }
            4 => match self.repr {
                2 => Epoch::from_tdb_seconds(epoch_delta.to_seconds()),
                5 => Epoch::from_jde_tdb(epoch_delta.to_unit(Unit::Day)),
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
            Err(e) => {
                return Err(NyxError::LoadingError {
                    msg: format!("{e}"),
                })
            }
            Ok(mut f) => {
                if f.read_to_end(&mut input_xb_buf).is_err() {
                    return Err(NyxError::LoadingError {
                        msg: "Could not read buffer".to_string(),
                    });
                }
            }
        };

        Self::from_buffer(&input_xb_buf)
    }

    /// Loads the provided input buffer as an XB
    pub fn from_buffer(input_xb_buf: &[u8]) -> Result<Self, NyxError> {
        use prost::Message;
        if input_xb_buf.is_empty() {
            return Err(NyxError::LoadingError {
                msg: "XB buffer is empty".to_string(),
            });
        }

        #[cfg(not(target_arch = "wasm32"))]
        let decode_start = Instant::now();

        match Self::decode(input_xb_buf) {
            Ok(xb) => {
                #[cfg(not(target_arch = "wasm32"))]
                debug!("Loaded XB in {} ms.", decode_start.elapsed().as_millis());
                Ok(xb)
            }
            Err(e) => Err(NyxError::LoadingError {
                msg: format!("Could not decode XB: {e}"),
            }),
        }
    }

    /// Finds the ephemeris provided the path as usize, e.g. [3,1] would return the Moon with any DE xb.
    pub fn ephemeris_from_path<'a>(&'a self, path: &[usize]) -> Result<&'a Ephemeris, NyxError> {
        match &self.ephemeris_root {
            None => Err(NyxError::ObjectNotFound {
                needle: "not ephemeris root".to_string(),
                haystack: self.ephemeris_get_names(),
            }),
            Some(root) => {
                if path.is_empty() {
                    return Ok(root);
                }
                for pos in path {
                    if root.children.get(*pos).is_none() {
                        let hpath = path.iter().fold(String::new(), |mut output, p| {
                            let _ = write!(output, "{p}");
                            output
                        });
                        return Err(NyxError::ObjectNotFound {
                            needle: hpath,
                            haystack: self.ephemeris_get_names(),
                        });
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
        cur_path: &[usize],
        e: &Ephemeris,
    ) -> Result<Vec<usize>, NyxError> {
        if e.name == name {
            Ok(cur_path.to_vec())
        } else if e.children.is_empty() {
            Err(NyxError::ObjectNotFound {
                needle: name.to_string(),
                haystack: e.children.iter().map(|c| c.name.clone()).collect(),
            })
        } else {
            for (cno, child) in e.children.iter().enumerate() {
                let mut this_path = cur_path.to_owned();
                this_path.push(cno);
                let child_attempt = Self::ephemeris_seek_by_name(name, &this_path, child);
                if let Ok(found_path) = child_attempt {
                    return Ok(found_path);
                }
            }
            // Could not find name in iteration, fail
            Err(NyxError::ObjectNotFound {
                needle: name.to_string(),
                haystack: e.children.iter().map(|c| c.name.clone()).collect(),
            })
        }
    }

    /// Returns the machine path of the requested ephemeris
    pub fn ephemeris_find_path(&self, name: String) -> Result<Vec<usize>, NyxError> {
        match &self.ephemeris_root {
            None => Err(NyxError::ObjectNotFound {
                needle: "No root!".to_string(),
                haystack: self.ephemeris_get_names(),
            }),
            Some(root) => {
                if root.name == name {
                    // Return an empty vector (but OK because we're asking for the root)
                    Ok(Vec::new())
                } else {
                    let path = Vec::new();
                    Self::ephemeris_seek_by_name(&name, &path, root)
                }
            }
        }
    }

    fn ephemeris_names(names: &mut Vec<String>, e: &Ephemeris) {
        names.push(e.name.clone());
        for child in &e.children {
            Self::ephemeris_names(names, child);
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
