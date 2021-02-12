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

/*
use super::rotations::*;
use crate::na::Matrix3;
use crate::time::Epoch;
*/
// pub use celestia::xb::Identifier as XbId;
use super::Bodies;
use std::cmp::PartialEq;
use std::convert::TryFrom;
use std::fmt;

#[allow(non_snake_case)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum Frame {
    /// Any celestial frame which only has a GM (e.g. 3 body frames)
    Celestial {
        gm: f64,
        ephem_path: [Option<usize>; 3],
        frame_path: [Option<usize>; 3],
    },
    /// Any Geoid which has a GM, flattening value, etc.
    Geoid {
        gm: f64,
        flattening: f64,
        equatorial_radius: f64,
        semi_major_radius: f64,
        ephem_path: [Option<usize>; 3],
        frame_path: [Option<usize>; 3],
    },
    /// Velocity, Normal, Cross
    VNC,
    /// Radial, Cross, Normal
    RCN,
    /// Radial, in-track, normal
    RIC,
    /// Used as a placeholder only
    Inertial,
}

impl Frame {
    pub fn is_geoid(&self) -> bool {
        matches!(self, Frame::Geoid { .. })
    }

    pub fn is_celestial(&self) -> bool {
        matches!(self, Frame::Celestial { .. })
    }

    pub fn ephem_path(&self) -> Vec<usize> {
        match self {
            Frame::Celestial { ephem_path, .. } | Frame::Geoid { ephem_path, .. } => {
                let mut path = Vec::with_capacity(3);
                for p in ephem_path {
                    if let Some(f) = p {
                        path.push(*f)
                    }
                }
                path
            }
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    pub fn frame_path(&self) -> Vec<usize> {
        match self {
            Frame::Celestial { frame_path, .. } | Frame::Geoid { frame_path, .. } => {
                let mut path = Vec::with_capacity(3);
                for p in frame_path {
                    if let Some(f) = p {
                        path.push(*f)
                    }
                }
                path
            }
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    pub fn gm(&self) -> f64 {
        match self {
            Frame::Celestial { gm, .. } | Frame::Geoid { gm, .. } => *gm,
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    /// Allows mutuating the GM for this frame
    pub fn gm_mut(&mut self, new_gm: f64) {
        match self {
            Self::Geoid { ref mut gm, .. } | Self::Celestial { ref mut gm, .. } => *gm = new_gm,
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    pub fn equatorial_radius(&self) -> f64 {
        match self {
            Frame::Geoid {
                equatorial_radius, ..
            } => *equatorial_radius,
            _ => panic!("Frame is not Geoid in kind"),
        }
    }

    pub fn flattening(&self) -> f64 {
        match self {
            Frame::Geoid { flattening, .. } => *flattening,
            _ => panic!("Frame is not Geoid in kind"),
        }
    }

    pub fn semi_major_radius(&self) -> f64 {
        match self {
            Frame::Geoid {
                semi_major_radius, ..
            } => *semi_major_radius,
            _ => panic!("Frame is not Geoid in kind"),
        }
    }
}

impl fmt::Display for Frame {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                write!(
                    f,
                    "{} {}",
                    Bodies::try_from(self.ephem_path()).unwrap().name(),
                    match self.frame_path().len() {
                        0 => "J2000".to_string(),
                        1 => "IAU Fixed".to_string(),
                        2 => "IAU Poles Fixed".to_string(),
                        _ => "Custom".to_string(),
                    }
                )
            }
            othframe => write!(f, "{:?}", othframe),
        }
    }
}
