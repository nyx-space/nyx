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

use super::Bodies;
use crate::time::{Duration, Unit};
use std::cmp::PartialEq;
use std::convert::TryFrom;
use std::f64::consts::PI;
use std::fmt;

#[allow(non_snake_case, clippy::upper_case_acronyms)]
#[derive(Copy, Clone, PartialEq)]
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
    /// Velocity, Normal, Cross (called VNB in GMAT)
    VNC,
    /// Radial, Cross, Normal
    RCN,
    /// Radial, in-track, normal
    RIC,
    /// SEZ or topocentric frame. The positive horizontal vector S is due south , the positive horizontal vector E is east, and the vector Z normal to the surface of the earth (up) is the third axis.
    SEZ,
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
                for p in ephem_path.iter().flatten() {
                    path.push(*p)
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
                for p in frame_path.iter().flatten() {
                    path.push(*p)
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

    pub fn flattening_mut(&mut self, new_flattening: f64) {
        match self {
            Self::Geoid {
                ref mut flattening, ..
            } => *flattening = new_flattening,
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

    /// Returns the angular velocity for _some_ planets and moons
    /// Source for Earth: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016 (confirmed by https://hpiers.obspm.fr/eop-pc/models/constants.html)
    /// Source for everything else: https://en.wikipedia.org/w/index.php?title=Day&oldid=1008298887
    #[allow(clippy::identity_op)]
    pub fn angular_velocity(&self) -> f64 {
        let period_to_mean_motion = |dur: Duration| -> f64 { 2.0 * PI / dur.in_seconds() };
        match Bodies::try_from(self.ephem_path()).unwrap() {
            Bodies::MercuryBarycenter | Bodies::Mercury => {
                period_to_mean_motion(58 * Unit::Day + 15 * Unit::Hour + 30 * Unit::Minute)
            }
            Bodies::VenusBarycenter | Bodies::Venus => period_to_mean_motion(243 * Unit::Day),
            Bodies::Earth => 7.292_115_146_706_4e-5,
            Bodies::Luna => {
                period_to_mean_motion(27 * Unit::Day + 7 * Unit::Hour + 12 * Unit::Minute)
            }
            Bodies::MarsBarycenter => period_to_mean_motion(1 * Unit::Day + 37 * Unit::Minute),
            Bodies::JupiterBarycenter => period_to_mean_motion(9 * Unit::Hour + 56 * Unit::Minute),
            Bodies::SaturnBarycenter => period_to_mean_motion(10 * Unit::Hour + 30 * Unit::Minute),
            Bodies::UranusBarycenter => period_to_mean_motion(17 * Unit::Hour + 14 * Unit::Minute),
            Bodies::NeptuneBarycenter => period_to_mean_motion(16 * Unit::Hour + 6 * Unit::Minute),
            _ => unimplemented!(),
        }
    }

    /// Returns whether this frame is body fixed or not
    pub fn is_body_fixed(&self) -> bool {
        self.frame_path().len() == 2 || self.frame_path().len() == 3
    }
}

impl fmt::Display for Frame {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                write!(
                    f,
                    "{} {}",
                    Bodies::try_from(self.ephem_path()).unwrap().name(),
                    match self.frame_path().len() {
                        0 | 1 => "J2000".to_string(),
                        2 => "IAU Fixed".to_string(),
                        3 => "IAU Poles Fixed".to_string(),
                        _ => "Custom".to_string(),
                    }
                )
            }
            othframe => write!(f, "{:?}", othframe),
        }
    }
}

#[allow(non_snake_case, clippy::upper_case_acronyms)]
impl fmt::Debug for Frame {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Frame::Celestial { gm, .. } => {
                write!(
                    f,
                    "{} {} (μ = {:.06} km^3/s^2)",
                    Bodies::try_from(self.ephem_path()).unwrap().name(),
                    match self.frame_path().len() {
                        0 | 1 => "J2000".to_string(),
                        2 => "IAU Fixed".to_string(),
                        3 => "IAU Poles Fixed".to_string(),
                        _ => "Custom".to_string(),
                    },
                    gm,
                )
            }
            Frame::Geoid {
                gm,
                equatorial_radius,
                flattening,
                ..
            } => {
                write!(
                    f,
                    "{} {} (μ = {:.06} km^3/s^2 , r = {:.06} km, f = {:.09})",
                    Bodies::try_from(self.ephem_path()).unwrap().name(),
                    match self.frame_path().len() {
                        0 | 1 => "J2000".to_string(),
                        2 => "IAU Fixed".to_string(),
                        3 => "IAU Poles Fixed".to_string(),
                        _ => "Custom".to_string(),
                    },
                    gm,
                    equatorial_radius,
                    flattening,
                )
            }
            Frame::VNC => write!(f, "VNC"),
            Frame::RCN => write!(f, "RCN"),
            Frame::RIC => write!(f, "RIC"),
            Frame::SEZ => write!(f, "SEZ"),
            Frame::Inertial => write!(f, "Inertial"),
        }
    }
}
