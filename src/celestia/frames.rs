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
        axb_id: i32,
        exb_id: i32,
        gm: f64,
        parent_axb_id: Option<i32>,
        parent_exb_id: Option<i32>,
        ephem_path: [Option<usize>; 3],
        frame_path: [Option<usize>; 3],
    },
    /// Any Geoid which has a GM, flattening value, etc.
    Geoid {
        axb_id: i32,
        exb_id: i32,
        gm: f64,
        parent_axb_id: Option<i32>,
        parent_exb_id: Option<i32>,
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

    pub fn axb_id(&self) -> i32 {
        match self {
            Frame::Geoid { axb_id, .. } | Frame::Celestial { axb_id, .. } => *axb_id,
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    pub fn exb_id(&self) -> i32 {
        match self {
            Frame::Geoid { exb_id, .. } | Frame::Celestial { exb_id, .. } => *exb_id,
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    pub fn parent_axb_id(&self) -> Option<i32> {
        match self {
            Frame::Geoid { parent_axb_id, .. } | Frame::Celestial { parent_axb_id, .. } => {
                *parent_axb_id
            }
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    pub fn parent_exb_id(&self) -> Option<i32> {
        match self {
            Frame::Geoid { parent_exb_id, .. } | Frame::Celestial { parent_exb_id, .. } => {
                *parent_exb_id
            }
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
            Frame::Celestial { axb_id, exb_id, .. } | Frame::Geoid { axb_id, exb_id, .. } => {
                write!(
                    f,
                    "{} {}",
                    Bodies::try_from(self.ephem_path()).unwrap().name(),
                    if exb_id - axb_id == 99 {
                        "IAU Fixed".to_string()
                    } else {
                        match axb_id / 100 {
                            0 => "J2000".to_string(),
                            10 => "IAU Sun".to_string(),
                            1 => "Mercury IAU Fixed".to_string(),
                            2 => "Venus IAU Fixed".to_string(),
                            3 => "Earth IAU Fixed".to_string(),
                            4 => "Mars IAU Fixed".to_string(),
                            5 => "Jupiter IAU Fixed".to_string(),
                            6 => "Saturn IAU Fixed".to_string(),
                            7 => "Uranus IAU Fixed".to_string(),
                            8 => "Neptune IAU Fixed".to_string(),
                            _ => format!("{:3}", axb_id),
                        }
                    }
                )
            }
            othframe => write!(f, "{:?}", othframe),
        }
    }
}
