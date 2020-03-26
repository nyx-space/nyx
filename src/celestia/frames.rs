/*
use super::rotations::*;
use crate::na::Matrix3;
use crate::time::Epoch;
*/
pub use celestia::xb::Identifier as XbId;
use std::cmp::PartialEq;
use std::fmt;

/// Defines a Frame of some FrameInfo, with an identifier, and an optional parent frame.
#[derive(Debug)]
pub struct Frame {
    pub id: XbId,
    pub info: FrameInfo,
    pub exb_id: Option<XbId>,
}

impl PartialEq for Frame {
    fn eq(&self, other: &Self) -> bool {
        let exb_eq = if let Some(my_exb_id) = &self.exb_id {
            if let Some(oth_exb_id) = &other.exb_id {
                my_exb_id == oth_exb_id
            } else {
                false
            }
        } else {
            // Does not have any exb_id, other shouldn't either
            other.exb_id.is_none()
        };
        self.id == other.id && self.info == other.info && exb_eq
    }
}

// impl<'a> Frame<'a> {
/*
pub fn try_dcm_to_parent(&self, datetime: Epoch) -> Option<Matrix3<f64>> {
    if let Some(parent) = self.parent.as_ref() {
        if let Some(dcm) = parent.1.dcm_to_parent(datetime) {
            return Some(dcm);
        }
    }
    None
}

pub fn dcm_to_parent(&self, datetime: Epoch) -> Matrix3<f64> {
    self.try_dcm_to_parent(datetime).unwrap()
}

pub fn try_dcm_from_parent(&self, datetime: Epoch) -> Option<Matrix3<f64>> {
    if let Some(dcm) = self.try_dcm_to_parent(datetime) {
        Some(dcm.transpose())
    } else {
        None
    }
}

pub fn dcm_from_parent(&self, datetime: Epoch) -> Matrix3<f64> {
    self.try_dcm_from_parent(datetime).unwrap()
}
*/
// }

// TODO: Rename to FrameInfo, and add an ID. Then then &Frame will be stored only in the Cosm
// and the state can be Clonable again. Also removes any lifetime problem.
// All transformations need to happen with the Cosm again, which isn't a bad thing!
// Should this also have a exb ID so that we know the center object of the frame?
// Celestial {id: i32, exb: i32, gm: f64}
// Think about printing the state. If [399] this gives only the center, not rotation.
// So maybe Celestial{fxb:i32, exb:i32, ...} since eventually everything here will be in an fxb?
#[allow(non_snake_case)]
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum FrameInfo {
    /// Any celestial frame which only has a GM (e.g. 3 body frames)
    Celestial { fxb_id: i32, exb_id: i32, gm: f64 },
    /// Any Geoid which has a GM, flattening value, etc.
    Geoid {
        fxb_id: i32,
        exb_id: i32,
        gm: f64,
        flattening: f64,
        equatorial_radius: f64,
        semi_major_radius: f64,
    },
    /// Velocity, Normal, Cross
    VNC,
    /// Radial, Cross, Normal
    RCN,
    /// Radial, in-track, normal
    RIC,
}

impl FrameInfo {
    pub fn is_geoid(&self) -> bool {
        match self {
            FrameInfo::Geoid { .. } => true,
            _ => false,
        }
    }

    pub fn is_celestial(&self) -> bool {
        match self {
            FrameInfo::Celestial { .. } => true,
            _ => false,
        }
    }

    pub fn gm(&self) -> f64 {
        match self {
            FrameInfo::Celestial { gm, .. } | FrameInfo::Geoid { gm, .. } => *gm,
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    pub fn fxb_id(&self) -> i32 {
        match self {
            FrameInfo::Geoid { fxb_id, .. } | FrameInfo::Celestial { fxb_id, .. } => *fxb_id,
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    pub fn exb_id(&self) -> i32 {
        match self {
            FrameInfo::Geoid { exb_id, .. } | FrameInfo::Celestial { exb_id, .. } => *exb_id,
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    pub fn equatorial_radius(&self) -> f64 {
        match self {
            FrameInfo::Geoid {
                equatorial_radius, ..
            } => *equatorial_radius,
            _ => panic!("Frame is not Geoid in kind"),
        }
    }

    pub fn flattening(&self) -> f64 {
        match self {
            FrameInfo::Geoid { flattening, .. } => *flattening,
            _ => panic!("Frame is not Geoid in kind"),
        }
    }

    pub fn semi_major_radius(&self) -> f64 {
        match self {
            FrameInfo::Geoid {
                semi_major_radius, ..
            } => *semi_major_radius,
            _ => panic!("Frame is not Geoid in kind"),
        }
    }
}

impl fmt::Display for FrameInfo {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            FrameInfo::Celestial { fxb_id, exb_id, .. }
            | FrameInfo::Geoid { fxb_id, exb_id, .. } => write!(f, "{} ({})", exb_id, fxb_id),
            othframe => write!(f, "{:?}", othframe),
        }
    }
}

#[test]
fn frame_def() {
    /*
    let ssb = Frame {
        id: XbId {
            number: 0,
            name: "Solar System Barycenter".to_owned(),
        },
        exb_id: None,
        info: FrameInfo::Celestial {
            gm: 1.0014 * 132_712_440_041.939_38,
        },
        parent: None,
    };

    let sun2ssb_ctx = HashMap::new();
    let right_asc: meval::Expr = "289.13".parse().unwrap();
    let declin: meval::Expr = "63.87".parse().unwrap();
    let w_expr: meval::Expr = "84.176 + 14.18440000*d".parse().unwrap();
    let sun2ssb_rot =
        Euler3AxisDt::from_ra_dec_w(right_asc, declin, w_expr, &sun2ssb_ctx, AngleUnit::Degrees);

    let sun_fixed = Frame {
        id: XbId {
            number: 10,
            name: "Sun body fixed".to_owned(),
        },
        exb_id: None,
        info: FrameInfo::Celestial { gm: 1.0 },
        parent: Some((&ssb, Box::new(sun2ssb_rot))),
    };

    let dt = Epoch::from_gregorian_tai_at_noon(2000, 1, 1);
    println!("{}", sun_fixed.dcm_to_parent(dt));
    println!("{}", sun_fixed.dcm_to_parent(dt + 1.0));
    println!(
        "{}",
        sun_fixed.dcm_to_parent(dt) * sun_fixed.dcm_from_parent(dt)
    );

    let mut map: HashMap<i32, Frame> = HashMap::new();
    map.insert(ssb.id.number, ssb);

    println!("{:?}", map);
    */
}
