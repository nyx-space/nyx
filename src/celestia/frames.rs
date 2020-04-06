/*
use super::rotations::*;
use crate::na::Matrix3;
use crate::time::Epoch;
*/
pub use celestia::xb::Identifier as XbId;
use std::cmp::PartialEq;
use std::fmt;

// TODO: Rename to Frame, and add an ID. Then then &Frame will be stored only in the Cosm
// and the state can be Clonable again. Also removes any lifetime problem.
// All transformations need to happen with the Cosm again, which isn't a bad thing!
// Should this also have a exb ID so that we know the center object of the frame?
// Celestial {id: i32, exb: i32, gm: f64}
// Think about printing the state. If [399] this gives only the center, not rotation.
// So maybe Celestial{fxb:i32, exb:i32, ...} since eventually everything here will be in an fxb?
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
    },
    /// Velocity, Normal, Cross
    VNC,
    /// Radial, Cross, Normal
    RCN,
    /// Radial, in-track, normal
    RIC,
}

impl Frame {
    pub fn is_geoid(&self) -> bool {
        match self {
            Frame::Geoid { .. } => true,
            _ => false,
        }
    }

    pub fn is_celestial(&self) -> bool {
        match self {
            Frame::Celestial { .. } => true,
            _ => false,
        }
    }

    pub fn gm(&self) -> f64 {
        match self {
            Frame::Celestial { gm, .. } | Frame::Geoid { gm, .. } => *gm,
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
                write!(f, "{} ({})", exb_id, axb_id)
            }
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
        info: Frame::Celestial {
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
        info: Frame::Celestial { gm: 1.0 },
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
