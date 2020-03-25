extern crate meval;
use self::meval::{Context, Expr};
use crate::log::error;
use crate::na::Matrix3;
use crate::time::{Epoch, J2000_OFFSET, MJD_OFFSET};
use crate::utils::{r1, r2, r3};
pub use celestia::xb::Identifier as XbId;
use std::cmp::PartialEq;
use std::collections::HashMap;
use std::fmt;

/// Defines a Frame of some FrameInfo, with an identifier, and an optional parent frame.
#[derive(Clone, Debug)]
pub struct Frame<'a> {
    pub id: XbId,
    pub kind: FrameInfo,
    pub exb_id: Option<XbId>,
    pub parent: Option<(&'a Self, Box<dyn ParentRotation>)>,
}

impl<'a> PartialEq for Frame<'a> {
    fn eq(&self, other: &Self) -> bool {
        let parent_eq = if let Some(my_parent) = self.parent {
            if let Some(oth_parent) = other.parent {
                my_parent.0 == oth_parent.0
            }
            false
        } else {
            // Does not have any parent, other shouldn't either
            other.parent.is_none()
        };
        let exb_eq = if let Some(my_exb_id) = self.exb_id {
            if let Some(oth_exb_id) = other.exb_id {
                my_exb_id == oth_exb_id
            }
            false
        } else {
            // Does not have any exb_id, other shouldn't either
            other.exb_id.is_none()
        };
        self.id == other.id && self.kind == other.kind && parent_eq && exb_eq
    }
}

impl<'a> Frame<'a> {
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
}

// TODO: Rename to FrameInfo, and add an ID. Then then &Frame will be stored only in the Cosm
// and the state can be Clonable again. Also removes any lifetime problem.
// All transformations need to happen with the Cosm again, which isn't a bad thing!
// Should this also have a exb ID so that we know the center object of the frame?
// Celestial {id: i32, exb: i32, gm: f64}
// Think about printing the state. If [399] this gives only the center, not rotation.
// So maybe Celestial{fxb:i32, exb:i32, ...} since eventually everything here will be in an fxb?
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

pub trait ParentRotation: fmt::Debug {
    fn dcm_to_parent(&self, datetime: Epoch) -> Option<Matrix3<f64>>;
}

/// Defines an Euler rotation, angle must be in radians
#[derive(Clone, Copy, Debug)]
pub enum EulerRotation {
    R1(f64),
    R2(f64),
    R3(f64),
}

impl EulerRotation {
    pub fn r1_from_degrees(angle_deg: f64) -> Self {
        Self::R1(angle_deg.to_radians())
    }
    pub fn r2_from_degrees(angle_deg: f64) -> Self {
        Self::R2(angle_deg.to_radians())
    }
    pub fn r3_from_degrees(angle_deg: f64) -> Self {
        Self::R3(angle_deg.to_radians())
    }
    /// Get the DCM from this Euler rotation
    pub fn dcm(&self) -> Matrix3<f64> {
        match *self {
            Self::R1(angle) => r1(angle),
            Self::R2(angle) => r2(angle),
            Self::R3(angle) => r3(angle),
        }
    }
}

impl ParentRotation for EulerRotation {
    fn dcm_to_parent(&self, _: Epoch) -> Option<Matrix3<f64>> {
        Some(self.dcm())
    }
}

/// A fixed three-axis Euler rotation
#[derive(Debug)]
pub struct Euler3Axis {
    /// The first rotation (e.g. R3)
    pub first: EulerRotation,
    /// The second rotation (e.g. R1)
    pub second: EulerRotation,
    /// The third and final rotation (e.g. R3, to complete a 3-1-1 rotation)
    pub third: EulerRotation,
}

impl ParentRotation for Euler3Axis {
    fn dcm_to_parent(&self, _: Epoch) -> Option<Matrix3<f64>> {
        Some(self.first.dcm() * self.second.dcm() * self.third.dcm())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AngleUnit {
    Degrees,
    Radians,
}

/// A time varying three-axis Euler rotation
pub struct Euler3AxisDt<'a> {
    pub base_context: Context<'a>,
    pub rot_order: [(EulerRotation, Expr); 3],
    pub unit: AngleUnit,
    pub is_ra_dec_w: bool,
}

impl<'a> Euler3AxisDt<'a> {
    /// Initialize a new time varying transformation.
    /// Reserved keywords in the context are "T" for centuries past 2000 Jan 1 12h TBD
    /// epoch (JDE 2451545.0), and "d" for days since that epoch.
    fn new(
        first_rot: (EulerRotation, Expr),
        second_rot: (EulerRotation, Expr),
        third_rot: (EulerRotation, Expr),
        context: &HashMap<String, f64>,
        unit: AngleUnit,
        is_ra_dec_w: bool,
    ) -> Self {
        let mut ctx = Context::default();
        for (var, value) in context {
            ctx.var(var, *value);
        }
        let rot_order = [first_rot, second_rot, third_rot];
        Self {
            base_context: ctx,
            rot_order,
            unit,
            is_ra_dec_w,
        }
    }

    /// Specify how to compute this frame from the provided Euler angles and their time varying expressions.
    /// Note that these angles define how to go from THIS frame TO the PARENT frame (e.g. Sun fixed to ICRF).
    pub fn from_euler_angles(
        first_rot: (EulerRotation, Expr),
        second_rot: (EulerRotation, Expr),
        third_rot: (EulerRotation, Expr),
        context: &HashMap<String, f64>,
        unit: AngleUnit,
    ) -> Self {
        Self::new(first_rot, second_rot, third_rot, context, unit, false)
    }

    /// A time varying Right ascension, Declination, and W frame
    /// Conversion TO parent frame (e.g. Sun body to ICRF) defined as:
    /// R3(-(alpha-90 deg)) * R1(delta - 90 deg) * R3(-W)
    /// Where alpha is the right ascension and delta the declination
    pub fn from_ra_dec_w(
        right_asc: Expr,
        declin: Expr,
        w: Expr,
        context: &HashMap<String, f64>,
        unit: AngleUnit,
    ) -> Self {
        Self::new(
            (EulerRotation::R3(0.0), right_asc),
            (EulerRotation::R1(0.0), declin),
            (EulerRotation::R3(0.0), w),
            context,
            unit,
            true,
        )
    }
}

impl<'a> fmt::Debug for Euler3AxisDt<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:?}-{:?}-{:?}",
            self.rot_order[0].0, self.rot_order[1].0, self.rot_order[2].0
        )
    }
}

impl<'a> ParentRotation for Euler3AxisDt<'a> {
    fn dcm_to_parent(&self, datetime: Epoch) -> Option<Matrix3<f64>> {
        let days_d = datetime.as_jde_et_days() - MJD_OFFSET - J2000_OFFSET;
        let centuries_t = days_d / 36_525.0;
        // Now let's clone the context, and add the time variables.
        let mut ctx = self.base_context.clone();
        ctx.var("d", days_d);
        ctx.var("T", centuries_t);
        let mut dcm = Matrix3::identity();
        for (rot, expr) in &self.rot_order {
            // Compute the correct angle
            match expr.eval_with_context(&ctx) {
                Ok(eval_angle) => {
                    // Convert the angle to radians if needed
                    let angle = if self.unit == AngleUnit::Degrees {
                        eval_angle.to_radians()
                    } else {
                        eval_angle
                    };
                    let rot_with_angl = match rot {
                        EulerRotation::R1(_) => EulerRotation::R1(angle),
                        EulerRotation::R2(_) => EulerRotation::R2(angle),
                        EulerRotation::R3(_) => EulerRotation::R3(angle),
                    };
                    dcm *= rot_with_angl.dcm();
                }
                Err(e) => {
                    error!("{}", e);
                    // Stop here if something when wrong
                    return None;
                }
            }
        }
        Some(dcm)
    }
}

/*
pub struct LocalStateFrame {
    pub state: State,
}

*/

/*
pub trait Frame: Any + Copy + Clone + fmt::Debug + fmt::Display
where
    Self: Sized,
{
    fn id(&self) -> i32; // Returns the ID of this frame
    fn center_id(&self) -> i32; // Returns the integer of the center of this frame
    fn orientation_id(&self) -> i32; // Returns the orientation (1: J2000)
}

#[derive(Clone, Copy, Debug)]
pub struct Geoid {
    pub id: i32,
    pub center_id: i32,
    pub orientation_id: i32,
    pub gm: f64,
    pub flattening: f64,
    pub equatorial_radius: f64,
    pub semi_major_radius: f64,
}

impl Geoid {
    pub fn perfect_sphere(id: i32, center_id: i32, orientation_id: i32, gm: f64) -> Geoid {
        Geoid {
            id,
            center_id,
            orientation_id,
            gm,
            flattening: 0.0,
            equatorial_radius: 0.0,
            semi_major_radius: 0.0,
        }
    }
}

impl fmt::Display for Geoid {
    // Prints the same frame ID as in the FXB itself.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Geoid {}", self.id)
    }
}

impl Frame for Geoid {
    fn id(&self) -> i32 {
        self.id
    }

    fn center_id(&self) -> i32 {
        self.center_id
    }

    fn orientation_id(&self) -> i32 {
        self.orientation_id
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Spacecraft {
    pub id: i32,
    pub center_id: i32,
    pub orientation_id: i32,
}

impl Frame for Spacecraft {
    fn id(&self) -> i32 {
        self.id
    }

    fn center_id(&self) -> i32 {
        self.center_id
    }

    fn orientation_id(&self) -> i32 {
        self.orientation_id
    }
}

impl fmt::Display for Spacecraft {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "SC {}", self.id)
    }
}

#[derive(Copy, Clone, Debug)]
pub enum LocalFrame {
    /// Radial, In-track, Cross-track
    RIC,
    /// Velocity, Normal, Cross
    VNC,
    /// Radial, Cross, Normal
    RCN,
}

impl fmt::Display for LocalFrame {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Frame for LocalFrame {
    fn id(&self) -> i32 {
        match *self {
            LocalFrame::RIC => 1,
            LocalFrame::VNC => 2,
            LocalFrame::RCN => 3,
        }
    }

    fn center_id(&self) -> i32 {
        0
    }

    fn orientation_id(&self) -> i32 {
        0
    }
}
*/
#[test]
fn frame_def() {
    let ssb = Frame {
        id: XbId {
            number: 0,
            name: "Solar System Barycenter".to_owned(),
        },
        exb_id: None,
        kind: FrameInfo::Celestial {
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
        kind: FrameInfo::Celestial { gm: 1.0 },
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
}
