extern crate meval;
extern crate serde_derive;
extern crate toml;

use self::meval::Expr;
use self::serde_derive::Deserialize;
use crate::celestia::{AngleUnit, Euler3AxisDt, Frame};
use std::collections::HashMap;
use std::str::FromStr;

#[derive(Deserialize)]
pub struct FramesToml {
    pub frames: HashMap<String, FrameToml>,
}

// TODO: Support copying info from the J2000 frames somehow by specifying the name, e.g. "Venus barycenter J2000"

#[derive(Clone, Deserialize)]
pub struct FrameToml {
    axb: i32,
    exb: i32,
    // Set this in the TOML to clone the negative values of this frame from another frame
    pub clonefrom: Option<String>,
    gm: f64,
    parent_axb: Option<i32>,
    parent_exb: Option<i32>,
    flattening: f64,
    equatorial_radius: f64,
    semi_major_radius: f64,
    pub rotation: RotationToml,
}

impl FrameToml {
    pub fn to_frame(&self) -> Frame {
        Frame::Geoid {
            axb_id: self.axb,
            exb_id: self.exb,
            gm: self.gm,
            parent_axb_id: self.parent_axb,
            parent_exb_id: self.parent_exb,
            flattening: self.flattening,
            equatorial_radius: self.equatorial_radius,
            semi_major_radius: self.semi_major_radius,
        }
    }

    pub fn update_from(&mut self, src: &Frame) {
        if self.gm < 0.0 {
            self.gm = src.gm();
        }
        if self.flattening < 0.0 {
            self.flattening = src.flattening();
        }
        if self.equatorial_radius < 0.0 {
            self.equatorial_radius = src.equatorial_radius();
        }
        if self.semi_major_radius < 0.0 {
            self.semi_major_radius = src.semi_major_radius();
        }
        if let Some(p_axb) = self.parent_axb {
            if p_axb == -1 {
                self.parent_axb = src.parent_axb_id();
            }
        }
        if let Some(p_exb) = self.parent_exb {
            if p_exb == -1 {
                self.parent_exb = src.parent_exb_id();
            }
        }
    }

    pub fn as_frame(&self) -> Frame {
        Frame::Geoid {
            axb_id: self.axb,
            exb_id: self.exb,
            gm: self.gm,
            parent_axb_id: self.parent_axb,
            parent_exb_id: self.parent_exb,
            flattening: self.flattening,
            equatorial_radius: self.equatorial_radius,
            semi_major_radius: self.semi_major_radius,
        }
    }
}

#[derive(Clone, Deserialize)]
pub struct RotationToml {
    right_asc: String,
    declin: String,
    w: String,
    angle_unit: Option<String>,
    context: Option<HashMap<String, f64>>,
}

impl RotationToml {
    pub fn to_euler3_axis_dt(&self) -> Euler3AxisDt {
        let right_asc: Expr = self.right_asc.parse().unwrap();
        let declin: Expr = self.declin.parse().unwrap();
        let w_expr: Expr = self.w.parse().unwrap();

        Euler3AxisDt::from_ra_dec_w(
            right_asc,
            declin,
            w_expr,
            match &self.context {
                Some(ctx) => ctx.clone(),
                None => HashMap::new(),
            },
            match &self.angle_unit {
                Some(val) => AngleUnit::from_str(val.as_str()).unwrap(),
                None => AngleUnit::Degrees,
            },
        )
    }
}

#[test]
fn test_deser_frame_toml() {
    let frames: FramesToml = toml::from_str(
        r#"
        [frames.iau_sun]
        axb = 10
        exb = 10
        gm = 132712440041.93938
        parent_axb = 0
        parent_exb = 0
        flattening = 0
        equatorial_radius = 696342.0
        semi_major_radius = 696342.0
        [frames.iau_sun.rotation]
        right_asc = "289.13"
        declin = "63.87"
        w = "84.176 + 14.18440000*d"
        angle_unit = "degrees"
        [frames.iau_sun.rotation.context]
        t_prime = 1.0  # For example

        [frames.iau_sun2]
        axb = 10
        exb = 10
        clonefrom = "Sun J2000"
        gm = -1
        parent_axb = -1
        parent_exb = -1
        flattening = -1
        equatorial_radius = -1
        semi_major_radius = -1
        [frames.iau_sun2.rotation]
        right_asc = "289.13"
        declin = "63.87"
        w = "84.176 + 14.18440000*d"
        angle_unit = "degrees"
        [frames.iau_sun2.rotation.context]
        t_prime = 1.0  # For example
    "#,
    )
    .unwrap();

    let iau_sun = &frames.frames["iau_sun"];

    assert_eq!(iau_sun.axb, 10);
    assert_eq!(iau_sun.exb, 10);
    assert!((iau_sun.gm - 132_712_440_041.939_38).abs() < std::f64::EPSILON);
    assert_eq!(iau_sun.parent_axb.unwrap(), 0);
    assert_eq!(iau_sun.parent_exb.unwrap(), 0);
    assert!((iau_sun.equatorial_radius - 696_342.0).abs() < std::f64::EPSILON);
    assert!((iau_sun.semi_major_radius - 696_342.0).abs() < std::f64::EPSILON);

    let iau_sun_rot = &iau_sun.rotation;
    assert_eq!(iau_sun_rot.right_asc, "289.13");
    assert_eq!(iau_sun_rot.declin, "63.87");
    assert_eq!(iau_sun_rot.w, "84.176 + 14.18440000*d");
    assert_eq!(iau_sun_rot.angle_unit.as_ref().unwrap(), "degrees");

    assert!((iau_sun_rot.context.as_ref().unwrap()["t_prime"] - 1.0).abs() < std::f64::EPSILON);

    // And test the clonable frame
    let iau_sun = &frames.frames["iau_sun2"];

    assert_eq!(iau_sun.axb, 10);
    assert_eq!(iau_sun.exb, 10);
    assert_eq!(iau_sun.clonefrom.as_ref().unwrap(), "Sun J2000");
    assert!((iau_sun.gm - -1.0).abs() < std::f64::EPSILON);
    assert_eq!(iau_sun.parent_axb.unwrap(), -1);
    assert_eq!(iau_sun.parent_exb.unwrap(), -1);
    assert!((iau_sun.equatorial_radius - -1.0).abs() < std::f64::EPSILON);
    assert!((iau_sun.semi_major_radius - -1.0).abs() < std::f64::EPSILON);

    let iau_sun_rot = &iau_sun.rotation;
    assert_eq!(iau_sun_rot.right_asc, "289.13");
    assert_eq!(iau_sun_rot.declin, "63.87");
    assert_eq!(iau_sun_rot.w, "84.176 + 14.18440000*d");
    assert_eq!(iau_sun_rot.angle_unit.as_ref().unwrap(), "degrees");

    assert!((iau_sun_rot.context.as_ref().unwrap()["t_prime"] - 1.0).abs() < std::f64::EPSILON);
}
