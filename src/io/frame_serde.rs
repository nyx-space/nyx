extern crate toml;

use super::serde_derive::Deserialize;
use crate::celestia::Frame;
use std::collections::HashMap;

#[derive(Clone, Deserialize)]
pub struct FramesSerde {
    pub frames: HashMap<String, FrameSerde>,
}

/// A structure specifying the format of a frame defined in TOML, or some other serialization.
#[derive(Clone, Deserialize)]
pub struct FrameSerde {
    /// Refers to, or create, a unique identifier of the orientation defined by this frame.
    orientation: i32,
    /// Refers to, or create, a unique identifier of the center object defined by this frame.
    center: i32,
    // Set this in the TOML to clone the negative values of this frame from another frame
    pub inherit: Option<String>,
    gm: f64,
    /// Refers to the a unique identifier of the parent orientation of this frame (e.g. J2000 Ecliptic)
    parent_orientation: Option<i32>,
    /// Refers to the a unique identifier of the parent center object of this frame (e.g. Solar System Barycenter)
    parent_center: Option<i32>,
    flattening: f64,
    equatorial_radius: f64,
    semi_major_radius: f64,
    pub rotation: RotationToml,
}

impl FrameSerde {
    pub fn to_frame(&self) -> Frame {
        Frame::Geoid {
            axb_id: self.orientation,
            exb_id: self.center,
            gm: self.gm,
            parent_axb_id: self.parent_orientation,
            parent_exb_id: self.parent_center,
            flattening: self.flattening,
            equatorial_radius: self.equatorial_radius,
            semi_major_radius: self.semi_major_radius,
            frame_path: [None, None, None],
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
        if let Some(p_axb) = self.parent_orientation {
            if p_axb == -1 {
                self.parent_orientation = src.parent_axb_id();
            }
        }
        if let Some(p_exb) = self.parent_center {
            if p_exb == -1 {
                self.parent_center = src.parent_exb_id();
            }
        }
    }

    pub fn as_frame(&self) -> Frame {
        Frame::Geoid {
            axb_id: self.orientation,
            exb_id: self.center,
            gm: self.gm,
            parent_axb_id: self.parent_orientation,
            parent_exb_id: self.parent_center,
            flattening: self.flattening,
            equatorial_radius: self.equatorial_radius,
            semi_major_radius: self.semi_major_radius,
            frame_path: [None, None, None],
        }
    }
}

#[derive(Clone, Deserialize)]
pub struct RotationToml {
    pub right_asc: String,
    pub declin: String,
    pub w: String,
    pub angle_unit: Option<String>,
    pub context: Option<HashMap<String, String>>,
}

#[test]
fn test_deser_frame_toml() {
    let frames: FramesSerde = toml::from_str(
        r#"
        [frames.iau_sun]
        orientation = 10
        center = 10
        gm = 132712440041.93938
        parent_orientation = 0
        parent_center = 0
        flattening = 0
        equatorial_radius = 696342.0
        semi_major_radius = 696342.0
        [frames.iau_sun.rotation]
        right_asc = "289.13"
        declin = "63.87"
        w = "84.176 + 14.18440000*d"
        angle_unit = "degrees"
        [frames.iau_sun.rotation.context]
        t_prime = "1.0"  # Must be encasted in quote even if just a floating point value

        [frames.iau_sun2]
        orientation = 10
        center = 10
        inherit = "Sun J2000"
        gm = -1
        parent_orientation = -1
        parent_center = -1
        flattening = -1
        equatorial_radius = -1
        semi_major_radius = -1
        [frames.iau_sun2.rotation]
        right_asc = "289.13"
        declin = "63.87"
        w = "84.176 + 14.18440000*d"
        angle_unit = "degrees"
    "#,
    )
    .unwrap();

    let iau_sun = &frames.frames["iau_sun"];

    assert_eq!(iau_sun.orientation, 10);
    assert_eq!(iau_sun.center, 10);
    assert!((iau_sun.gm - 132_712_440_041.939_38).abs() < std::f64::EPSILON);
    assert_eq!(iau_sun.parent_orientation.unwrap(), 0);
    assert_eq!(iau_sun.parent_center.unwrap(), 0);
    assert!((iau_sun.equatorial_radius - 696_342.0).abs() < std::f64::EPSILON);
    assert!((iau_sun.semi_major_radius - 696_342.0).abs() < std::f64::EPSILON);

    let iau_sun_rot = &iau_sun.rotation;
    assert_eq!(iau_sun_rot.right_asc, "289.13");
    assert_eq!(iau_sun_rot.declin, "63.87");
    assert_eq!(iau_sun_rot.w, "84.176 + 14.18440000*d");
    assert_eq!(iau_sun_rot.angle_unit.as_ref().unwrap(), "degrees");

    // And test the clonable frame
    let iau_sun = &frames.frames["iau_sun2"];

    assert_eq!(iau_sun.orientation, 10);
    assert_eq!(iau_sun.center, 10);
    assert_eq!(iau_sun.inherit.as_ref().unwrap(), "Sun J2000");
    assert!((iau_sun.gm - -1.0).abs() < std::f64::EPSILON);
    assert_eq!(iau_sun.parent_orientation.unwrap(), -1);
    assert_eq!(iau_sun.parent_center.unwrap(), -1);
    assert!((iau_sun.equatorial_radius - -1.0).abs() < std::f64::EPSILON);
    assert!((iau_sun.semi_major_radius - -1.0).abs() < std::f64::EPSILON);

    let iau_sun_rot = &iau_sun.rotation;
    assert_eq!(iau_sun_rot.right_asc, "289.13");
    assert_eq!(iau_sun_rot.declin, "63.87");
    assert_eq!(iau_sun_rot.w, "84.176 + 14.18440000*d");
    assert_eq!(iau_sun_rot.angle_unit.as_ref().unwrap(), "degrees");
}
