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

extern crate toml;

use super::serde_derive::Deserialize;
use crate::cosmic::Frame;
use std::collections::HashMap;

#[derive(Clone, Deserialize)]
pub struct FramesSerde {
    pub frames: HashMap<String, FrameSerde>,
}

/// A structure specifying the format of a frame defined in TOML, or some other serialization.
#[derive(Clone, Deserialize)]
pub struct FrameSerde {
    // Set this in the TOML to clone the negative values of this frame from another frame
    pub inherit: Option<String>,
    gm: f64,
    flattening: f64,
    equatorial_radius: f64,
    semi_major_radius: f64,
    pub rotation: RotationToml,
}

impl FrameSerde {
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
    }

    pub fn as_frame(&self) -> Frame {
        Frame::Geoid {
            gm: self.gm,
            flattening: self.flattening,
            equatorial_radius: self.equatorial_radius,
            semi_major_radius: self.semi_major_radius,
            ephem_path: [None, None, None],
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
        gm = 132712440041.93938
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
        inherit = "Sun J2000"
        gm = -1
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

    assert!((iau_sun.gm - 132_712_440_041.939_38).abs() < std::f64::EPSILON);
    assert!((iau_sun.equatorial_radius - 696_342.0).abs() < std::f64::EPSILON);
    assert!((iau_sun.semi_major_radius - 696_342.0).abs() < std::f64::EPSILON);

    let iau_sun_rot = &iau_sun.rotation;
    assert_eq!(iau_sun_rot.right_asc, "289.13");
    assert_eq!(iau_sun_rot.declin, "63.87");
    assert_eq!(iau_sun_rot.w, "84.176 + 14.18440000*d");
    assert_eq!(iau_sun_rot.angle_unit.as_ref().unwrap(), "degrees");

    // And test the clonable frame
    let iau_sun = &frames.frames["iau_sun2"];

    assert_eq!(iau_sun.inherit.as_ref().unwrap(), "Sun J2000");
    assert!((iau_sun.gm - -1.0).abs() < std::f64::EPSILON);
    assert!((iau_sun.equatorial_radius - -1.0).abs() < std::f64::EPSILON);
    assert!((iau_sun.semi_major_radius - -1.0).abs() < std::f64::EPSILON);

    let iau_sun_rot = &iau_sun.rotation;
    assert_eq!(iau_sun_rot.right_asc, "289.13");
    assert_eq!(iau_sun_rot.declin, "63.87");
    assert_eq!(iau_sun_rot.w, "84.176 + 14.18440000*d");
    assert_eq!(iau_sun_rot.angle_unit.as_ref().unwrap(), "degrees");
}
