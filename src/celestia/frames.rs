pub use celestia::exb::ephemeris::Identifier as EXBID;
pub use celestia::fxb::frame::Identifier as FrameID;
use std::fmt;

pub trait Frame: Clone + fmt::Debug + fmt::Display {
    fn center_id(&self) -> &i32; // Returns the integer of the center of this frame
    fn orientation_id(&self) -> &i32; // Returns the orientation (1: J2000)
}

#[derive(Clone, Copy, Debug)]
pub struct Geoid {
    pub center_id: i32,
    pub orientation_id: i32,
    pub gm: f64,
    pub flattening: f64,
    pub equatorial_radius: f64,
    pub semi_major_radius: f64,
}

impl Geoid {
    pub fn perfect_sphere(center_id: i32, orientation_id: i32, gm: f64) -> Geoid {
        Geoid {
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
        write!(f, "Geoid {}{:05}", self.orientation_id, self.center_id)
    }
}

impl Frame for Geoid {
    fn center_id(&self) -> &i32 {
        &self.center_id
    }

    fn orientation_id(&self) -> &i32 {
        &self.orientation_id
    }
}

#[derive(Clone, Debug)]
pub struct Spacecraft {
    pub center_id: i32,
    pub orientation_id: i32,
}

impl Frame for Spacecraft {
    fn center_id(&self) -> &i32 {
        &self.center_id
    }

    fn orientation_id(&self) -> &i32 {
        &self.orientation_id
    }
}

impl fmt::Display for Spacecraft {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "SC {}{:05}", self.orientation_id, self.center_id)
    }
}
