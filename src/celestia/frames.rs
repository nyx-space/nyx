pub use celestia::xb::Identifier;
use std::fmt;

pub trait Frame: Clone + fmt::Debug + fmt::Display {
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

#[derive(Clone, Debug)]
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
