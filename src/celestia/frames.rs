pub use celestia::exb::ephemeris::Identifier as EXBID;
pub use celestia::fxb::frame::Identifier as FrameID;

pub trait Body {
    fn frame_id(&self) -> &FrameID;
}

#[derive(Clone, Debug)]
pub struct Geoid {
    pub frame_id: FrameID,
    pub gm: f64,
    pub flattening: f64,
    pub equatorial_radius: f64,
    pub semi_major_radius: f64,
}

impl Geoid {
    pub fn perfect_sphere(frame_id: FrameID, gm: f64) -> Geoid {
        Geoid {
            frame_id,
            gm,
            flattening: 0.0,
            equatorial_radius: 0.0,
            semi_major_radius: 0.0,
        }
    }
}

impl Body for Geoid {
    fn frame_id(&self) -> &FrameID {
        &self.frame_id
    }
}

pub struct Spacecraft {
    pub frame_id: FrameID,
}

impl Body for Spacecraft {
    fn frame_id(&self) -> &FrameID {
        &self.frame_id
    }
}
