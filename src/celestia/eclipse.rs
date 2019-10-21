use super::hifitime::Epoch;
use super::na::Vector3;
use super::{Cosm, Geoid, State};

/// Stores the eclipse state
#[derive(Clone, Copy, Debug)]
pub enum EclipseState {
    Umbra,
    Penumbra,
    Visibilis,
}

/// Computes the state of an eclipse at the provided time.
/// Warning: Direction must be in the same frame as the Geoid
///
/// For example, if computing whether an Earth-orbiting spacecraft can see the Sun despite the Earth,
/// the direction should be the direction of the spaceraft to the Sun, and the Geoid should be the Earth.
pub fn eclipse_state(dir: Vector3, dir_frame: i32, dt: Epoch, geoid: Geoid, cosm: &Cosm) -> EclipseState {
    let unit_dir = direction / direction.norm();
    let eclipse_body = self.cosm.unwrap().geoid_from_id(geoid.id());
    // State of the eclipsing body as seen from primary body
    let body_state = self.cosm.unwrap().celestial_state(*exb_id, dt.as_jde_et_days(), dir_frame);
    // let discriminant =
}
