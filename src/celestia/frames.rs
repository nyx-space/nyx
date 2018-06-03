extern crate hifitime;
extern crate nalgebra as na;

use self::hifitime::instant::*;
use self::na::Matrix3;
use super::{EARTH, NAIF};

/// Defines a coordinate frame trait around the body B which implements the trait NAIF
pub trait CoordinateFrame {
    /// Defines the center of this frame.
    type Center: NAIF;
    /// Returns the rotation matrix (or DCM) which converts this coordinate frame to an inertial coordinate frame with the same center object.
    fn to_inertial(at: Instant) -> Matrix3<f64>;
}

pub struct ECI {}
impl CoordinateFrame for ECI {
    type Center = EARTH;
    fn to_inertial(_: Instant) -> Matrix3<f64> {
        Matrix3::identity()
    }
}
