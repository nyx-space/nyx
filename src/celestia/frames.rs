extern crate hifitime;
extern crate nalgebra as na;
use self::hifitime::instant::{Era, Instant};
use self::na::Matrix3;
use super::{EARTH, NAIF, SSB};
use std::fmt;

/// Defines a coordinate frame trait around the body B which implements the trait NAIF
pub trait CoordinateFrame: fmt::Display + Copy + fmt::Debug {
    /// Defines the center of this frame.
    type Center: NAIF;
    /// Returns the rotation matrix (or DCM) which converts this coordinate frame to an inertial coordinate frame with the same center object.
    fn to_inertial(at: Instant) -> Matrix3<f64>;
}

#[derive(Copy, Clone, Debug)]
pub struct ICRF {}
impl CoordinateFrame for ICRF {
    type Center = SSB;
    fn to_inertial(_: Instant) -> Matrix3<f64> {
        Matrix3::identity()
    }
}

impl fmt::Display for ICRF {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ICRF")
    }
}

#[derive(Copy, Clone, Debug)]
pub struct ECI {}
impl CoordinateFrame for ECI {
    type Center = EARTH;
    fn to_inertial(_: Instant) -> Matrix3<f64> {
        Matrix3::identity()
    }
}

impl fmt::Display for ECI {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ECI")
    }
}

/// DummyECEF is a significantly simplified version of the ECEF frame. It uses only the mean Earth rotation rate
/// and may be initialized at an initial date time.
#[derive(Copy, Clone, Debug)]
pub struct DummyECEF {}

impl CoordinateFrame for DummyECEF {
    type Center = EARTH;
    fn to_inertial(at: Instant) -> Matrix3<f64> {
        let earth_rot_rate = 7.292115900231276e-5; // in radians per second.
        let delta_secs = at - Instant::new(0, 0, Era::Present);
        let (s, c) = (earth_rot_rate * delta_secs).sin_cos(); // Compute theta GST
        Matrix3::new(c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0)
    }
}

impl fmt::Display for DummyECEF {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "dECEF")
    }
}
