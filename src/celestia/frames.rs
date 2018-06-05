extern crate hifitime;
extern crate nalgebra as na;
use self::hifitime::instant::{Era, Instant};
use self::hifitime::julian::ModifiedJulian;
use self::hifitime::TimeSystem;
use self::na::Matrix3;
use super::{EARTH, NAIF, SSB};
use std::fmt;

/// Defines a coordinate frame trait around the body B which implements the trait NAIF
pub trait CoordinateFrame: fmt::Display + Copy + fmt::Debug {
    /// Defines the center of this frame.
    type Center: NAIF;
    /// Returns the rotation matrix (or DCM) which converts this coordinate frame to an inertial coordinate frame with the same center object.
    fn to_inertial(at: Instant) -> Matrix3<f64>;

    /// Returns the rotation matrix (or DCM) which converts this coordinate frame back from an inertial coordinate frame with the same center object.
    fn from_inertial(at: Instant) -> Matrix3<f64> {
        Self::to_inertial(at).transpose()
    }
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

/// ECEF is the Earth Centered Earth Fixed frame. It's a decent approximation of a body fixed frame.
#[derive(Copy, Clone, Debug)]
pub struct ECEF {}

impl ECEF {
    /// Computes the Greenwich Mean Sideral Time in degrees (between 0 360) at the provided time.
    ///
    /// The computation here uses the Theta GMST computation of Algo 15 from Vallado, 4th Ed., page 188.
    pub fn theta_gmst(at: Instant) -> f64 {
        use utils::between_0_360;

        let t_ut1 = (ModifiedJulian::from_instant(at).julian_days() - 2451545.0) / 36_525.0;
        let theta_gmst_secs = (67_310.548_41 + (876_600.0 * 3_600.0 + 8_640_184.812_866) * t_ut1 + 0.093_104 * t_ut1.powi(2)
            - 6.2e-6 * t_ut1.powi(3)) % 86_400.0;
        between_0_360(theta_gmst_secs / 240.0)
    }
}

impl CoordinateFrame for ECEF {
    type Center = EARTH;
    fn to_inertial(at: Instant) -> Matrix3<f64> {
        let theta_gmst = Self::theta_gmst(at);
        let delta_secs = at - Instant::new(0, 0, Era::Present);
        // Compute theta GST
        let (s, c) = (EARTH::rotation_rate() * delta_secs).sin_cos();
        // It's an R3 rotation to go from ECI to ECEF, so let's take the transpose to go from ECEF to ECI.
        Matrix3::new(c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0).transpose()
    }
}

impl fmt::Display for ECEF {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "dECEF")
    }
}
