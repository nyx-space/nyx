extern crate hifitime;
extern crate nalgebra as na;
use self::hifitime::instant::Instant;
use self::hifitime::julian::ModifiedJulian;
use self::hifitime::TimeSystem;
use self::na::Matrix3;
use super::{EARTH, NAIF, SSB};
use std::fmt;
use utils::r3;

/// Defines a coordinate frame trait around the body B which implements the trait NAIF.
///
/// A CoordinateFrame only needs to define either to_inertial or from_inertial.
pub trait CoordinateFrame: fmt::Display + Copy + fmt::Debug {
    /// Defines the center of this frame.
    type Center: NAIF;
    /// Returns the rotation matrix (or DCM) which converts this coordinate frame to an inertial coordinate frame with the same center object.
    fn to_inertial(at: Instant) -> Matrix3<f64> {
        Self::from_inertial(at).transpose()
    }

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

/// ECEF is the Earth Centered Earth Fixed frame. It's an **approximation** of a body fixed frame which currently only account for the rotation of Earth.
#[derive(Copy, Clone, Debug)]
pub struct ECEF {}

impl ECEF {
    /// Computes the Greenwich Mean Sideral Time in degrees (between 0 360) at the provided time.
    ///
    /// The computation here uses the Theta GMST computation of Algo 15 from Vallado, 4th Ed., page 188.
    pub fn gmst(at: Instant) -> f64 {
        use utils::between_0_360;

        let t_ut1 = (ModifiedJulian::from_instant(at).julian_days() - 2_451_545.0) / 36_525.0;
        let theta_gmst_secs = (67_310.548_41 + (876_600.0 * 3_600.0 + 8_640_184.812_866) * t_ut1 + 0.093_104 * t_ut1.powi(2)
            - 6.2e-6 * t_ut1.powi(3))
            % 86_400.0;
        between_0_360(theta_gmst_secs / 240.0)
    }

    /// Computes the (approximate) Greenwich Apparent Sideral Time as per IAU2000.
    ///
    /// NOTE: This is an approximation valid to within 0.9 seconds in absolute value.
    /// In fact, hifitime does not support UT1, but according to the [IERS](https://www.iers.org/IERS/EN/Science/EarthRotation/UTC.html;jsessionid=A6E88EB4CF0FC2E1A3C10D807F51B829.live2?nn=12932)
    /// UTC with leap seconds is always within 0.9 seconds to UT1, and hifitime inherently supports leap seconds.
    /// Reference: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016
    pub fn gast(at: Instant) -> f64 {
        use std::f64::consts::PI;
        let tu = ModifiedJulian::from_instant(at).days - 51_544.5;
        2.0 * PI * (0.779_057_273_264_0 + 1.002_737_811_911_354_48 * tu)
    }
}

impl CoordinateFrame for ECEF {
    type Center = EARTH;
    /// WARNING: The conversion to and from ECEF currently ignores precession, nutation, and polar motion.
    fn from_inertial(at: Instant) -> Matrix3<f64> {
        r3(ECEF::gast(at))
    }
}

impl fmt::Display for ECEF {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ECEF")
    }
}
