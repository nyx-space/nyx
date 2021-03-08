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

extern crate approx;
extern crate hifitime;
extern crate serde;

use super::hyperdual::linalg::norm;
use super::hyperdual::{Float, Hyperdual};
use super::na::{Vector3, U7};
use super::{Frame, Orbit};
use crate::time::Epoch;
use crate::{NyxError, TimeTagged};
use std::f64::consts::PI;
use std::f64::EPSILON;

/// If an orbit has an eccentricity below the following value, it is considered circular (only affects warning messages)
pub const ECC_EPSILON: f64 = 1e-11;

/// Orbit defines an orbital state
///
/// Unless noted otherwise, algorithms are from GMAT 2016a [StateConversionUtil.cpp](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/StateConversionUtil.cpp).
/// Regardless of the constructor used, this struct stores all the state information in Cartesian coordinates
/// as these are always non singular.
/// _Note:_ although not yet supported, this struct may change once True of Date or other nutation frames
/// are added to the toolkit.
#[derive(Copy, Clone, Debug)]
pub struct OrbitDual {
    /// in km
    pub x: Hyperdual<f64, U7>,
    /// in km
    pub y: Hyperdual<f64, U7>,
    /// in km
    pub z: Hyperdual<f64, U7>,
    /// in km/s
    pub vx: Hyperdual<f64, U7>,
    /// in km/s
    pub vy: Hyperdual<f64, U7>,
    /// in km/s
    pub vz: Hyperdual<f64, U7>,
    pub dt: Epoch,
    /// Frame contains everything we need to compute state information
    pub frame: Frame,
}

impl From<Orbit> for OrbitDual {
    /// Initialize a new OrbitDual from an orbit, no other initializers
    fn from(orbit: Orbit) -> Self {
        Self {
            x: Hyperdual::from_slice(&[orbit.x, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            y: Hyperdual::from_slice(&[orbit.y, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
            z: Hyperdual::from_slice(&[orbit.z, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]),
            vx: Hyperdual::from_slice(&[orbit.vx, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
            vy: Hyperdual::from_slice(&[orbit.vy, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
            vz: Hyperdual::from_slice(&[orbit.vz, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
            dt: orbit.dt,
            frame: orbit.frame,
        }
    }
}

impl OrbitDual {
    /// Returns the magnitude of the radius vector in km
    pub fn rmag(&self) -> Hyperdual<f64, U7> {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }

    /// Returns the magnitude of the velocity vector in km/s
    pub fn vmag(&self) -> Hyperdual<f64, U7> {
        (self.vx.powi(2) + self.vy.powi(2) + self.vz.powi(2)).sqrt()
    }

    /// Returns the radius vector of this Orbit in [km, km, km]
    pub fn radius(&self) -> Vector3<Hyperdual<f64, U7>> {
        Vector3::new(self.x, self.y, self.z)
    }

    /// Returns the velocity vector of this Orbit in [km/s, km/s, km/s]
    pub fn velocity(&self) -> Vector3<Hyperdual<f64, U7>> {
        Vector3::new(self.vx, self.vy, self.vz)
    }

    /// Returns the orbital momentum vector
    pub fn hvec(&self) -> Vector3<Hyperdual<f64, U7>> {
        self.radius().cross(&self.velocity())
    }

    /// Returns the orbital momentum value on the X axis
    pub fn hx(&self) -> Hyperdual<f64, U7> {
        self.hvec()[0]
    }

    /// Returns the orbital momentum value on the Y axis
    pub fn hy(&self) -> Hyperdual<f64, U7> {
        self.hvec()[1]
    }

    /// Returns the orbital momentum value on the Z axis
    pub fn hz(&self) -> Hyperdual<f64, U7> {
        self.hvec()[2]
    }

    /// Returns the norm of the orbital momentum
    pub fn hmag(&self) -> Hyperdual<f64, U7> {
        norm(&self.hvec())
    }

    /// Returns the specific mechanical energy
    pub fn energy(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Geoid { gm, .. } | Frame::Celestial { gm, .. } => {
                self.vmag().powi(2) / Hyperdual::from(2.0) - Hyperdual::from(gm) / self.rmag()
            }
            _ => panic!("orbital energy not defined in this frame"),
        }
    }

    /// Returns the semi-major axis in km
    pub fn sma(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Geoid { gm, .. } | Frame::Celestial { gm, .. } => {
                -Hyperdual::from(gm) / (Hyperdual::from(2.0) * self.energy())
            }
            _ => panic!("sma not defined in this frame"),
        }
    }

    /// Returns the eccentricity vector (no unit)
    pub fn evec(&self) -> Vector3<Hyperdual<f64, U7>> {
        match self.frame {
            Frame::Geoid { gm, .. } | Frame::Celestial { gm, .. } => {
                let r = self.radius();
                let v = self.velocity();
                let hgm = Hyperdual::from(gm);
                // Split up this operation because it doesn't seem to be implemented
                // ((norm(&v).powi(2) - hgm / norm(&r)) * r - (r.dot(&v)) * v) / hgm
                Vector3::new(
                    ((norm(&v).powi(2) - hgm / norm(&r)) * r[0] - (r.dot(&v)) * v[0]) / hgm,
                    ((norm(&v).powi(2) - hgm / norm(&r)) * r[1] - (r.dot(&v)) * v[1]) / hgm,
                    ((norm(&v).powi(2) - hgm / norm(&r)) * r[2] - (r.dot(&v)) * v[2]) / hgm,
                )
            }
            _ => panic!("eccentricity not defined in this frame"),
        }
    }

    /// Returns the eccentricity (no unit)
    pub fn ecc(&self) -> Hyperdual<f64, U7> {
        norm(&self.evec())
    }

    /// Returns the inclination in degrees
    pub fn inc(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                (self.hvec()[(2, 0)] / self.hmag()).acos().to_degrees()
            }
            _ => panic!("inclination not defined in this frame"),
        }
    }

    /// Returns the argument of periapsis in degrees
    pub fn aop(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                let n = Vector3::new(
                    Hyperdual::from(0.0),
                    Hyperdual::from(0.0),
                    Hyperdual::from(1.0),
                )
                .cross(&self.hvec());
                let aop = (n.dot(&self.evec()) / (norm(&n) * self.ecc())).acos();
                if aop.is_nan() {
                    error!("AoP is NaN");
                    Hyperdual::from(0.0)
                } else if self.evec()[2].real() < 0.0 {
                    (Hyperdual::from(2.0 * PI) - aop).to_degrees()
                } else {
                    aop.to_degrees()
                }
            }
            _ => panic!("aop not defined in this frame"),
        }
    }

    /// Returns the right ascension of ther ascending node in degrees
    pub fn raan(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                let n = Vector3::new(
                    Hyperdual::from(0.0),
                    Hyperdual::from(0.0),
                    Hyperdual::from(1.0),
                )
                .cross(&self.hvec());
                let raan = (n[(0, 0)] / norm(&n)).acos();
                if raan.is_nan() {
                    warn!("RAAN is NaN");
                    Hyperdual::from(0.0)
                } else if n[(1, 0)] < 0.0 {
                    (Hyperdual::from(2.0 * PI) - raan).to_degrees()
                } else {
                    raan.to_degrees()
                }
            }
            _ => panic!("RAAN not defined in this frame"),
        }
    }

    /// Returns the true anomaly in degrees between 0 and 360.0
    ///
    /// NOTE: This function will emit a warning stating that the TA should be avoided if in a very near circular orbit
    /// Code from https://github.com/ChristopherRabotin/GMAT/blob/80bde040e12946a61dae90d9fc3538f16df34190/src/gmatutil/util/StateConversionUtil.cpp#L6835
    pub fn ta(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                if self.ecc() < ECC_EPSILON {
                    warn!(
                        "true anomaly ill-defined for circular orbit (e = {})",
                        self.ecc()
                    );
                }
                let cos_nu = self.evec().dot(&self.radius()) / (self.ecc() * self.rmag());
                if (cos_nu.real().abs() - 1.0).abs() < EPSILON {
                    // This bug drove me crazy when writing SMD in Go in 2017.
                    if cos_nu > 1.0 {
                        Hyperdual::from(180.0)
                    } else {
                        Hyperdual::from(0.0)
                    }
                } else {
                    let ta = cos_nu.acos();
                    if ta.is_nan() {
                        warn!("TA is NaN");
                        Hyperdual::from(0.0)
                    } else if self.radius().dot(&self.velocity()) < 0.0 {
                        (Hyperdual::from(2.0 * PI) - ta).to_degrees()
                    } else {
                        ta.to_degrees()
                    }
                }
            }
            _ => panic!("true anomaly not defined in this frame"),
        }
    }

    /// Returns the true longitude in degrees
    pub fn tlong(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                // Angles already in degrees
                self.aop() + self.raan() + self.ta()
            }
            _ => panic!("true longitude not defined in this frame"),
        }
    }

    /// Returns the argument of latitude in degrees
    ///
    /// NOTE: If the orbit is near circular, the AoL will be computed from the true longitude
    /// instead of relying on the ill-defined true anomaly.
    pub fn aol(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                if self.ecc() < ECC_EPSILON {
                    self.tlong() - self.raan()
                } else {
                    self.aop() + self.ta()
                }
            }
            _ => panic!("argument of latitude not defined in this frame"),
        }
    }

    /// Returns the radius of periapsis (or perigee around Earth), in kilometers.
    pub fn periapsis(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                self.sma() * (Hyperdual::from(1.0) - self.ecc())
            }
            _ => panic!("periapsis not defined in this frame"),
        }
    }

    /// Returns the radius of apoapsis (or apogee around Earth), in kilometers.
    pub fn apoapsis(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                self.sma() * (Hyperdual::from(1.0) + self.ecc())
            }
            _ => panic!("apoapsis not defined in this frame"),
        }
    }

    /// Returns the eccentric anomaly in degrees
    ///
    /// This is a conversion from GMAT's StateConversionUtil::TrueToEccentricAnomaly
    pub fn ea(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                let (sin_ta, cos_ta) = self.ta().to_radians().sin_cos();
                let ecc_cos_ta = self.ecc() * cos_ta;
                let sin_ea = ((Hyperdual::from(1.0) - self.ecc().powi(2)).sqrt() * sin_ta)
                    / (Hyperdual::from(1.0) + ecc_cos_ta);
                let cos_ea = (self.ecc() + cos_ta) / (Hyperdual::from(1.0) + ecc_cos_ta);
                // The atan2 function is a bit confusing: https://doc.rust-lang.org/std/primitive.f64.html#method.atan2 .
                sin_ea.atan2(cos_ea).to_degrees()
            }
            _ => panic!("eccentric anomaly is not defined in this frame"),
        }
    }

    /// Returns the flight path angle in degrees
    pub fn fpa(&self) -> Hyperdual<f64, U7> {
        let nu = self.ta().to_radians();
        let ecc = self.ecc();
        let denom =
            (Hyperdual::from(1.0) + Hyperdual::from(2.0) * ecc * nu.cos() + ecc.powi(2)).sqrt();
        let sin_fpa = ecc * nu.sin() / denom;
        let cos_fpa = Hyperdual::from(1.0) + ecc * nu.cos() / denom;
        sin_fpa.atan2(cos_fpa).to_degrees()
    }

    /// Returns the mean anomaly in degrees
    ///
    /// This is a conversion from GMAT's StateConversionUtil::TrueToMeanAnomaly
    pub fn ma(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                if self.ecc().real() < 1.0 {
                    (self.ea().to_radians() - self.ecc() * self.ea().to_radians().sin())
                        .to_degrees()
                } else if self.ecc().real() > 1.0 {
                    info!("computing the hyperbolic anomaly");
                    // From GMAT's TrueToHyperbolicAnomaly
                    ((self.ta().to_radians().sin() * (self.ecc().powi(2) - Hyperdual::from(1.0)))
                        .sqrt()
                        / (Hyperdual::from(1.0) + self.ecc() * self.ta().to_radians().cos()))
                    .asinh()
                    .to_degrees()
                } else {
                    error!("parabolic orbit: setting mean anomaly to 0.0");
                    Hyperdual::from(0.0)
                }
            }
            _ => panic!("mean anomaly is not defined in this frame"),
        }
    }

    /// Returns the semi parameter (or semilatus rectum)
    pub fn semi_parameter(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                self.sma() * (Hyperdual::from(1.0) - self.ecc().powi(2))
            }
            _ => panic!("semi parameter is not defined in this frame"),
        }
    }

    /// Returns the geodetic longitude (λ) in degrees. Value is between 0 and 360 degrees.
    ///
    /// Although the reference is not Vallado, the math from Vallado proves to be equivalent.
    /// Reference: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016
    pub fn geodetic_longitude(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Geoid { .. } => self.y.atan2(self.x).to_degrees(),
            _ => panic!("geodetic elements only defined in a Geoid frame"),
        }
    }

    /// Returns the geodetic latitude (φ) in degrees. Value is between -180 and +180 degrees.
    ///
    /// Reference: Vallado, 4th Ed., Algorithm 12 page 172.
    pub fn geodetic_latitude(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Geoid {
                flattening,
                semi_major_radius,
                ..
            } => {
                let eps = 1e-12;
                let max_attempts = 20;
                let mut attempt_no = 0;
                let r_delta = (self.x.powi(2) + self.y.powi(2)).sqrt();
                let mut latitude = (self.z / self.rmag()).asin();
                let e2 = Hyperdual::from(flattening * (2.0 - flattening));
                loop {
                    attempt_no += 1;
                    let c_earth = Hyperdual::from(semi_major_radius)
                        / ((Hyperdual::from(1.0) - e2 * (latitude).sin().powi(2)).sqrt());
                    let new_latitude = (self.z + c_earth * e2 * (latitude).sin()).atan2(r_delta);
                    if (latitude - new_latitude).abs() < eps {
                        return new_latitude.to_degrees();
                    } else if attempt_no >= max_attempts {
                        error!(
                            "geodetic latitude failed to converge -- error = {}",
                            (latitude - new_latitude).abs()
                        );
                        return new_latitude.to_degrees();
                    }
                    latitude = new_latitude;
                }
            }
            _ => panic!("geodetic elements only defined in a Geoid frame"),
        }
    }

    /// Returns the geodetic height in km.
    ///
    /// Reference: Vallado, 4th Ed., Algorithm 12 page 172.
    pub fn geodetic_height(&self) -> Hyperdual<f64, U7> {
        match self.frame {
            Frame::Geoid {
                flattening,
                semi_major_radius,
                ..
            } => {
                let e2 = Hyperdual::from(flattening * (2.0 - flattening));
                let latitude = self.geodetic_latitude().to_radians();
                let sin_lat = latitude.sin();
                if (latitude - Hyperdual::from(1.0)).abs() < 0.1 {
                    // We are near poles, let's use another formulation.
                    let s_earth0: f64 = semi_major_radius * (1.0 - flattening).powi(2);
                    let s_earth = Hyperdual::from(s_earth0)
                        / ((Hyperdual::from(1.0) - e2 * sin_lat.powi(2)).sqrt());
                    self.z / latitude.sin() - s_earth
                } else {
                    let c_earth = Hyperdual::from(semi_major_radius)
                        / ((Hyperdual::from(1.0) - e2 * sin_lat.powi(2)).sqrt());
                    let r_delta = (self.x.powi(2) + self.y.powi(2)).sqrt();
                    r_delta / latitude.cos() - c_earth
                }
            }
            _ => panic!("geodetic elements only defined in a Geoid frame"),
        }
    }

    /// Returns the right ascension of this orbit in degrees
    pub fn right_ascension(&self) -> Hyperdual<f64, U7> {
        (self.y.atan2(self.x)).to_degrees()
    }

    /// Returns the declination of this orbit in degrees
    pub fn declination(&self) -> Hyperdual<f64, U7> {
        (self.z / self.rmag()).asin().to_degrees()
    }

    /// Returns the semi minor axis in km, includes code for a hyperbolic orbit
    pub fn semi_minor_axis(&self) -> Hyperdual<f64, U7> {
        if self.ecc() <= 1.0 {
            ((self.sma() * self.ecc()).powi(2) - self.sma().powi(2)).sqrt()
        } else {
            self.hmag().powi(2)
                / (Hyperdual::from(self.frame.gm())
                    * (self.ecc().powi(2) - Hyperdual::from(1.0)).sqrt())
        }
    }

    /// Returns the velocity declination of this orbit in degrees
    pub fn velocity_declination(&self) -> Hyperdual<f64, U7> {
        (self.vz / self.vmag()).asin().to_degrees()
    }

    /// Returns the $C_3$ of this orbit
    pub fn c3(&self) -> Hyperdual<f64, U7> {
        -Hyperdual::from(self.frame.gm()) / self.sma()
    }

    /// Returns the hyperbolic anomaly in degrees between 0 and 360.0
    pub fn hyperbolic_anomaly(&self) -> Result<Hyperdual<f64, U7>, NyxError> {
        if self.ecc() <= 1.0 {
            Err(NyxError::NotHyperbolic(
                "Orbit is not hyperbolic so there is no hyperbolic anomaly.".to_string(),
            ))
        } else {
            let (sin_ta, cos_ta) = self.ta().to_radians().sin_cos();
            let sinh_h = (sin_ta * (self.ecc().powi(2) - Hyperdual::from(1.0)).sqrt())
                / (Hyperdual::from(1.0) + self.ecc() * cos_ta);
            Ok(sinh_h.asinh().to_degrees())
        }
    }
}

impl TimeTagged for OrbitDual {
    fn epoch(&self) -> Epoch {
        self.dt
    }

    fn set_epoch(&mut self, epoch: Epoch) {
        self.dt = epoch
    }
}
