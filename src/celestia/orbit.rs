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

use self::approx::{abs_diff_eq, relative_eq};
use self::serde::ser::SerializeStruct;
use self::serde::{Serialize, Serializer};
use super::na::{Matrix3, Matrix6, Vector3, Vector6};
use super::{BPlane, Frame};
use crate::time::{Duration, Epoch, TimeUnit};
use crate::utils::{between_0_360, between_pm_180, perpv, r1, r3};
use crate::{NyxError, TimeTagged};
use std::f64::consts::PI;
use std::f64::EPSILON;
use std::fmt;
use std::ops::{Add, Neg, Sub};

/// If an orbit has an eccentricity below the following value, it is considered circular (only affects warning messages)
pub const ECC_EPSILON: f64 = 1e-11;

pub fn assert_orbit_eq_or_abs<'a>(left: &Orbit, right: &Orbit, epsilon: f64, msg: &'a str) {
    if !(left.to_cartesian_vec() == right.to_cartesian_vec())
        && !abs_diff_eq!(
            left.to_cartesian_vec(),
            right.to_cartesian_vec(),
            epsilon = epsilon
        )
        && left.epoch() != right.epoch()
    {
        panic!(
            r#"assertion failed: `(left == right)`
  left: `{:?}`,
 right: `{:?}`: {}"#,
            left.to_cartesian_vec(),
            right.to_cartesian_vec(),
            msg
        )
    }
}

pub fn assert_orbit_eq_or_rel<'a>(left: &Orbit, right: &Orbit, epsilon: f64, msg: &'a str) {
    if !(left.to_cartesian_vec() == right.to_cartesian_vec())
        && !relative_eq!(
            left.to_cartesian_vec(),
            right.to_cartesian_vec(),
            max_relative = epsilon
        )
        && left.epoch() != right.epoch()
    {
        panic!(
            r#"assertion failed: `(left == right)`
  left: `{:?}`,
 right: `{:?}`: {}"#,
            left.to_cartesian_vec(),
            right.to_cartesian_vec(),
            msg
        )
    }
}

/// Defines the kind of state transition matrix stored in the orbit.
#[derive(Copy, Clone, Debug)]
pub enum StmKind {
    /// For navigation (default): Corresponds to the linearization between the previous step of the propagator and the new step ($\Phi_k$)
    Step,
    /// For trajectory optimization: Corresponds to the linearization between from the first state until the current state ($\Phi_{k\to0}$)
    Traj,
    /// For propagation: no STM is set
    Unset,
}

/// Orbit defines an orbital state
///
/// Unless noted otherwise, algorithms are from GMAT 2016a [StateConversionUtil.cpp](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/StateConversionUtil.cpp).
/// Regardless of the constructor used, this struct stores all the state information in Cartesian coordinates
/// as these are always non singular.
/// _Note:_ although not yet supported, this struct may change once True of Date or other nutation frames
/// are added to the toolkit.
#[derive(Copy, Clone, Debug)]
pub struct Orbit {
    /// in km
    pub x: f64,
    /// in km
    pub y: f64,
    /// in km
    pub z: f64,
    /// in km/s
    pub vx: f64,
    /// in km/s
    pub vy: f64,
    /// in km/s
    pub vz: f64,
    pub dt: Epoch,
    /// Frame contains everything we need to compute state information
    pub frame: Frame,
    /// Optionally stores an STM
    pub stm: Option<Matrix6<f64>>,
    pub stm_kind: StmKind,
}

impl Orbit {
    /// Creates a new Orbit in the provided frame at the provided Epoch.
    ///
    /// **Units:** km, km, km, km/s, km/s, km/s
    pub fn cartesian(
        x: f64,
        y: f64,
        z: f64,
        vx: f64,
        vy: f64,
        vz: f64,
        dt: Epoch,
        frame: Frame,
    ) -> Self {
        Orbit {
            x,
            y,
            z,
            vx,
            vy,
            vz,
            dt,
            frame,
            stm: None,
            stm_kind: StmKind::Unset,
        }
    }

    /// Creates a new Orbit and initializes its STM.
    pub fn cartesian_stm(
        x: f64,
        y: f64,
        z: f64,
        vx: f64,
        vy: f64,
        vz: f64,
        dt: Epoch,
        frame: Frame,
    ) -> Self {
        let mut me = Self::cartesian(x, y, z, vx, vy, vz, dt, frame);
        me.enable_stm();
        me
    }

    /// Creates a new Orbit in the provided frame at the provided Epoch in time with 0.0 velocity.
    ///
    /// **Units:** km, km, km
    pub fn from_position(x: f64, y: f64, z: f64, dt: Epoch, frame: Frame) -> Self {
        Orbit {
            x,
            y,
            z,
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
            dt,
            frame,
            stm: None,
            stm_kind: StmKind::Unset,
        }
    }

    /// Creates a new Orbit around in the provided frame from the borrowed state vector
    ///
    /// The state vector **must** be x, y, z, vx, vy, vz. This function is a shortcut to `cartesian`
    /// and as such it has the same unit requirements.
    pub fn cartesian_vec(state: &Vector6<f64>, dt: Epoch, frame: Frame) -> Self {
        Orbit {
            x: state[0],
            y: state[1],
            z: state[2],
            vx: state[3],
            vy: state[4],
            vz: state[5],
            dt,
            frame,
            stm: None,
            stm_kind: StmKind::Unset,
        }
    }

    /// Returns the magnitude of the radius vector in km
    pub fn rmag(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }

    /// Returns the magnitude of the velocity vector in km/s
    pub fn vmag(&self) -> f64 {
        (self.vx.powi(2) + self.vy.powi(2) + self.vz.powi(2)).sqrt()
    }

    /// Returns the radius vector of this Orbit in [km, km, km]
    pub fn radius(&self) -> Vector3<f64> {
        Vector3::new(self.x, self.y, self.z)
    }

    /// Returns the velocity vector of this Orbit in [km/s, km/s, km/s]
    pub fn velocity(&self) -> Vector3<f64> {
        Vector3::new(self.vx, self.vy, self.vz)
    }

    /// Returns this state as a Cartesian Vector6 in [km, km, km, km/s, km/s, km/s]
    ///
    /// Note that the time is **not** returned in the vector.
    pub fn to_cartesian_vec(&self) -> Vector6<f64> {
        Vector6::new(self.x, self.y, self.z, self.vx, self.vy, self.vz)
    }

    /// Returns the distance in kilometers between this state and another state.
    /// Will **panic** is the frames are different
    pub fn distance_to(&self, other: &Orbit) -> f64 {
        assert_eq!(
            self.frame, other.frame,
            "cannot compute the distance between two states in different frames"
        );
        self.distance_to_point(&other.radius())
    }

    /// Returns the distance in kilometers between this state and a point assumed to be in the same frame.
    pub fn distance_to_point(&self, other: &Vector3<f64>) -> f64 {
        ((self.x - other.x).powi(2) + (self.y - other.y).powi(2) + (self.z - other.z).powi(2))
            .sqrt()
    }

    /// Returns the unit vector in the direction of the state radius
    pub fn r_hat(&self) -> Vector3<f64> {
        self.radius() / self.rmag()
    }

    /// Returns the unit vector in the direction of the state velocity
    pub fn v_hat(&self) -> Vector3<f64> {
        perpv(&self.velocity(), &self.r_hat()) / self.rmag()
    }

    /// Creates a new Orbit around the provided Celestial or Geoid frame from the Keplerian orbital elements.
    ///
    /// **Units:** km, none, degrees, degrees, degrees, degrees
    ///
    /// WARNING: This function will panic if the singularities in the conversion are expected.
    /// NOTE: The state is defined in Cartesian coordinates as they are non-singular. This causes rounding
    /// errors when creating a state from its Keplerian orbital elements (cf. the state tests).
    /// One should expect these errors to be on the order of 1e-12.
    pub fn keplerian(
        sma: f64,
        ecc: f64,
        inc: f64,
        raan: f64,
        aop: f64,
        ta: f64,
        dt: Epoch,
        frame: Frame,
    ) -> Self {
        match frame {
            Frame::Geoid { gm, .. } | Frame::Celestial { gm, .. } => {
                if gm.abs() < std::f64::EPSILON {
                    warn!(
                        "GM is near zero ({}): expect math errors in Keplerian to Cartesian conversion",
                        gm
                    );
                }
                // Algorithm from GMAT's StateConversionUtil::KeplerianToCartesian
                let ecc = if ecc < 0.0 {
                    warn!("eccentricity cannot be negative: sign of eccentricity changed");
                    ecc * -1.0
                } else {
                    ecc
                };
                let sma = if ecc > 1.0 && sma > 0.0 {
                    warn!("eccentricity > 1 (hyperbolic) BUT SMA > 0 (elliptical): sign of SMA changed");
                    sma * -1.0
                } else if ecc < 1.0 && sma < 0.0 {
                    warn!("eccentricity < 1 (elliptical) BUT SMA < 0 (hyperbolic): sign of SMA changed");
                    sma * -1.0
                } else {
                    sma
                };
                if (sma * (1.0 - ecc)).abs() < 1e-3 {
                    // GMAT errors below one meter. Let's warn for below that, but not panic, might be useful for landing scenarios?
                    warn!("radius of periapsis is less than one meter");
                }
                if (1.0 - ecc).abs() < EPSILON {
                    panic!("parabolic orbits have ill-defined Keplerian orbital elements");
                }
                if ecc > 1.0 {
                    let ta = between_0_360(ta);
                    if ta > (PI - (1.0 / ecc).acos()).to_degrees() {
                        panic!(
                            "true anomaly value ({}) physically impossible for a hyperbolic orbit",
                            ta
                        );
                    }
                }
                if (1.0 + ecc * ta.to_radians().cos()).is_infinite() {
                    panic!("radius of orbit is infinite");
                }
                // Done with all the warnings and errors supported by GMAT
                // The conversion algorithm itself comes from GMAT's StateConversionUtil::ComputeKeplToCart
                // NOTE: GMAT supports mean anomaly instead of true anomaly, but only for backward compatibility reasons
                // so it isn't supported here.
                let inc = inc.to_radians();
                let raan = raan.to_radians();
                let aop = aop.to_radians();
                let ta = ta.to_radians();
                let p = sma * (1.0 - ecc.powi(2));
                if p.abs() < EPSILON {
                    panic!("Semilatus rectum ~= 0.0: parabolic orbit");
                }
                // NOTE: At this point GMAT computes 1+ecc**2 and checks whether it's very small.
                // It then reports that the radius may be too large. We've effectively already done
                // this check above (and panicked if needed), so it isn't repeated here.
                let radius = p / (1.0 + ecc * ta.cos());
                let (sin_aop_ta, cos_aop_ta) = (aop + ta).sin_cos();
                let (sin_inc, cos_inc) = inc.sin_cos();
                let (sin_raan, cos_raan) = raan.sin_cos();
                let (sin_aop, cos_aop) = aop.sin_cos();
                let x = radius * (cos_aop_ta * cos_raan - cos_inc * sin_aop_ta * sin_raan);
                let y = radius * (cos_aop_ta * sin_raan + cos_inc * sin_aop_ta * cos_raan);
                let z = radius * sin_aop_ta * sin_inc;
                let sqrt_gm_p = (gm / p).sqrt();
                let cos_ta_ecc = ta.cos() + ecc;
                let sin_ta = ta.sin();

                let vx =
                    sqrt_gm_p * cos_ta_ecc * (-sin_aop * cos_raan - cos_inc * sin_raan * cos_aop)
                        - sqrt_gm_p * sin_ta * (cos_aop * cos_raan - cos_inc * sin_raan * sin_aop);
                let vy =
                    sqrt_gm_p * cos_ta_ecc * (-sin_aop * sin_raan + cos_inc * cos_raan * cos_aop)
                        - sqrt_gm_p * sin_ta * (cos_aop * sin_raan + cos_inc * cos_raan * sin_aop);
                let vz = sqrt_gm_p * (cos_ta_ecc * sin_inc * cos_aop - sin_ta * sin_inc * sin_aop);
                Orbit {
                    x,
                    y,
                    z,
                    vx,
                    vy,
                    vz,
                    dt,
                    frame,
                    stm: None,
                    stm_kind: StmKind::Unset,
                }
            }
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    /// Creates a new Orbit from the provided semi-major axis altitude in kilometers
    pub fn keplerian_alt(
        sma_altitude: f64,
        ecc: f64,
        inc: f64,
        raan: f64,
        aop: f64,
        ta: f64,
        dt: Epoch,
        frame: Frame,
    ) -> Self {
        Self::keplerian(
            sma_altitude + frame.equatorial_radius(),
            ecc,
            inc,
            raan,
            aop,
            ta,
            dt,
            frame,
        )
    }

    /// Creates a new Orbit around the provided frame from the borrowed state vector
    ///
    /// The state vector **must** be sma, ecc, inc, raan, aop, ta. This function is a shortcut to `cartesian`
    /// and as such it has the same unit requirements.
    pub fn keplerian_vec(state: &Vector6<f64>, dt: Epoch, frame: Frame) -> Self {
        match frame {
            Frame::Geoid { .. } | Frame::Celestial { .. } => Self::keplerian(
                state[0], state[1], state[2], state[3], state[4], state[5], dt, frame,
            ),
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    /// Creates a new Orbit from the geodetic latitude (φ), longitude (λ) and height with respect to the ellipsoid of the frame.
    ///
    /// **Units:** degrees, degrees, km
    /// NOTE: This computation differs from the spherical coordinates because we consider the flattening of body.
    /// Reference: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016
    /// **WARNING:** This uses the rotational rates known to Nyx. For other objects, use `from_altlatlong` for other celestial bodies.
    pub fn from_geodesic(
        latitude: f64,
        longitude: f64,
        height: f64,
        dt: Epoch,
        frame: Frame,
    ) -> Self {
        Self::from_altlatlong(
            latitude,
            longitude,
            height,
            frame.angular_velocity(),
            dt,
            frame,
        )
    }

    /// Creates a new Orbit from the latitude (φ), longitude (λ) and height with respect to the frame's ellipsoid.
    ///
    /// **Units:** degrees, degrees, km, rad/s
    /// NOTE: This computation differs from the spherical coordinates because we consider the flattening of body.
    /// Reference: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016
    pub fn from_altlatlong(
        latitude: f64,
        longitude: f64,
        height: f64,
        angular_velocity: f64,
        dt: Epoch,
        frame: Frame,
    ) -> Self {
        match frame {
            Frame::Geoid {
                flattening,
                semi_major_radius,
                ..
            } => {
                let e2 = 2.0 * flattening - flattening.powi(2);
                let (sin_long, cos_long) = longitude.to_radians().sin_cos();
                let (sin_lat, cos_lat) = latitude.to_radians().sin_cos();
                // page 144
                let c_body = semi_major_radius / ((1.0 - e2 * sin_lat.powi(2)).sqrt());
                let s_body = (semi_major_radius * (1.0 - flattening).powi(2))
                    / ((1.0 - e2 * sin_lat.powi(2)).sqrt());
                let ri = (c_body + height) * cos_lat * cos_long;
                let rj = (c_body + height) * cos_lat * sin_long;
                let rk = (s_body + height) * sin_lat;
                let radius = Vector3::new(ri, rj, rk);
                let velocity = Vector3::new(0.0, 0.0, angular_velocity).cross(&radius);
                Orbit::cartesian(
                    radius[0],
                    radius[1],
                    radius[2],
                    velocity[0],
                    velocity[1],
                    velocity[2],
                    dt,
                    frame,
                )
            }
            _ => panic!("Frame is not Geoid in kind"),
        }
    }

    /// Returns this state as a Keplerian Vector6 in [km, none, degrees, degrees, degrees, degrees]
    ///
    /// Note that the time is **not** returned in the vector.
    pub fn to_keplerian_vec(&self) -> Vector6<f64> {
        Vector6::new(
            self.sma(),
            self.ecc(),
            self.inc(),
            self.raan(),
            self.aop(),
            self.ta(),
        )
    }

    /// Returns the orbital momentum vector
    pub fn hvec(&self) -> Vector3<f64> {
        self.radius().cross(&self.velocity())
    }

    /// Returns the orbital momentum value on the X axis
    pub fn hx(&self) -> f64 {
        self.hvec()[0]
    }

    /// Returns the orbital momentum value on the Y axis
    pub fn hy(&self) -> f64 {
        self.hvec()[1]
    }

    /// Returns the orbital momentum value on the Z axis
    pub fn hz(&self) -> f64 {
        self.hvec()[2]
    }

    /// Returns the norm of the orbital momentum
    pub fn hmag(&self) -> f64 {
        self.hvec().norm()
    }

    /// Returns the specific mechanical energy
    pub fn energy(&self) -> f64 {
        match self.frame {
            Frame::Geoid { gm, .. } | Frame::Celestial { gm, .. } => {
                self.vmag().powi(2) / 2.0 - gm / self.rmag()
            }
            _ => panic!("orbital energy not defined in this frame"),
        }
    }

    /// Returns the semi-major axis in km
    pub fn sma(&self) -> f64 {
        match self.frame {
            Frame::Geoid { gm, .. } | Frame::Celestial { gm, .. } => -gm / (2.0 * self.energy()),
            _ => panic!("sma not defined in this frame"),
        }
    }

    /// Mutates this orbit to change the SMA
    pub fn set_sma(&mut self, new_sma_km: f64) {
        let me = Self::keplerian(
            new_sma_km,
            self.ecc(),
            self.inc(),
            self.raan(),
            self.aop(),
            self.ta(),
            self.dt,
            self.frame,
        );

        self.x = me.x;
        self.y = me.y;
        self.z = me.z;
        self.vx = me.vx;
        self.vy = me.vy;
        self.vz = me.vz;
    }

    /// Returns a copy of the state with a new SMA
    pub fn with_sma(self, new_sma_km: f64) -> Self {
        let mut me = self;
        me.set_sma(new_sma_km);
        me
    }

    /// Returns a copy of the state with a provided SMA added to the current one
    pub fn add_sma(self, delta_sma: f64) -> Self {
        let mut me = self;
        me.set_sma(me.sma() + delta_sma);
        me
    }

    /// Returns the period in seconds
    pub fn period(&self) -> Duration {
        match self.frame {
            Frame::Geoid { gm, .. } | Frame::Celestial { gm, .. } => {
                2.0 * PI * (self.sma().powi(3) / gm).sqrt() * TimeUnit::Second
            }
            _ => panic!("orbital period not defined in this frame"),
        }
    }

    /// Returns the eccentricity vector (no unit)
    pub fn evec(&self) -> Vector3<f64> {
        match self.frame {
            Frame::Geoid { gm, .. } | Frame::Celestial { gm, .. } => {
                let r = self.radius();
                let v = self.velocity();
                ((v.norm().powi(2) - gm / r.norm()) * r - (r.dot(&v)) * v) / gm
            }
            _ => panic!("eccentricity not defined in this frame"),
        }
    }

    /// Returns the eccentricity (no unit)
    pub fn ecc(&self) -> f64 {
        self.evec().norm()
    }

    /// Mutates this orbit to change the ECC
    pub fn set_ecc(&mut self, new_ecc: f64) {
        let me = Self::keplerian(
            self.sma(),
            new_ecc,
            self.inc(),
            self.raan(),
            self.aop(),
            self.ta(),
            self.dt,
            self.frame,
        );

        self.x = me.x;
        self.y = me.y;
        self.z = me.z;
        self.vx = me.vx;
        self.vy = me.vy;
        self.vz = me.vz;
    }

    /// Returns a copy of the state with a new ECC
    pub fn with_ecc(self, new_ecc: f64) -> Self {
        let mut me = self;
        me.set_ecc(new_ecc);
        me
    }

    /// Returns a copy of the state with a provided ECC added to the current one
    pub fn add_ecc(self, delta_ecc: f64) -> Self {
        let mut me = self;
        me.set_ecc(me.ecc() + delta_ecc);
        me
    }

    /// Returns the inclination in degrees
    pub fn inc(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                (self.hvec()[2] / self.hmag()).acos().to_degrees()
            }
            _ => panic!("inclination not defined in this frame"),
        }
    }

    /// Mutates this orbit to change the INC
    pub fn set_inc(&mut self, new_inc: f64) {
        let me = Self::keplerian(
            self.sma(),
            self.ecc(),
            new_inc,
            self.raan(),
            self.aop(),
            self.ta(),
            self.dt,
            self.frame,
        );

        self.x = me.x;
        self.y = me.y;
        self.z = me.z;
        self.vx = me.vx;
        self.vy = me.vy;
        self.vz = me.vz;
    }

    /// Returns a copy of the state with a new INC
    pub fn with_inc(self, new_inc: f64) -> Self {
        let mut me = self;
        me.set_inc(new_inc);
        me
    }

    /// Returns a copy of the state with a provided INC added to the current one
    pub fn add_inc(self, delta_inc: f64) -> Self {
        let mut me = self;
        me.set_inc(me.inc() + delta_inc);
        me
    }

    /// Returns the argument of periapsis in degrees
    pub fn aop(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                let n = Vector3::new(0.0, 0.0, 1.0).cross(&self.hvec());
                let aop = (n.dot(&self.evec()) / (n.norm() * self.ecc())).acos();
                if aop.is_nan() {
                    error!("AoP is NaN");
                    0.0
                } else if self.evec()[2] < 0.0 {
                    (2.0 * PI - aop).to_degrees()
                } else {
                    aop.to_degrees()
                }
            }
            _ => panic!("aop not defined in this frame"),
        }
    }

    /// Mutates this orbit to change the AOP
    pub fn set_aop(&mut self, new_aop: f64) {
        let me = Self::keplerian(
            self.sma(),
            self.ecc(),
            self.inc(),
            self.raan(),
            new_aop,
            self.ta(),
            self.dt,
            self.frame,
        );

        self.x = me.x;
        self.y = me.y;
        self.z = me.z;
        self.vx = me.vx;
        self.vy = me.vy;
        self.vz = me.vz;
    }

    /// Returns a copy of the state with a new AOP
    pub fn with_aop(self, new_aop: f64) -> Self {
        let mut me = self;
        me.set_aop(new_aop);
        me
    }

    /// Returns a copy of the state with a provided AOP added to the current one
    pub fn add_aop(self, delta_aop: f64) -> Self {
        let mut me = self;
        me.set_aop(me.aop() + delta_aop);
        me
    }

    /// Returns the right ascension of ther ascending node in degrees
    pub fn raan(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                let n = Vector3::new(0.0, 0.0, 1.0).cross(&self.hvec());
                let raan = (n[0] / n.norm()).acos();
                if raan.is_nan() {
                    warn!("RAAN is NaN");
                    0.0
                } else if n[1] < 0.0 {
                    (2.0 * PI - raan).to_degrees()
                } else {
                    raan.to_degrees()
                }
            }
            _ => panic!("RAAN not defined in this frame"),
        }
    }

    /// Mutates this orbit to change the RAAN
    pub fn set_raan(&mut self, new_raan: f64) {
        let me = Self::keplerian(
            self.sma(),
            self.ecc(),
            self.inc(),
            new_raan,
            self.aop(),
            self.ta(),
            self.dt,
            self.frame,
        );

        self.x = me.x;
        self.y = me.y;
        self.z = me.z;
        self.vx = me.vx;
        self.vy = me.vy;
        self.vz = me.vz;
    }

    /// Returns a copy of the state with a new RAAN
    pub fn with_raan(self, new_raan: f64) -> Self {
        let mut me = self;
        me.set_raan(new_raan);
        me
    }

    /// Returns a copy of the state with a provided RAAN added to the current one
    pub fn add_raan(self, delta_raan: f64) -> Self {
        let mut me = self;
        me.set_raan(me.raan() + delta_raan);
        me
    }

    /// Returns the true anomaly in degrees between 0 and 360.0
    ///
    /// NOTE: This function will emit a warning stating that the TA should be avoided if in a very near circular orbit
    /// Code from https://github.com/ChristopherRabotin/GMAT/blob/80bde040e12946a61dae90d9fc3538f16df34190/src/gmatutil/util/StateConversionUtil.cpp#L6835
    pub fn ta(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                if self.ecc() < ECC_EPSILON {
                    warn!(
                        "true anomaly ill-defined for circular orbit (e = {})",
                        self.ecc()
                    );
                }
                let cos_nu = self.evec().dot(&self.radius()) / (self.ecc() * self.rmag());
                if (cos_nu.abs() - 1.0).abs() < EPSILON {
                    // This bug drove me crazy when writing SMD in Go in 2017.
                    if cos_nu > 1.0 {
                        180.0
                    } else {
                        0.0
                    }
                } else {
                    let ta = cos_nu.acos();
                    if ta.is_nan() {
                        warn!("TA is NaN");
                        0.0
                    } else if self.radius().dot(&self.velocity()) < 0.0 {
                        (2.0 * PI - ta).to_degrees()
                    } else {
                        ta.to_degrees()
                    }
                }
            }
            _ => panic!("true anomaly not defined in this frame"),
        }
    }

    /// Mutates this orbit to change the TA
    pub fn set_ta(&mut self, new_ta: f64) {
        let me = Self::keplerian(
            self.sma(),
            self.ecc(),
            self.inc(),
            self.raan(),
            self.aop(),
            new_ta,
            self.dt,
            self.frame,
        );

        self.x = me.x;
        self.y = me.y;
        self.z = me.z;
        self.vx = me.vx;
        self.vy = me.vy;
        self.vz = me.vz;
    }

    /// Returns a copy of the state with a new TA
    pub fn with_ta(self, new_ta: f64) -> Self {
        let mut me = self;
        me.set_ta(new_ta);
        me
    }

    /// Returns a copy of the state with a provided TA added to the current one
    pub fn add_ta(self, delta_ta: f64) -> Self {
        let mut me = self;
        me.set_ta(me.ta() + delta_ta);
        me
    }

    /// Returns the true longitude in degrees
    pub fn tlong(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                // Angles already in degrees
                between_0_360(self.aop() + self.raan() + self.ta())
            }
            _ => panic!("true longitude not defined in this frame"),
        }
    }

    /// Returns the argument of latitude in degrees
    ///
    /// NOTE: If the orbit is near circular, the AoL will be computed from the true longitude
    /// instead of relying on the ill-defined true anomaly.
    pub fn aol(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                between_0_360(if self.ecc() < ECC_EPSILON {
                    self.tlong() - self.raan()
                } else {
                    self.aop() + self.ta()
                })
            }
            _ => panic!("argument of latitude not defined in this frame"),
        }
    }

    /// Returns the radius of periapsis (or perigee around Earth), in kilometers.
    pub fn periapsis(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => self.sma() * (1.0 - self.ecc()),
            _ => panic!("periapsis not defined in this frame"),
        }
    }

    /// Returns the radius of apoapsis (or apogee around Earth), in kilometers.
    pub fn apoapsis(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => self.sma() * (1.0 + self.ecc()),
            _ => panic!("apoapsis not defined in this frame"),
        }
    }

    /// Returns the eccentric anomaly in degrees
    ///
    /// This is a conversion from GMAT's StateConversionUtil::TrueToEccentricAnomaly
    pub fn ea(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                let (sin_ta, cos_ta) = self.ta().to_radians().sin_cos();
                let ecc_cos_ta = self.ecc() * cos_ta;
                let sin_ea = ((1.0 - self.ecc().powi(2)).sqrt() * sin_ta) / (1.0 + ecc_cos_ta);
                let cos_ea = (self.ecc() + cos_ta) / (1.0 + ecc_cos_ta);
                // The atan2 function is a bit confusing: https://doc.rust-lang.org/std/primitive.f64.html#method.atan2 .
                sin_ea.atan2(cos_ea).to_degrees()
            }
            _ => panic!("eccentric anomaly is not defined in this frame"),
        }
    }

    /// Returns the flight path angle in degrees
    pub fn fpa(&self) -> f64 {
        let nu = self.ta().to_radians();
        let ecc = self.ecc();
        let denom = (1.0 + 2.0 * ecc * nu.cos() + ecc.powi(2)).sqrt();
        let sin_fpa = ecc * nu.sin() / denom;
        let cos_fpa = 1.0 + ecc * nu.cos() / denom;
        sin_fpa.atan2(cos_fpa).to_degrees()
    }

    /// Returns the mean anomaly in degrees
    ///
    /// This is a conversion from GMAT's StateConversionUtil::TrueToMeanAnomaly
    pub fn ma(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                if self.ecc() < 1.0 {
                    between_0_360(
                        (self.ea().to_radians() - self.ecc() * self.ea().to_radians().sin())
                            .to_degrees(),
                    )
                } else if self.ecc() > 1.0 {
                    info!("computing the hyperbolic anomaly");
                    // From GMAT's TrueToHyperbolicAnomaly
                    ((self.ta().to_radians().sin() * (self.ecc().powi(2) - 1.0)).sqrt()
                        / (1.0 + self.ecc() * self.ta().to_radians().cos()))
                    .asinh()
                    .to_degrees()
                } else {
                    error!("parabolic orbit: setting mean anomaly to 0.0");
                    0.0
                }
            }
            _ => panic!("mean anomaly is not defined in this frame"),
        }
    }

    /// Returns the semi parameter (or semilatus rectum)
    pub fn semi_parameter(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                self.sma() * (1.0 - self.ecc().powi(2))
            }
            _ => panic!("semi parameter is not defined in this frame"),
        }
    }

    /// Returns whether this state satisfies the requirement to compute the Mean Brouwer Short orbital
    /// element set.
    ///
    /// This is a conversion from GMAT's StateConversionUtil::CartesianToBrouwerMeanShort.
    /// The details are at the log level `info`.
    /// NOTE: Mean Brouwer Short are only defined around Earth. However, `nyx` does *not* check the
    /// main celestial body around which the state is defined (GMAT does perform this verification).
    pub fn is_brouwer_short_valid(&self) -> bool {
        if self.inc() > 180.0 {
            info!("Brouwer Mean Short only applicable for inclinations less than 180.0");
            false
        } else if self.ecc() >= 1.0 || self.ecc() < 0.0 {
            info!("Brouwer Mean Short only applicable for elliptical orbits");
            false
        } else if self.periapsis() < 3000.0 {
            // NOTE: GMAT emits a warning if the periagee is less than the Earth radius, but we do not do that here.
            info!("Brouwer Mean Short only applicable for if perigee is greater than 3000 km");
            false
        } else {
            true
        }
    }

    /// Returns the geodetic longitude (λ) in degrees. Value is between 0 and 360 degrees.
    ///
    /// Although the reference is not Vallado, the math from Vallado proves to be equivalent.
    /// Reference: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016
    pub fn geodetic_longitude(&self) -> f64 {
        match self.frame {
            Frame::Geoid { .. } => between_0_360(self.y.atan2(self.x).to_degrees()),
            _ => panic!("geodetic elements only defined in a Geoid frame"),
        }
    }

    /// Returns the geodetic latitude (φ) in degrees. Value is between -180 and +180 degrees.
    ///
    /// Reference: Vallado, 4th Ed., Algorithm 12 page 172.
    pub fn geodetic_latitude(&self) -> f64 {
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
                let e2 = flattening * (2.0 - flattening);
                loop {
                    attempt_no += 1;
                    let c_earth =
                        semi_major_radius / ((1.0 - e2 * (latitude).sin().powi(2)).sqrt());
                    let new_latitude = (self.z + c_earth * e2 * (latitude).sin()).atan2(r_delta);
                    if (latitude - new_latitude).abs() < eps {
                        return between_pm_180(new_latitude.to_degrees());
                    } else if attempt_no >= max_attempts {
                        error!(
                            "geodetic latitude failed to converge -- error = {}",
                            (latitude - new_latitude).abs()
                        );
                        return between_pm_180(new_latitude.to_degrees());
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
    pub fn geodetic_height(&self) -> f64 {
        match self.frame {
            Frame::Geoid {
                flattening,
                semi_major_radius,
                ..
            } => {
                let e2 = flattening * (2.0 - flattening);
                let latitude = self.geodetic_latitude().to_radians();
                let sin_lat = latitude.sin();
                if (latitude - 1.0).abs() < 0.1 {
                    // We are near poles, let's use another formulation.
                    let s_earth = (semi_major_radius * (1.0 - flattening).powi(2))
                        / ((1.0 - e2 * sin_lat.powi(2)).sqrt());
                    self.z / latitude.sin() - s_earth
                } else {
                    let c_earth = semi_major_radius / ((1.0 - e2 * sin_lat.powi(2)).sqrt());
                    let r_delta = (self.x.powi(2) + self.y.powi(2)).sqrt();
                    r_delta / latitude.cos() - c_earth
                }
            }
            _ => panic!("geodetic elements only defined in a Geoid frame"),
        }
    }

    /// Returns the right ascension of this orbit in degrees
    pub fn right_ascension(&self) -> f64 {
        between_0_360((self.y.atan2(self.x)).to_degrees())
    }

    /// Returns the declination of this orbit in degrees
    pub fn declination(&self) -> f64 {
        between_pm_180((self.z / self.rmag()).asin().to_degrees())
    }

    /// Returns the semi minor axis in km, includes code for a hyperbolic orbit
    pub fn semi_minor_axis(&self) -> f64 {
        if self.ecc() <= 1.0 {
            ((self.sma() * self.ecc()).powi(2) - self.sma().powi(2)).sqrt()
        } else {
            self.hmag().powi(2) / (self.frame.gm() * (self.ecc().powi(2) - 1.0).sqrt())
        }
    }

    /// Returns the velocity declination of this orbit in degrees
    pub fn velocity_declination(&self) -> f64 {
        between_pm_180((self.vz / self.vmag()).asin().to_degrees())
    }

    pub fn b_plane(&self) -> Result<BPlane, NyxError> {
        BPlane::new(*self)
    }

    /// Returns the $C_3$ of this orbit
    pub fn c3(&self) -> f64 {
        -self.frame.gm() / self.sma()
    }

    /// Returns the radius of periapse in kilometers for the provided turn angle of this hyperbolic orbit.
    pub fn vinf_periapsis(&self, turn_angle_degrees: f64) -> Result<f64, NyxError> {
        if self.ecc() <= 1.0 {
            Err(NyxError::NotHyperbolic(
                "Orbit is not hyperbolic. Convert to target object first".to_string(),
            ))
        } else {
            let cos_rho = (0.5 * (PI - turn_angle_degrees.to_radians())).cos();

            Ok((1.0 / cos_rho - 1.0) * self.frame.gm() / self.vmag().powi(2))
        }
    }

    /// Returns the turn angle in degrees for the provided radius of periapse passage of this hyperbolic orbit
    pub fn vinf_turn_angle(&self, periapsis_km: f64) -> Result<f64, NyxError> {
        if self.ecc() <= 1.0 {
            Err(NyxError::NotHyperbolic(
                "Orbit is not hyperbolic. Convert to target object first".to_string(),
            ))
        } else {
            let rho = (1.0 / (1.0 + self.vmag().powi(2) * (periapsis_km / self.frame.gm()))).acos();
            Ok(between_0_360((PI - 2.0 * rho).to_degrees()))
        }
    }

    /// Returns the hyperbolic anomaly in degrees between 0 and 360.0
    pub fn hyperbolic_anomaly(&self) -> Result<f64, NyxError> {
        if self.ecc() <= 1.0 {
            Err(NyxError::NotHyperbolic(
                "Orbit is not hyperbolic so there is no hyperbolic anomaly.".to_string(),
            ))
        } else {
            let (sin_ta, cos_ta) = self.ta().to_radians().sin_cos();
            let sinh_h = (sin_ta * (self.ecc().powi(2) - 1.0).sqrt()) / (1.0 + self.ecc() * cos_ta);
            Ok(between_0_360(sinh_h.asinh().to_degrees()))
        }
    }

    /// Returns the direct cosine rotation matrix to convert to this inertial state.
    pub fn dcm_to_inertial(&self, from: Frame) -> Matrix3<f64> {
        match from {
            Frame::RIC => {
                r3(-self.raan().to_radians())
                    * r1(-self.inc().to_radians())
                    * r3(-self.aol().to_radians())
            }
            Frame::VNC => {
                let v = self.velocity() / self.vmag();
                let n = self.hvec() / self.hmag();
                let c = v.cross(&n);
                Matrix3::new(v[0], v[1], v[2], n[0], n[1], n[2], c[0], c[1], c[2]).transpose()
            }
            Frame::RCN => {
                let r = self.radius() / self.rmag();
                let n = self.hvec() / self.hmag();
                let c = n.cross(&r);
                Matrix3::new(r[0], r[1], r[2], c[0], c[1], c[2], n[0], n[1], n[2]).transpose()
            }
            _ => panic!("did not provide a local frame"),
        }
    }

    /// Apply the provided delta-v (in km/s)
    pub fn apply_dv(&mut self, dv: Vector3<f64>) {
        self.vx += dv[0];
        self.vy += dv[1];
        self.vz += dv[2];
    }

    /// Copies this orbit after applying the provided delta-v (in km/s)
    pub fn with_dv(self, dv: Vector3<f64>) -> Self {
        let mut me = self;
        me.apply_dv(dv);
        me
    }

    /// Rotate this state provided a direct cosine matrix
    /// **Bug:** This does not account for the transport theorem!
    pub fn apply_dcm(&mut self, dcm: Matrix3<f64>) {
        let new_r = dcm * self.radius();
        let new_v = dcm * self.velocity();

        self.x = new_r[0];
        self.y = new_r[1];
        self.z = new_r[2];

        self.vx = new_v[0];
        self.vy = new_v[1];
        self.vz = new_v[2];
    }

    /// Sets the STM of this state of identity, which also enables computation of the STM for spacecraft navigation
    pub fn enable_stm(&mut self) {
        self.stm = Some(Matrix6::identity());
        self.stm_kind = StmKind::Step;
    }

    /// Sets the STM of this state of identity, which also enables computation of the STM for trajectory optimization
    pub fn enable_traj_stm(&mut self) {
        self.stm = Some(Matrix6::identity());
        self.stm_kind = StmKind::Traj;
    }

    /// Copies the current state but sets the STM to identity
    pub fn with_stm(self) -> Self {
        let mut me = self;
        me.enable_stm();
        me
    }

    /// Sets the STM of this state of identity
    pub fn stm_identity(&mut self) {
        self.stm = Some(Matrix6::identity());
    }

    /// Unwraps this STM, or panics if unset.
    pub fn stm(&self) -> Matrix6<f64> {
        self.stm.unwrap()
    }
}

impl TimeTagged for Orbit {
    fn epoch(&self) -> Epoch {
        self.dt
    }

    fn set_epoch(&mut self, epoch: Epoch) {
        self.dt = epoch
    }
}

impl PartialEq for Orbit {
    /// Two states are equal if their position are equal within one centimeter and their velocities within one centimeter per second.
    fn eq(&self, other: &Orbit) -> bool {
        let distance_tol = 1e-5; // centimeter
        let velocity_tol = 1e-5; // centimeter per second
        self.dt == other.dt
            && (self.x - other.x).abs() < distance_tol
            && (self.y - other.y).abs() < distance_tol
            && (self.z - other.z).abs() < distance_tol
            && (self.vx - other.vx).abs() < velocity_tol
            && (self.vy - other.vy).abs() < velocity_tol
            && (self.vz - other.vz).abs() < velocity_tol
            && self.frame == other.frame
    }
}

impl Add for Orbit {
    type Output = Orbit;

    /// Add one state from another. Frame must be manually changed if needed. STM will be copied from &self.
    fn add(self, other: Orbit) -> Orbit {
        Orbit {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
            vx: self.vx + other.vx,
            vy: self.vy + other.vy,
            vz: self.vz + other.vz,
            dt: self.dt,
            frame: self.frame,
            stm: self.stm,
            stm_kind: self.stm_kind,
        }
    }
}

impl Sub for Orbit {
    type Output = Orbit;

    /// Subtract one state from another. STM will be copied from &self.
    fn sub(self, other: Orbit) -> Orbit {
        Orbit {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
            vx: self.vx - other.vx,
            vy: self.vy - other.vy,
            vz: self.vz - other.vz,
            dt: self.dt,
            frame: self.frame,
            stm: self.stm,
            stm_kind: self.stm_kind,
        }
    }
}

impl Neg for Orbit {
    type Output = Orbit;

    /// Subtract one state from another. STM will be copied from &self.
    fn neg(self) -> Self::Output {
        Orbit {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            vx: -self.vx,
            vy: -self.vy,
            vz: -self.vz,
            dt: self.dt,
            frame: self.frame,
            stm: self.stm,
            stm_kind: self.stm_kind,
        }
    }
}

impl Add for &Orbit {
    type Output = Orbit;

    /// Add one state from another. Frame must be manually changed if needed. STM will be copied from &self.
    fn add(self, other: &Orbit) -> Orbit {
        Orbit {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
            vx: self.vx + other.vx,
            vy: self.vy + other.vy,
            vz: self.vz + other.vz,
            dt: self.dt,
            frame: self.frame,
            stm: self.stm,
            stm_kind: self.stm_kind,
        }
    }
}

impl Sub for &Orbit {
    type Output = Orbit;

    /// Subtract one state from another. STM will be copied from &self.
    fn sub(self, other: &Orbit) -> Orbit {
        Orbit {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
            vx: self.vx - other.vx,
            vy: self.vy - other.vy,
            vz: self.vz - other.vz,
            dt: self.dt,
            frame: self.frame,
            stm: self.stm,
            stm_kind: self.stm_kind,
        }
    }
}

impl Neg for &Orbit {
    type Output = Orbit;

    /// Subtract one state from another. STM will be copied form &self.
    fn neg(self) -> Self::Output {
        Orbit {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            vx: -self.vx,
            vy: -self.vy,
            vz: -self.vz,
            dt: self.dt,
            frame: self.frame,
            stm: self.stm,
            stm_kind: self.stm_kind,
        }
    }
}

impl Serialize for Orbit {
    /// NOTE: This is not part of unit testing because there is no deseralization of Orbit (yet)
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("Orbit", 7)?;
        state.serialize_field("dt", &self.dt.as_jde_et_days())?;
        state.serialize_field("x", &self.x)?;
        state.serialize_field("y", &self.y)?;
        state.serialize_field("z", &self.z)?;
        state.serialize_field("vx", &self.vx)?;
        state.serialize_field("vy", &self.vy)?;
        state.serialize_field("vz", &self.vz)?;
        state.end()
    }
}

impl fmt::Display for Orbit {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "[{}] {}\tposition = [{:.6}, {:.6}, {:.6}] km\tvelocity = [{:.6}, {:.6}, {:.6}] km/s",
            self.frame, self.dt, self.x, self.y, self.z, self.vx, self.vy, self.vz
        )
    }
}

impl fmt::LowerExp for Orbit {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "[{}] {}\tposition = [{:e}, {:e}, {:e}] km\tvelocity = [{:e}, {:e}, {:e}] km/s",
            self.frame, self.dt, self.x, self.y, self.z, self.vx, self.vy, self.vz
        )
    }
}

impl fmt::Octal for Orbit {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "[{}] {}\tsma = {:.6} km\tecc = {:.6}\tinc = {:.6} deg\traan = {:.6} deg\taop = {:.6} deg\tta = {:.6} deg",
            self.frame,
            self.dt,
            self.sma(),
            self.ecc(),
            self.inc(),
            self.raan(),
            self.aop(),
            self.ta()
        )
    }
}
