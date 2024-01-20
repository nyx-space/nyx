/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::AstroError;
use super::Cosm;
use super::State;
use super::{BPlane, Frame};
use crate::dynamics::DynamicsError;
use crate::io::orbit::OrbitSerde;
use crate::io::{
    epoch_from_str, epoch_to_str, frame_from_str, frame_to_str, ConfigRepr, Configurable,
};
use crate::linalg::{Const, Matrix3, Matrix6, OVector, Vector3, Vector6};
use crate::mc::MultivariateNormal;
use crate::md::prelude::Objective;
use crate::md::StateParameter;

use crate::time::{Duration, Epoch, Unit};
use crate::utils::{
    between_0_360, between_pm_180, cartesian_to_spherical, perpv, r1, r3, rss_orbit_errors,
    spherical_to_cartesian,
};
use crate::NyxError;
use approx::{abs_diff_eq, relative_eq};
use serde::{Deserialize, Serialize};
use std::f64::consts::PI;
use std::f64::consts::TAU;
use std::f64::EPSILON;
use std::fmt;
use std::ops::{Add, AddAssign, Neg, Sub, SubAssign};
use std::sync::Arc;

#[cfg(feature = "python")]
use crate::io::ConfigError;
#[cfg(feature = "python")]
use crate::python::cosmic::Frame as PyFrame;
#[cfg(feature = "python")]
use pyo3::class::basic::CompareOp;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyType;
#[cfg(feature = "python")]
use std::collections::BTreeMap;

/// If an orbit has an eccentricity below the following value, it is considered circular (only affects warning messages)
pub const ECC_EPSILON: f64 = 1e-11;
pub const MA_EPSILON: f64 = 1e-16;

pub fn assert_orbit_eq_or_abs(left: &Orbit, right: &Orbit, epsilon: f64, msg: &str) {
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

pub fn assert_orbit_eq_or_rel(left: &Orbit, right: &Orbit, epsilon: f64, msg: &str) {
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

/// Orbit defines an orbital state
///
/// Unless noted otherwise, algorithms are from GMAT 2016a [StateConversionUtil.cpp](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/StateConversionUtil.cpp).
/// Regardless of the constructor used, this struct stores all the state information in Cartesian coordinates
/// as these are always non singular.
/// _Note:_ although not yet supported, this struct may change once True of Date or other nutation frames
/// are added to the toolkit.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "python", pyclass)]
pub struct Orbit {
    /// in km
    pub x_km: f64,
    /// in km
    pub y_km: f64,
    /// in km
    pub z_km: f64,
    /// in km/s
    pub vx_km_s: f64,
    /// in km/s
    pub vy_km_s: f64,
    /// in km/s
    pub vz_km_s: f64,
    #[serde(serialize_with = "epoch_to_str", deserialize_with = "epoch_from_str")]
    pub epoch: Epoch,
    /// Frame contains everything we need to compute state information
    #[serde(serialize_with = "frame_to_str", deserialize_with = "frame_from_str")]
    pub frame: Frame,
    /// Optionally stores the state transition matrix from the start of the propagation until the current time (i.e. trajectory STM, not step-size STM)
    #[serde(skip)]
    pub stm: Option<Matrix6<f64>>,
}

impl Orbit {
    /// Creates a new Orbit in the provided frame at the provided Epoch.
    ///
    /// **Units:** km, km, km, km/s, km/s, km/s
    pub fn cartesian(
        x_km: f64,
        y_km: f64,
        z_km: f64,
        vx_km_s: f64,
        vy_km_s: f64,
        vz_km_s: f64,
        epoch: Epoch,
        frame: Frame,
    ) -> Self {
        Self {
            x_km,
            y_km,
            z_km,
            vx_km_s,
            vy_km_s,
            vz_km_s,
            epoch,
            frame,
            stm: None,
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
        epoch: Epoch,
        frame: Frame,
    ) -> Self {
        let mut me = Self::cartesian(x, y, z, vx, vy, vz, epoch, frame);
        me.enable_stm();
        me
    }

    /// Creates a new Orbit in the provided frame at the provided Epoch in time with 0.0 velocity.
    ///
    /// **Units:** km, km, km
    pub fn from_position(x: f64, y: f64, z: f64, epoch: Epoch, frame: Frame) -> Self {
        Orbit {
            x_km: x,
            y_km: y,
            z_km: z,
            vx_km_s: 0.0,
            vy_km_s: 0.0,
            vz_km_s: 0.0,
            epoch,
            frame,
            stm: None,
        }
    }

    /// Creates a new Orbit around in the provided frame from the borrowed state vector
    ///
    /// The state vector **must** be x, y, z, vx, vy, vz. This function is a shortcut to `cartesian`
    /// and as such it has the same unit requirements.
    pub fn cartesian_vec(state: &Vector6<f64>, epoch: Epoch, frame: Frame) -> Self {
        Orbit {
            x_km: state[0],
            y_km: state[1],
            z_km: state[2],
            vx_km_s: state[3],
            vy_km_s: state[4],
            vz_km_s: state[5],
            epoch,
            frame,
            stm: None,
        }
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
        sma_km: f64,
        ecc: f64,
        inc_deg: f64,
        raan_deg: f64,
        aop_deg: f64,
        ta_deg: f64,
        epoch: Epoch,
        frame: Frame,
    ) -> Self {
        match frame {
            Frame::Geoid { gm, .. } | Frame::Celestial { gm, .. } => {
                if gm.abs() < EPSILON {
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
                let sma = if ecc > 1.0 && sma_km > 0.0 {
                    warn!("eccentricity > 1 (hyperbolic) BUT SMA > 0 (elliptical): sign of SMA changed");
                    sma_km * -1.0
                } else if ecc < 1.0 && sma_km < 0.0 {
                    warn!("eccentricity < 1 (elliptical) BUT SMA < 0 (hyperbolic): sign of SMA changed");
                    sma_km * -1.0
                } else {
                    sma_km
                };
                if (sma * (1.0 - ecc)).abs() < 1e-3 {
                    // GMAT errors below one meter. Let's warn for below that, but not panic, might be useful for landing scenarios?
                    warn!("radius of periapsis is less than one meter");
                }
                if (1.0 - ecc).abs() < EPSILON {
                    panic!("parabolic orbits have ill-defined Keplerian orbital elements");
                }
                if ecc > 1.0 {
                    let ta = between_0_360(ta_deg);
                    if ta > (PI - (1.0 / ecc).acos()).to_degrees() {
                        panic!(
                            "true anomaly value ({ta}) physically impossible for a hyperbolic orbit",
                        );
                    }
                }
                if (1.0 + ecc * ta_deg.to_radians().cos()).is_infinite() {
                    panic!("radius of orbit is infinite");
                }
                // Done with all the warnings and errors supported by GMAT
                // The conversion algorithm itself comes from GMAT's StateConversionUtil::ComputeKeplToCart
                // NOTE: GMAT supports mean anomaly instead of true anomaly, but only for backward compatibility reasons
                // so it isn't supported here.
                let inc = inc_deg.to_radians();
                let raan = raan_deg.to_radians();
                let aop = aop_deg.to_radians();
                let ta = ta_deg.to_radians();
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
                    x_km: x,
                    y_km: y,
                    z_km: z,
                    vx_km_s: vx,
                    vy_km_s: vy,
                    vz_km_s: vz,
                    epoch,
                    frame,
                    stm: None,
                }
            }
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    /// Creates a new Orbit from the provided semi-major axis altitude in kilometers
    pub fn keplerian_altitude(
        sma_altitude: f64,
        ecc: f64,
        inc: f64,
        raan: f64,
        aop: f64,
        ta: f64,
        epoch: Epoch,
        frame: Frame,
    ) -> Self {
        Self::keplerian(
            sma_altitude + frame.equatorial_radius(),
            ecc,
            inc,
            raan,
            aop,
            ta,
            epoch,
            frame,
        )
    }

    /// Creates a new Orbit from the provided radii of apoapsis and periapsis, in kilometers
    pub fn keplerian_apsis_radii(
        r_apo_km: f64,
        r_peri_km: f64,
        inc_deg: f64,
        raan_deg: f64,
        aop_deg: f64,
        ta_deg: f64,
        epoch: Epoch,
        frame: Frame,
    ) -> Self {
        let sma_km = (r_apo_km + r_peri_km) / 2.0;
        let ecc = r_apo_km / sma_km - 1.0;
        Self::keplerian(
            sma_km, ecc, inc_deg, raan_deg, aop_deg, ta_deg, epoch, frame,
        )
    }

    /// Creates a new Orbit from the provided altitudes of apoapsis and periapsis, in kilometers
    pub fn keplerian_apsis_altitude(
        apo_alt_km: f64,
        peri_alt_km: f64,
        inc_deg: f64,
        raan_deg: f64,
        aop_deg: f64,
        ta_deg: f64,
        epoch: Epoch,
        frame: Frame,
    ) -> Self {
        Self::keplerian_apsis_radii(
            apo_alt_km + frame.equatorial_radius(),
            peri_alt_km + frame.equatorial_radius(),
            inc_deg,
            raan_deg,
            aop_deg,
            ta_deg,
            epoch,
            frame,
        )
    }

    /// Initializes a new orbit from the Keplerian orbital elements using the mean anomaly instead of the true anomaly.
    ///
    /// # Implementation notes
    /// This function starts by converting the mean anomaly to true anomaly, and then it initializes the orbit
    /// using the keplerian(..) method.
    /// The conversion is from GMAT's MeanToTrueAnomaly function, transliterated originally by Claude and GPT4 with human adjustments.
    pub fn keplerian_mean_anomaly(
        sma_km: f64,
        ecc: f64,
        inc_deg: f64,
        raan_deg: f64,
        aop_deg: f64,
        ma_deg: f64,
        epoch: Epoch,
        frame: Frame,
    ) -> Result<Self, NyxError> {
        // Start by computing the true anomaly
        let ta_rad = compute_mean_to_true_anomaly(ma_deg.to_radians(), ecc, MA_EPSILON)?;

        Ok(Self::keplerian(
            sma_km,
            ecc,
            inc_deg,
            raan_deg,
            aop_deg,
            ta_rad.to_degrees(),
            epoch,
            frame,
        ))
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
        epoch: Epoch,
        frame: Frame,
    ) -> Self {
        Self::from_altlatlong(
            latitude,
            longitude,
            height,
            frame.angular_velocity(),
            epoch,
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
        epoch: Epoch,
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
                    epoch,
                    frame,
                )
            }
            _ => panic!("Frame is not Geoid in kind"),
        }
    }

    /// Returns the radius vector of this Orbit in [km, km, km]
    pub fn radius(&self) -> Vector3<f64> {
        Vector3::new(self.x_km, self.y_km, self.z_km)
    }

    /// Returns the velocity vector of this Orbit in [km/s, km/s, km/s]
    pub fn velocity(&self) -> Vector3<f64> {
        Vector3::new(self.vx_km_s, self.vy_km_s, self.vz_km_s)
    }

    /// Returns the unit vector in the direction of the state radius
    pub fn r_hat(&self) -> Vector3<f64> {
        self.radius() / self.rmag_km()
    }

    /// Returns the unit vector in the direction of the state velocity
    pub fn v_hat(&self) -> Vector3<f64> {
        perpv(&self.velocity(), &self.r_hat()) / self.rmag_km()
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

    /// Returns this state as a Cartesian Vector6 in [km, km, km, km/s, km/s, km/s]
    ///
    /// Note that the time is **not** returned in the vector.
    pub fn to_cartesian_vec(self) -> Vector6<f64> {
        Vector6::new(
            self.x_km,
            self.y_km,
            self.z_km,
            self.vx_km_s,
            self.vy_km_s,
            self.vz_km_s,
        )
    }

    /// Returns this state as a Keplerian Vector6 in [km, none, degrees, degrees, degrees, degrees]
    ///
    /// Note that the time is **not** returned in the vector.
    pub fn to_keplerian_vec(self) -> Vector6<f64> {
        Vector6::new(
            self.sma_km(),
            self.ecc(),
            self.inc_deg(),
            self.raan_deg(),
            self.aop_deg(),
            self.ta_deg(),
        )
    }

    /// Creates a new Orbit around the provided frame from the borrowed state vector
    ///
    /// The state vector **must** be sma, ecc, inc, raan, aop, ta. This function is a shortcut to `cartesian`
    /// and as such it has the same unit requirements.
    pub fn keplerian_vec(state: &Vector6<f64>, epoch: Epoch, frame: Frame) -> Self {
        match frame {
            Frame::Geoid { .. } | Frame::Celestial { .. } => Self::keplerian(
                state[0], state[1], state[2], state[3], state[4], state[5], epoch, frame,
            ),
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    /// Returns the orbital momentum vector
    pub fn hvec(&self) -> Vector3<f64> {
        self.radius().cross(&self.velocity())
    }

    /// Returns a copy of the state with a new radius
    pub fn with_radius(self, new_radius: &Vector3<f64>) -> Self {
        let mut me = self;
        me.x_km = new_radius[0];
        me.y_km = new_radius[1];
        me.z_km = new_radius[2];
        me
    }

    /// Returns a copy of the state with a new radius
    pub fn with_velocity(self, new_velocity: &Vector3<f64>) -> Self {
        let mut me = self;
        me.vx_km_s = new_velocity[0];
        me.vy_km_s = new_velocity[1];
        me.vz_km_s = new_velocity[2];
        me
    }

    /// Returns a copy of the state with a new SMA
    pub fn with_sma(self, new_sma_km: f64) -> Self {
        let mut me = self;
        me.set_sma_km(new_sma_km);
        me
    }

    /// Returns a copy of the state with a provided SMA added to the current one
    pub fn add_sma(self, delta_sma: f64) -> Self {
        let mut me = self;
        me.set_sma_km(me.sma_km() + delta_sma);
        me
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

    /// Returns a copy of the state with a new INC
    pub fn with_inc(self, new_inc: f64) -> Self {
        let mut me = self;
        me.set_inc_deg(new_inc);
        me
    }

    /// Returns a copy of the state with a provided INC added to the current one
    pub fn add_inc(self, delta_inc: f64) -> Self {
        let mut me = self;
        me.set_inc_deg(me.inc_deg() + delta_inc);
        me
    }

    /// Returns a copy of the state with a new AOP
    pub fn with_aop(self, new_aop: f64) -> Self {
        let mut me = self;
        me.set_aop_deg(new_aop);
        me
    }

    /// Returns a copy of the state with a provided AOP added to the current one
    pub fn add_aop(self, delta_aop: f64) -> Self {
        let mut me = self;
        me.set_aop_deg(me.aop_deg() + delta_aop);
        me
    }

    /// Returns a copy of the state with a new RAAN
    pub fn with_raan(self, new_raan: f64) -> Self {
        let mut me = self;
        me.set_raan_deg(new_raan);
        me
    }

    /// Returns a copy of the state with a provided RAAN added to the current one
    pub fn add_raan(self, delta_raan: f64) -> Self {
        let mut me = self;
        me.set_raan_deg(me.raan_deg() + delta_raan);
        me
    }

    /// Returns a copy of the state with a new TA
    pub fn with_ta(self, new_ta: f64) -> Self {
        let mut me = self;
        me.set_ta_deg(new_ta);
        me
    }

    /// Returns a copy of the state with a provided TA added to the current one
    pub fn add_ta(self, delta_ta: f64) -> Self {
        let mut me = self;
        me.set_ta_deg(me.ta_deg() + delta_ta);
        me
    }

    /// Returns a copy of this state with the provided apoasis and periapse
    pub fn with_apoapsis_periapsis(self, new_ra: f64, new_rp: f64) -> Self {
        Self::keplerian_apsis_radii(
            new_ra,
            new_rp,
            self.inc_deg(),
            self.raan_deg(),
            self.aop_deg(),
            self.ta_deg(),
            self.epoch,
            self.frame,
        )
    }

    /// Returns a copy of this state with the provided apoasis and periapse added to the current values
    pub fn add_apoapsis_periapsis(self, delta_ra: f64, delta_rp: f64) -> Self {
        Self::keplerian_apsis_radii(
            self.apoapsis_km() + delta_ra,
            self.periapsis_km() + delta_rp,
            self.inc_deg(),
            self.raan_deg(),
            self.aop_deg(),
            self.ta_deg(),
            self.epoch,
            self.frame,
        )
    }

    /// Returns the direct cosine rotation matrix to convert to this state's frame (inertial or otherwise).
    /// ## Example
    /// let dcm_vnc2inertial = orbit.dcm_from_traj_frame(Frame::VNC)?;
    /// let vector_inertial = dcm_vnc2inertial * vector_vnc;
    pub fn dcm_from_traj_frame(&self, from: Frame) -> Result<Matrix3<f64>, AstroError> {
        match from {
            Frame::RIC => Ok(r3(-self.raan_deg().to_radians())
                * r1(-self.inc_deg().to_radians())
                * r3(-self.aol_deg().to_radians())),
            Frame::VNC => {
                let v = self.velocity() / self.vmag_km_s();
                let n = self.hvec() / self.hmag_km2_s();
                let c = v.cross(&n);
                Ok(Matrix3::new(v[0], v[1], v[2], n[0], n[1], n[2], c[0], c[1], c[2]).transpose())
            }
            Frame::RCN => {
                let r = self.radius() / self.rmag_km();
                let n = self.hvec() / self.hmag_km2_s();
                let c = n.cross(&r);
                Ok(Matrix3::new(r[0], r[1], r[2], c[0], c[1], c[2], n[0], n[1], n[2]).transpose())
            }
            Frame::SEZ => {
                // From the GMAT MathSpec, page 30 section 2.6.9 and from `Calculate_RFT` in `TopocentricAxes.cpp`, this returns the
                // rotation matrix from the topocentric frame (SEZ) to body fixed frame.
                // In the GMAT MathSpec notation, R_{IF} is the DCM from body fixed to inertial. Similarly, R{FT} is from topocentric
                // to body fixed.
                if !self.frame.is_body_fixed() {
                    warn!("Computation of SEZ rotation matrix must be done in a body fixed frame and {} is not one!", self.frame);
                }
                if (self.x_km.powi(2) + self.y_km.powi(2)).sqrt() < 1e-3 {
                    warn!("SEZ frame ill-defined when close to the poles");
                }
                let phi = self.geodetic_latitude_deg().to_radians();
                let lambda = self.geodetic_longitude_deg().to_radians();
                let z_hat = Vector3::new(
                    phi.cos() * lambda.cos(),
                    phi.cos() * lambda.sin(),
                    phi.sin(),
                );
                // y_hat MUST be renormalized otherwise it's about 0.76 and therefore the rotation looses the norms conservation property.
                let mut y_hat = Vector3::new(0.0, 0.0, 1.0).cross(&z_hat);
                y_hat /= y_hat.norm();
                let x_hat = y_hat.cross(&z_hat);
                Ok(Matrix3::new(
                    x_hat[0], y_hat[0], z_hat[0], x_hat[1], y_hat[1], z_hat[1], x_hat[2], y_hat[2],
                    z_hat[2],
                ))
            }
            _ => Err(AstroError::NotLocalFrame),
        }
    }

    /// Returns a 6x6 DCM to convert to this inertial state.
    /// WARNING: This DCM does NOT contain the corrections needed for the transport theorem, and therefore the velocity rotation is wrong.
    pub fn dcm6x6_from_traj_frame(&self, from: Frame) -> Result<Matrix6<f64>, AstroError> {
        let dcm3x3 = self.dcm_from_traj_frame(from)?;

        let mut dcm = Matrix6::zeros();
        for i in 0..3 {
            for j in 0..3 {
                dcm[(i, j)] = dcm3x3[(i, j)];
                dcm[(i + 3, j + 3)] = dcm3x3[(i, j)];
            }
        }

        Ok(dcm)
    }

    /// Apply the provided delta-v (in km/s)
    pub fn apply_dv(&mut self, dv: Vector3<f64>) {
        self.vx_km_s += dv[0];
        self.vy_km_s += dv[1];
        self.vz_km_s += dv[2];
    }

    /// Copies this orbit after applying the provided delta-v (in km/s)
    pub fn with_dv(self, dv: Vector3<f64>) -> Self {
        let mut me = self;
        me.apply_dv(dv);
        me
    }

    /// Rotate this state provided a direct cosine matrix of position and velocity
    pub fn rotate_by(&mut self, dcm: Matrix6<f64>) {
        let new_orbit = dcm * self.to_cartesian_vec();
        self.x_km = new_orbit[0];
        self.y_km = new_orbit[1];
        self.z_km = new_orbit[2];

        self.vx_km_s = new_orbit[3];
        self.vy_km_s = new_orbit[4];
        self.vz_km_s = new_orbit[5];
    }

    /// Rotate this state provided a direct cosine matrix of position and velocity
    pub fn with_rotation_by(&self, dcm: Matrix6<f64>) -> Self {
        let mut me = *self;
        me.rotate_by(dcm);
        me
    }

    /// Rotate the position and the velocity of this state provided a direct cosine matrix of position and velocity
    /// WARNING: You only want to use this if you'll only be using the position components of the rotated state.
    /// This does not account for the transport theorem and therefore is physically WRONG.
    pub fn position_rotated_by(&mut self, dcm: Matrix3<f64>) {
        let new_radius = dcm * self.radius();
        self.x_km = new_radius[0];
        self.y_km = new_radius[1];
        self.z_km = new_radius[2];

        let new_velocity = dcm * self.velocity();
        self.vx_km_s = new_velocity[0];
        self.vy_km_s = new_velocity[1];
        self.vz_km_s = new_velocity[2];
    }

    /// Rotate the position of this state provided a direct cosine matrix of position and velocity
    /// WARNING: You only want to use this if you'll only be using the position components of the rotated state.
    /// This does not account for the transport theorem and therefore is physically WRONG.
    pub fn with_position_rotated_by(&self, dcm: Matrix3<f64>) -> Self {
        let mut me = *self;
        me.position_rotated_by(dcm);
        me
    }

    /// Copies the current state but sets the STM to identity
    pub fn with_stm(self) -> Self {
        let mut me = self;
        me.enable_stm();
        me
    }

    /// Copies the current state but disables the STM
    pub fn without_stm(self) -> Self {
        let mut me = self;
        me.disable_stm();
        me
    }

    /// Returns the distance in kilometers between this state and a point assumed to be in the same frame.
    pub fn distance_to_point(&self, other: &Vector3<f64>) -> f64 {
        ((self.x_km - other.x).powi(2)
            + (self.y_km - other.y).powi(2)
            + (self.z_km - other.z).powi(2))
        .sqrt()
    }

    /// Use the current orbit as a template to generate mission design objectives.
    /// Note: this sets the objective tolerances to be quite tight, so consider modifying them.
    pub fn to_objectives(&self, params: &[StateParameter]) -> Result<Vec<Objective>, NyxError> {
        let mut rtn = Vec::with_capacity(params.len());
        for parameter in params {
            rtn.push(Objective::new(*parameter, self.value(*parameter)?));
        }
        Ok(rtn)
    }

    /// Create a multivariate normal dispersion structure from this orbit with the provided mean and covariance,
    /// specified as {X, Y, Z, VX, VY, VZ} in km and km/s
    pub fn disperse(
        &self,
        mean: Vector6<f64>,
        cov: Matrix6<f64>,
    ) -> Result<MultivariateNormal<Orbit>, NyxError> {
        MultivariateNormal::new(
            *self,
            vec![
                StateParameter::X,
                StateParameter::Y,
                StateParameter::Z,
                StateParameter::VX,
                StateParameter::VY,
                StateParameter::VZ,
            ],
            mean,
            cov,
        )
    }

    /// Create a multivariate normal dispersion structure from this orbit with the provided covariance,
    /// specified as {X, Y, Z, VX, VY, VZ} in km and km/s
    pub fn disperse_zero_mean(
        &self,
        cov: Matrix6<f64>,
    ) -> Result<MultivariateNormal<Orbit>, NyxError> {
        MultivariateNormal::zero_mean(
            *self,
            vec![
                StateParameter::X,
                StateParameter::Y,
                StateParameter::Z,
                StateParameter::VX,
                StateParameter::VY,
                StateParameter::VZ,
            ],
            cov,
        )
    }
}

#[cfg_attr(feature = "python", pymethods)]
impl Orbit {
    #[cfg(feature = "python")]
    #[getter]
    fn get_epoch(&self) -> Epoch {
        self.epoch
    }
    #[cfg(feature = "python")]
    #[getter]
    fn x_km(&self) -> f64 {
        self.x_km
    }
    #[cfg(feature = "python")]
    #[getter]
    fn y_km(&self) -> f64 {
        self.y_km
    }
    #[cfg(feature = "python")]
    #[getter]
    fn z_km(&self) -> f64 {
        self.z_km
    }
    #[cfg(feature = "python")]
    #[getter]
    fn vx_km_s(&self) -> f64 {
        self.vx_km_s
    }
    #[cfg(feature = "python")]
    #[getter]
    fn vy_km_s(&self) -> f64 {
        self.vy_km_s
    }
    #[cfg(feature = "python")]
    #[getter]
    fn vz_km_s(&self) -> f64 {
        self.vz_km_s
    }
    /// Returns the magnitude of the radius vector in km
    pub fn rmag_km(&self) -> f64 {
        (self.x_km.powi(2) + self.y_km.powi(2) + self.z_km.powi(2)).sqrt()
    }

    /// Returns the magnitude of the velocity vector in km/s
    pub fn vmag_km_s(&self) -> f64 {
        (self.vx_km_s.powi(2) + self.vy_km_s.powi(2) + self.vz_km_s.powi(2)).sqrt()
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

    /// Returns the orbital momentum value on the X axis
    pub fn hx_km2_s(&self) -> f64 {
        self.hvec()[0]
    }

    /// Returns the orbital momentum value on the Y axis
    pub fn hy_km2_s(&self) -> f64 {
        self.hvec()[1]
    }

    /// Returns the orbital momentum value on the Z axis
    pub fn hz_km2_s(&self) -> f64 {
        self.hvec()[2]
    }

    /// Returns the norm of the orbital momentum
    pub fn hmag_km2_s(&self) -> f64 {
        self.hvec().norm()
    }

    /// Returns the specific mechanical energy in km^2/s^2
    pub fn energy_km2_s2(&self) -> f64 {
        match self.frame {
            Frame::Geoid { gm, .. } | Frame::Celestial { gm, .. } => {
                self.vmag_km_s().powi(2) / 2.0 - gm / self.rmag_km()
            }
            _ => panic!("orbital energy not defined in this frame"),
        }
    }

    /// Returns the semi-major axis in km
    pub fn sma_km(&self) -> f64 {
        match self.frame {
            Frame::Geoid { gm, .. } | Frame::Celestial { gm, .. } => {
                -gm / (2.0 * self.energy_km2_s2())
            }
            _ => panic!("sma not defined in this frame"),
        }
    }

    /// Mutates this orbit to change the SMA
    pub fn set_sma_km(&mut self, new_sma_km: f64) {
        let me = Self::keplerian(
            new_sma_km,
            self.ecc(),
            self.inc_deg(),
            self.raan_deg(),
            self.aop_deg(),
            self.ta_deg(),
            self.epoch,
            self.frame,
        );

        self.x_km = me.x_km;
        self.y_km = me.y_km;
        self.z_km = me.z_km;
        self.vx_km_s = me.vx_km_s;
        self.vy_km_s = me.vy_km_s;
        self.vz_km_s = me.vz_km_s;
    }

    /// Returns the SMA altitude in km
    pub fn sma_altitude_km(&self) -> f64 {
        self.sma_km() - self.frame.equatorial_radius()
    }

    /// Returns the period in seconds
    pub fn period(&self) -> Duration {
        match self.frame {
            Frame::Geoid { gm, .. } | Frame::Celestial { gm, .. } => {
                2.0 * PI * (self.sma_km().powi(3) / gm).sqrt() * Unit::Second
            }
            _ => panic!("orbital period not defined in this frame"),
        }
    }

    /// Returns the eccentricity (no unit)
    pub fn ecc(&self) -> f64 {
        self.evec().norm()
    }

    /// Mutates this orbit to change the ECC
    pub fn set_ecc(&mut self, new_ecc: f64) {
        let me = Self::keplerian(
            self.sma_km(),
            new_ecc,
            self.inc_deg(),
            self.raan_deg(),
            self.aop_deg(),
            self.ta_deg(),
            self.epoch,
            self.frame,
        );

        self.x_km = me.x_km;
        self.y_km = me.y_km;
        self.z_km = me.z_km;
        self.vx_km_s = me.vx_km_s;
        self.vy_km_s = me.vy_km_s;
        self.vz_km_s = me.vz_km_s;
    }

    /// Returns the inclination in degrees
    pub fn inc_deg(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                (self.hvec()[2] / self.hmag_km2_s()).acos().to_degrees()
            }
            _ => panic!("inclination not defined in this frame"),
        }
    }

    /// Mutates this orbit to change the INC
    pub fn set_inc_deg(&mut self, new_inc: f64) {
        let me = Self::keplerian(
            self.sma_km(),
            self.ecc(),
            new_inc,
            self.raan_deg(),
            self.aop_deg(),
            self.ta_deg(),
            self.epoch,
            self.frame,
        );

        self.x_km = me.x_km;
        self.y_km = me.y_km;
        self.z_km = me.z_km;
        self.vx_km_s = me.vx_km_s;
        self.vy_km_s = me.vy_km_s;
        self.vz_km_s = me.vz_km_s;
    }

    /// Returns the argument of periapsis in degrees
    pub fn aop_deg(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                let n = Vector3::new(0.0, 0.0, 1.0).cross(&self.hvec());
                let cos_aop = n.dot(&self.evec()) / (n.norm() * self.ecc());
                let aop = cos_aop.acos();
                if aop.is_nan() {
                    if cos_aop > 1.0 {
                        180.0
                    } else {
                        0.0
                    }
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
    pub fn set_aop_deg(&mut self, new_aop: f64) {
        let me = Self::keplerian(
            self.sma_km(),
            self.ecc(),
            self.inc_deg(),
            self.raan_deg(),
            new_aop,
            self.ta_deg(),
            self.epoch,
            self.frame,
        );

        self.x_km = me.x_km;
        self.y_km = me.y_km;
        self.z_km = me.z_km;
        self.vx_km_s = me.vx_km_s;
        self.vy_km_s = me.vy_km_s;
        self.vz_km_s = me.vz_km_s;
    }

    /// Returns the right ascension of ther ascending node in degrees
    pub fn raan_deg(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                let n = Vector3::new(0.0, 0.0, 1.0).cross(&self.hvec());
                let cos_raan = n[0] / n.norm();
                let raan = cos_raan.acos();
                if raan.is_nan() {
                    if cos_raan > 1.0 {
                        180.0
                    } else {
                        0.0
                    }
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
    pub fn set_raan_deg(&mut self, new_raan: f64) {
        let me = Self::keplerian(
            self.sma_km(),
            self.ecc(),
            self.inc_deg(),
            new_raan,
            self.aop_deg(),
            self.ta_deg(),
            self.epoch,
            self.frame,
        );

        self.x_km = me.x_km;
        self.y_km = me.y_km;
        self.z_km = me.z_km;
        self.vx_km_s = me.vx_km_s;
        self.vy_km_s = me.vy_km_s;
        self.vz_km_s = me.vz_km_s;
    }

    /// Returns the true anomaly in degrees between 0 and 360.0
    ///
    /// NOTE: This function will emit a warning stating that the TA should be avoided if in a very near circular orbit
    /// Code from https://github.com/ChristopherRabotin/GMAT/blob/80bde040e12946a61dae90d9fc3538f16df34190/src/gmatutil/util/StateConversionUtil.cpp#L6835
    ///
    /// LIMITATION: For an orbit whose true anomaly is (very nearly) 0.0 or 180.0, this function may return either 0.0 or 180.0 with a very small time increment.
    /// This is due to the precision of the cosine calculation: if the arccosine calculation is out of bounds, the sign of the cosine of the true anomaly is used
    /// to determine whether the true anomaly should be 0.0 or 180.0. **In other words**, there is an ambiguity in the computation in the true anomaly exactly at 180.0 and 0.0.
    pub fn ta_deg(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                if self.ecc() < ECC_EPSILON {
                    warn!(
                        "true anomaly ill-defined for circular orbit (e = {})",
                        self.ecc()
                    );
                }
                let cos_nu = self.evec().dot(&self.radius()) / (self.ecc() * self.rmag_km());
                // If we're close the valid bounds, let's just do a sign check and return the true anomaly
                let ta = cos_nu.acos();
                if ta.is_nan() {
                    if cos_nu > 1.0 {
                        180.0
                    } else {
                        0.0
                    }
                } else if self.radius().dot(&self.velocity()) < 0.0 {
                    (2.0 * PI - ta).to_degrees()
                } else {
                    ta.to_degrees()
                }
            }
            _ => panic!("true anomaly not defined in this frame"),
        }
    }

    /// Mutates this orbit to change the TA
    pub fn set_ta_deg(&mut self, new_ta: f64) {
        let me = Self::keplerian(
            self.sma_km(),
            self.ecc(),
            self.inc_deg(),
            self.raan_deg(),
            self.aop_deg(),
            new_ta,
            self.epoch,
            self.frame,
        );

        self.x_km = me.x_km;
        self.y_km = me.y_km;
        self.z_km = me.z_km;
        self.vx_km_s = me.vx_km_s;
        self.vy_km_s = me.vy_km_s;
        self.vz_km_s = me.vz_km_s;
    }

    /// Returns the true longitude in degrees
    pub fn tlong_deg(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                // Angles already in degrees
                between_0_360(self.aop_deg() + self.raan_deg() + self.ta_deg())
            }
            _ => panic!("true longitude not defined in this frame"),
        }
    }

    /// Returns the argument of latitude in degrees
    ///
    /// NOTE: If the orbit is near circular, the AoL will be computed from the true longitude
    /// instead of relying on the ill-defined true anomaly.
    pub fn aol_deg(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                between_0_360(if self.ecc() < ECC_EPSILON {
                    self.tlong_deg() - self.raan_deg()
                } else {
                    self.aop_deg() + self.ta_deg()
                })
            }
            _ => panic!("argument of latitude not defined in this frame"),
        }
    }

    /// Returns the radius of periapsis (or perigee around Earth), in kilometers.
    pub fn periapsis_km(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => self.sma_km() * (1.0 - self.ecc()),
            _ => panic!("periapsis not defined in this frame"),
        }
    }

    /// Returns the radius of apoapsis (or apogee around Earth), in kilometers.
    pub fn apoapsis_km(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => self.sma_km() * (1.0 + self.ecc()),
            _ => panic!("apoapsis not defined in this frame"),
        }
    }

    /// Returns the altitude of periapsis (or perigee around Earth), in kilometers.
    pub fn periapsis_altitude_km(&self) -> f64 {
        self.periapsis_km() - self.frame.equatorial_radius()
    }

    /// Returns the altitude of apoapsis (or apogee around Earth), in kilometers.
    pub fn apoapsis_altitude_km(&self) -> f64 {
        self.apoapsis_km() - self.frame.equatorial_radius()
    }

    /// Returns the eccentric anomaly in degrees
    ///
    /// This is a conversion from GMAT's StateConversionUtil::TrueToEccentricAnomaly
    pub fn ea_deg(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                let (sin_ta, cos_ta) = self.ta_deg().to_radians().sin_cos();
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
    pub fn fpa_deg(&self) -> f64 {
        let nu = self.ta_deg().to_radians();
        let ecc = self.ecc();
        let denom = (1.0 + 2.0 * ecc * nu.cos() + ecc.powi(2)).sqrt();
        let sin_fpa = ecc * nu.sin() / denom;
        let cos_fpa = 1.0 + ecc * nu.cos() / denom;
        sin_fpa.atan2(cos_fpa).to_degrees()
    }

    /// Returns the mean anomaly in degrees
    ///
    /// This is a conversion from GMAT's StateConversionUtil::TrueToMeanAnomaly
    pub fn ma_deg(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                if self.ecc() < 1.0 {
                    between_0_360(
                        (self.ea_deg().to_radians()
                            - self.ecc() * self.ea_deg().to_radians().sin())
                        .to_degrees(),
                    )
                } else if self.ecc() > 1.0 {
                    info!(
                        "computing the hyperbolic anomaly (ecc = {:.6} @ {})",
                        self.ecc(),
                        self.epoch
                    );
                    // From GMAT's TrueToHyperbolicAnomaly
                    ((self.ta_deg().to_radians().sin() * (self.ecc().powi(2) - 1.0)).sqrt()
                        / (1.0 + self.ecc() * self.ta_deg().to_radians().cos()))
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
    pub fn semi_parameter_km(&self) -> f64 {
        match self.frame {
            Frame::Celestial { .. } | Frame::Geoid { .. } => {
                self.sma_km() * (1.0 - self.ecc().powi(2))
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
        if self.inc_deg() > 180.0 {
            info!("Brouwer Mean Short only applicable for inclinations less than 180.0");
            false
        } else if self.ecc() >= 1.0 || self.ecc() < 0.0 {
            info!("Brouwer Mean Short only applicable for elliptical orbits");
            false
        } else if self.periapsis_km() < 3000.0 {
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
    pub fn geodetic_longitude_deg(&self) -> f64 {
        match self.frame {
            Frame::Geoid { .. } => between_0_360(self.y_km.atan2(self.x_km).to_degrees()),
            _ => panic!("geodetic elements only defined in a Geoid frame"),
        }
    }

    /// Returns the geodetic latitude (φ) in degrees. Value is between -180 and +180 degrees.
    ///
    /// Reference: Vallado, 4th Ed., Algorithm 12 page 172.
    pub fn geodetic_latitude_deg(&self) -> f64 {
        match self.frame {
            Frame::Geoid {
                flattening,
                semi_major_radius,
                ..
            } => {
                if !self.frame.is_body_fixed() {
                    warn!("Computation of geodetic latitude must be done in a body fixed frame and {} is not one!", self.frame);
                }
                let eps = 1e-12;
                let max_attempts = 20;
                let mut attempt_no = 0;
                let r_delta = (self.x_km.powi(2) + self.y_km.powi(2)).sqrt();
                let mut latitude = (self.z_km / self.rmag_km()).asin();
                let e2 = flattening * (2.0 - flattening);
                loop {
                    attempt_no += 1;
                    let c_earth =
                        semi_major_radius / ((1.0 - e2 * (latitude).sin().powi(2)).sqrt());
                    let new_latitude = (self.z_km + c_earth * e2 * (latitude).sin()).atan2(r_delta);
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
    pub fn geodetic_height_km(&self) -> f64 {
        match self.frame {
            Frame::Geoid {
                flattening,
                semi_major_radius,
                ..
            } => {
                if !self.frame.is_body_fixed() {
                    warn!("Computation of geodetic height must be done in a body fixed frame and {} is not one!", self.frame);
                }
                let e2 = flattening * (2.0 - flattening);
                let latitude = self.geodetic_latitude_deg().to_radians();
                let sin_lat = latitude.sin();
                if (latitude - 1.0).abs() < 0.1 {
                    // We are near poles, let's use another formulation.
                    let s_earth = (semi_major_radius * (1.0 - flattening).powi(2))
                        / ((1.0 - e2 * sin_lat.powi(2)).sqrt());
                    self.z_km / latitude.sin() - s_earth
                } else {
                    let c_earth = semi_major_radius / ((1.0 - e2 * sin_lat.powi(2)).sqrt());
                    let r_delta = (self.x_km.powi(2) + self.y_km.powi(2)).sqrt();
                    r_delta / latitude.cos() - c_earth
                }
            }
            _ => panic!("geodetic elements only defined in a Geoid frame"),
        }
    }

    /// Returns the right ascension of this orbit in degrees
    pub fn right_ascension_deg(&self) -> f64 {
        between_0_360((self.y_km.atan2(self.x_km)).to_degrees())
    }

    /// Returns the declination of this orbit in degrees
    pub fn declination_deg(&self) -> f64 {
        between_pm_180((self.z_km / self.rmag_km()).asin().to_degrees())
    }

    /// Returns the semi minor axis in km, includes code for a hyperbolic orbit
    pub fn semi_minor_axis_km(&self) -> f64 {
        if self.ecc() <= 1.0 {
            self.sma_km() * (1.0 - self.ecc().powi(2)).sqrt()
        } else {
            self.hmag_km2_s().powi(2) / (self.frame.gm() * (self.ecc().powi(2) - 1.0).sqrt())
        }
    }

    /// Returns the velocity declination of this orbit in degrees
    pub fn velocity_declination_deg(&self) -> f64 {
        between_pm_180((self.vz_km_s / self.vmag_km_s()).asin().to_degrees())
    }

    pub fn b_plane(&self) -> Result<BPlane, AstroError> {
        BPlane::new(*self)
    }

    /// Returns the $C_3$ of this orbit in km^2/s^2
    pub fn c3_km2_s2(&self) -> f64 {
        -self.frame.gm() / self.sma_km()
    }

    /// Returns the radius of periapse in kilometers for the provided turn angle of this hyperbolic orbit.
    pub fn vinf_periapsis_km(&self, turn_angle_degrees: f64) -> Result<f64, NyxError> {
        if self.ecc() <= 1.0 {
            Err(NyxError::NotHyperbolic {
                msg: "Orbit is not hyperbolic. Convert to target object first".to_string(),
            })
        } else {
            let cos_rho = (0.5 * (PI - turn_angle_degrees.to_radians())).cos();

            Ok((1.0 / cos_rho - 1.0) * self.frame.gm() / self.vmag_km_s().powi(2))
        }
    }

    /// Returns the turn angle in degrees for the provided radius of periapse passage of this hyperbolic orbit
    pub fn vinf_turn_angle_deg(&self, periapsis_km: f64) -> Result<f64, NyxError> {
        if self.ecc() <= 1.0 {
            Err(NyxError::NotHyperbolic {
                msg: "Orbit is not hyperbolic. Convert to target object first".to_string(),
            })
        } else {
            let rho =
                (1.0 / (1.0 + self.vmag_km_s().powi(2) * (periapsis_km / self.frame.gm()))).acos();
            Ok(between_0_360((PI - 2.0 * rho).to_degrees()))
        }
    }

    /// Returns the hyperbolic anomaly in degrees between 0 and 360.0
    pub fn hyperbolic_anomaly_deg(&self) -> Result<f64, NyxError> {
        if self.ecc() <= 1.0 {
            Err(NyxError::NotHyperbolic {
                msg: "Orbit is not hyperbolic so there is no hyperbolic anomaly.".to_string(),
            })
        } else {
            let (sin_ta, cos_ta) = self.ta_deg().to_radians().sin_cos();
            let sinh_h = (sin_ta * (self.ecc().powi(2) - 1.0).sqrt()) / (1.0 + self.ecc() * cos_ta);
            Ok(between_0_360(sinh_h.asinh().to_degrees()))
        }
    }

    /// Sets the STM of this state of identity, which also enables computation of the STM for spacecraft navigation
    pub fn enable_stm(&mut self) {
        self.stm = Some(Matrix6::identity());
    }

    /// Disable the STM of this state
    pub fn disable_stm(&mut self) {
        self.stm = None;
    }

    /// Returns the root sum square error between this state and the other, in kilometers for the position and kilometers per second in velocity
    pub fn rss(&self, other: &Self) -> (f64, f64) {
        rss_orbit_errors(self, other)
    }

    /// Returns whether this orbit and another are equal within the specified radial and velocity absolute tolerances
    pub fn eq_within(&self, other: &Self, radial_tol: f64, velocity_tol: f64) -> bool {
        self.epoch == other.epoch
            && (self.x_km - other.x_km).abs() < radial_tol
            && (self.y_km - other.y_km).abs() < radial_tol
            && (self.z_km - other.z_km).abs() < radial_tol
            && (self.vx_km_s - other.vx_km_s).abs() < velocity_tol
            && (self.vy_km_s - other.vy_km_s).abs() < velocity_tol
            && (self.vz_km_s - other.vz_km_s).abs() < velocity_tol
            && self.frame == other.frame
            && self.stm.is_some() == other.stm.is_some()
            && if self.stm.is_some() {
                self.stm.unwrap() == other.stm.unwrap()
            } else {
                true
            }
    }

    /// Adjusts the true anomaly of this orbit using the mean anomaly.
    ///
    /// # Astrodynamics note
    /// This is akin to a two body propagation.
    pub fn at_epoch(&self, new_epoch: Epoch) -> Result<Self, NyxError> {
        let m0_rad = self.ma_deg().to_radians();
        let mt_rad = m0_rad
            + (self.frame.gm() / self.sma_km().powi(3)).sqrt()
                * (new_epoch - self.epoch).to_seconds();

        Self::keplerian_mean_anomaly(
            self.sma_km(),
            self.ecc(),
            self.inc_deg(),
            self.raan_deg(),
            self.aop_deg(),
            mt_rad.to_degrees(),
            new_epoch,
            self.frame,
        )
    }

    /// Prints this orbit in Cartesian form
    #[cfg(feature = "python")]
    fn __repr__(&self) -> String {
        format!("{self}")
    }

    /// Prints this orbit in Keplerian form
    #[cfg(feature = "python")]
    fn __str__(&self) -> String {
        format!("{self:x}")
    }

    #[cfg(feature = "python")]
    fn __richcmp__(&self, other: &Self, op: CompareOp) -> Result<bool, NyxError> {
        match op {
            CompareOp::Eq => Ok(self == other),
            CompareOp::Ne => Ok(self != other),
            _ => Err(NyxError::CustomError(format!("{op:?} not available"))),
        }
    }

    /// Creates a new Orbit in the provided frame at the provided Epoch given the position and velocity components.
    ///
    /// **Units:** km, km, km, km/s, km/s, km/s
    #[cfg(feature = "python")]
    #[new]
    fn py_new(
        x_km: f64,
        y_km: f64,
        z_km: f64,
        vx_km_s: f64,
        vy_km_s: f64,
        vz_km_s: f64,
        epoch: Epoch,
        frame: PyRef<PyFrame>,
    ) -> Self {
        Self::cartesian(
            x_km,
            y_km,
            z_km,
            vx_km_s,
            vy_km_s,
            vz_km_s,
            epoch,
            frame.inner,
        )
    }

    #[cfg(feature = "python")]
    #[classmethod]
    fn load(_cls: &PyType, path: &str) -> Result<Self, ConfigError> {
        let serde = OrbitSerde::load(path)?;

        let cosm = Cosm::de438();

        Self::from_config(serde, cosm)
    }

    #[cfg(feature = "python")]
    #[classmethod]
    fn load_many(_cls: &PyType, path: &str) -> Result<Vec<Self>, ConfigError> {
        let orbits = OrbitSerde::load_many(path)?;

        let cosm = Cosm::de438();

        let mut selves = Vec::with_capacity(orbits.len());

        for serde in orbits {
            selves.push(Self::from_config(serde, cosm.clone())?);
        }

        Ok(selves)
    }

    #[cfg(feature = "python")]
    #[classmethod]
    fn load_named(_cls: &PyType, path: &str) -> Result<BTreeMap<String, Self>, ConfigError> {
        let orbits = OrbitSerde::load_named(path)?;

        let cosm = Cosm::de438();

        let mut selves = BTreeMap::new();

        for (k, v) in orbits {
            selves.insert(k, Self::from_config(v, cosm.clone())?);
        }

        Ok(selves)
    }

    /// Creates a new Orbit around the provided Celestial or Geoid frame from the Keplerian orbital elements.
    ///
    /// **Units:** km, none, degrees, degrees, degrees, degrees
    ///
    /// WARNING: This function will panic if the singularities in the conversion are expected.
    /// NOTE: The state is defined in Cartesian coordinates as they are non-singular. This causes rounding
    /// errors when creating a state from its Keplerian orbital elements (cf. the state tests).
    /// One should expect these errors to be on the order of 1e-12.
    #[cfg(feature = "python")]
    #[classmethod]
    fn from_keplerian(
        _cls: &PyType,
        sma_km: f64,
        ecc: f64,
        inc_deg: f64,
        raan_deg: f64,
        aop_deg: f64,
        ta_deg: f64,
        epoch: Epoch,
        frame: PyRef<PyFrame>,
    ) -> Self {
        Self::keplerian(
            sma_km,
            ecc,
            inc_deg,
            raan_deg,
            aop_deg,
            ta_deg,
            epoch,
            frame.inner,
        )
    }

    #[cfg(feature = "python")]
    #[classmethod]
    fn from_keplerian_mean_anomaly(
        _cls: &PyType,
        sma_km: f64,
        ecc: f64,
        inc_deg: f64,
        raan_deg: f64,
        aop_deg: f64,
        ma_deg: f64,
        epoch: Epoch,
        frame: PyRef<PyFrame>,
    ) -> Result<Self, NyxError> {
        Self::keplerian_mean_anomaly(
            sma_km,
            ecc,
            inc_deg,
            raan_deg,
            aop_deg,
            ma_deg,
            epoch,
            frame.inner,
        )
    }

    /// Creates a new Orbit from the provided semi-major axis altitude in kilometers
    #[cfg(feature = "python")]
    #[classmethod]
    fn from_keplerian_altitude(
        _cls: &PyType,
        sma_altitude_km: f64,
        ecc: f64,
        inc_deg: f64,
        raan_deg: f64,
        aop_deg: f64,
        ta_deg: f64,
        epoch: Epoch,
        frame: PyRef<PyFrame>,
    ) -> Self {
        Self::keplerian_altitude(
            sma_altitude_km,
            ecc,
            inc_deg,
            raan_deg,
            aop_deg,
            ta_deg,
            epoch,
            frame.inner,
        )
    }

    /// Creates a new Orbit from the provided altitudes of apoapsis and periapsis, in kilometers
    #[cfg(feature = "python")]
    #[classmethod]
    fn from_keplerian_apsis_altitude(
        _cls: &PyType,
        apo_alt_km: f64,
        peri_alt_km: f64,
        inc_deg: f64,
        raan_deg: f64,
        aop_deg: f64,
        ta_deg: f64,
        epoch: Epoch,
        frame: PyRef<PyFrame>,
    ) -> Self {
        Self::keplerian_apsis_altitude(
            apo_alt_km,
            peri_alt_km,
            inc_deg,
            raan_deg,
            aop_deg,
            ta_deg,
            epoch,
            frame.inner,
        )
    }

    #[cfg(feature = "python")]
    #[setter(sma_km)]
    fn py_set_sma(&mut self, new_sma_km: f64) -> PyResult<()> {
        self.set_sma_km(new_sma_km);
        Ok(())
    }

    #[cfg(feature = "python")]
    #[setter(ecc)]
    fn py_set_ecc(&mut self, new_ecc: f64) -> PyResult<()> {
        self.set_ecc(new_ecc);
        Ok(())
    }

    #[cfg(feature = "python")]
    #[setter(inc_deg)]
    fn py_set_inc(&mut self, new_inc_deg: f64) -> PyResult<()> {
        self.set_inc_deg(new_inc_deg);
        Ok(())
    }

    #[cfg(feature = "python")]
    #[setter(inc_deg)]
    fn py_set_raan(&mut self, new_raan_deg: f64) -> PyResult<()> {
        self.set_raan_deg(new_raan_deg);
        Ok(())
    }

    #[cfg(feature = "python")]
    #[setter(aop_deg)]
    fn py_set_aop(&mut self, new_aop_deg: f64) -> PyResult<()> {
        self.set_aop_deg(new_aop_deg);
        Ok(())
    }

    #[cfg(feature = "python")]
    #[setter(ta_deg)]
    fn py_set_ta(&mut self, new_ta_deg: f64) -> PyResult<()> {
        self.set_ta_deg(new_ta_deg);
        Ok(())
    }

    /// Returns the value of the provided state parameter if available
    #[cfg(feature = "python")]
    fn value_of(&self, param: StateParameter) -> Result<f64, NyxError> {
        self.value(param)
    }
}

impl PartialEq for Orbit {
    /// Two states are equal if their position are equal within one centimeter and their velocities within one centimeter per second.
    fn eq(&self, other: &Orbit) -> bool {
        let radial_tol = 1e-5; // centimeter
        let velocity_tol = 1e-5; // centimeter per second
        self.eq_within(other, radial_tol, velocity_tol)
    }
}

impl Add for Orbit {
    type Output = Orbit;

    /// Add one state from another. Frame must be manually changed if needed. STM will be copied from &self.
    fn add(self, other: Orbit) -> Orbit {
        Orbit {
            x_km: self.x_km + other.x_km,
            y_km: self.y_km + other.y_km,
            z_km: self.z_km + other.z_km,
            vx_km_s: self.vx_km_s + other.vx_km_s,
            vy_km_s: self.vy_km_s + other.vy_km_s,
            vz_km_s: self.vz_km_s + other.vz_km_s,
            epoch: self.epoch,
            frame: self.frame,
            stm: self.stm,
        }
    }
}

impl Sub for Orbit {
    type Output = Orbit;

    /// Subtract one state from another. STM will be copied from &self.
    fn sub(self, other: Orbit) -> Orbit {
        Orbit {
            x_km: self.x_km - other.x_km,
            y_km: self.y_km - other.y_km,
            z_km: self.z_km - other.z_km,
            vx_km_s: self.vx_km_s - other.vx_km_s,
            vy_km_s: self.vy_km_s - other.vy_km_s,
            vz_km_s: self.vz_km_s - other.vz_km_s,
            epoch: self.epoch,
            frame: self.frame,
            stm: self.stm,
        }
    }
}

impl Neg for Orbit {
    type Output = Orbit;

    /// Subtract one state from another. STM will be copied from &self.
    fn neg(self) -> Self::Output {
        Orbit {
            x_km: -self.x_km,
            y_km: -self.y_km,
            z_km: -self.z_km,
            vx_km_s: -self.vx_km_s,
            vy_km_s: -self.vy_km_s,
            vz_km_s: -self.vz_km_s,
            epoch: self.epoch,
            frame: self.frame,
            stm: self.stm,
        }
    }
}

impl Add for &Orbit {
    type Output = Orbit;

    /// Add one state from another. Frame must be manually changed if needed. STM will be copied from &self.
    fn add(self, other: &Orbit) -> Orbit {
        Orbit {
            x_km: self.x_km + other.x_km,
            y_km: self.y_km + other.y_km,
            z_km: self.z_km + other.z_km,
            vx_km_s: self.vx_km_s + other.vx_km_s,
            vy_km_s: self.vy_km_s + other.vy_km_s,
            vz_km_s: self.vz_km_s + other.vz_km_s,
            epoch: self.epoch,
            frame: self.frame,
            stm: self.stm,
        }
    }
}

impl AddAssign for Orbit {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl Sub for &Orbit {
    type Output = Orbit;

    /// Subtract one state from another. STM will be copied from &self.
    fn sub(self, other: &Orbit) -> Orbit {
        Orbit {
            x_km: self.x_km - other.x_km,
            y_km: self.y_km - other.y_km,
            z_km: self.z_km - other.z_km,
            vx_km_s: self.vx_km_s - other.vx_km_s,
            vy_km_s: self.vy_km_s - other.vy_km_s,
            vz_km_s: self.vz_km_s - other.vz_km_s,
            epoch: self.epoch,
            frame: self.frame,
            stm: self.stm,
        }
    }
}

impl SubAssign for Orbit {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl Neg for &Orbit {
    type Output = Orbit;

    /// Subtract one state from another. STM will be copied form &self.
    fn neg(self) -> Self::Output {
        Orbit {
            x_km: -self.x_km,
            y_km: -self.y_km,
            z_km: -self.z_km,
            vx_km_s: -self.vx_km_s,
            vy_km_s: -self.vy_km_s,
            vz_km_s: -self.vz_km_s,
            epoch: self.epoch,
            frame: self.frame,
            stm: self.stm,
        }
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::Display for Orbit {
    // Prints as Cartesian in floating point with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "[{}] {}\tposition = [{}, {}, {}] km\tvelocity = [{}, {}, {}] km/s",
            self.frame,
            self.epoch,
            format!("{:.*}", decimals, self.x_km),
            format!("{:.*}", decimals, self.y_km),
            format!("{:.*}", decimals, self.z_km),
            format!("{:.*}", decimals, self.vx_km_s),
            format!("{:.*}", decimals, self.vy_km_s),
            format!("{:.*}", decimals, self.vz_km_s)
        )
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::LowerExp for Orbit {
    // Prints as Cartesian in scientific notation with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "[{}] {}\tposition = [{}, {}, {}] km\tvelocity = [{}, {}, {}] km/s",
            self.frame,
            self.epoch,
            format!("{:.*e}", decimals, self.x_km),
            format!("{:.*e}", decimals, self.y_km),
            format!("{:.*e}", decimals, self.z_km),
            format!("{:.*e}", decimals, self.vx_km_s),
            format!("{:.*e}", decimals, self.vy_km_s),
            format!("{:.*e}", decimals, self.vz_km_s)
        )
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::UpperExp for Orbit {
    // Prints as Cartesian in scientific notation with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "[{}] {}\tposition = [{}, {}, {}] km\tvelocity = [{}, {}, {}] km/s",
            self.frame,
            self.epoch,
            format!("{:.*E}", decimals, self.x_km),
            format!("{:.*E}", decimals, self.y_km),
            format!("{:.*E}", decimals, self.z_km),
            format!("{:.*E}", decimals, self.vx_km_s),
            format!("{:.*E}", decimals, self.vy_km_s),
            format!("{:.*E}", decimals, self.vz_km_s)
        )
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::LowerHex for Orbit {
    // Prints the Keplerian orbital elements in floating point with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "[{}] {}\tsma = {} km\tecc = {}\tinc = {} deg\traan = {} deg\taop = {} deg\tta = {} deg",
            self.frame,
            self.epoch,
            format!("{:.*}", decimals, self.sma_km()),
            format!("{:.*}", decimals, self.ecc()),
            format!("{:.*}", decimals, self.inc_deg()),
            format!("{:.*}", decimals, self.raan_deg()),
            format!("{:.*}", decimals, self.aop_deg()),
            format!("{:.*}", decimals, self.ta_deg()),
        )
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::UpperHex for Orbit {
    // Prints the Keplerian orbital elements in scientific notation with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "[{}] {}\tsma = {} km\tecc = {}\tinc = {} deg\traan = {} deg\taop = {} deg\tta = {} deg",
            self.frame,
            self.epoch,
            format!("{:.*e}", decimals, self.sma_km()),
            format!("{:.*e}", decimals, self.ecc()),
            format!("{:.*e}", decimals, self.inc_deg()),
            format!("{:.*e}", decimals, self.raan_deg()),
            format!("{:.*e}", decimals, self.aop_deg()),
            format!("{:.*e}", decimals, self.ta_deg()),
        )
    }
}

/// Implementation of Orbit as a State for orbital dynamics with STM
impl State for Orbit {
    type Size = Const<6>;
    type VecLength = Const<42>;

    fn reset_stm(&mut self) {
        self.stm = Some(Matrix6::identity());
    }

    /// Returns a state whose position, velocity and frame are zero, and STM is I_{6x6}.
    fn zeros() -> Self {
        let frame = Frame::Celestial {
            gm: 1.0,
            ephem_path: [None, None, None],
            frame_path: [None, None, None],
        };

        Self {
            x_km: 0.0,
            y_km: 0.0,
            z_km: 0.0,
            vx_km_s: 0.0,
            vy_km_s: 0.0,
            vz_km_s: 0.0,
            epoch: Epoch::from_tai_seconds(0.0),
            frame,
            stm: Some(Matrix6::identity()),
        }
    }

    fn as_vector(&self) -> OVector<f64, Const<42>> {
        let mut as_vec = OVector::<f64, Const<42>>::zeros();
        as_vec[0] = self.x_km;
        as_vec[1] = self.y_km;
        as_vec[2] = self.z_km;
        as_vec[3] = self.vx_km_s;
        as_vec[4] = self.vy_km_s;
        as_vec[5] = self.vz_km_s;
        if let Some(stm) = self.stm {
            for (idx, stm_val) in stm.as_slice().iter().enumerate() {
                as_vec[idx + 6] = *stm_val;
            }
        }
        as_vec
    }

    fn set(&mut self, epoch: Epoch, vector: &OVector<f64, Const<42>>) {
        self.set_epoch(epoch);
        self.x_km = vector[0];
        self.y_km = vector[1];
        self.z_km = vector[2];
        self.vx_km_s = vector[3];
        self.vy_km_s = vector[4];
        self.vz_km_s = vector[5];
        // And update the STM if applicable
        if self.stm.is_some() {
            let stm_k_to_0 = Matrix6::from_column_slice(&vector.as_slice()[6..]);
            self.stm = Some(stm_k_to_0);
        }
    }

    fn stm(&self) -> Result<Matrix6<f64>, DynamicsError> {
        match self.stm {
            Some(stm) => Ok(stm),
            None => Err(DynamicsError::StateTransitionMatrixUnset),
        }
    }

    fn epoch(&self) -> Epoch {
        self.epoch
    }

    fn set_epoch(&mut self, epoch: Epoch) {
        self.epoch = epoch
    }

    fn add(self, other: OVector<f64, Self::Size>) -> Self {
        self + other
    }

    fn value(&self, param: StateParameter) -> Result<f64, NyxError> {
        match param {
            StateParameter::ApoapsisRadius => Ok(self.apoapsis_km()),
            StateParameter::AoL => Ok(self.aol_deg()),
            StateParameter::AoP => Ok(self.aop_deg()),
            StateParameter::BdotR => Ok(BPlane::new(*self).unwrap().b_r.real()),
            StateParameter::BdotT => Ok(BPlane::new(*self).unwrap().b_t.real()),
            StateParameter::BLTOF => Ok(BPlane::new(*self).unwrap().ltof_s.real()),
            StateParameter::C3 => Ok(self.c3_km2_s2()),
            StateParameter::Declination => Ok(self.declination_deg()),
            StateParameter::EccentricAnomaly => Ok(self.ea_deg()),
            StateParameter::Eccentricity => Ok(self.ecc()),
            StateParameter::Energy => Ok(self.energy_km2_s2()),
            StateParameter::FlightPathAngle => Ok(self.fpa_deg()),
            StateParameter::GeodeticHeight => Ok(self.geodetic_height_km()),
            StateParameter::GeodeticLatitude => Ok(self.geodetic_latitude_deg()),
            StateParameter::GeodeticLongitude => Ok(self.geodetic_longitude_deg()),
            StateParameter::Hmag => Ok(self.hmag_km2_s()),
            StateParameter::HX => Ok(self.hx_km2_s()),
            StateParameter::HY => Ok(self.hy_km2_s()),
            StateParameter::HZ => Ok(self.hz_km2_s()),
            StateParameter::HyperbolicAnomaly => self.hyperbolic_anomaly_deg(),
            StateParameter::Inclination => Ok(self.inc_deg()),
            StateParameter::MeanAnomaly => Ok(self.ma_deg()),
            StateParameter::PeriapsisRadius => Ok(self.periapsis_km()),
            StateParameter::Period => Ok(self.period().to_seconds()),
            StateParameter::RightAscension => Ok(self.right_ascension_deg()),
            StateParameter::RAAN => Ok(self.raan_deg()),
            StateParameter::Rmag => Ok(self.rmag_km()),
            StateParameter::SemiMinorAxis => Ok(self.semi_minor_axis_km()),
            StateParameter::SemiParameter => Ok(self.semi_parameter_km()),
            StateParameter::SMA => Ok(self.sma_km()),
            StateParameter::TrueAnomaly => Ok(self.ta_deg()),
            StateParameter::TrueLongitude => Ok(self.tlong_deg()),
            StateParameter::VelocityDeclination => Ok(self.velocity_declination_deg()),
            StateParameter::Vmag => Ok(self.vmag_km_s()),
            StateParameter::X => Ok(self.x_km),
            StateParameter::Y => Ok(self.y_km),
            StateParameter::Z => Ok(self.z_km),
            StateParameter::VX => Ok(self.vx_km_s),
            StateParameter::VY => Ok(self.vy_km_s),
            StateParameter::VZ => Ok(self.vz_km_s),
            _ => Err(NyxError::StateParameterUnavailable {
                param,
                msg: "no such parameter for orbit structure".to_string(),
            }),
        }
    }

    fn set_value(&mut self, param: StateParameter, val: f64) -> Result<(), NyxError> {
        match param {
            StateParameter::AoP => self.set_aop_deg(val),
            StateParameter::Eccentricity => self.set_ecc(val),
            StateParameter::Inclination => self.set_inc_deg(val),
            StateParameter::RAAN => self.set_raan_deg(val),
            StateParameter::SMA => self.set_sma_km(val),
            StateParameter::TrueAnomaly => self.set_ta_deg(val),
            StateParameter::X => self.x_km = val,
            StateParameter::Y => self.y_km = val,
            StateParameter::Z => self.z_km = val,
            StateParameter::Rmag => {
                // Convert the position to spherical coordinates
                let (_, θ, φ) = cartesian_to_spherical(&self.radius());
                // Convert back to cartesian after setting the new range value
                let new_radius = spherical_to_cartesian(val, θ, φ);
                self.x_km = new_radius.x;
                self.y_km = new_radius.y;
                self.z_km = new_radius.z;
            }
            StateParameter::VX => self.vx_km_s = val,
            StateParameter::VY => self.vy_km_s = val,
            StateParameter::VZ => self.vz_km_s = val,
            StateParameter::Vmag => {
                // Convert the velocity to spherical coordinates
                let (_, θ, φ) = cartesian_to_spherical(&self.velocity());
                // Convert back to cartesian after setting the new range value
                let new_radius = spherical_to_cartesian(val, θ, φ);
                self.vx_km_s = new_radius.x;
                self.vy_km_s = new_radius.y;
                self.vz_km_s = new_radius.z;
            }
            _ => {
                return Err(NyxError::StateParameterUnavailable {
                    param,
                    msg: "not settable for orbit structure with set_value".to_string(),
                })
            }
        }
        Ok(())
    }

    fn unset_stm(&mut self) {
        self.stm = None;
    }
}

impl Add<OVector<f64, Const<6>>> for Orbit {
    type Output = Self;

    /// Adds the provided state deviation to this orbit
    fn add(self, other: OVector<f64, Const<6>>) -> Self {
        let mut me = self;
        me.x_km += other[0];
        me.y_km += other[1];
        me.z_km += other[2];
        me.vx_km_s += other[3];
        me.vy_km_s += other[4];
        me.vz_km_s += other[5];

        me
    }
}

impl Default for Orbit {
    fn default() -> Self {
        Self::zeros()
    }
}

impl ConfigRepr for OrbitSerde {}

impl Configurable for Orbit {
    type IntermediateRepr = OrbitSerde;

    fn from_config(
        cfg: Self::IntermediateRepr,
        _cosm: Arc<Cosm>,
    ) -> Result<Self, crate::io::ConfigError>
    where
        Self: Sized,
    {
        let orbit: Orbit = cfg.into();
        Ok(orbit)
    }

    fn to_config(&self) -> Result<Self::IntermediateRepr, crate::io::ConfigError> {
        let serded: Self::IntermediateRepr = (*self).into();
        Ok(serded)
    }
}

/// Computes the true anomaly from the given mean anomaly for an orbit.
///
/// The computation process varies depending on whether the orbit is elliptical (eccentricity less than or equal to 1)
/// or hyperbolic (eccentricity greater than 1). In each case, the method uses an iterative algorithm to find a
/// sufficiently accurate approximation of the true anomaly.
///
/// # Arguments
///
/// * `ma_radians` - The mean anomaly in radians.
/// * `ecc` - The eccentricity of the orbit.
/// * `tol` - The tolerance for accuracy in the resulting true anomaly.
///
/// # Returns
///
/// * Ok(f64) - The calculated true anomaly in radians, if the computation was successful.
/// * Err(NyxError) - An error that occurred during the computation. Possible errors include:
///     - MaxIterReached: The iterative process did not converge within 1000 iterations.
///     - MathDomain: A mathematical error occurred during computation, such as division by zero.
///
/// # Remarks
///
/// This function uses GTDS MathSpec Equations 3-180, 3-181, and 3-186 for the iterative computation process.
///
/// If a numerical error occurs during computation, the function may return a MathDomain error. In the case of a
/// non-converging iterative process, the function will return a MaxIterReached error after 1000 iterations.
fn compute_mean_to_true_anomaly(ma_radians: f64, ecc: f64, tol: f64) -> Result<f64, NyxError> {
    let rm = ma_radians;
    if ecc <= 1.0 {
        // Elliptical orbit
        let mut e2 = rm + ecc * rm.sin(); // GTDS MathSpec Equation 3-182

        let mut iter = 0;

        loop {
            iter += 1;
            if iter > 1000 {
                return Err(NyxError::MaxIterReached {
                    msg: format!("{iter}"),
                });
            }

            // GTDS MathSpec Equation 3-180  Note: a little difference here is that it uses Cos(E) instead of Cos(E-0.5*f)
            let normalized_anomaly = 1.0 - ecc * e2.cos();

            if normalized_anomaly.abs() < MA_EPSILON {
                return Err(NyxError::MathDomain {
                    msg: format!("normalizer too small {normalized_anomaly}"),
                });
            }

            // GTDS MathSpec Equation 3-181
            let e1 = e2 - (e2 - ecc * e2.sin() - rm) / normalized_anomaly;

            if (e2 - e1).abs() < tol {
                break;
            }

            e2 = e1;
        }

        let mut e = e2;

        if e < 0.0 {
            e += TAU;
        }

        let c = (e - PI).abs();

        let mut ta = if c >= 1.0e-08 {
            let normalized_anomaly = 1.0 - ecc;

            if (normalized_anomaly).abs() < MA_EPSILON {
                return Err(NyxError::MathDomain {
                    msg: format!("normalized anomaly too small {normalized_anomaly}"),
                });
            }

            let eccentricity_ratio = (1.0 + ecc) / normalized_anomaly; // temp2 = (1+ecc)/(1-ecc)

            if eccentricity_ratio < 0.0 {
                return Err(NyxError::MathDomain {
                    msg: format!("eccentric ratio too small {eccentricity_ratio}"),
                });
            }

            let f = eccentricity_ratio.sqrt();
            let g = (e / 2.0).tan();
            // tan(TA/2) = Sqrt[(1+ecc)/(1-ecc)] * tan(E/2)
            2.0 * (f * g).atan()
        } else {
            e
        };

        if ta < 0.0 {
            ta += TAU;
        }
        Ok(ta)
    } else {
        //---------------------------------------------------------
        // hyperbolic orbit
        //---------------------------------------------------------

        // For hyperbolic orbit, anomaly is nolonger to be an angle so we cannot use mod of 2*PI to mean anomaly.
        // We need to keep its original value for calculation.
        //if (rm > PI)                       // incorrect
        //   rm = rm - TWO_PI;               // incorrect

        //f2 = ecc * Sinh(rm) - rm;          // incorrect
        //f2 = rm / 2;                       // incorrect  // GTDS MathSpec Equation 3-186
        let mut f2: f64 = 0.0; // This is the correct initial value for hyperbolic eccentric anomaly.
        let mut iter = 0;

        loop {
            iter += 1;
            if iter > 1000 {
                return Err(NyxError::MaxIterReached {
                    msg: format!("{iter}"),
                });
            }

            let normalizer = ecc * f2.cosh() - 1.0;

            if normalizer.abs() < MA_EPSILON {
                return Err(NyxError::MathDomain {
                    msg: format!("normalizer too small {normalizer}"),
                });
            }

            let f1 = f2 - (ecc * f2.sinh() - f2 - rm) / normalizer; // GTDS MathSpec Equation 3-186
            if (f2 - f1).abs() < tol {
                break;
            }
            f2 = f1;
        }

        let f = f2;
        let normalized_anomaly = ecc - 1.0;

        if normalized_anomaly.abs() < MA_EPSILON {
            return Err(NyxError::MathDomain {
                msg: format!("eccentric ratio too small {normalized_anomaly}"),
            });
        }

        let eccentricity_ratio = (ecc + 1.0) / normalized_anomaly; // temp2 = (ecc+1)/(ecc-1)

        if eccentricity_ratio < 0.0 {
            return Err(NyxError::MathDomain {
                msg: format!("eccentric ratio too small {eccentricity_ratio}"),
            });
        }

        let e = eccentricity_ratio.sqrt();
        let g = (f / 2.0).tanh();
        let mut ta = 2.0 * (e * g).atan(); // tan(TA/2) = Sqrt[(ecc+1)/(ecc-1)] * Tanh(F/2)    where: F is hyperbolic centric anomaly

        if ta < 0.0 {
            ta += TAU;
        }
        Ok(ta)
    }
}

#[test]
fn test_serde() {
    use super::Cosm;
    use crate::io::orbit::OrbitSerde;
    use crate::utils::rss_orbit_errors;
    use serde_yaml;
    use std::env;
    use std::path::PathBuf;
    use std::str::FromStr;

    // Get the path to the root directory of the current Cargo project
    let s = r#"
x_km: -9042.862234
y_km: 18536.333069
z_km: 6999.957069
vx_km_s: -3.288789
vy_km_s: -2.226285
vz_km_s: 1.646738
frame: EME2000
epoch: 2018-09-15T00:15:53.098 UTC
    "#;

    let orbit: Orbit = serde_yaml::from_str(s).unwrap();

    dbg!(&orbit);

    let cosm = Cosm::de438();
    let exp = Orbit::cartesian(
        -9042.862234,
        18536.333069,
        6999.957069,
        -3.288789,
        -2.226285,
        1.646738,
        Epoch::from_str("2018-09-15T00:15:53.098 UTC").unwrap(),
        cosm.frame("EME2000"),
    );

    assert_eq!(exp, orbit);

    let cosm = Cosm::de438();
    let orbit = Orbit::cartesian(
        -9042.862234,
        18536.333069,
        6999.957069,
        -3.288789,
        -2.226285,
        1.646738,
        Epoch::from_str("2018-09-15T00:15:53.098 UTC").unwrap(),
        cosm.frame("EME2000"),
    );

    dbg!(serde_yaml::to_string(&orbit).unwrap());

    let as_serde: OrbitSerde = orbit.into();
    println!("{}", serde_yaml::to_string(&as_serde).unwrap());

    // Try to deserialize from Keplerian
    let s = format!("sma_km: {}\necc: {}\ninc_deg: {}\nraan_deg: {}\naop_deg: {}\nta_deg: {}\nepoch: {}\nframe: {}", orbit.sma_km(), orbit.ecc(), orbit.inc_deg(), orbit.raan_deg(), orbit.aop_deg(), orbit.ta_deg(), orbit.epoch, orbit.frame);
    println!("{s}");
    let deserd: OrbitSerde = serde_yaml::from_str(&s).unwrap();
    let deser_orbit: Orbit = deserd.into();
    dbg!(deser_orbit);

    // Check that the orbits (mostly) match -- there will be rounding errors
    assert_eq!(deser_orbit.frame, orbit.frame);
    assert_eq!(deser_orbit.epoch, orbit.epoch);
    let (pos_err, vel_err) = rss_orbit_errors(&deser_orbit, &orbit);
    assert!(pos_err < 1e-5);
    assert!(vel_err < 1e-9);

    // Check that we can deserialize from a file
    let test_data: PathBuf = [
        env::var("CARGO_MANIFEST_DIR").unwrap(),
        "data".to_string(),
        "tests".to_string(),
        "config".to_string(),
        "orbit.yaml".to_string(),
    ]
    .iter()
    .collect();

    let orbit = Orbit::from_config(OrbitSerde::load(test_data).unwrap(), cosm).unwrap();
    assert_eq!(exp, orbit);

    let test_data: PathBuf = [
        env::var("CARGO_MANIFEST_DIR").unwrap(),
        "data".to_string(),
        "tests".to_string(),
        "config".to_string(),
        "orbits.yaml".to_string(),
    ]
    .iter()
    .collect();

    let serded_orbits = OrbitSerde::load_many(test_data).unwrap();
    for orbit_s in serded_orbits {
        let orbit: Orbit = orbit_s.into();
        // Check that the orbits (mostly) match -- there will be rounding errors
        assert_eq!(exp.frame, orbit.frame);
        assert_eq!(exp.epoch, orbit.epoch);
        let (pos_err, vel_err) = rss_orbit_errors(&exp, &orbit);
        assert!(pos_err < 1e-5);
        assert!(vel_err < 1e-9);
    }
}
