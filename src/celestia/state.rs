extern crate hifitime;
extern crate serde;

use self::hifitime::Epoch;
use self::serde::ser::SerializeStruct;
use self::serde::{Serialize, Serializer};
use super::na::{Matrix3, Vector3, Vector6};
use super::{Frame, FrameKind};
use celestia::xb::ephem_registry::State as XBState;
use celestia::xb::Epoch as XBEpoch;
use celestia::xb::Vector as XBVector;
use celestia::xb::{TimeRepr, TimeSystem, Unit};
use std::collections::HashMap;
use std::f64::consts::PI;
use std::f64::EPSILON;
use std::fmt;
use std::ops::{Add, Neg, Sub};
use utils::{between_0_360, between_pm_180, perpv, r1, r3};

/// If an orbit has an eccentricity below the following value, it is considered circular (only affects warning messages)
pub const ECC_EPSILON: f64 = 1e-11;

/// State defines an orbital state parameterized  by a `CelestialBody`.
///
/// Unless noted otherwise, algorithms are from GMAT 2016a [StateConversionUtil.cpp](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/StateConversionUtil.cpp).
/// Regardless of the constructor used, this struct stores all the state information in Cartesian coordinates
/// as these are always non singular.
/// _Note:_ although not yet supported, this struct may change once True of Date or other nutation frames
/// are added to the toolkit.
#[derive(Copy, Clone, Debug)]
pub struct State<'a> {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    pub dt: Epoch,
    pub frame: &'a Frame<'a>,
}

impl<'a> State<'a> {
    /// Creates a new State in the provided frame at the provided Epoch.
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
        frame: &'a Frame<'a>,
    ) -> Self {
        State {
            x,
            y,
            z,
            vx,
            vy,
            vz,
            dt,
            frame,
        }
    }

    /// Creates a new State in the provided frame at the provided Epoch in time with 0.0 velocity.
    ///
    /// **Units:** km, km, km
    pub fn from_position(x: f64, y: f64, z: f64, dt: Epoch, frame: &'a Frame<'a>) -> Self {
        State {
            x,
            y,
            z,
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
            dt,
            frame,
        }
    }

    /// Creates a new State around in the provided frame from the borrowed state vector
    ///
    /// The state vector **must** be x, y, z, vx, vy, vz. This function is a shortcut to `cartesian`
    /// and as such it has the same unit requirements.
    pub fn cartesian_vec(state: &Vector6<f64>, dt: Epoch, frame: &'a Frame<'a>) -> Self {
        State {
            x: state[(0, 0)],
            y: state[(1, 0)],
            z: state[(2, 0)],
            vx: state[(3, 0)],
            vy: state[(4, 0)],
            vz: state[(5, 0)],
            dt,
            frame,
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

    /// Returns the radius vector of this State in [km, km, km]
    pub fn radius(&self) -> Vector3<f64> {
        Vector3::new(self.x, self.y, self.z)
    }

    /// Returns the radius vector of this State in [km, km, km]
    pub fn velocity(&self) -> Vector3<f64> {
        Vector3::new(self.vx, self.vy, self.vz)
    }

    /// Returns this state as a Cartesian Vector6 in [km, km, km, km/s, km/s, km/s]
    ///
    /// Note that the time is **not** returned in the vector.
    pub fn to_cartesian_vec(&self) -> Vector6<f64> {
        Vector6::new(self.x, self.y, self.z, self.vx, self.vy, self.vz)
    }

    /// Returns the distancein kilometers between this state and another state.
    /// Will **panic** is the frames are different
    pub fn distance_to(&self, other: &'a State) -> f64 {
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

    pub fn to_exb_state(&self) -> XBState {
        XBState {
            epoch: Some(XBEpoch {
                ts: i32::from(TimeSystem::Tai),
                repr: i32::from(TimeRepr::DaysJ1900),
                value: self.dt.as_mjd_tai_days(),
            }),
            position: Some(XBVector {
                x: self.x,
                y: self.y,
                z: self.z,
                is_zero: false,
                unit: i32::from(Unit::Km),
            }),
            velocity: Some(XBVector {
                x: self.vx,
                y: self.vy,
                z: self.vz,
                is_zero: false,
                unit: i32::from(Unit::KmS),
            }),
            covariance: None,
            covariance_exponent: 0.0,
            parameters: HashMap::new(),
        }
    }

    /// Returns the unit vector in the direction of the state radius
    pub fn r_hat(&self) -> Vector3<f64> {
        self.radius() / self.rmag()
    }

    /// Returns the unit vector in the direction of the state velocity
    pub fn v_hat(&self) -> Vector3<f64> {
        perpv(&self.velocity(), &self.r_hat()) / self.rmag()
    }
}

impl<'a> PartialEq for State<'a> {
    /// Two states are equal if their position are equal within one centimeter and their velocities within one centimeter per second.
    fn eq(&self, other: &'a State<'a>) -> bool {
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

impl<'a> Add for State<'a> {
    type Output = State<'a>;

    /// Add one state from another. Frame must be manually changed if needed.
    fn add(self, other: State<'a>) -> State<'a> {
        State {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
            vx: self.vx + other.vx,
            vy: self.vy + other.vy,
            vz: self.vz + other.vz,
            dt: self.dt,
            frame: self.frame,
        }
    }
}

impl<'a> Sub for State<'a> {
    type Output = State<'a>;

    /// Subtract one state from another
    fn sub(self, other: State<'a>) -> State<'a> {
        State {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
            vx: self.vx - other.vx,
            vy: self.vy - other.vy,
            vz: self.vz - other.vz,
            dt: self.dt,
            frame: self.frame,
        }
    }
}

impl<'a> Neg for State<'a> {
    type Output = State<'a>;

    /// Subtract one state from another
    fn neg(self) -> Self::Output {
        State {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            vx: -self.vx,
            vy: -self.vy,
            vz: -self.vz,
            dt: self.dt,
            frame: self.frame,
        }
    }
}

impl<'a> Add for &State<'a> {
    type Output = State<'a>;

    /// Add one state from another. Frame must be manually changed if needed.
    fn add(self, other: &State<'a>) -> State<'a> {
        State {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
            vx: self.vx + other.vx,
            vy: self.vy + other.vy,
            vz: self.vz + other.vz,
            dt: self.dt,
            frame: self.frame,
        }
    }
}

impl<'a> Sub for &State<'a> {
    type Output = State<'a>;

    /// Subtract one state from another
    fn sub(self, other: &State<'a>) -> State<'a> {
        State {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
            vx: self.vx - other.vx,
            vy: self.vy - other.vy,
            vz: self.vz - other.vz,
            dt: self.dt,
            frame: self.frame,
        }
    }
}

impl<'a> Neg for &State<'a> {
    type Output = State<'a>;

    /// Subtract one state from another
    fn neg(self) -> Self::Output {
        State {
            x: -self.x,
            y: -self.y,
            z: -self.z,
            vx: -self.vx,
            vy: -self.vy,
            vz: -self.vz,
            dt: self.dt,
            frame: self.frame,
        }
    }
}

impl<'a> Serialize for State<'a> {
    /// NOTE: This is not part of unit testing because there is no deseralization of State (yet)
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut state = serializer.serialize_struct("State", 7)?;
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

impl<'a> State<'a> {
    /// Creates a new State around the provided Celestial or Geoid frame from the Keplerian orbital elements.
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
            FrameKind::Geoid { gm, .. } | FrameKind::Celestial { gm } => {
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
                State {
                    x,
                    y,
                    z,
                    vx,
                    vy,
                    vz,
                    dt,
                    frame,
                }
            }
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    /// Creates a new State around the provided CelestialBody from the borrowed state vector
    ///
    /// The state vector **must** be sma, ecc, inc, raan, aop, ta. This function is a shortcut to `cartesian`
    /// and as such it has the same unit requirements.
    pub fn keplerian_vec(state: &Vector6<f64>, dt: Epoch, frame: Frame) -> Self {
        match frame {
            FrameKind::Geoid { .. } | FrameKind::Celestial { .. } => Self::keplerian(
                state[(0, 0)],
                state[(1, 0)],
                state[(2, 0)],
                state[(3, 0)],
                state[(4, 0)],
                state[(5, 0)],
                dt,
                frame,
            ),
            _ => panic!("Frame is not Celestial or Geoid in kind"),
        }
    }

    /// Creates a new State from the geodetic latitude (φ), longitude (λ) and height with respect to Earth's ellipsoid.
    ///
    /// **Units:** degrees, degrees, km
    /// NOTE: This computation differs from the spherical coordinates because we consider the flattening of Earth.
    /// Reference: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016
    pub fn from_geodesic(
        latitude: f64,
        longitude: f64,
        height: f64,
        dt: Epoch,
        frame: Frame,
    ) -> Self {
        match frame {
            FrameKind::Geoid {
                _gm,
                flattening,
                _eq_rad,
                semi_major_radius,
            } => {
                let e2 = 2.0 * flattening - flattening.powi(2);
                let (sin_long, cos_long) = longitude.to_radians().sin_cos();
                let (sin_lat, cos_lat) = latitude.to_radians().sin_cos();
                // page 144
                let c_earth = semi_major_radius / ((1.0 - e2 * sin_lat.powi(2)).sqrt());
                let s_earth = (semi_major_radius * (1.0 - flattening).powi(2))
                    / ((1.0 - e2 * sin_lat.powi(2)).sqrt());
                let ri = (c_earth + height) * cos_lat * cos_long;
                let rj = (c_earth + height) * cos_lat * sin_long;
                let rk = (s_earth + height) * sin_lat;
                let radius = Vector3::new(ri, rj, rk);
                let velocity = Vector3::new(0.0, 0.0, 7.292_115_146_706_4e-5).cross(&radius);
                OrbitState::cartesian(
                    radius[(0, 0)],
                    radius[(1, 0)],
                    radius[(2, 0)],
                    velocity[(0, 0)],
                    velocity[(1, 0)],
                    velocity[(2, 0)],
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
        self.hvec()[(0, 0)]
    }

    /// Returns the orbital momentum value on the Y axis
    pub fn hy(&self) -> f64 {
        self.hvec()[(1, 0)]
    }

    /// Returns the orbital momentum value on the Z axis
    pub fn hz(&self) -> f64 {
        self.hvec()[(2, 0)]
    }

    /// Returns the norm of the orbital momentum
    pub fn hmag(&self) -> f64 {
        self.hvec().norm()
    }

    /// Returns the specific mechanical energy
    pub fn energy(&self) -> f64 {
        self.vmag().powi(2) / 2.0 - self.frame.gm / self.rmag()
    }

    /// Returns the semi-major axis in km
    pub fn sma(&self) -> f64 {
        -self.frame.gm / (2.0 * self.energy())
    }

    /// Returns the period in seconds
    pub fn period(&self) -> f64 {
        2.0 * PI * (self.sma().powi(3) / self.frame.gm).sqrt()
    }

    /// Returns the eccentricity vector (no unit)
    pub fn evec(&self) -> Vector3<f64> {
        let r = self.radius();
        let v = self.velocity();
        ((v.norm().powi(2) - self.frame.gm / r.norm()) * r - (r.dot(&v)) * v) / self.frame.gm
    }

    /// Returns the eccentricity (no unit)
    pub fn ecc(&self) -> f64 {
        self.evec().norm()
    }

    /// Returns the inclination in degrees
    pub fn inc(&self) -> f64 {
        (self.hvec()[(2, 0)] / self.hmag()).acos().to_degrees()
    }

    /// Returns the argument of periapsis in degrees
    pub fn aop(&self) -> f64 {
        let n = Vector3::new(0.0, 0.0, 1.0).cross(&self.hvec());
        let aop = (n.dot(&self.evec()) / (n.norm() * self.ecc())).acos();
        if aop.is_nan() {
            warn!("AoP is NaN");
            0.0
        } else if self.evec()[(2, 0)] < 0.0 {
            (2.0 * PI - aop).to_degrees()
        } else {
            aop.to_degrees()
        }
    }

    /// Returns the right ascension of ther ascending node in degrees
    pub fn raan(&self) -> f64 {
        let n = Vector3::new(0.0, 0.0, 1.0).cross(&self.hvec());
        let raan = (n[(0, 0)] / n.norm()).acos();
        if raan.is_nan() {
            warn!("RAAN is NaN");
            0.0
        } else if n[(1, 0)] < 0.0 {
            (2.0 * PI - raan).to_degrees()
        } else {
            raan.to_degrees()
        }
    }

    /// Returns the true anomaly in degrees between 0 and 360.0
    ///
    /// NOTE: This function will emit a warning stating that the TA should be avoided if in a very near circular orbit
    /// Code from https://github.com/ChristopherRabotin/GMAT/blob/80bde040e12946a61dae90d9fc3538f16df34190/src/gmatutil/util/StateConversionUtil.cpp#L6835
    pub fn ta(&self) -> f64 {
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

    /// Returns the true longitude in degrees
    pub fn tlong(&self) -> f64 {
        // Angles already in degrees
        between_0_360(self.aop() + self.raan() + self.ta())
    }

    /// Returns the argument of latitude in degrees
    ///
    /// NOTE: If the orbit is near circular, the AoL will be computed from the true longitude
    /// instead of relying on the ill-defined true anomaly.
    pub fn aol(&self) -> f64 {
        between_0_360(if self.ecc() < ECC_EPSILON {
            self.tlong() - self.raan()
        } else {
            self.aop() + self.ta()
        })
    }

    /// Returns the radius of periapsis (or perigee around Earth), in kilometers.
    pub fn periapsis(&self) -> f64 {
        self.sma() * (1.0 - self.ecc())
    }

    /// Returns the radius of apoapsis (or apogee around Earth), in kilometers.
    pub fn apoapsis(&self) -> f64 {
        self.sma() * (1.0 + self.ecc())
    }

    /// Returns the eccentric anomaly in degrees
    ///
    /// This is a conversion from GMAT's StateConversionUtil::TrueToEccentricAnomaly
    pub fn ea(&self) -> f64 {
        let (sin_ta, cos_ta) = self.ta().to_radians().sin_cos();
        let ecc_cos_ta = self.ecc() * cos_ta;
        let sin_ea = ((1.0 - self.ecc().powi(2)).sqrt() * sin_ta) / (1.0 + ecc_cos_ta);
        let cos_ea = (self.ecc() + cos_ta) / (1.0 + ecc_cos_ta);
        // The atan2 function is a bit confusing: https://doc.rust-lang.org/std/primitive.f64.html#method.atan2 .
        sin_ea.atan2(cos_ea).to_degrees()
    }

    /// Returns the mean anomaly in degrees
    ///
    /// This is a conversion from GMAT's StateConversionUtil::TrueToMeanAnomaly
    pub fn ma(&self) -> f64 {
        if self.ecc() < 1.0 {
            between_0_360(
                (self.ea().to_radians() - self.ecc() * self.ea().to_radians().sin()).to_degrees(),
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

    /// Returns the semi parameter (or semilatus rectum)
    pub fn semi_parameter(&self) -> f64 {
        self.sma() * (1.0 - self.ecc().powi(2))
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
        between_0_360(self.y.atan2(self.x).to_degrees())
    }

    /// Returns the geodetic latitude (φ) in degrees. Value is between -180 and +180 degrees.
    ///
    /// Reference: Vallado, 4th Ed., Algorithm 12 page 172.
    pub fn geodetic_latitude(&self) -> f64 {
        let eps = 1e-12;
        let max_attempts = 20;
        let mut attempt_no = 0;
        let r_delta = (self.x.powi(2) + self.y.powi(2)).sqrt();
        let mut latitude = (self.z / self.rmag()).asin();
        let e2 = self.frame.flattening * (2.0 - self.frame.flattening);
        loop {
            attempt_no += 1;
            let c_earth =
                self.frame.semi_major_radius / ((1.0 - e2 * (latitude).sin().powi(2)).sqrt());
            let new_latitude = (self.z + c_earth * e2 * (latitude).sin()).atan2(r_delta);
            if (latitude - new_latitude).abs() < eps {
                return between_pm_180(new_latitude.to_degrees());
            } else if attempt_no >= max_attempts {
                println!(
                    "geodetic latitude failed to converge -- error = {}",
                    (latitude - new_latitude).abs()
                );
                return between_pm_180(new_latitude.to_degrees());
            }
            latitude = new_latitude;
        }
    }

    /// Returns the geodetic height in km.
    ///
    /// Reference: Vallado, 4th Ed., Algorithm 12 page 172.
    pub fn geodetic_height(&self) -> f64 {
        let e2 = self.frame.flattening * (2.0 - self.frame.flattening);
        let latitude = self.geodetic_latitude().to_radians();
        let sin_lat = latitude.sin();
        if (latitude - 1.0).abs() < 0.1 {
            // We are near poles, let's use another formulation.
            let s_earth = (self.frame.semi_major_radius * (1.0 - self.frame.flattening).powi(2))
                / ((1.0 - e2 * sin_lat.powi(2)).sqrt());
            self.z / latitude.sin() - s_earth
        } else {
            let c_earth = self.frame.semi_major_radius / ((1.0 - e2 * sin_lat.powi(2)).sqrt());
            let r_delta = (self.x.powi(2) + self.y.powi(2)).sqrt();
            r_delta / latitude.cos() - c_earth
        }
    }

    /// Returns the direct cosine rotation matrix to convert to this inertial state.
    pub fn dcm_to_inertial(&self, from: FrameKind) -> Matrix3<f64> {
        match from {
            FrameKind::RIC => {
                r3(-self.raan().to_radians())
                    * r1(-self.inc().to_radians())
                    * r3(-self.aol().to_radians())
            }
            FrameKind::VNC => {
                let v = self.velocity() / self.vmag();
                let n = self.hvec() / self.hmag();
                let c = v.cross(&n);
                Matrix3::new(v[0], v[1], v[2], n[0], n[1], n[2], c[0], c[1], c[2]).transpose()
            }
            FrameKind::RCN => {
                let r = self.radius() / self.rmag();
                let n = self.hvec() / self.hmag();
                let c = n.cross(&r);
                Matrix3::new(r[0], r[1], r[2], c[0], c[1], c[2], n[0], n[1], n[2]).transpose()
            }
            _ => panic!("did not provide a local frame"),
        }
    }
}

impl<'a> fmt::Display for State<'a> {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "[{}] {}\tposition = [{:.6}, {:.6}, {:.6}] km\tvelocity = [{:.6}, {:.6}, {:.6}] km/s",
            self.frame,
            self.dt.as_gregorian_utc_tai(),
            self.x,
            self.y,
            self.z,
            self.vx,
            self.vy,
            self.vz
        )
    }
}

impl<'a> fmt::LowerExp for State<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "[{}] {}\tposition = [{:e}, {:e}, {:e}] km\tvelocity = [{:e}, {:e}, {:e}] km/s",
            self.frame,
            self.dt.as_gregorian_utc_str(),
            self.x,
            self.y,
            self.z,
            self.vx,
            self.vy,
            self.vz
        )
    }
}

impl<'a> fmt::Octal for State<'a> {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "[{}] {}\tsma = {:.6} km\tecc = {:.6}\tinc = {:.6} deg\traan = {:.6} deg\taop = {:.6} deg\tta = {:.6} deg",
            self.frame,
            self.dt.as_gregorian_utc_tai(),
            self.sma(),
            self.ecc(),
            self.inc(),
            self.raan(),
            self.aop(),
            self.ta()
        )
    }
}

// TODO: REMOVE ME
/// An orbit is simply a State typed around a Geoid.
pub type OrbitState<'a> = State<'a>;
