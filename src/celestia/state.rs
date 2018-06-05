extern crate hifitime;
use self::hifitime::instant::{Duration, Instant};
use self::hifitime::TimeSystem;
use self::hifitime::{datetime, julian};

use super::na::{Vector3, Vector6};
use super::{CelestialBody, CoordinateFrame, EARTH, ECEF, ECI};
use utils::between_0_360;

use std::f64::consts::PI;
use std::f64::EPSILON;
use std::fmt;

/// If an orbit has an eccentricity below the following value, it is considered circular (only affects warning messages)
pub const ECC_EPSILON: f64 = 1e-4;

// A warning will be logged if a division operation is planned with a value smaller than this following value.
const ZERO_DIV_TOL: f64 = 1e-15;

/// State defines an orbital state parameterized  by a `CelestialBody`.
///
/// Unless noted otherwise, algorithms are from GMAT 2016a [StateConversionUtil.cpp](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/StateConversionUtil.cpp).
/// Regardless of the constructor used, this struct stores all the state information in Cartesian coordinates
/// as these are always non singular.
/// _Note:_ although not yet supported, this struct may change once True of Date or other nutation frames
/// are added to the toolkit.
#[derive(Copy, Clone, Debug)]
pub struct State<F>
where
    F: CoordinateFrame,
{
    gm: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    /// The date time is stored as a a hifitime::instant::Instant, which is the NTP representation of time.
    pub dt: Instant,
    /// The frame will later allow for coordinate frame transformations.
    pub frame: F,
}

impl<F: CoordinateFrame> PartialEq for State<F> {
    /// Two states are equal if their position are equal within one centimeter and their velocities within one centimeter per second.
    /// For time equality, we're relying on the high fidelity time computation of `hifitime` provided through the `Instant` representation.
    fn eq(&self, other: &State<F>) -> bool
    where
        F: CoordinateFrame,
    {
        let distance_tol = 1e-5; // centimeter
        let velocity_tol = 1e-5; // centimeter per second
        self.dt == other.dt
            && (self.gm - other.gm).abs() < 1e-4
            && (self.x - other.x).abs() < distance_tol
            && (self.y - other.y).abs() < distance_tol
            && (self.z - other.z).abs() < distance_tol
            && (self.vx - other.vx).abs() < velocity_tol
            && (self.vy - other.vy).abs() < velocity_tol
            && (self.vz - other.vz).abs() < velocity_tol
    }
}

impl<F> State<F>
where
    F: CoordinateFrame,
{
    /// Creates a new State around the provided CelestialBody
    ///
    /// **Units:** km, km, km, km/s, km/s, km/s
    pub fn from_cartesian<B: CelestialBody, T: TimeSystem>(
        x: f64,
        y: f64,
        z: f64,
        vx: f64,
        vy: f64,
        vz: f64,
        dt: T,
        frame: F,
    ) -> State<F>
    where
        F: CoordinateFrame,
    {
        State {
            gm: B::gm(),
            x,
            y,
            z,
            vx,
            vy,
            vz,
            dt: dt.into_instant(),
            frame,
        }
    }

    /// Creates a new State around the provided CelestialBody from the borrowed state vector
    ///
    /// The state vector **must** be x, y, z, dx, dy, dz. This function is a shortcut to `from_cartesian`
    /// and as such it has the same unit requirements.
    pub fn from_cartesian_vec<B: CelestialBody, T: TimeSystem>(state: &Vector6<f64>, dt: T, frame: F) -> State<F>
    where
        F: CoordinateFrame,
    {
        State {
            gm: B::gm(),
            x: state[(0, 0)],
            y: state[(1, 0)],
            z: state[(2, 0)],
            vx: state[(3, 0)],
            vy: state[(4, 0)],
            vz: state[(5, 0)],
            dt: dt.into_instant(),
            frame,
        }
    }

    /// Creates a new State around the provided CelestialBody from the Keplerian orbital elements.
    ///
    /// **Units:** km, none, degrees, degrees, degrees, degrees
    ///
    /// WARNING: This function will panic if the singularities in the conversion are expected.
    /// NOTE: The state is defined in Cartesian coordinates as they are non-singular. This causes rounding
    /// errors when creating a state from its Keplerian orbital elements (cf. the state tests).
    /// One should expect these errors to be on the order of 1e-12.
    pub fn from_keplerian<B: CelestialBody, T: TimeSystem>(
        sma: f64,
        ecc: f64,
        inc: f64,
        raan: f64,
        aop: f64,
        ta: f64,
        dt: T,
        frame: F,
    ) -> State<F>
    where
        F: CoordinateFrame,
    {
        if B::gm().abs() < ZERO_DIV_TOL {
            warn!(
                "GM very low ({}): expect math errors in Keplerian to Cartesian conversion",
                B::gm()
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
                panic!("true anomaly value ({}) physically impossible for a hyperbolic orbit", ta);
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

        let sqrt_gm_p = (B::gm() / p).sqrt();
        let cos_ta_ecc = ta.cos() + ecc;
        let sin_ta = ta.sin();

        let vx = sqrt_gm_p * cos_ta_ecc * (-sin_aop * cos_raan - cos_inc * sin_raan * cos_aop)
            - sqrt_gm_p * sin_ta * (cos_aop * cos_raan - cos_inc * sin_raan * sin_aop);
        let vy = sqrt_gm_p * cos_ta_ecc * (-sin_aop * sin_raan + cos_inc * cos_raan * cos_aop)
            - sqrt_gm_p * sin_ta * (cos_aop * sin_raan + cos_inc * cos_raan * sin_aop);

        let vz = sqrt_gm_p * (cos_ta_ecc * sin_inc * cos_aop - sin_ta * sin_inc * sin_aop);

        State {
            gm: B::gm(),
            x,
            y,
            z,
            vx,
            vy,
            vz,
            dt: dt.into_instant(),
            frame,
        }
    }

    /// Creates a new State around the provided CelestialBody from the borrowed state vector
    ///
    /// The state vector **must** be sma, ecc, inc, raan, aop, ta. This function is a shortcut to `from_cartesian`
    /// and as such it has the same unit requirements.
    pub fn from_keplerian_vec<B: CelestialBody, T: TimeSystem>(state: &Vector6<f64>, dt: T, frame: F) -> State<F>
    where
        F: CoordinateFrame,
    {
        Self::from_keplerian::<B, T>(
            state[(0, 0)],
            state[(1, 0)],
            state[(2, 0)],
            state[(3, 0)],
            state[(4, 0)],
            state[(5, 0)],
            dt,
            frame,
        )
    }

    /// Returns the date time as its Modified Julian time representation
    pub fn dt_as_modified_julian(self) -> julian::ModifiedJulian {
        julian::ModifiedJulian::from_instant(self.dt)
    }

    /// Returns the date time as its UTC representation
    pub fn dt_as_utc(self) -> datetime::Datetime {
        datetime::Datetime::from_instant(self.dt)
    }

    /// Converts this state to a state in the destination frame (`other`).
    ///
    /// Reference: Vallado, 4th Ed., page 167, equation 3-27.
    /// Note that we compute the derivative of the DCM by taking the difference between
    /// said DCMs at a 0.1 second interval.
    pub fn in_frame<G: CoordinateFrame>(self, other: G) -> State<G>
    where
        F: CoordinateFrame,
    {
        // if self.frame.Center != other.Center {
        //     unimplemented!("transformations between differently centered frames depends on #4");
        // }
        let radius = Vector3::new(self.x, self.y, self.z);
        let velocity = Vector3::new(self.vx, self.vy, self.vz);
        // Convert from our frame to the destination frame via their inertial frames.
        let dest_radius = G::from_inertial(self.dt) * F::to_inertial(self.dt) * radius;
        // We need to compute the derivate of the DCM for the velocity.
        let h = Duration::from_millis(100);
        let frame_dt = (F::from_inertial(self.dt + h) - F::from_inertial(self.dt)) / 0.1;
        let other_dt = (G::from_inertial(self.dt + h) - G::from_inertial(self.dt)) / 0.1;
        let dest_velocity = G::from_inertial(self.dt) * F::to_inertial(self.dt) * velocity + other_dt * frame_dt * radius;
        State {
            gm: self.gm,
            x: dest_radius[(0, 0)],
            y: dest_radius[(1, 0)],
            z: dest_radius[(2, 0)],
            vx: dest_velocity[(0, 0)],
            vy: dest_velocity[(1, 0)],
            vz: dest_velocity[(2, 0)],
            dt: self.dt,
            frame: other,
        }
    }

    /// Returns this state as a Cartesian Vector6 in [km, km, km, km/s, km/s, km/s]
    ///
    /// Note that the time is **not** returned in the vector.
    pub fn to_cartesian_vec(self) -> Vector6<f64> {
        Vector6::new(self.x, self.y, self.z, self.vx, self.vy, self.vz)
    }

    /// Returns this state as a Keplerian Vector6 in [km, none, degrees, degrees, degrees, degrees]
    ///
    /// Note that the time is **not** returned in the vector.
    pub fn to_keplerian_vec(self) -> Vector6<f64> {
        Vector6::new(self.sma(), self.ecc(), self.inc(), self.raan(), self.aop(), self.ta())
    }

    /// Returns the magnitude of the radius vector in km
    pub fn rmag(self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }

    /// Returns the magnitude of the velocity vector in km/s
    pub fn vmag(self) -> f64 {
        (self.vx.powi(2) + self.vy.powi(2) + self.vz.powi(2)).sqrt()
    }

    /// Returns the orbital momentum vector
    pub fn hvec(self) -> Vector3<f64> {
        let r = Vector3::new(self.x, self.y, self.z);
        let v = Vector3::new(self.vx, self.vy, self.vz);
        r.cross(&v)
    }

    /// Returns the orbital momentum value on the X axis
    pub fn hx(self) -> f64 {
        self.hvec()[(0, 0)]
    }

    /// Returns the orbital momentum value on the Y axis
    pub fn hy(self) -> f64 {
        self.hvec()[(1, 0)]
    }

    /// Returns the orbital momentum value on the Z axis
    pub fn hz(self) -> f64 {
        self.hvec()[(2, 0)]
    }

    /// Returns the norm of the orbital momentum
    pub fn hmag(self) -> f64 {
        self.hvec().norm()
    }

    /// Returns the specific mechanical energy
    pub fn energy(self) -> f64 {
        self.vmag().powi(2) / 2.0 - self.gm / self.rmag()
    }

    /// Returns the semi-major axis in km
    pub fn sma(self) -> f64 {
        -self.gm / (2.0 * self.energy())
    }

    /// Returns the period in seconds
    pub fn period(self) -> f64 {
        2.0 * PI * (self.sma().powi(3) / self.gm).sqrt()
    }

    /// Returns the eccentricity vector (no unit)
    pub fn evec(self) -> Vector3<f64> {
        let r = Vector3::new(self.x, self.y, self.z);
        let v = Vector3::new(self.vx, self.vy, self.vz);
        ((v.norm().powi(2) - self.gm / r.norm()) * r - (r.dot(&v)) * v) / self.gm
    }

    /// Returns the eccentricity (no unit)
    pub fn ecc(self) -> f64 {
        self.evec().norm()
    }

    /// Returns the inclination in degrees
    pub fn inc(self) -> f64 {
        (self.hvec()[(2, 0)] / self.hmag()).acos().to_degrees()
    }

    /// Returns the argument of periapsis in degrees
    pub fn aop(self) -> f64 {
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
    pub fn raan(self) -> f64 {
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

    /// Returns the true anomaly in degrees
    ///
    /// NOTE: This function will emit a warning stating that the TA should be avoided if in a very near circular orbit
    pub fn ta(self) -> f64 {
        if self.ecc() < ECC_EPSILON {
            warn!("true anomaly ill-defined (eccentricity too low, e = {})", self.ecc());
        }
        let cos_nu = self.evec().dot(&Vector3::new(self.x, self.y, self.z)) / (self.ecc() * self.rmag());
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
            } else if self.hmag() < 0.0 {
                (2.0 * PI - ta).to_degrees()
            } else {
                ta.to_degrees()
            }
        }
    }

    /// Returns the true longitude in degrees
    pub fn tlong(self) -> f64 {
        // Angles already in degrees
        between_0_360(self.aop() + self.raan() + self.ta())
    }

    /// Returns the argument of latitude in degrees
    ///
    /// NOTE: If the orbit is near circular, the AoL will be computed from the true longitude
    /// instead of relying on the ill-defined true anomaly.
    pub fn aol(self) -> f64 {
        between_0_360(if self.ecc() < ECC_EPSILON {
            self.tlong() - self.raan()
        } else {
            self.aop() + self.ta()
        })
    }

    /// Returns the radius of periapsis (or perigee around Earth), in kilometers.
    pub fn periapsis(self) -> f64 {
        self.sma() * (1.0 - self.ecc())
    }

    /// Returns the radius of apoapsis (or apogee around Earth), in kilometers.
    pub fn apoapsis(self) -> f64 {
        self.sma() * (1.0 + self.ecc())
    }

    /// Returns the eccentric anomaly in degrees
    ///
    /// This is a conversion from GMAT's StateConversionUtil::TrueToEccentricAnomaly
    pub fn ea(self) -> f64 {
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
    pub fn ma(self) -> f64 {
        if self.ecc() < 1.0 {
            between_0_360((self.ea().to_radians() - self.ecc() * self.ea().to_radians().sin()).to_degrees())
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
    pub fn semi_parameter(self) -> f64 {
        self.sma() * (1.0 - self.ecc().powi(2))
    }

    /// Returns whether this state satisfies the requirement to compute the Mean Brouwer Short orbital
    /// element set.
    ///
    /// This is a conversion from GMAT's StateConversionUtil::CartesianToBrouwerMeanShort.
    /// The details are at the log level `info`.
    /// NOTE: Mean Brouwer Short are only defined around Earth. However, `nyx` does *not* check the
    /// main celestial body around which the state is defined (GMAT does perform this verification).
    pub fn is_brouwer_short_valid(self) -> bool {
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
}

impl<F> fmt::Display for State<F>
where
    F: CoordinateFrame,
{
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "[{}] {}\tposition = [{:.6}, {:.6}, {:.6}] km\tvelocity = [{:.6}, {:.6}, {:.6}] km/s",
            self.frame,
            self.dt_as_utc(),
            self.x,
            self.y,
            self.z,
            self.vx,
            self.vy,
            self.vz
        )
    }
}

impl<F> fmt::Octal for State<F>
where
    F: CoordinateFrame,
{
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "[{}] {}\tsma = {:.6} km\tecc = {:.6}\tinc = {:.6} deg\traan = {:.6} deg\taop = {:.6} deg\tta = {:.6} deg",
            self.frame,
            self.dt_as_utc(),
            self.sma(),
            self.ecc(),
            self.inc(),
            self.raan(),
            self.aop(),
            self.ta()
        )
    }
}

impl State<ECI> {
    /// Creates a new State around Earth in the ECI frame
    ///
    /// **Units:** km, km, km, km/s, km/s, km/s
    pub fn from_cartesian_eci<T: TimeSystem>(x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64, dt: T) -> State<ECI> {
        State::from_cartesian::<EARTH, T>(x, y, z, vx, vy, vz, dt, ECI {})
    }

    /// Creates a new State around Earth in the ECI frame from the Keplerian orbital elements.
    ///
    /// **Units:** km, none, degrees, degrees, degrees, degrees
    pub fn from_keplerian_eci<T: TimeSystem>(x: f64, y: f64, z: f64, vx: f64, vy: f64, vz: f64, dt: T) -> State<ECI> {
        State::from_keplerian::<EARTH, T>(x, y, z, vx, vy, vz, dt, ECI {})
    }
}

impl State<ECEF> {
    /// Creates a new State from the geodetic latitude (φ), longitude (λ) and height with respect to Earth's ellipsoid.
    ///
    /// **Units:** degrees, degrees, km
    /// NOTE: This computation differs from the spherical coordinates because we consider the flatenning of Earth.
    /// Reference: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016
    pub fn from_geodesic<T: TimeSystem>(latitude: f64, longitude: f64, height: f64, dt: T) -> State<ECEF> {
        let (sin_long, cos_long) = longitude.to_radians().sin_cos();
        let (sin_lat, cos_lat) = latitude.to_radians().sin_cos();
        let ri = (EARTH::semi_major_radius() + height) * cos_lat * cos_long;
        let rj = (EARTH::semi_major_radius() + height) * cos_lat * sin_long;
        let e2 = 2.0 * EARTH::flatenning() - EARTH::flatenning().powi(2);
        let rk = (EARTH::semi_major_radius() * (1.0 - e2) + height) * sin_lat;
        let radius = Vector3::new(ri, rj, rk);
        let velocity = Vector3::new(0.0, 0.0, EARTH::rotation_rate()).cross(&radius);
        State::from_cartesian::<EARTH, T>(
            radius[(0, 0)],
            radius[(1, 0)],
            radius[(2, 0)],
            velocity[(0, 0)],
            velocity[(1, 0)],
            velocity[(2, 0)],
            dt,
            ECEF {},
        )
    }

    /// Creates a new ECEF state at the provided position.
    ///
    /// NOTE: This has the same container as the normal State. Hence, we set the velocity at zero.
    pub fn from_position<T: TimeSystem>(i: f64, j: f64, k: f64, dt: T) -> State<ECEF> {
        State::from_cartesian::<EARTH, T>(i, j, k, 0.0, 0.0, 0.0, dt, ECEF {})
    }

    /// Returns the I component of this ECEF frame
    pub fn ri(&self) -> f64 {
        self.x
    }

    /// Returns the J component of this ECEF frame
    pub fn rj(&self) -> f64 {
        self.y
    }

    /// Returns the K component of this ECEF frame
    pub fn rk(&self) -> f64 {
        self.z
    }

    /// Returns the geodetic longitude (λ) in degrees
    ///
    /// Although the reference is not Vallado, the math from Vallado proves to be equivalent.
    /// Reference: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016
    pub fn geodetic_longitude(&self) -> f64 {
        self.y.atan2(self.x).to_degrees()
    }

    /// Returns the geodetic latitude (φ) in degrees
    ///
    /// Reference: Vallado, 4th Ed., Algorithm 12 page 172.
    pub fn geodetic_latitude(&self) -> f64 {
        let eps = 1e-12;
        let max_attempts = 20;
        let mut attempt_no = 0;
        let r_delta = (self.x.powi(2) + self.y.powi(2)).sqrt();
        let mut latitude = (self.z / self.rmag()).asin();
        let e2 = EARTH::flatenning() * (2.0 - EARTH::flatenning());
        loop {
            attempt_no += 1;
            let c_earth = EARTH::semi_major_radius() / ((1.0 - e2 * (latitude).sin().powi(2)).sqrt());
            let new_latitude = (self.z + c_earth * e2 * (latitude).sin()).atan2(r_delta);
            if (latitude - new_latitude).abs() < eps {
                return new_latitude.to_degrees();
            } else if attempt_no >= max_attempts {
                warn!(
                    "geodetic latitude failed to converge -- error = {}",
                    (latitude - new_latitude).abs()
                );
                return new_latitude.to_degrees();
            }
            latitude = new_latitude;
        }
    }

    /// Returns the geodetic height in km.
    ///
    /// Reference: Vallado, 4th Ed., Algorithm 12 page 172.
    pub fn geodetic_height(&self) -> f64 {
        let e2 = EARTH::flatenning() * (2.0 - EARTH::flatenning());
        let latitude = self.geodetic_latitude().to_radians();
        if (latitude - 1.0).abs() < 0.1 {
            // We are near poles, let's use another formulation.
            let s_earth = (EARTH::semi_major_radius() * (1.0 - EARTH::flatenning()).powi(2))
                / ((1.0 - e2 * (latitude).sin().powi(2)).sqrt());
            self.z / latitude.sin() - s_earth
        } else {
            let c_earth = EARTH::semi_major_radius() / ((1.0 - e2 * (latitude).sin().powi(2)).sqrt());
            let r_delta = (self.x.powi(2) + self.y.powi(2)).sqrt();
            r_delta / latitude.cos() - c_earth
        }
    }
}
