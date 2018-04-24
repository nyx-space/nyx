use super::na::{Vector3, Vector6};
use super::CelestialBody;

use std::f64::consts::PI;
use std::f64::EPSILON;

/// If an orbit has an eccentricity below the following value, it is considered circular (only affects warning messages)
pub const ECC_EPSILON: f64 = 1e-4;

/// State defines an orbital state in the Celestial Reference frame of the parameterized CelestialBody.
///
/// Regardless of the constructor used, this struct stores all the state information in Cartesian coordinates
/// as these are always non singular.
/// _Note:_ although not yet supported, this struct may change once True of Date or other nutation frames
/// are added to the toolkit.
#[derive(Copy, Clone, Debug)]
pub struct State {
    gm: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
}

impl State {
    /// Creates a new State around the provided CelestialBody
    pub fn from_cartesian<B: CelestialBody>(
        x: f64,
        y: f64,
        z: f64,
        vx: f64,
        vy: f64,
        vz: f64,
    ) -> State {
        State {
            gm: B::gm(),
            x,
            y,
            z,
            vx,
            vy,
            vz,
        }
    }

    /// Creates a new State around the provided CelestialBody from the borrowed state vector
    pub fn from_cartesian_vec<B: CelestialBody>(state: &Vector6<f64>) -> State {
        State {
            gm: B::gm(),
            x: state[(0, 0)],
            y: state[(1, 0)],
            z: state[(2, 0)],
            vx: state[(3, 0)],
            vy: state[(4, 0)],
            vz: state[(5, 0)],
        }
    }

    /// Returns this state as a Cartesian Vector6 in [km, km, km, km/s, km/s, km/s]
    pub fn to_cartsian_vec(self) -> Vector6<f64> {
        Vector6::new(self.x, self.y, self.z, self.vx, self.vy, self.vz)
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
            warn!(
                "true anomaly improperly defined (eccentricity too low, e = {})",
                self.ecc()
            );
        }
        let cos_nu =
            self.evec().dot(&Vector3::new(self.x, self.y, self.z)) / (self.ecc() * self.rmag());
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
        let mut tlong = self.aop() + self.raan() + self.ta(); // Already in degrees
        while tlong < 0.0 {
            tlong += 360.0;
        }
        while tlong > 360.0 {
            tlong -= 360.0;
        }
        tlong
    }

    /// Returns the argument of latitude in degrees
    ///
    /// NOTE: If the orbit is near circular, the AoL will be computed from the true longitude
    /// instead of relying on the improperly defined true anomaly.
    pub fn aol(self) -> f64 {
        let mut aol = if self.ecc() < ECC_EPSILON {
            self.tlong() - self.raan()
        } else {
            self.aop() + self.ta()
        };
        while aol < 0.0 {
            aol += 360.0;
        }
        while aol > 360.0 {
            aol -= 360.0;
        }
        aol
    }
}
