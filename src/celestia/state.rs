use super::na::{Vector3, Vector6};
use super::CelestialBody;
use super::pretty_env_logger;

use std::f64::consts::PI;

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
        pretty_env_logger::init_custom_env("NYX_LOG");
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
        pretty_env_logger::init_custom_env("NYX_LOG");
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
        ((v.norm().powi(2) - self.gm / r.norm()) * v - (r.dot(&v)) * v) / self.gm
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
        let mut aop = (n.dot(&self.evec()) / (n.norm() * self.ecc())).acos();
        if aop.is_nan() {
            warn!("AoP is NaN");
            return 0.0;
        }
        if self.evec()[(2, 0)] < 0.0 {
            aop = 2.0 * PI - aop
        }
        aop.to_degrees()
    }
}
