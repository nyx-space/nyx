use super::na::{Vector3, Vector6};
use super::CelestialBody;

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

    /// Returns this state as a Cartesian Vector6
    pub fn to_cartsian_vec(self) -> Vector6<f64> {
        Vector6::new(self.x, self.y, self.z, self.vx, self.vy, self.vz)
    }

    /// Returns the specific mechanical energy
    pub fn energy(self) -> f64 {
        let r = (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt();
        let v2 = self.vx.powi(2) + self.vy.powi(2) + self.vz.powi(2);
        v2 / 2.0 - self.gm / r
    }

    /// Returns the semi-major axis in kilometers
    pub fn sma(self) -> f64 {
        -self.gm / (2.0 * self.energy())
    }

    /// Returns the period in seconds
    pub fn period(self) -> f64 {
        use std::f64::consts::PI;
        println!("a = {} ==> {}", self.sma(), self.sma().powi(3));
        2.0 * PI * (self.sma().powi(3) / self.gm).sqrt()
    }

    pub fn ecc(self) -> f64 {
        // Ugh, this is a bit of a heavy computation, so I should probably find a way to cache this? Or does Rust handle that by itself
        // if the state is immutable?
        let r = Vector3::new(self.x, self.y, self.z);
        let v = Vector3::new(self.vx, self.vy, self.vz);
        let e = ((v.norm().powi(2) - self.gm / r.norm()) * v - (r.dot(&v)) * v) / self.gm;
        e.norm()
    }
}
