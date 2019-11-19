use super::na::Vector3;
use super::thrustctrl::ThrustControl;
use celestia::{Geoid, State};

// TODO: Change to a trait to enable variable thrust and isp
#[derive(Copy, Clone, Debug)]
pub struct Thruster {
    pub thrust: f64,
    pub isp: f64,
}

const NORM_ERR: f64 = 1e-12;
const STD_GRAVITY: f64 = 9.80665; // From NIST special publication 330, 2008 edition

pub struct Propulsion<'a, T: ThrustControl> {
    /// in kg, set in Spacecraft constructors.
    pub dry_mass: f64,
    /// in kg
    pub fuel_mass: f64,
    /// Will eventually support tank depletion
    pub thrusters: Vec<Thruster>,
    /// Set to true to decrement the fuel mass
    pub decrement_mass: bool,
    pub ctrl: &'a mut T,
}

impl<'a, T: ThrustControl> Propulsion<'a, T> {
    pub fn new(
        ctrl: &'a mut T,
        fuel_mass: f64,
        thrusters: Vec<Thruster>,
        decrement_mass: bool,
    ) -> Self {
        Self {
            dry_mass: 0.0,
            fuel_mass,
            thrusters,
            decrement_mass,
            ctrl,
        }
    }

    pub fn set_state(&mut self, new_orbit: &State<Geoid>, new_fuel: f64) {
        self.ctrl.next(new_orbit);
        self.fuel_mass = new_fuel;
        if self.decrement_mass {
            assert!(
                self.fuel_mass >= 0.0,
                "negative fuel mass at {:?}",
                new_orbit.dt
            );
        }
    }

    /// Returns the thrust and the fuel usage
    pub fn eom(&self, osc: &State<Geoid>, mass: f64) -> (Vector3<f64>, f64) {
        let thrust_power = self.ctrl.throttle(osc);
        if thrust_power > 0.0 {
            // Thrust arc
            let thrust_inertial = self.ctrl.direction(osc);
            if (thrust_inertial.norm() - 1.0).abs() > NORM_ERR {
                panic!(
                    "thrust orientation is not a unit vector: norm = {}\n{}",
                    thrust_inertial.norm(),
                    thrust_inertial
                );
            }

            // Compute the thrust in Newtons and Isp
            let mut total_thrust = 0.0;
            let mut fuel_usage = 0.0;
            // TODO: Implement thrusting from different thrusters independently
            for thruster in &self.thrusters {
                total_thrust += thrust_power * thruster.thrust;
                fuel_usage += thrust_power * thruster.thrust / (thruster.isp * STD_GRAVITY);
            }
            // TODO: Try to remove the need for mass here so this can be closer to a ForceModel
            total_thrust /= self.dry_mass + mass;
            total_thrust *= 1e-3; // Convert m/s^-2 to km/s^-2
            (
                thrust_inertial * total_thrust,
                if self.decrement_mass {
                    -fuel_usage
                } else {
                    0.0
                },
            )
        } else {
            (Vector3::zeros(), 0.0)
        }
    }
}
