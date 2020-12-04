use super::thrustctrl::ThrustControl;
use crate::celestia::Orbit;
use crate::dimensions::Vector3;
use std::fmt;
use std::sync::Arc;

// TODO: Change to a trait to enable variable thrust and isp
#[derive(Copy, Clone, Debug)]
pub struct Thruster {
    pub thrust: f64,
    pub isp: f64,
}

const NORM_ERR: f64 = 1e-12;
const STD_GRAVITY: f64 = 9.80665; // From NIST special publication 330, 2008 edition

/// Propulsion allows returning a propulsion FORCE (not acceleration).
#[derive(Clone)]
pub struct Propulsion {
    /// Will eventually support tank depletion. TODO: Implement thrusting from different thrusters independently
    pub thruster: Thruster,
    /// Set to true to decrement the fuel mass
    pub decrement_mass: bool,
    pub ctrl: Arc<dyn ThrustControl>,
}

impl fmt::Debug for Propulsion {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Prop with {:?} thruster", self.thruster)
    }
}

impl Propulsion {
    pub fn new(ctrl: Arc<dyn ThrustControl>, thruster: Thruster, decrement_mass: bool) -> Self {
        Self {
            thruster,
            decrement_mass,
            ctrl,
        }
    }

    pub fn set_state(&mut self, new_orbit: &Orbit) {
        // TODO: Fix this.
        // self.ctrl.next(new_orbit);
    }

    /// Returns the thrust and the fuel usage
    pub fn eom(&self, osc: &Orbit) -> (Vector3<f64>, f64) {
        let thrust_power = self.ctrl.throttle(osc);
        if thrust_power > 0.0 {
            // Thrust arc
            let thrust_inertial = self.ctrl.direction(osc);
            if (thrust_inertial.norm() - 1.0).abs() > NORM_ERR {
                panic!(
                    "Invalid control low: thrust orientation is not a unit vector: norm = {}\n{}",
                    thrust_inertial.norm(),
                    thrust_inertial
                );
            }

            // Compute the thrust in Newtons and Isp
            let mut total_thrust = 0.0;
            let mut fuel_usage = 0.0;

            total_thrust += thrust_power * self.thruster.thrust;
            fuel_usage += thrust_power * self.thruster.thrust / (self.thruster.isp * STD_GRAVITY);
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
