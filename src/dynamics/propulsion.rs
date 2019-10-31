use super::hifitime::Epoch;
use super::na::{VectorN, U6, U7};
use super::thrustctrl::ThrustControl;
use super::Dynamics;

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
    /// Maneuver execution timing
    pub abs_time: Epoch,
    /// Will eventually support tank depletion
    pub thrusters: Vec<Thruster>,
    /// Set to true to decrement the fuel mass
    pub decrement_mass: bool,
    pub ctrl: &'a mut T,
    relative_time: f64,
    init_tai_secs: f64,
}

impl<'a, T: ThrustControl> Propulsion<'a, T> {
    pub fn new(
        ctrl: &'a mut T,
        fuel_mass: f64,
        abs_time: Epoch,
        thrusters: Vec<Thruster>,
        decrement_mass: bool,
    ) -> Self {
        Self {
            dry_mass: 0.0,
            fuel_mass,
            abs_time,
            thrusters,
            decrement_mass,
            ctrl,
            relative_time: 0.0,
            init_tai_secs: abs_time.as_tai_seconds(),
        }
    }
}

impl<'a, T: ThrustControl> Dynamics for Propulsion<'a, T> {
    type StateSize = U7;
    type StateType = f64;

    /// Returns the relative time
    fn time(&self) -> f64 {
        self.relative_time
    }

    /// State of the propulsion is only the fuel mass
    fn state(&self) -> Self::StateType {
        self.fuel_mass
    }

    /// Propulsion state is a vector of six zeros followed by the fuel mass
    fn state_vector(&self) -> VectorN<f64, Self::StateSize> {
        let mut s = VectorN::<f64, U7>::zeros();
        s[6] = self.fuel_mass;
        s
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        self.relative_time = new_t;
        self.abs_time = Epoch::from_tai_seconds(self.init_tai_secs + new_t);
        self.ctrl
            .next(self.abs_time, &new_state.fixed_rows::<U6>(0).into_owned());
        self.fuel_mass = new_state[6];
        if self.decrement_mass {
            assert!(
                self.fuel_mass >= 0.0,
                "negative fuel mass after {:?} seconds",
                self.abs_time.as_tai_seconds() - self.init_tai_secs
            );
        }
    }

    fn eom(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let celestial_state = state.fixed_rows::<U6>(0).into_owned();
        let dt = Epoch::from_tai_seconds(self.init_tai_secs + t);
        let thrust_power = self.ctrl.throttle(dt, &celestial_state);
        if thrust_power > 0.0 {
            // Thrust arc
            let thrust_inertial = self.ctrl.direction(dt, &celestial_state);
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
            total_thrust /= self.dry_mass + state[6];
            total_thrust *= 1e-3; // Convert m/s^-2 to km/s^-2
            let mut rtn = VectorN::<f64, U7>::zeros();
            for i in 0..3 {
                rtn[i + 3] = thrust_inertial[i] * total_thrust;
            }
            if self.decrement_mass {
                rtn[6] = -fuel_usage;
            } else {
                rtn[6] = 0.0;
            }
            rtn
        } else {
            let mut rtn = VectorN::<f64, U7>::zeros();
            rtn[6] = 0.0;
            rtn
        }
    }
}
