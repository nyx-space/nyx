use super::celestial::CelestialDynamics;
use super::deltavctrl::DeltaVctrl;
use super::na::{VectorN, U6};
use super::Dynamics;
use celestia::{Geoid, State};

pub struct MissionArc<'a, D: DeltaVctrl> {
    pub celestial: &'a mut CelestialDynamics<'a>,
    pub ctrl: &'a mut D,
}

impl<'a, D: DeltaVctrl> MissionArc<'a, D> {
    pub fn new(celestial: &'a mut CelestialDynamics<'a>, ctrl: &'a mut D) -> Self {
        Self { celestial, ctrl }
    }
}

impl<'a, D: DeltaVctrl> Dynamics for MissionArc<'a, D> {
    type StateSize = U6;
    type StateType = State<Geoid>;

    /// Returns the relative time
    fn time(&self) -> f64 {
        self.celestial.time()
    }

    /// State of the mission arc is always the celestial state
    fn state(&self) -> Self::StateType {
        self.celestial.state()
    }

    fn build_state(&self, t: f64, in_state: &VectorN<f64, Self::StateSize>) -> Self::StateType {
        self.celestial.build_state(t, in_state)
    }

    /// Mission arc state is a vector of six zeros followed by the fuel mass
    fn state_vector(&self) -> VectorN<f64, Self::StateSize> {
        self.celestial.state_vector()
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        self.celestial.set_state(new_t, new_state);
        self.ctrl.next(&self.celestial.state());
    }

    fn eom(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let mut d_x = self.celestial.eom(t, &state);
        let eom_state = self.celestial.state_ctor(t, &state);
        let dv_inertial = self.ctrl.ctrl_vector(&eom_state);

        for i in 0..3 {
            d_x[i + 3] = dv_inertial[i];
        }

        d_x
    }
}
