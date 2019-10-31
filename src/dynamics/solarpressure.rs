use super::hifitime::Epoch;
use super::na::{Vector1, VectorN, U1, U6, U7};
use super::Dynamics;
use super::{Cosm, Geoid, State};

pub struct SolarPressure<'a> {
    /// in kg, set in the Spacecraft's eom.
    pub sc_mass: f64,
    pub shadow_bodies: Vec<Geoid>,
    pub cosm: &'a Cosm,
}

impl<'a> SolarPressure<'a> {
    pub fn new(shadow_bodies: Vec<Geoid>, cosm: &'a Cosm) -> Self {
        Self {
            sc_mass: 0.0,
            shadow_bodies,
            cosm,
        }
    }
}

impl<'a> Dynamics for SolarPressure<'a> {
    type StateSize = U1;
    type StateType = Vector1<f64>;

    /// Time of SolarPressure will always return zero!
    fn time(&self) -> f64 {
        0.0
    }

    /// SolarPressure state is a vector of six zeros followed by the fuel mass
    fn state_vector(&self) -> VectorN<f64, Self::StateSize> {
        Vector1::zeros()
    }

    /// There is no state to set for SolarPressure, so this function does nothing.
    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {}
}
