use super::orbital::OrbitalDynamics;
use super::propulsion::Propulsion;
use super::thrustctrl::ThrustControl;
use super::{Dynamics, ForceModel};
use crate::dimensions::{Vector1, VectorN, U6, U7};
use celestia::State;
use std::fmt;
use std::marker::PhantomData;

pub use super::solarpressure::SolarPressure;

pub struct Spacecraft<'a, T: ThrustControl> {
    pub celestial: OrbitalDynamics<'a>,
    pub force_models: Vec<Box<dyn ForceModel + 'a>>,
    pub prop: Option<Propulsion<T>>,
    /// in kg
    pub dry_mass: f64,
    /// in kg
    pub fuel_mass: f64,
    _marker: PhantomData<T>,
}

impl<'a, T: ThrustControl> Spacecraft<'a, T> {
    /// Initialize a Spacecraft with a set of celestial dynamics and a propulsion subsystem.
    pub fn with_prop(
        celestial: OrbitalDynamics<'a>,
        prop: Propulsion<T>,
        dry_mass: f64,
        fuel_mass: f64,
    ) -> Self {
        Self {
            celestial,
            prop: Some(prop),
            force_models: Vec::new(),
            dry_mass,
            fuel_mass,
            _marker: PhantomData,
        }
    }

    /// Initialize a Spacecraft with a set of celestial dynamics and with SRP enabled.
    pub fn new(celestial: OrbitalDynamics<'a>, dry_mass: f64) -> Self {
        // Set the dry mass of the propulsion system
        Self {
            celestial,
            prop: None,
            force_models: Vec::new(),
            dry_mass,
            fuel_mass: 0.0,
            _marker: PhantomData,
        }
    }

    pub fn add_model(&mut self, force_model: Box<dyn ForceModel + 'a>) {
        self.force_models.push(force_model);
    }
}

impl<'a, T: ThrustControl> Dynamics for Spacecraft<'a, T> {
    type StateSize = U7;
    type StateType = SpacecraftState;

    fn time(&self) -> f64 {
        self.celestial.time()
    }

    fn state(&self) -> Self::StateType {
        SpacecraftState {
            orbit: self.celestial.state(),
            dry_mass: self.dry_mass,
            fuel_mass: self.fuel_mass,
        }
    }

    fn state_vector(&self) -> VectorN<f64, Self::StateSize> {
        VectorN::<f64, U7>::from_iterator(
            self.celestial
                .state_vector()
                .iter()
                .chain(Vector1::new(self.fuel_mass).iter())
                .cloned(),
        )
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        let celestial_state = new_state.fixed_rows::<U6>(0).into_owned();
        self.celestial.set_state(new_t, &celestial_state);
        self.fuel_mass = new_state[6];
        if let Some(prop) = &self.prop {
            if prop.decrement_mass {
                assert!(
                    self.fuel_mass >= 0.0,
                    "negative fuel mass at {:?}",
                    self.celestial.state().dt
                );
            }
        }
        if let Some(prop) = self.prop.as_mut() {
            prop.set_state(&self.celestial.state());
        }
    }

    fn eom(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        // Compute the celestial dynamics
        let celestial_vec = state.fixed_rows::<U6>(0).into_owned();
        let d_x_celestial = self.celestial.eom(t, &celestial_vec);
        let mut d_x = VectorN::<f64, U7>::from_iterator(
            d_x_celestial
                .iter()
                .chain(Vector1::new(0.0).iter())
                .cloned(),
        );

        let celestial_state = self.celestial.state_ctor(t, &celestial_vec);
        let mut total_mass = self.dry_mass;

        // Now compute the other dynamics as needed.
        if let Some(prop) = &self.prop {
            let (thrust_force, fuel_usage) = prop.eom(&celestial_state);
            // Add the fuel mass to the total mass, minus the change in fuel
            total_mass += self.fuel_mass + fuel_usage;
            for i in 0..3 {
                d_x[i + 3] += thrust_force[i] / (self.dry_mass + state[6]);
            }
            d_x[6] += fuel_usage;
        }

        // Compute additional force models as needed
        for model in &self.force_models {
            let model_frc = model.eom(&celestial_state) / total_mass;
            for i in 0..3 {
                d_x[i + 3] += model_frc[i];
            }
        }

        d_x
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SpacecraftState {
    pub orbit: State,
    pub dry_mass: f64,
    pub fuel_mass: f64,
}

impl fmt::Display for SpacecraftState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:o}\t{} kg", self.orbit, self.dry_mass + self.fuel_mass)
    }
}
