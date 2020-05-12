use super::orbital::{OrbitalDynamics, OrbitalDynamicsStm};
use super::propulsion::Propulsion;
use super::{Dynamics, ForceModel};
use crate::dimensions::{DimName, Vector1, VectorN, U42, U43, U6, U7};
use celestia::State;
use std::fmt;

pub use super::solarpressure::SolarPressure;

pub struct Spacecraft<'a, D: Dynamics> {
    pub orbital_dyn: D,
    pub force_models: Vec<Box<dyn ForceModel + 'a>>,
    pub prop: Option<Propulsion>,
    /// in kg
    pub dry_mass: f64,
    /// in kg
    pub fuel_mass: f64,
}

impl<'a> Spacecraft<'a, OrbitalDynamics<'a>> {
    /// Initialize a Spacecraft with a set of orbital dynamics and a propulsion subsystem.
    pub fn with_prop(
        orbital_dyn: OrbitalDynamics<'a>,
        prop: Propulsion,
        dry_mass: f64,
        fuel_mass: f64,
    ) -> Self {
        Self {
            orbital_dyn,
            prop: Some(prop),
            force_models: Vec::new(),
            dry_mass,
            fuel_mass,
        }
    }

    /// Initialize a Spacecraft with a set of orbital dynamics and with SRP enabled.
    pub fn new(orbital_dyn: OrbitalDynamics<'a>, dry_mass: f64) -> Self {
        // Set the dry mass of the propulsion system
        Self {
            orbital_dyn,
            prop: None,
            force_models: Vec::new(),
            dry_mass,
            fuel_mass: 0.0,
        }
    }

    pub fn add_model(&mut self, force_model: Box<dyn ForceModel + 'a>) {
        self.force_models.push(force_model);
    }
}

impl<'a> Dynamics for Spacecraft<'a, OrbitalDynamics<'a>> {
    type StateSize = U7;
    type StateType = SpacecraftState;

    fn time(&self) -> f64 {
        self.orbital_dyn.time()
    }

    fn state(&self) -> Self::StateType {
        SpacecraftState {
            orbit: self.orbital_dyn.state,
            dry_mass: self.dry_mass,
            fuel_mass: self.fuel_mass,
        }
    }

    fn state_vector(&self) -> VectorN<f64, Self::StateSize> {
        VectorN::<f64, Self::StateSize>::from_iterator(
            self.orbital_dyn
                .state_vector()
                .iter()
                .chain(Vector1::new(self.fuel_mass).iter())
                .cloned(),
        )
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        let orbital_dyn_state = new_state.fixed_rows::<U6>(0).into_owned();
        self.orbital_dyn.set_state(new_t, &orbital_dyn_state);
        self.fuel_mass = new_state[Self::StateSize::dim() - 1];
        if let Some(prop) = &self.prop {
            if prop.decrement_mass {
                assert!(
                    self.fuel_mass >= 0.0,
                    "negative fuel mass at {:?}",
                    self.orbital_dyn.state().dt
                );
            }
        }
        if let Some(prop) = self.prop.as_mut() {
            prop.set_state(&self.orbital_dyn.state());
        }
    }

    fn eom(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        // Compute the orbital dynamics
        let orbital_dyn_vec = state.fixed_rows::<U6>(0).into_owned();
        let d_x_orbital_dyn = self.orbital_dyn.eom(t, &orbital_dyn_vec);
        let mut d_x = VectorN::<f64, Self::StateSize>::from_iterator(
            d_x_orbital_dyn
                .iter()
                .chain(Vector1::new(0.0).iter())
                .cloned(),
        );

        let orbital_dyn_state = self.orbital_dyn.state_ctor(t, &orbital_dyn_vec);
        let mut total_mass = self.dry_mass;

        // Now compute the other dynamics as needed.
        if let Some(prop) = &self.prop {
            let (thrust_force, fuel_usage) = prop.eom(&orbital_dyn_state);
            // Add the fuel mass to the total mass, minus the change in fuel
            total_mass += self.fuel_mass + fuel_usage;
            for i in 0..3 {
                d_x[i + 3] += thrust_force[i] / (self.dry_mass + state[Self::StateSize::dim() - 1]);
            }
            d_x[Self::StateSize::dim() - 1] += fuel_usage;
        }

        // Compute additional force models as needed
        for model in &self.force_models {
            let model_frc = model.eom(&orbital_dyn_state) / total_mass;
            for i in 0..3 {
                d_x[i + 3] += model_frc[i];
            }
        }

        d_x
    }
}

impl<'a> Spacecraft<'a, OrbitalDynamicsStm<'a>> {
    /// Initialize a Spacecraft with a set of orbital dynamics, with SRP enabled, and the STM computation
    pub fn with_stm(orbital_dyn: OrbitalDynamicsStm<'a>, dry_mass: f64) -> Self {
        // Set the dry mass of the propulsion system
        Self {
            orbital_dyn,
            prop: None,
            force_models: Vec::new(),
            dry_mass,
            fuel_mass: 0.0,
        }
    }

    pub fn add_model(&mut self, force_model: Box<dyn ForceModel + 'a>) {
        self.force_models.push(force_model);
    }
}

impl<'a> Dynamics for Spacecraft<'a, OrbitalDynamicsStm<'a>> {
    type StateSize = U43;
    type StateType = SpacecraftState;

    fn time(&self) -> f64 {
        self.orbital_dyn.time()
    }

    fn state(&self) -> Self::StateType {
        SpacecraftState {
            orbit: self.orbital_dyn.state,
            dry_mass: self.dry_mass,
            fuel_mass: self.fuel_mass,
        }
    }

    fn state_vector(&self) -> VectorN<f64, Self::StateSize> {
        VectorN::<f64, Self::StateSize>::from_iterator(
            self.orbital_dyn
                .state_vector()
                .iter()
                .chain(Vector1::new(self.fuel_mass).iter())
                .cloned(),
        )
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
        let orbital_dyn_state = new_state.fixed_rows::<U42>(0).into_owned();
        self.orbital_dyn.set_state(new_t, &orbital_dyn_state);
        self.fuel_mass = new_state[Self::StateSize::dim() - 1];
        if let Some(prop) = &self.prop {
            if prop.decrement_mass {
                assert!(
                    self.fuel_mass >= 0.0,
                    "negative fuel mass at {:?}",
                    self.orbital_dyn.state.dt
                );
            }
        }
        if let Some(prop) = self.prop.as_mut() {
            prop.set_state(&self.orbital_dyn.state);
        }
    }

    fn eom(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        // Compute the orbital dynamics
        let orbital_dyn_vec = state.fixed_rows::<U42>(0).into_owned();
        let d_x_orbital_dyn = self.orbital_dyn.eom(t, &orbital_dyn_vec);
        let mut d_x = VectorN::<f64, Self::StateSize>::from_iterator(
            d_x_orbital_dyn
                .iter()
                .chain(Vector1::new(0.0).iter())
                .cloned(),
        );

        let orbital_dyn_state = self.orbital_dyn.state_ctor(t, &orbital_dyn_vec).0;
        let mut total_mass = self.dry_mass;

        // Now compute the other dynamics as needed.
        if let Some(prop) = &self.prop {
            let (thrust_force, fuel_usage) = prop.eom(&orbital_dyn_state);
            // Add the fuel mass to the total mass, minus the change in fuel
            total_mass += self.fuel_mass + fuel_usage;
            for i in 0..3 {
                d_x[i + 3] += thrust_force[i] / (self.dry_mass + state[Self::StateSize::dim() - 1]);
            }
            d_x[Self::StateSize::dim() - 1] += fuel_usage;
        }

        // Compute additional force models as needed
        for model in &self.force_models {
            let model_frc = model.eom(&orbital_dyn_state) / total_mass;
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
