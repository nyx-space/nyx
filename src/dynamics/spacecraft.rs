use super::orbital::{OrbitalDynamics, OrbitalDynamicsStm, OrbitalDynamicsT};
use super::propulsion::Propulsion;
use super::{Dynamics, ForceModel};
use crate::dimensions::allocator::Allocator;
use crate::dimensions::dimension::{DimNameAdd, DimNameSum};
use crate::dimensions::{DefaultAllocator, DimName, Matrix6, Vector1, VectorN, U1, U6};
use crate::od::Estimable;
use crate::time::Epoch;
use celestia::State;
use std::fmt;

pub use super::solarpressure::SolarPressure;

pub struct Spacecraft<'a, D: OrbitalDynamicsT> {
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

impl<'a, D: OrbitalDynamicsT> Dynamics for Spacecraft<'a, D>
where
    D::StateSize: DimNameAdd<U1>,
    DefaultAllocator: Allocator<f64, DimNameSum<D::StateSize, U1>>
        + Allocator<f64, D::StateSize>
        + Allocator<f64, U1>
        + Allocator<f64, U1, D::StateSize>
        + Allocator<f64, D::StateSize, U1>,
{
    type StateSize = DimNameSum<D::StateSize, U1>;
    type StateType = SpacecraftState;

    fn time(&self) -> f64 {
        self.orbital_dyn.time()
    }

    fn state(&self) -> Self::StateType {
        SpacecraftState {
            orbit: self.orbital_dyn.orbital_state(),
            dry_mass: self.dry_mass,
            fuel_mass: self.fuel_mass,
            stm: self.orbital_dyn.stm(),
        }
    }

    fn state_vector(&self) -> VectorN<f64, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize> + Allocator<f64, D::StateSize>,
    {
        VectorN::<f64, Self::StateSize>::from_iterator(
            self.orbital_dyn
                .state_vector()
                .iter()
                .chain(Vector1::new(self.fuel_mass).iter())
                .cloned(),
        )
    }

    fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>)
    where
        DefaultAllocator: Allocator<f64, Self::StateSize> + Allocator<f64, D::StateSize>,
    {
        let orbital_dyn_state = new_state.fixed_rows::<D::StateSize>(0).into_owned();
        self.orbital_dyn.set_state(new_t, &orbital_dyn_state);
        self.fuel_mass = new_state[Self::StateSize::dim() - 1];
        if let Some(prop) = &self.prop {
            if prop.decrement_mass {
                assert!(
                    self.fuel_mass >= 0.0,
                    "negative fuel mass at {:?}",
                    self.orbital_dyn.orbital_state().dt
                );
            }
        }
        if let Some(prop) = self.prop.as_mut() {
            prop.set_state(&self.orbital_dyn.orbital_state());
        }
    }

    fn eom(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize>,
    {
        // Compute the orbital dynamics
        let orbital_dyn_vec = state.fixed_rows::<D::StateSize>(0).into_owned();
        let d_x_orbital_dyn = self.orbital_dyn.eom(t, &orbital_dyn_vec);
        let mut d_x = VectorN::<f64, Self::StateSize>::from_iterator(
            d_x_orbital_dyn
                .iter()
                .chain(Vector1::new(0.0).iter())
                .cloned(),
        );

        let orbital_dyn_state = self.orbital_dyn.orbital_state_ctor(t, &orbital_dyn_vec);
        let mut total_mass = self.dry_mass;

        // Now compute the other dynamics as needed.
        if let Some(prop) = &self.prop {
            let (thrust_force, fuel_rate) = prop.eom(&orbital_dyn_state);
            // Add the fuel mass to the total mass, minus the change in fuel
            total_mass += self.fuel_mass + fuel_rate;
            for i in 0..3 {
                d_x[i + 3] += thrust_force[i] / (self.dry_mass + state[Self::StateSize::dim() - 1]);
            }
            d_x[Self::StateSize::dim() - 1] += fuel_rate;
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

impl<'a> Estimable<SpacecraftState> for Spacecraft<'a, OrbitalDynamicsStm<'a>> {
    type LinStateSize = U6;

    fn to_measurement(&self, prop_state: &Self::StateType) -> (Epoch, SpacecraftState) {
        (prop_state.orbit.dt, *prop_state)
    }

    fn extract_stm(&self, prop_state: &Self::StateType) -> Matrix6<f64> {
        prop_state.stm.unwrap()
    }

    fn extract_estimated_state(
        &self,
        prop_state: &Self::StateType,
    ) -> VectorN<f64, Self::LinStateSize> {
        prop_state.orbit.to_cartesian_vec()
    }

    /// Returns the estimated state
    fn set_estimated_state(&mut self, new_state: VectorN<f64, Self::LinStateSize>) {
        self.orbital_dyn.state.x = new_state[0];
        self.orbital_dyn.state.y = new_state[1];
        self.orbital_dyn.state.z = new_state[2];
        self.orbital_dyn.state.vx = new_state[3];
        self.orbital_dyn.state.vy = new_state[4];
        self.orbital_dyn.state.vz = new_state[5];
    }
}

#[derive(Clone, Copy, Debug)]
pub struct SpacecraftState {
    pub orbit: State,
    pub dry_mass: f64,
    pub fuel_mass: f64,
    pub stm: Option<Matrix6<f64>>,
}

impl fmt::Display for SpacecraftState {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:o}\t{} kg", self.orbit, self.dry_mass + self.fuel_mass)
    }
}
