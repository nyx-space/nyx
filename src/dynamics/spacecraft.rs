use super::orbital::{OrbitalDynamics, OrbitalDynamicsStm, OrbitalDynamicsT};
use super::propulsion::Propulsion;
use super::thrustctrl::ThrustControl;
use super::{Dynamics, ForceModel};
use crate::dimensions::allocator::Allocator;
use crate::dimensions::dimension::{DimNameAdd, DimNameSum};
use crate::dimensions::{DefaultAllocator, DimName, Matrix6, Vector1, Vector3, VectorN, U1, U6};
// use crate::od::Estimable;
use celestia::SpacecraftState;
use std::sync::Arc;
use State;

pub use super::solarpressure::SolarPressure;

const NORM_ERR: f64 = 1e-12;
const STD_GRAVITY: f64 = 9.80665; // From NIST special publication 330, 2008 edition

#[derive(Clone)]
pub struct Spacecraft<D: OrbitalDynamicsT> {
    pub orbital_dyn: D,
    pub force_models: Vec<Arc<dyn ForceModel<CtxType = SpacecraftState>>>,
    pub ctrl: Option<Arc<dyn ThrustControl>>,
    pub decrement_mass: bool,
}

impl<'a> Spacecraft<OrbitalDynamics<'a>> {
    /// Initialize a Spacecraft with a set of orbital dynamics and a propulsion subsystem.
    pub fn with_ctrl(orbital_dyn: OrbitalDynamics<'a>, ctrl: Arc<dyn ThrustControl>) -> Self {
        Self {
            orbital_dyn,
            ctrl: Some(ctrl),
            force_models: Vec::new(),
            decrement_mass: true,
        }
    }

    /// Initialize a Spacecraft with a set of orbital dynamics and with SRP enabled.
    pub fn new(orbital_dyn: OrbitalDynamics<'a>) -> Self {
        // Set the dry mass of the propulsion system
        Self {
            orbital_dyn,
            ctrl: None,
            force_models: Vec::new(),
            decrement_mass: true,
        }
    }

    pub fn add_model(&mut self, force_model: Arc<dyn ForceModel<CtxType = SpacecraftState>>) {
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
    // fn state_vector(&self) -> VectorN<f64, Self::StateSize>
    // where
    //     DefaultAllocator: Allocator<f64, Self::StateSize> + Allocator<f64, D::StateSize>,
    // {
    //     VectorN::<f64, Self::StateSize>::from_iterator(
    //         self.orbital_dyn
    //             .state_vector()
    //             .iter()
    //             .chain(Vector1::new(self.fuel_mass).iter())
    //             .cloned(),
    //     )
    // }

    // fn set_state(
    //     &mut self,
    //     new_t: f64,
    //     new_state: &VectorN<f64, Self::StateSize>,
    // ) -> Result<(), NyxError>
    // where
    //     DefaultAllocator: Allocator<f64, Self::StateSize> + Allocator<f64, D::StateSize>,
    // {
    //     let orbital_dyn_state = new_state.fixed_rows::<D::StateSize>(0).into_owned();
    //     self.orbital_dyn.set_state(new_t, &orbital_dyn_state)?;
    //     self.fuel_mass = new_state[Self::StateSize::dim() - 1];
    //     if let Some(prop) = &self.prop {
    //         if prop.decrement_mass && self.fuel_mass < 0.0 {
    //             error!(
    //                 "negative fuel mass at {:?}",
    //                 self.orbital_dyn.orbital_state().dt
    //             );
    //             return Err(NyxError::FuelExhausted);
    //         }
    //     }
    //     if let Some(prop) = self.prop.as_mut() {
    //         prop.set_state(&self.orbital_dyn.orbital_state());
    //     }

    //     Ok(())
    // }

    fn eom(
        &self,
        delta_t: f64,
        state: &VectorN<f64, Self::StateSize>,
        ctx: &SpacecraftState,
    ) -> VectorN<f64, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize>,
    {
        // Compute the orbital dynamics
        let orbital_dyn_vec = state.fixed_rows::<U6>(0).into_owned();
        let d_x_orbital_dyn = self.orbital_dyn.eom(delta_t, &orbital_dyn_vec, &ctx.orbit);
        let mut d_x = VectorN::<f64, Self::StateSize>::from_iterator(
            d_x_orbital_dyn
                .iter()
                .chain(Vector1::new(0.0).iter())
                .cloned(),
        );

        // let orbital_dyn_state = self.orbital_dyn.orbital_state_ctor(t, &orbital_dyn_vec);
        let orbital_dyn_state = ctx.orbit.ctor_from(delta_t, &orbital_dyn_vec);
        let mut total_mass = ctx.dry_mass;

        // Now compute the other dynamics as needed.
        if let Some(ctrl) = &self.ctrl {
            let (thrust_force, fuel_rate) = {
                let osc = ctx.ctor_from(delta_t, state);
                let thruster = osc.thruster.unwrap();
                let thrust_power = ctrl.throttle(&osc.orbit);
                if thrust_power > 0.0 {
                    // Thrust arc
                    let thrust_inertial = ctrl.direction(&osc.orbit);
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

                    total_thrust += thrust_power * thruster.thrust;
                    fuel_usage += thrust_power * thruster.thrust / (thruster.isp * STD_GRAVITY);
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
            };
            // Add the fuel mass to the total mass, minus the change in fuel
            total_mass += ctx.fuel_mass + fuel_rate;
            for i in 0..3 {
                d_x[i + 3] += thrust_force[i] / (ctx.dry_mass + state[Self::StateSize::dim() - 1]);
            }
            d_x[Self::StateSize::dim() - 1] += fuel_rate;
        }

        // Compute additional force models as needed
        for model in &self.force_models {
            let model_frc = model.eom(&orbital_dyn_state, &ctx.orbit) / total_mass;
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

    pub fn add_model(&mut self, force_model: Box<dyn ForceModel<CtxType = SpacecraftState> + 'a>) {
        self.force_models.push(force_model);
    }
}

// impl<'a> Estimable<SpacecraftState> for Spacecraft<'a, OrbitalDynamicsStm<'a>> {
//     type LinStateSize = U6;

//     fn to_measurement(&self, prop_state: &Self::StateType) -> SpacecraftState {
//         *prop_state
//     }

//     fn extract_stm(&self, prop_state: &Self::StateType) -> Matrix6<f64> {
//         prop_state.stm.unwrap()
//     }

//     fn extract_estimated_state(
//         &self,
//         prop_state: &Self::StateType,
//     ) -> VectorN<f64, Self::LinStateSize> {
//         prop_state.orbit.to_cartesian_vec()
//     }

//     /// Returns the estimated state
//     fn set_estimated_state(&mut self, new_state: VectorN<f64, Self::LinStateSize>) {
//         self.orbital_dyn.state.x = new_state[0];
//         self.orbital_dyn.state.y = new_state[1];
//         self.orbital_dyn.state.z = new_state[2];
//         self.orbital_dyn.state.vx = new_state[3];
//         self.orbital_dyn.state.vy = new_state[4];
//         self.orbital_dyn.state.vz = new_state[5];
//     }
// }
