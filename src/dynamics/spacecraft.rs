use super::orbital::OrbitalDynamics;
use super::thrustctrl::ThrustControl;
use super::{Dynamics, ForceModel};
use crate::dimensions::{DimName, MatrixN, Vector1, Vector3, VectorN, U3, U4, U42, U43, U6, U7};
use dynamics::Hyperdual;
// use crate::od::Estimable;
use crate::time::TimeUnit;
use crate::TimeTagged;
use celestia::SpacecraftState;
use errors::NyxError;
use std::sync::Arc;
use State;

pub use super::solarpressure::SolarPressure;

const NORM_ERR: f64 = 1e-12;
const STD_GRAVITY: f64 = 9.80665; // From NIST special publication 330, 2008 edition

#[derive(Clone)]
pub struct Spacecraft<'a> {
    pub orbital_dyn: OrbitalDynamics<'a>,
    pub force_models: Vec<Arc<dyn ForceModel + 'a>>,
    pub ctrl: Option<Arc<dyn ThrustControl + 'a>>,
    pub decrement_mass: bool,
}

impl<'a> Spacecraft<'a> {
    /// Initialize a Spacecraft with a set of orbital dynamics and a propulsion subsystem.
    pub fn with_ctrl(orbital_dyn: OrbitalDynamics<'a>, ctrl: Arc<dyn ThrustControl + 'a>) -> Self {
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

    pub fn add_model(&mut self, force_model: Arc<dyn ForceModel + 'a>) {
        self.force_models.push(force_model);
    }
}

impl<'a> Dynamics for Spacecraft<'a>
// where
//     // D::StateSize: DimNameAdd<U1>,
//     DefaultAllocator: Allocator<f64, DimNameSum<D::StateSize, U1>>
//         + Allocator<f64, D::StateSize>
//         + Allocator<f64, U1>
//         + Allocator<f64, U1, D::StateSize>
//         + Allocator<f64, D::StateSize, U1>
//         + Allocator<f64, DimNameSum<D::StateSize, U1>, U1>
//         + Allocator<f64, U1, DimNameSum<D::StateSize, U1>>,
{
    // type StateSize = U7;
    // type PropVecSize = U43; // This implies that we are NOT computing the STM with the mass
    type HyperdualSize = U7; // This implies that we are NOT computing the STM with the mass
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
        state: &VectorN<f64, U43>,
        ctx: &SpacecraftState,
    ) -> Result<VectorN<f64, U43>, NyxError> {
        // Compute the orbital dynamics
        let orbital_dyn_vec = state.fixed_rows::<U42>(0).into_owned();
        let d_x_orbital_dyn = self
            .orbital_dyn
            .eom(delta_t, &orbital_dyn_vec, &ctx.orbit)?;
        let mut d_x = VectorN::<f64, U43>::from_iterator(
            d_x_orbital_dyn
                .iter()
                .chain(Vector1::new(0.0).iter())
                .cloned(),
        );

        // let orbital_dyn_state = self.orbital_dyn.orbital_state_ctor(t, &orbital_dyn_vec);
        // let orbital_dyn_state = ctx.orbit.ctor_from(delta_t, &orbital_dyn_vec);
        let mut total_mass = ctx.dry_mass;

        // Rebuild the osculating state for the EOM context.
        let osc_sc = ctx.ctor_from(delta_t, &d_x);

        // Now compute the other dynamics as needed.
        if let Some(ctrl) = &self.ctrl {
            let (thrust_force, fuel_rate) = {
                let thruster = osc_sc.thruster.unwrap();
                let thrust_power = ctrl.throttle(&osc_sc.orbit);
                if thrust_power > 0.0 {
                    // Thrust arc
                    let thrust_inertial = ctrl.direction(&osc_sc.orbit);
                    if (thrust_inertial.norm() - 1.0).abs() > NORM_ERR {
                        return Err(NyxError::InvalidThrusterCtrlNorm(thrust_inertial.norm()));
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
                d_x[i + 3] += thrust_force[i] / (ctx.dry_mass + state[U7::dim() - 1]);
            }
            d_x[U7::dim() - 1] += fuel_rate;
        }

        // Compute additional force models as needed
        for model in &self.force_models {
            let model_frc = model.eom(&osc_sc)? / total_mass;
            for i in 0..3 {
                d_x[i + 3] += model_frc[i];
            }
        }

        Ok(d_x)
    }

    fn dual_eom(
        &self,
        delta_t_s: f64,
        state_vec: &VectorN<Hyperdual<f64, Self::HyperdualSize>, U7>,
        ctx: &Self::StateType,
    ) -> Result<(VectorN<f64, U7>, MatrixN<f64, U7>), NyxError> {
        let pos_vel = state_vec.fixed_rows::<U6>(0).into_owned();
        let (orb_state, orb_grad) = self.orbital_dyn.dual_eom(delta_t_s, &pos_vel, &ctx.orbit)?;
        // Rebuild the appropriately sized state and STM.
        let mut d_x = VectorN::<f64, U7>::from_iterator(
            orb_state.iter().chain(Vector1::new(0.0).iter()).cloned(),
        );
        let mut grad = MatrixN::<f64, U7>::zeros();
        for i in 0..U6::dim() {
            for j in 0..U6::dim() {
                grad[(i, j)] = orb_grad[(i, j)];
            }
        }

        // Call the EOMs
        let mut total_mass = ctx.dry_mass;
        let radius = state_vec.fixed_rows::<U3>(0).into_owned();

        // Recreate the osculating state.
        let mut osc_sc = ctx;
        osc_sc.set_epoch(ctx.epoch() + delta_t_s * TimeUnit::Second);
        osc_sc.orbit.x = orb_state[0];
        osc_sc.orbit.y = orb_state[1];
        osc_sc.orbit.z = orb_state[2];
        osc_sc.orbit.vx = orb_state[3];
        osc_sc.orbit.vy = orb_state[4];
        osc_sc.orbit.vz = orb_state[5];

        for model in &self.force_models {
            // let model_frc = model.dual_eom(delta_t, &radius, &osc_sc)? / total_mass;
            let (model_frc, model_grad) = model.dual_eom(&radius, &osc_sc)?;
            for i in 0..U3::dim() {
                d_x[i + 3] += model_frc[i] / total_mass;
                for j in 1..U4::dim() {
                    grad[(i + 3, j - 1)] += model_grad[(i, j - 1)] / total_mass;
                }
            }
        }

        Ok((d_x, grad))
    }
}

// impl<'a> Spacecraft<'a, OrbitalDynamicsStm<'a>> {
//     /// Initialize a Spacecraft with a set of orbital dynamics, with SRP enabled, and the STM computation
//     pub fn with_stm(orbital_dyn: OrbitalDynamicsStm<'a>, dry_mass: f64) -> Self {
//         // Set the dry mass of the propulsion system
//         Self {
//             orbital_dyn,
//             prop: None,
//             force_models: Vec::new(),
//             dry_mass,
//             fuel_mass: 0.0,
//         }
//     }

//     pub fn add_model(&mut self, force_model: Box<dyn ForceModel<CtxType = SpacecraftState> + 'a>) {
//         self.force_models.push(force_model);
//     }
// }

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
