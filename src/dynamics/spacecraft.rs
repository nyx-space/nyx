/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

use super::orbital::OrbitalDynamics;
use super::thrustctrl::ThrustControl;
use super::{Dynamics, ForceModel};
use crate::dimensions::{DimName, MatrixN, Vector1, Vector3, VectorN, U3, U4, U42, U43, U6, U7};
use crate::dynamics::Hyperdual;
// use crate::od::Estimable;
use crate::celestia::Spacecraft;
use crate::errors::NyxError;
use crate::state::State;
use crate::time::TimeUnit;
use crate::TimeTagged;
use std::sync::Arc;

pub use super::solarpressure::SolarPressure;

const NORM_ERR: f64 = 1e-12;
const STD_GRAVITY: f64 = 9.80665; // From NIST special publication 330, 2008 edition

#[derive(Clone)]
pub struct SpacecraftDynamics<'a> {
    pub orbital_dyn: Arc<OrbitalDynamics<'a>>,
    pub force_models: Vec<Arc<dyn ForceModel + 'a>>,
    pub ctrl: Option<Arc<dyn ThrustControl + 'a>>,
    pub decrement_mass: bool,
}

impl<'a> SpacecraftDynamics<'a> {
    /// Initialize a Spacecraft with a set of orbital dynamics and a propulsion subsystem.
    /// By default, the mass of the vehicle will be decremented as propellant is consummed.
    pub fn with_ctrl(
        orbital_dyn: Arc<OrbitalDynamics<'a>>,
        ctrl: Arc<dyn ThrustControl + 'a>,
    ) -> Arc<Self> {
        Arc::new(Self {
            orbital_dyn,
            ctrl: Some(ctrl),
            force_models: Vec::new(),
            decrement_mass: true,
        })
    }

    /// Initialize a Spacecraft with a set of orbital dynamics and a propulsion subsystem.
    /// Will _not_ decrement the fuel mass as propellant is consummed.
    pub fn with_ctrl_no_decr(
        orbital_dyn: Arc<OrbitalDynamics<'a>>,
        ctrl: Arc<dyn ThrustControl + 'a>,
    ) -> Arc<Self> {
        Arc::new(Self {
            orbital_dyn,
            ctrl: Some(ctrl),
            force_models: Vec::new(),
            decrement_mass: false,
        })
    }

    /// Initialize a Spacecraft with a set of orbital dynamics and with SRP enabled.
    pub fn new(orbital_dyn: Arc<OrbitalDynamics<'a>>) -> Arc<Self> {
        Arc::new(Self::new_raw(orbital_dyn))
    }

    /// Initialize a Spacecraft with a set of orbital dynamics and with SRP enabled.
    pub fn new_raw(orbital_dyn: Arc<OrbitalDynamics<'a>>) -> Self {
        // Set the dry mass of the propulsion system
        Self {
            orbital_dyn,
            ctrl: None,
            force_models: Vec::new(),
            decrement_mass: true,
        }
    }

    /// Initialize new spacecraft dynamics with the provided orbital mechanics and with the provided force model.
    pub fn with_model(
        orbital_dyn: Arc<OrbitalDynamics<'a>>,
        force_model: Arc<dyn ForceModel + 'a>,
    ) -> Arc<Self> {
        let mut me = Self::new_raw(orbital_dyn);
        me.add_model(force_model);
        Arc::new(me)
    }

    /// Initialize new spacecraft dynamics with a vector of force models.
    pub fn with_models(
        orbital_dyn: Arc<OrbitalDynamics<'a>>,
        force_models: Vec<Arc<dyn ForceModel + 'a>>,
    ) -> Arc<Self> {
        let mut me = Self::new_raw(orbital_dyn);
        me.force_models = force_models;
        Arc::new(me)
    }

    pub fn add_model(&mut self, force_model: Arc<dyn ForceModel + 'a>) {
        self.force_models.push(force_model);
    }

    /// A shortcut to spacecraft.ctrl if the control is defined
    pub fn ctrl_achieved(&self, state: &Spacecraft) -> Result<bool, NyxError> {
        match &self.ctrl {
            Some(ctrl) => ctrl.achieved(state),
            None => Err(NyxError::NoObjectiveDefined),
        }
    }
}

impl<'a> Dynamics for SpacecraftDynamics<'a> {
    type HyperdualSize = U7; // This implies that we are NOT computing the STM with the mass
    type StateType = Spacecraft;

    fn finally(&self, next_state: Self::StateType) -> Result<Self::StateType, NyxError> {
        if next_state.fuel_mass_kg < 0.0 {
            error!("negative fuel mass at {}", next_state.epoch());
            return Err(NyxError::FuelExhausted);
        }

        if let Some(ctrl) = &self.ctrl {
            let mut state = next_state;
            // Update the control mode
            state.mode = ctrl.next(&state);
            Ok(state)
        } else {
            Ok(next_state)
        }
    }

    fn eom(
        &self,
        delta_t: f64,
        state: &VectorN<f64, U43>,
        ctx: &Spacecraft,
    ) -> Result<VectorN<f64, U43>, NyxError> {
        // Compute the orbital dynamics
        let orbital_dyn_vec = state.fixed_rows::<U42>(0).into_owned();
        let d_x_orbital_dyn = self
            .orbital_dyn
            .eom(delta_t, &orbital_dyn_vec, &ctx.orbit)?;
        // Note: 0.0 is the current fuel usage at this point.
        let mut d_x = VectorN::<f64, U43>::from_iterator(
            d_x_orbital_dyn
                .iter()
                .chain(Vector1::new(0.0).iter())
                .cloned(),
        );

        let mut total_mass = ctx.dry_mass_kg;

        // Rebuild the osculating state for the EOM context.
        let osc_sc = ctx.ctor_from(delta_t, state);

        // Now include the control as needed.
        if let Some(ctrl) = &self.ctrl {
            let (thrust_force, fuel_rate) = {
                if osc_sc.thruster.is_none() {
                    return Err(NyxError::CtrlExistsButNoThrusterAvail);
                }
                let thruster = osc_sc.thruster.unwrap();
                let thrust_power = ctrl.throttle(&osc_sc);
                if !(0.0..=1.0).contains(&thrust_power) {
                    return Err(NyxError::CtrlThrottleRangeErr(thrust_power));
                } else if thrust_power > 0.0 {
                    // Thrust arc
                    let thrust_inertial = ctrl.direction(&osc_sc);
                    if (thrust_inertial.norm() - 1.0).abs() > NORM_ERR {
                        return Err(NyxError::CtrlNotAUnitVector(thrust_inertial.norm()));
                    }
                    // Compute the thrust in Newtons and Isp
                    let total_thrust = (thrust_power * thruster.thrust) * 1e-3; // Convert m/s^-2 to km/s^-2
                    (
                        thrust_inertial * total_thrust,
                        if self.decrement_mass {
                            let fuel_usage =
                                thrust_power * thruster.thrust / (thruster.isp * STD_GRAVITY);
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
            total_mass += ctx.fuel_mass_kg + fuel_rate;
            for i in 0..3 {
                d_x[i + 3] += thrust_force[i] / (ctx.dry_mass_kg + state[U43::dim() - 1]);
            }
            d_x[U43::dim() - 1] += fuel_rate;
        }

        // Compute additional force models as needed.
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

        if self.ctrl.is_some() {
            return Err(NyxError::PartialsUndefined);
        }

        // Call the EOMs
        let total_mass = ctx.dry_mass_kg;
        let radius = state_vec.fixed_rows::<U3>(0).into_owned();

        // Recreate the osculating state.
        let mut osc_sc = *ctx;
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
