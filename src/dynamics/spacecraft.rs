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

use super::hyperdual::{hyperspace_from_vector, Hyperdual};
use super::orbital::OrbitalDynamics;
use super::thrustctrl::ThrustControl;
use super::{Dynamics, ForceModel};
use crate::cosmic::Spacecraft;
use crate::dimensions::{Const, DimName, OMatrix, OVector, Vector3};
use crate::errors::NyxError;

use crate::{State, TimeTagged};
use std::fmt;
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
    pub fn from_model(
        orbital_dyn: Arc<OrbitalDynamics<'a>>,
        force_model: Arc<dyn ForceModel + 'a>,
    ) -> Arc<Self> {
        let mut me = Self::new_raw(orbital_dyn);
        me.add_model(force_model);
        Arc::new(me)
    }

    /// Initialize new spacecraft dynamics with a vector of force models.
    pub fn from_models(
        orbital_dyn: Arc<OrbitalDynamics<'a>>,
        force_models: Vec<Arc<dyn ForceModel + 'a>>,
    ) -> Arc<Self> {
        let mut me = Self::new_raw(orbital_dyn);
        me.force_models = force_models;
        Arc::new(me)
    }

    /// Add a model to the currently defined spacecraft dynamics
    pub fn add_model(&mut self, force_model: Arc<dyn ForceModel + 'a>) {
        self.force_models.push(force_model);
    }

    /// Clone these dynamics and add a model to the currently defined orbital dynamics
    pub fn with_model(self, force_model: Arc<dyn ForceModel + 'a>) -> Arc<Self> {
        let mut me = self.clone();
        me.add_model(force_model);
        Arc::new(me)
    }

    /// A shortcut to spacecraft.ctrl if the control is defined
    pub fn ctrl_achieved(&self, state: &Spacecraft) -> Result<bool, NyxError> {
        match &self.ctrl {
            Some(ctrl) => ctrl.achieved(state),
            None => Err(NyxError::NoObjectiveDefined),
        }
    }
}

impl<'a> fmt::Display for SpacecraftDynamics<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let force_models: String = self
            .force_models
            .iter()
            .map(|x| format!("{}; ", x))
            .collect();
        write!(
            f,
            "Spacecraft dynamics (with ctrl = {}): {}\t{}",
            self.ctrl.is_some(),
            force_models,
            self.orbital_dyn
        )
    }
}

impl<'a> Dynamics for SpacecraftDynamics<'a> {
    type HyperdualSize = Const<9>; // The STM includes Cr and Cd but not the mass (which is assumed constant)
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
        state: &OVector<f64, Const<73>>,
        ctx: &Spacecraft,
    ) -> Result<OVector<f64, Const<73>>, NyxError> {
        // Rebuild the osculating state for the EOM context.
        let osc_sc = ctx.ctor_from(delta_t, state);
        let mut d_x = OVector::<f64, Const<73>>::zeros();

        if ctx.orbit.stm.is_some() {
            let (state, grad) = self.dual_eom(
                delta_t,
                &hyperspace_from_vector(&osc_sc.as_vector()?.fixed_rows::<8>(0).into_owned()),
                &osc_sc,
            )?;

            for (i, val) in state.iter().enumerate() {
                d_x[i] = *val;
            }

            for (i, val) in grad.iter().enumerate() {
                d_x[i + <Spacecraft as State>::Size::dim()] = *val;
            }
        } else {
            // Compute the orbital dynamics
            let orbital_dyn_vec = state.fixed_rows::<42>(0).into_owned();
            // Copy the d orbit dt data
            for (i, val) in self
                .orbital_dyn
                .eom(delta_t, &orbital_dyn_vec, &ctx.orbit)?
                .iter()
                .enumerate()
            {
                d_x[i] = *val;
            }

            // Apply the force models for non STM propagation
            for model in &self.force_models {
                let model_frc = model.eom(&osc_sc)? / osc_sc.mass_kg();
                for i in 0..3 {
                    d_x[i + 3] += model_frc[i];
                }
            }
        }

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

            for i in 0..3 {
                d_x[i + 3] += thrust_force[i] / osc_sc.mass_kg();
            }
            d_x[72] += fuel_rate;
        }
        Ok(d_x)
    }

    fn dual_eom(
        &self,
        delta_t_s: f64,
        state_vec: &OVector<Hyperdual<f64, Self::HyperdualSize>, Const<8>>,
        ctx: &Self::StateType,
    ) -> Result<(OVector<f64, Const<8>>, OMatrix<f64, Const<8>, Const<8>>), NyxError> {
        // Rebuild the appropriately sized state and STM.
        // This is the orbital state followed by Cr and Cd
        let mut d_x = OVector::<f64, Const<8>>::zeros();
        let mut grad = OMatrix::<f64, Const<8>, Const<8>>::zeros();

        // This also means that the STM for the Spacecraft must also be copied into the spacecraft structure
        let one = Const::<1> {};
        let six = Const::<6> {};
        let pos_vel: OVector<Hyperdual<f64, Const<7>>, Const<6>> = OVector::from_iterator_generic(
            six,
            one,
            state_vec
                .fixed_rows::<6>(0)
                .into_iter()
                .map(|x| Hyperdual::from_slice(x.fixed_rows::<7>(0).as_slice())),
        );

        let (orb_state, orb_grad) = self.orbital_dyn.dual_eom(delta_t_s, &pos_vel, &ctx.orbit)?;

        // Copy the d orbit dt data
        for (i, val) in orb_state.iter().enumerate() {
            d_x[i] = *val;
        }

        for i in 0..6 {
            for j in 0..6 {
                grad[(i, j)] = orb_grad[(i, j)];
            }
        }

        if self.ctrl.is_some() {
            return Err(NyxError::PartialsUndefined);
        }

        // Call the EOMs
        let total_mass = ctx.mass_kg();
        let radius = state_vec.fixed_rows::<3>(0).into_owned();

        for model in &self.force_models {
            let (model_frc, model_grad) = model.dual_eom(&radius, ctx)?;
            for i in 0..3 {
                // Add the velocity changes
                d_x[i + 3] += model_frc[i] / total_mass;
                // Add the velocity partials
                for j in 1..4 {
                    grad[(i + 3, j - 1)] += model_grad[(i, j - 1)] / total_mass;
                }
            }
        }

        Ok((d_x, grad))
    }
}
