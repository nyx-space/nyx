/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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

use anise::prelude::Almanac;
use snafu::ResultExt;

use super::guidance::{ra_dec_from_unit_vector, GuidanceError, GuidanceLaw};
use super::orbital::OrbitalDynamics;
use super::{Dynamics, DynamicsGuidanceSnafu, ForceModel};
pub use crate::cosmic::{GuidanceMode, Spacecraft, STD_GRAVITY};
use crate::dynamics::DynamicsError;

use crate::linalg::{Const, DimName, OMatrix, OVector, Vector3};
pub use crate::md::prelude::SolarPressure;
use crate::State;

use std::fmt::{self, Write};
use std::sync::Arc;

use crate::cosmic::AstroError;
#[cfg(feature = "python")]
use crate::io::ConfigRepr;
#[cfg(feature = "python")]
use crate::python::PythonError;
#[cfg(feature = "python")]
use pyo3::class::basic::CompareOp;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyType;
#[cfg(feature = "python")]
use pythonize::depythonize;
#[cfg(feature = "python")]
use std::collections::BTreeMap;

const NORM_ERR: f64 = 1e-4;

/// A generic spacecraft dynamics with associated force models, guidance law, and flag specifying whether to decrement the fuel mass or not.
/// Note: when developing new guidance laws, it is recommended to _not_ enable fuel decrement until the guidance law seems to work without proper physics.
/// Note: if the spacecraft runs out of fuel, the propagation segment will return an error.
#[derive(Clone)]
#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(module = "nyx_space.mission_design"))]
pub struct SpacecraftDynamics {
    pub orbital_dyn: OrbitalDynamics,
    // TODO: https://github.com/nyx-space/nyx/issues/214
    pub force_models: Vec<Arc<dyn ForceModel>>,
    pub guid_law: Option<Arc<dyn GuidanceLaw>>,
    pub decrement_mass: bool,
}

impl SpacecraftDynamics {
    /// Initialize a Spacecraft with a set of orbital dynamics and a propulsion subsystem.
    /// By default, the mass of the vehicle will be decremented as propellant is consumed.
    pub fn from_guidance_law(orbital_dyn: OrbitalDynamics, guid_law: Arc<dyn GuidanceLaw>) -> Self {
        Self {
            orbital_dyn,
            guid_law: Some(guid_law),
            force_models: Vec::new(),
            decrement_mass: true,
        }
    }

    /// Initialize a Spacecraft with a set of orbital dynamics and a propulsion subsystem.
    /// Will _not_ decrement the fuel mass as propellant is consumed.
    pub fn from_guidance_law_no_decr(
        orbital_dyn: OrbitalDynamics,
        guid_law: Arc<dyn GuidanceLaw>,
    ) -> Self {
        Self {
            orbital_dyn,
            guid_law: Some(guid_law),
            force_models: Vec::new(),
            decrement_mass: false,
        }
    }

    /// Initialize a Spacecraft with a set of orbital dynamics and with SRP enabled.
    pub fn new(orbital_dyn: OrbitalDynamics) -> Self {
        Self {
            orbital_dyn,
            guid_law: None,
            force_models: Vec::new(),
            decrement_mass: true,
        }
    }

    /// Initialize new spacecraft dynamics with the provided orbital mechanics and with the provided force model.
    pub fn from_model(orbital_dyn: OrbitalDynamics, force_model: Arc<dyn ForceModel>) -> Self {
        Self {
            orbital_dyn,
            guid_law: None,
            force_models: vec![force_model],
            decrement_mass: true,
        }
    }

    /// Initialize new spacecraft dynamics with a vector of force models.
    pub fn from_models(
        orbital_dyn: OrbitalDynamics,
        force_models: Vec<Arc<dyn ForceModel>>,
    ) -> Self {
        let mut me = Self::new(orbital_dyn);
        me.force_models = force_models;
        me
    }

    /// A shortcut to spacecraft.guid_law if a guidance law is defined for these dynamics
    pub fn guidance_achieved(&self, state: &Spacecraft) -> Result<bool, GuidanceError> {
        match &self.guid_law {
            Some(guid_law) => guid_law.achieved(state),
            None => Err(GuidanceError::NoGuidanceObjectiveDefined),
        }
    }

    /// Clone these spacecraft dynamics and update the control to the one provided.
    pub fn with_guidance_law(&self, guid_law: Arc<dyn GuidanceLaw>) -> Self {
        Self {
            orbital_dyn: self.orbital_dyn.clone(),
            guid_law: Some(guid_law),
            force_models: self.force_models.clone(),
            decrement_mass: self.decrement_mass,
        }
    }
}

#[cfg_attr(feature = "python", pymethods)]
impl SpacecraftDynamics {
    #[cfg(feature = "python")]
    #[classmethod]
    fn load(_cls: &PyType, path: &str) -> Result<Self, ConfigError> {
        let serde = DynamicsSerde::load(path)?;

        let cosm = Cosm::de438();

        Self::from_config(serde, cosm)
    }

    #[cfg(feature = "python")]
    #[classmethod]
    fn load_many(_cls: &PyType, path: &str) -> Result<Vec<Self>, ConfigError> {
        let orbits = DynamicsSerde::load_many(path)?;

        let cosm = Cosm::de438();

        let mut selves = Vec::with_capacity(orbits.len());

        for serde in orbits {
            selves.push(Self::from_config(serde, cosm.clone())?);
        }

        Ok(selves)
    }

    #[cfg(feature = "python")]
    #[classmethod]
    fn load_named(_cls: &PyType, path: &str) -> Result<BTreeMap<String, Self>, ConfigError> {
        let orbits = DynamicsSerde::load_named(path)?;

        let cosm = Cosm::de438();

        let mut selves = BTreeMap::new();

        for (k, v) in orbits {
            selves.insert(k, Self::from_config(v, cosm.clone())?);
        }

        Ok(selves)
    }

    #[cfg(feature = "python")]
    fn __repr__(&self) -> String {
        format!("{self}")
    }

    #[cfg(feature = "python")]
    fn __richcmp__(&self, other: &Self, op: CompareOp) -> Result<bool, PythonError> {
        match op {
            CompareOp::Eq => Ok(self.__repr__() == other.__repr__()),
            CompareOp::Ne => Ok(self.__repr__() != other.__repr__()),
            _ => Err(PythonError::OperationError { op }),
        }
    }

    #[cfg(feature = "python")]
    #[classmethod]
    /// Loads the SpacecraftDynamics from its YAML representation
    fn loads(_cls: &PyType, state: &PyAny) -> Result<Self, ConfigError> {
        <Self as Configurable>::from_config(
            depythonize(state).map_err(|e| ConfigError::InvalidConfig { msg: e.to_string() })?,
            Cosm::de438(),
        )
    }
}

impl fmt::Display for SpacecraftDynamics {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let force_models: String = if self.force_models.is_empty() {
            "No force models;".to_string()
        } else {
            self.force_models
                .iter()
                .fold(String::new(), |mut output, x| {
                    let _ = write!(output, "{x}; ");
                    output
                })
        };
        write!(
            f,
            "Spacecraft dynamics (with guidance = {}): {} {}",
            self.guid_law.is_some(),
            force_models,
            self.orbital_dyn
        )
    }
}

impl Dynamics for SpacecraftDynamics {
    type HyperdualSize = Const<9>;
    type StateType = Spacecraft;

    fn finally(
        &self,
        next_state: Self::StateType,
        almanac: Arc<Almanac>,
    ) -> Result<Self::StateType, DynamicsError> {
        if next_state.fuel_mass_kg < 0.0 {
            error!("negative fuel mass at {}", next_state.epoch());
            return Err(DynamicsError::FuelExhausted {
                sc: Box::new(next_state),
            });
        }

        if let Some(guid_law) = &self.guid_law {
            let mut state = next_state;
            // Update the control mode
            guid_law.next(&mut state, almanac.clone());
            Ok(state)
        } else {
            Ok(next_state)
        }
    }

    fn eom(
        &self,
        delta_t_s: f64,
        state: &OVector<f64, Const<90>>,
        ctx: &Self::StateType,
        almanac: Arc<Almanac>,
    ) -> Result<OVector<f64, Const<90>>, DynamicsError> {
        // Rebuild the osculating state for the EOM context.
        let osc_sc = ctx.set_with_delta_seconds(delta_t_s, state);
        let mut d_x = OVector::<f64, Const<90>>::zeros();

        // Maybe I use this only when estimating the orbit state from a spacecraft, but that functionality will soon disappear.
        match ctx.stm {
            Some(stm) => {
                // Call the gradient (also called the dual EOM function of the force models)
                let (state, grad) = self.dual_eom(delta_t_s, &osc_sc, almanac)?;

                // Apply the gradient to the STM
                let stm_dt = stm * grad;

                // Rebuild the state vector
                for (i, val) in state.iter().enumerate() {
                    d_x[i] = *val;
                }

                for (i, val) in stm_dt.iter().copied().enumerate() {
                    d_x[i + <Spacecraft as State>::Size::dim()] = val;
                }
            }
            None => {
                // Compute the orbital dynamics
                for (i, val) in self
                    .orbital_dyn
                    .eom(&osc_sc.orbit, almanac.clone())?
                    .iter()
                    .copied()
                    .enumerate()
                {
                    d_x[i] = val;
                }

                // Apply the force models for non STM propagation
                for model in &self.force_models {
                    let model_frc = model.eom(&osc_sc, almanac.clone())? / osc_sc.mass_kg();
                    for i in 0..3 {
                        d_x[i + 3] += model_frc[i];
                    }
                }
            }
        };

        // Now include the control as needed.
        if let Some(guid_law) = &self.guid_law {
            let (thrust_force, fuel_rate) = {
                if osc_sc.thruster.is_none() {
                    return Err(DynamicsError::DynamicsGuidance {
                        source: GuidanceError::NoThrustersDefined,
                    });
                }
                let thruster = osc_sc.thruster.unwrap();
                let thrust_throttle_lvl =
                    guid_law.throttle(&osc_sc).context(DynamicsGuidanceSnafu)?;
                if !(0.0..=1.0).contains(&thrust_throttle_lvl) {
                    return Err(DynamicsError::DynamicsGuidance {
                        source: GuidanceError::ThrottleRatio {
                            ratio: thrust_throttle_lvl,
                        },
                    });
                } else if thrust_throttle_lvl > 0.0 {
                    // Thrust arc
                    let thrust_inertial =
                        guid_law.direction(&osc_sc).context(DynamicsGuidanceSnafu)?;
                    if (thrust_inertial.norm() - 1.0).abs() > NORM_ERR {
                        let (alpha, delta) = ra_dec_from_unit_vector(thrust_inertial);
                        return Err(DynamicsError::DynamicsGuidance {
                            source: GuidanceError::InvalidDirection {
                                x: thrust_inertial[0],
                                y: thrust_inertial[1],
                                z: thrust_inertial[2],
                                in_plane_deg: alpha.to_degrees(),
                                out_of_plane_deg: delta.to_degrees(),
                            },
                        });
                    } else if thrust_inertial.norm().is_normal() {
                        // Compute the thrust in Newtons and Isp
                        let total_thrust = (thrust_throttle_lvl * thruster.thrust_N) * 1e-3; // Convert m/s^-2 to km/s^-2
                        (
                            thrust_inertial * total_thrust,
                            if self.decrement_mass {
                                let fuel_usage = thrust_throttle_lvl * thruster.thrust_N
                                    / (thruster.isp_s * STD_GRAVITY);
                                -fuel_usage
                            } else {
                                0.0
                            },
                        )
                    } else {
                        warn!(
                            "Abnormal thrust direction vector\t|u| = {}",
                            thrust_inertial.norm()
                        );
                        (Vector3::zeros(), 0.0)
                    }
                } else {
                    (Vector3::zeros(), 0.0)
                }
            };

            for i in 0..3 {
                d_x[i + 3] += thrust_force[i] / osc_sc.mass_kg();
            }
            d_x[8] += fuel_rate;
        }
        Ok(d_x)
    }

    fn dual_eom(
        &self,
        delta_t_s: f64,
        ctx: &Self::StateType,
        almanac: Arc<Almanac>,
    ) -> Result<(OVector<f64, Const<9>>, OMatrix<f64, Const<9>, Const<9>>), DynamicsError> {
        // Rebuild the appropriately sized state and STM.
        // This is the orbital state followed by Cr and Cd
        let mut d_x = OVector::<f64, Const<9>>::zeros();
        let mut grad = OMatrix::<f64, Const<9>, Const<9>>::zeros();

        let (orb_state, orb_grad) =
            self.orbital_dyn
                .dual_eom(delta_t_s, &ctx.orbit, almanac.clone())?;

        // Copy the d orbit dt data
        for (i, val) in orb_state.iter().enumerate() {
            d_x[i] = *val;
        }

        for i in 0..6 {
            for j in 0..6 {
                grad[(i, j)] = orb_grad[(i, j)];
            }
        }

        if self.guid_law.is_some() {
            return Err(DynamicsError::DynamicsAstro {
                source: AstroError::PartialsUndefined,
            });
        }

        // Call the EOMs
        let total_mass = ctx.mass_kg();
        for model in &self.force_models {
            let (model_frc, model_grad) = model.dual_eom(ctx, almanac.clone())?;
            for i in 0..3 {
                // Add the velocity changes
                d_x[i + 3] += model_frc[i] / total_mass;
                // Add the velocity partials
                for j in 1..4 {
                    grad[(i + 3, j - 1)] += model_grad[(i, j - 1)] / total_mass;
                }
            }
            // Add this force model's estimation if applicable.
            if let Some(idx) = model.estimation_index() {
                for j in 0..3 {
                    grad[(idx, j)] += model_grad[(3, j)] / total_mass;
                }
                // d_x[idx] = 1.0;
            }
        }

        Ok((d_x, grad))
    }
}
