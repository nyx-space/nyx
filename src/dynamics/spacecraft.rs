/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::guidance::GuidanceLaw;
use super::orbital::OrbitalDynamics;
use super::{AccelModel, Dynamics, ForceModel};
pub use crate::cosmic::{GuidanceMode, Spacecraft, STD_GRAVITY};
use crate::errors::NyxError;
use crate::io::dynamics::DynamicsSerde;
use crate::io::gravity::HarmonicsMem;

use crate::io::{ConfigError, Configurable};
use crate::linalg::{Const, DimName, OMatrix, OVector, Vector3};
pub use crate::md::prelude::SolarPressure;
use crate::md::prelude::{Harmonics, PointMasses};
use crate::State;

use std::fmt::{self, Write};
use std::sync::Arc;

use crate::cosmic::Cosm;
#[cfg(feature = "python")]
use crate::io::ConfigRepr;
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

    /// Add a model to the currently defined spacecraft dynamics
    pub fn add_model(&mut self, force_model: Arc<dyn ForceModel>) {
        self.force_models.push(force_model);
    }

    /// Clone these dynamics and add a model to the currently defined orbital dynamics
    pub fn with_model(self, force_model: Arc<dyn ForceModel>) -> Self {
        let mut me = self;
        me.add_model(force_model);
        me
    }

    /// A shortcut to spacecraft.guid_law if a guidance law is defined for these dynamics
    pub fn guidance_achieved(&self, state: &Spacecraft) -> Result<bool, NyxError> {
        match &self.guid_law {
            Some(guid_law) => guid_law.achieved(state),
            None => Err(NyxError::NoObjectiveDefined),
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

    /// Clone these spacecraft dynamics and update the control to the one provided.
    pub fn with_guidance_law_no_decr(&self, guid_law: Arc<dyn GuidanceLaw>) -> Self {
        Self {
            orbital_dyn: self.orbital_dyn.clone(),
            guid_law: Some(guid_law),
            force_models: self.force_models.clone(),
            decrement_mass: false,
        }
    }

    /// Clone these spacecraft dynamics and remove any control model
    pub fn without_guidance_law(&self) -> Self {
        Self {
            orbital_dyn: self.orbital_dyn.clone(),
            guid_law: None,
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
    fn __richcmp__(&self, other: &Self, op: CompareOp) -> Result<bool, NyxError> {
        match op {
            CompareOp::Eq => Ok(self.__repr__() == other.__repr__()),
            CompareOp::Ne => Ok(self.__repr__() != other.__repr__()),
            _ => Err(NyxError::CustomError(format!("{op:?} not available"))),
        }
    }

    #[cfg(feature = "python")]
    #[classmethod]
    /// Loads the SpacecraftDynamics from its YAML representation
    fn loads(_cls: &PyType, state: &PyAny) -> Result<Self, ConfigError> {
        <Self as Configurable>::from_config(
            depythonize(state).map_err(|e| ConfigError::InvalidConfig(e.to_string()))?,
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

    fn finally(&self, next_state: Self::StateType) -> Result<Self::StateType, NyxError> {
        if next_state.fuel_mass_kg < 0.0 {
            error!("negative fuel mass at {}", next_state.epoch());
            return Err(NyxError::FuelExhausted(Box::new(next_state)));
        }

        if let Some(guid_law) = &self.guid_law {
            let mut state = next_state;
            // Update the control mode
            guid_law.next(&mut state);
            Ok(state)
        } else {
            Ok(next_state)
        }
    }

    fn eom(
        &self,
        delta_t: f64,
        state: &OVector<f64, Const<90>>,
        ctx: &Self::StateType,
    ) -> Result<OVector<f64, Const<90>>, NyxError> {
        // Rebuild the osculating state for the EOM context.
        let osc_sc = ctx.set_with_delta_seconds(delta_t, state);
        let mut d_x = OVector::<f64, Const<90>>::zeros();

        if ctx.orbit.stm.is_some() {
            // Call the gradient (also called the dual EOM function of the force models)
            let (state, grad) = self.dual_eom(delta_t, &osc_sc)?;

            // Apply the gradient to the STM
            let stm_dt = ctx.stm()? * grad;

            // Rebuild the state vectors
            for (i, val) in state.iter().enumerate() {
                d_x[i] = *val;
            }

            for (i, val) in stm_dt.iter().enumerate() {
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
        if let Some(guid_law) = &self.guid_law {
            let (thrust_force, fuel_rate) = {
                if osc_sc.thruster.is_none() {
                    return Err(NyxError::NoThrusterAvail);
                }
                let thruster = osc_sc.thruster.unwrap();
                let thrust_throttle_lvl = guid_law.throttle(&osc_sc);
                if !(0.0..=1.0).contains(&thrust_throttle_lvl) {
                    return Err(NyxError::CtrlThrottleRangeErr(thrust_throttle_lvl));
                } else if thrust_throttle_lvl > 0.0 {
                    // Thrust arc
                    let thrust_inertial = guid_law.direction(&osc_sc);
                    if (thrust_inertial.norm() - 1.0).abs() > NORM_ERR {
                        return Err(NyxError::CtrlNotAUnitVector(thrust_inertial.norm()));
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
    ) -> Result<(OVector<f64, Const<9>>, OMatrix<f64, Const<9>, Const<9>>), NyxError> {
        // Rebuild the appropriately sized state and STM.
        // This is the orbital state followed by Cr and Cd
        let mut d_x = OVector::<f64, Const<9>>::zeros();
        let mut grad = OMatrix::<f64, Const<9>, Const<9>>::zeros();

        let (orb_state, orb_grad) = self.orbital_dyn.dual_eom(delta_t_s, &ctx.orbit)?;

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
            return Err(NyxError::PartialsUndefined);
        }

        // Call the EOMs
        let total_mass = ctx.mass_kg();
        for model in &self.force_models {
            let (model_frc, model_grad) = model.dual_eom(ctx)?;
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

impl Configurable for SpacecraftDynamics {
    type IntermediateRepr = DynamicsSerde;

    fn from_config(cfg: Self::IntermediateRepr, cosm: Arc<Cosm>) -> Result<Self, ConfigError>
    where
        Self: Sized,
    {
        // Builds the orbital dynamics
        let mut accel_models: Vec<Arc<dyn AccelModel + Sync>> = Vec::new();

        accel_models.push(PointMasses::new(&cfg.point_masses, cosm.clone()));

        if let Some(harmonics_serde) = cfg.harmonics {
            for hh in harmonics_serde {
                let gunzipped = hh.coeffs.ends_with(".gz");
                let stor = if hh.coeffs.contains("cof") {
                    HarmonicsMem::from_cof(&hh.coeffs, hh.degree, hh.order, gunzipped)
                        .map_err(|e| ConfigError::InvalidConfig(e.to_string()))?
                } else if hh.coeffs.contains("sha") {
                    HarmonicsMem::from_shadr(&hh.coeffs, hh.degree, hh.order, gunzipped)
                        .map_err(|e| ConfigError::InvalidConfig(e.to_string()))?
                } else if hh.coeffs.contains("EGM") {
                    HarmonicsMem::from_egm(&hh.coeffs, hh.degree, hh.order, gunzipped)
                        .map_err(|e| ConfigError::InvalidConfig(e.to_string()))?
                } else {
                    return Err(ConfigError::InvalidConfig(
                        "Unknown coefficients file type".to_string(),
                    ));
                };

                // Grab the frame
                let frame = cosm
                    .try_frame(&hh.frame)
                    .map_err(|e| ConfigError::InvalidConfig(e.to_string()))?;

                accel_models.push(Harmonics::from_stor(frame, stor, cosm.clone()));
            }
        }

        let orbital_dyn = OrbitalDynamics::new(accel_models);

        let mut force_models: Vec<Arc<dyn ForceModel>> = Vec::new();

        // SRP
        if let Some(srp) = cfg.srp {
            force_models.push(SolarPressure::with_flux(
                srp.phi.map_or(1367.0, |v| v),
                srp.shadows,
                cosm,
            ));
        }

        // TODO: Drag -- https://github.com/nyx-space/nyx/issues/86

        Ok(SpacecraftDynamics::from_models(orbital_dyn, force_models))
    }

    fn to_config(&self) -> Result<Self::IntermediateRepr, ConfigError> {
        todo!()
    }
}
