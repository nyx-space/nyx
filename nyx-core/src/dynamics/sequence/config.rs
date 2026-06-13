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
use anise::frames::FrameUid;
use anise::prelude::Almanac;
use serde::{Deserialize, Serialize};
use serde_dhall::{SimpleType, StaticType};
use std::collections::HashMap;
use std::sync::Arc;

use crate::dynamics::SpacecraftDynamics;
use crate::dynamics::{GravityField, OrbitalDynamics};
use crate::io::gravity::GravityFieldData;
use crate::propagators::Propagator;
use crate::{
    dynamics::{
        Drag, PointMasses, SolarPressure,
        guidance::{Maneuver, ObjectiveEfficiency, ObjectiveWeight},
    },
    io::gravity::GravityFieldConfig,
    propagators::{IntegratorMethod, IntegratorOptions},
};

use crate::dynamics::sequence::discrete_event::DiscreteEvent;

#[cfg(feature = "python")]
use pyo3::prelude::*;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum Phase {
    Terminate,
    Activity {
        name: String,
        propagator: String,
        guidance: Option<Box<GuidanceConfig>>,
        /// The discrete event will be applied ONCE before the equation of motions are integrated.
        on_entry: Option<Box<DiscreteEvent>>,
        /// Allows disabling a phase without removing it
        disabled: bool,
    },
}

impl StaticType for Phase {
    fn static_type() -> SimpleType {
        let mut variants = HashMap::new();

        // Handle the Vector(Vector3<f64>) variant
        // Most math crates serialize Vector3 as a list of 3 doubles
        variants.insert("Terminate".to_string(), None);

        //  Activity variant (Record variant)
        let mut activity_fields = HashMap::new();

        activity_fields.insert("name".to_string(), String::static_type());
        activity_fields.insert("propagator".to_string(), String::static_type());

        // Use the StaticType impl of the boxed inner types
        activity_fields.insert(
            "guidance".to_string(),
            <Option<GuidanceConfig> as StaticType>::static_type(),
        );

        activity_fields.insert(
            "on_entry".to_string(),
            <Option<DiscreteEvent> as StaticType>::static_type(),
        );

        activity_fields.insert("disabled".to_string(), bool::static_type());

        variants.insert(
            "Activity".to_string(),
            Some(SimpleType::Record(activity_fields)),
        );

        SimpleType::Union(variants)
    }
}

/// Propagator config includes the method, options, and all dynamics
#[derive(Clone, Debug, Serialize, Deserialize, StaticType)]
#[cfg_attr(feature = "python", pyclass(from_py_object))]
pub struct PropagatorConfig {
    pub accel_models: AccelModels,
    pub force_models: ForceModels,
    pub method: IntegratorMethod,
    pub options: IntegratorOptions,
}

impl PropagatorConfig {
    pub fn build(&self, almanac: Arc<Almanac>) -> Result<Propagator<SpacecraftDynamics>, String> {
        // Build the orbital dynamics
        let mut orbital_dyn = OrbitalDynamics::two_body();
        if let Some(point_masses) = &self.accel_models.point_masses {
            orbital_dyn.accel_models.push(point_masses.clone());
        }
        if let Some((gravity_cfg, frame_uid)) = &self.accel_models.gravity_field {
            let grav_data =
                GravityFieldData::from_config(gravity_cfg.clone()).map_err(|e| e.to_string())?;
            let compute_frame = almanac.frame_info(*frame_uid).map_err(|e| e.to_string())?;
            let gravity_field = GravityField::from_stor(compute_frame, grav_data);
            orbital_dyn.accel_models.push(gravity_field);
        }
        // Build the spacecraft dynamics
        let mut sc_dyn = SpacecraftDynamics::new(orbital_dyn);

        if let Some(srp) = &self.force_models.solar_pressure {
            sc_dyn.force_models.push(srp.clone());
        }

        if let Some(drag) = &self.force_models.drag {
            sc_dyn.force_models.push(drag.clone());
        }

        // And set it all up!
        Ok(Propagator::new(sc_dyn, self.method, self.options))
    }
}

/// Acceleration models alter the orbital dynamics
#[derive(Clone, Serialize, Deserialize, Debug)]
#[cfg_attr(feature = "python", pyclass(from_py_object))]
pub struct AccelModels {
    pub point_masses: Option<Arc<PointMasses>>,
    pub gravity_field: Option<(GravityFieldConfig, FrameUid)>,
}

/// Force models alter the spacecraft dynamics (they need a mass).
#[derive(Clone, Serialize, Deserialize, Debug)]
#[cfg_attr(feature = "python", pyclass(from_py_object))]
pub struct ForceModels {
    pub solar_pressure: Option<Arc<SolarPressure>>,
    pub drag: Option<Arc<Drag>>,
}

#[derive(Clone, Debug, Serialize, Deserialize, StaticType)]
pub struct GuidanceConfig {
    pub thruster_model: String,
    pub disable_prop_mass: bool,
    pub law: SteeringLaw,
}

// NOTE: Steering laws are not yet available in Python =(
#[derive(Clone, Debug, Serialize, Deserialize, StaticType)]
pub enum SteeringLaw {
    FiniteBurn(Maneuver),
    Kluever {
        /// Stores the objectives, and their associated weights (set to zero to disable).
        objectives: Vec<ObjectiveWeight>,
        /// If defined, coast until vehicle is out of the provided eclipse state.
        max_eclipse_prct: Option<f64>,
    },
    Ruggiero {
        /// Stores the objectives, and their associated efficiency threshold (set to zero if not minimum efficiency).
        objectives: Vec<ObjectiveEfficiency>,
        /// If defined, coast until vehicle is out of the provided eclipse state.
        max_eclipse_prct: Option<f64>,
    },
}

impl StaticType for AccelModels {
    fn static_type() -> serde_dhall::SimpleType {
        let mut fields = HashMap::new();

        fields.insert(
            "point_masses".to_string(),
            SimpleType::Optional(Box::new(PointMasses::static_type())),
        );

        #[allow(dead_code)]
        #[derive(StaticType)]
        struct GravityFieldDhall(GravityFieldConfig, FrameUid);

        fields.insert(
            "gravity_field".to_string(),
            SimpleType::Optional(Box::new(GravityFieldDhall::static_type())),
        );

        SimpleType::Record(fields)
    }
}

impl StaticType for ForceModels {
    fn static_type() -> serde_dhall::SimpleType {
        let mut fields = HashMap::new();

        fields.insert(
            "solar_pressure".to_string(),
            SimpleType::Optional(Box::new(SolarPressure::static_type())),
        );

        fields.insert(
            "drag".to_string(),
            SimpleType::Optional(Box::new(Drag::static_type())),
        );

        SimpleType::Record(fields)
    }
}
