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
use serde::{Deserialize, Serialize};
use serde_dhall::{SimpleType, StaticType};
use std::collections::HashMap;
use std::sync::Arc;

use crate::{
    dynamics::{guidance::Maneuver, Drag, PointMasses, SolarPressure},
    io::gravity::GravityFieldConfig,
    md::objective::Objective,
    propagators::{IntegratorMethod, IntegratorOptions},
};

use crate::dynamics::sequence::discrete_event::DiscreteEvent;

#[derive(Clone, Debug)]
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

/// Propagator config includes the method, options, and all dynamics
#[derive(Clone, Debug, Serialize, Deserialize, StaticType)]
pub struct PropagatorConfig {
    pub method: IntegratorMethod,
    pub options: IntegratorOptions,
    pub accel_models: AccelModels,
    pub force_models: ForceModels,
}

/// Acceleration models alter the orbital dynamics
#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct AccelModels {
    pub point_masses: Option<Arc<PointMasses>>,
    pub gravity_field: Option<(GravityFieldConfig, FrameUid)>,
}

/// Force models alter the spacecraft dynamics (they need a mass).
#[derive(Clone, Serialize, Deserialize, Debug)]
pub struct ForceModels {
    pub solar_pressure: Option<Arc<SolarPressure>>,
    pub drag: Option<Arc<Drag>>,
}

#[derive(Clone, Debug)]
pub enum GuidanceConfig {
    FiniteBurn {
        thruster_model: String,
        disable_prop_mass: bool,
        maneuver: Maneuver,
    },
    Ruggiero {
        thruster_model: String,
        disable_prop_mass: bool,
        /// Stores the objectives, and their associated efficiency threshold (set to zero if not minimum efficiency).
        objectives: Vec<(Objective, f64)>,
        /// If defined, coast until vehicle is out of the provided eclipse state.
        max_eclipse_prct: Option<f64>,
    },
    Kluever {
        thruster_model: String,
        disable_prop_mass: bool,
        /// Stores the objectives, and their associated weights (set to zero to disable).
        objectives: Vec<(Objective, f64)>,
        /// If defined, coast until vehicle is out of the provided eclipse state.
        max_eclipse_prct: Option<f64>,
    },
}

impl GuidanceConfig {
    pub fn thruster_model(&self) -> &str {
        match self {
            Self::FiniteBurn { thruster_model, .. } => thruster_model,
            Self::Ruggiero { thruster_model, .. } => thruster_model,
            Self::Kluever { thruster_model, .. } => thruster_model,
        }
    }
    pub fn disable_prop_mass(&self) -> bool {
        match self {
            Self::FiniteBurn {
                disable_prop_mass, ..
            } => *disable_prop_mass,
            Self::Ruggiero {
                disable_prop_mass, ..
            } => *disable_prop_mass,
            Self::Kluever {
                disable_prop_mass, ..
            } => *disable_prop_mass,
        }
    }
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
