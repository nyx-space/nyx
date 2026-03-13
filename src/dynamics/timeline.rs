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

use std::collections::BTreeMap;

use anise::{
    frames::Frame,
    structure::spacecraft::{DragData, Inertia, Mass, SRPData},
};
use hifitime::Epoch;
use indexmap::IndexMap;

use crate::{
    dynamics::{
        guidance::{mnvr::ImpulsiveManeuver, Maneuver, Thruster},
        Drag, PointMasses, SolarPressure,
    },
    io::gravity::GravityFieldConfig,
    propagators::{IntegratorMethod, IntegratorOptions},
};

#[derive(Clone, Debug)]
pub struct SpacecraftTimeline {
    pub timeline: BTreeMap<Epoch, TimelinePhase>,
    pub thruster_sets: IndexMap<String, Thruster>,
    pub propagators: IndexMap<String, PropagatorConfig>,
}

impl SpacecraftTimeline {
    pub fn validate(&self) -> Result<(), String> {
        // Check that the last statement is a terminate
        if let Some((_, TimelinePhase::Phase { .. })) = self.timeline.iter().last() {
            return Err("final phase must be a Terminate".into());
        }

        // Check that all of the thruster set indexes reference an available thruster
        for (epoch, phase) in &self.timeline {
            if let TimelinePhase::Phase {
                name: _,
                propagator,
                guidance,
                on_entry: _,
            } = phase
            {
                // Check that the propagator exists
                if self.propagators.get(propagator).is_none() {
                    return Err(format!("{epoch}: no propagator named `{propagator}`"));
                }
                if let Some(guidance) = guidance {
                    let thruster = guidance.thruster_model();
                    if self.thruster_sets.get(thruster).is_none() {
                        return Err(format!("{epoch}: no thruster set named {thruster}"));
                    }
                }
            }
        }

        Ok(())
    }
}

#[derive(Clone, Debug)]
pub enum TimelinePhase {
    Terminate,
    Phase {
        name: String,
        propagator: String,
        guidance: Option<Box<GuidanceConfig>>,
        /// The discrete event will be applied ONCE before the equation of motions are integrated.
        on_entry: Option<Box<DiscreteEvent>>,
    },
}

/// Propagator config includes the method, options, and all dynamics
#[derive(Clone, Debug)]
pub struct PropagatorConfig {
    pub method: IntegratorMethod,
    pub options: IntegratorOptions,
    pub accel_models: AccelModels,
    pub force_models: ForceModels,
}

#[derive(Clone, Debug)]
pub struct AccelModels {
    pub point_masses: Option<PointMasses>,
    pub gravity_field: Option<GravityFieldConfig>,
}

#[derive(Clone, Debug)]
pub struct ForceModels {
    pub solar_pressure: Option<SolarPressure>,
    pub drag: Option<Drag>,
}

#[derive(Clone, Debug)]
pub enum GuidanceConfig {
    FiniteBurn {
        maneuver: Maneuver,
        thruster_model: String,
    },
    // TODO: Enable config
    Ruggiero {
        thruster_model: String,
    },
    Kluever {
        thruster_model: String,
    },
}

impl GuidanceConfig {
    pub fn thruster_model(&self) -> &str {
        match self {
            Self::FiniteBurn { thruster_model, .. } => thruster_model,
            Self::Ruggiero { thruster_model } => thruster_model,
            Self::Kluever { thruster_model } => thruster_model,
        }
    }
}

#[derive(Clone, Debug)]
pub enum DiscreteEvent {
    Staging {
        impulsive_maneuver: Option<ImpulsiveManeuver>,
        decrement_properties: Option<PhysicalProperties>,
    },
    Docking {
        impulsive_maneuver: Option<ImpulsiveManeuver>,
        increment_properties: Option<PhysicalProperties>,
    },
    CentralBodySwap {
        new_central_body: Frame,
    },
}

#[derive(Clone, Debug)]
pub struct PhysicalProperties {
    pub mass: Option<Mass>,
    pub srp: Option<SRPData>,
    pub drag: Option<DragData>,
    pub inertia: Option<Inertia>,
}
