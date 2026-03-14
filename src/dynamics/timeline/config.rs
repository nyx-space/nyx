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

use crate::{
    dynamics::{guidance::Maneuver, Drag, PointMasses, SolarPressure},
    io::gravity::GravityFieldConfig,
    propagators::{IntegratorMethod, IntegratorOptions},
};

use crate::dynamics::timeline::DiscreteEvent;

#[derive(Clone, Debug)]
pub enum TimelinePhase {
    Terminate,
    Phase {
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
#[derive(Clone, Debug)]
pub struct PropagatorConfig {
    pub method: IntegratorMethod,
    pub options: IntegratorOptions,
    pub accel_models: AccelModels,
    pub force_models: ForceModels,
}

/// Acceleration models alter the orbital dynamics
#[derive(Clone, Debug)]
pub struct AccelModels {
    pub point_masses: Option<PointMasses>,
    pub gravity_field: Option<(GravityFieldConfig, FrameUid)>,
}

/// Force models alter the spacecraft dynamics (they need a mass).
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
