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

use anise::{
    frames::Frame,
    structure::spacecraft::{DragData, Inertia, Mass, SRPData},
};

use crate::{
    dynamics::{guidance::mnvr::ImpulsiveManeuver, Drag, PointMasses, SolarPressure},
    io::gravity::GravityFieldConfig,
    propagators::{IntegratorMethod, IntegratorOptions},
};

#[derive(Clone, Debug)]
pub enum TimelinePhase {
    Terminate,
    Phase {
        name: String,
        propagator: Box<PropagatorConfig>,
        guidance: Option<Box<GuidanceConfig>>,
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
pub enum GuidanceConfig {}

#[derive(Clone, Debug)]
pub enum DiscreteEvent {
    Staging {
        impulsive_maneuver: Option<ImpulsiveManeuver>,
        decrement_properties: Option<PhysicalProperties>,
    },
    Docking {
        impulsive_maneuver: Option<ImpulsiveManeuver>,
        decrement_properties: Option<PhysicalProperties>,
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
