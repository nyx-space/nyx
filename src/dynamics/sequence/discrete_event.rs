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

use crate::dynamics::guidance::mnvr::ImpulsiveManeuver;

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
    FrameSwap {
        new_frame: Frame,
    },
}

#[derive(Clone, Debug)]
pub struct PhysicalProperties {
    pub mass: Option<Mass>,
    pub srp: Option<SRPData>,
    pub drag: Option<DragData>,
    pub inertia: Option<Inertia>,
}
