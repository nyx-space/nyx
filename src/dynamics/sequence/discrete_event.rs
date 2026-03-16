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
    structure::spacecraft::{DragData, Mass, SRPData},
};
use serde::{Deserialize, Serialize};
use serde_dhall::{SimpleType, StaticType};
use std::collections::HashMap;

use crate::dynamics::guidance::mnvr::ImpulsiveManeuver;

#[derive(Clone, Debug, Serialize, Deserialize, StaticType)]
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

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PhysicalProperties {
    pub mass: Option<Mass>,
    pub srp: Option<SRPData>,
    pub drag: Option<DragData>,
}

impl StaticType for PhysicalProperties {
    fn static_type() -> serde_dhall::SimpleType {
        let mut fields = HashMap::new();
        // TODO: Switch to ANISE's type once released.

        #[allow(dead_code)]
        #[derive(StaticType)]
        struct MassClone {
            dry_mass_kg: f64,
            prop_mass_kg: f64,
            extra_mass_kg: f64,
        }

        #[allow(dead_code)]
        #[derive(StaticType)]
        struct Srp {
            area_m2: f64,
            coeff_reflectivity: f64,
        }

        #[allow(dead_code)]
        #[derive(StaticType)]
        struct Drag {
            area_m2: f64,
            coeff_drag: f64,
        }

        fields.insert("mass".to_string(), MassClone::static_type());
        fields.insert("srp".to_string(), Srp::static_type());
        fields.insert("drag".to_string(), Drag::static_type());

        SimpleType::Record(fields)
    }
}
