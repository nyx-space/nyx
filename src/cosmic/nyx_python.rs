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

use crate::cosmic::GuidanceMode;
use crate::Spacecraft;
use anise::prelude::Orbit;
// use anise::structure::spacecraft::{DragData, Mass, SRPData};
use pyo3::prelude::*;

#[pymethods]
impl Spacecraft {
    #[new]
    fn py_new(
        orbit: Orbit,
        // TODO: Implement Python classes for all following elements
        // mass: Option<Mass>,
        // srp: Option<SRPData>,
        // drag: Option<DragData>,
        // thruster: Option<Thruster>,
        mode: Option<GuidanceMode>,
    ) -> Self {
        let builder = Spacecraft::builder().orbit(orbit);
        if let Some(mode) = mode {
            builder.mode(mode).build()
        } else {
            builder.build()
        }
    }
}
