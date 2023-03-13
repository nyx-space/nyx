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
use crate::{
    cosmic::{DragConfig, SrpConfig},
    dynamics::guidance::Thruster,
    io::{ConfigError, ConfigRepr},
    md::ui::GuidanceMode,
    Orbit, Spacecraft, State,
};
use crate::{md::StateParameter, NyxError};
use std::collections::HashMap;

use hifitime::Epoch;
use pyo3::prelude::*;

#[pymethods]
impl Spacecraft {
    /// Initialize a new Spacecraft with optional thruster, mode, SRP, and Drag parameters.
    #[new]
    pub fn py_new(
        orbit: Orbit,
        dry_mass_kg: f64,
        fuel_mass_kg: f64,
        srp_area_m2: f64,
        drag_area_m2: f64,
        cr: f64,
        cd: f64,
        thruster: Option<Thruster>,
        mode: Option<GuidanceMode>,
    ) -> Self {
        Self {
            orbit,
            dry_mass_kg,
            fuel_mass_kg,
            thruster,
            mode: mode.or_else(|| Some(GuidanceMode::Coast)).unwrap(),
            stm: None,
            srp: SrpConfig {
                area_m2: srp_area_m2,
                cr,
            },
            drag: DragConfig {
                area_m2: drag_area_m2,
                cd,
            },
        }
    }

    #[staticmethod]
    fn load_yaml(path: &str) -> Result<Self, ConfigError> {
        <Self as ConfigRepr>::load_yaml(path)
    }

    #[staticmethod]
    fn load_many_yaml(path: &str) -> Result<Vec<Self>, ConfigError> {
        <Self as ConfigRepr>::load_many_yaml(path)
    }

    #[staticmethod]
    fn load_named_yaml(path: &str) -> Result<HashMap<String, Self>, ConfigError> {
        <Self as ConfigRepr>::load_named_yaml(path)
    }

    fn __repr__(&self) -> String {
        format!("{self:?}")
    }

    fn __str__(&self) -> String {
        format!("{self}\n{self:x}")
    }

    /// Note: this returns a COPY of the orbit, not a mutable reference to it!
    #[getter]
    fn get_orbit(&self) -> Orbit {
        self.orbit
    }

    #[setter(orbit)]
    fn py_set_orbit(&mut self, orbit: Orbit) -> PyResult<()> {
        self.orbit = orbit;
        Ok(())
    }

    /// Note: this returns a COPY of the epoch, not a mutable reference to it!
    #[getter]
    fn get_epoch(&self) -> Epoch {
        self.orbit.epoch
    }

    /// Returns the value of the provided state parameter if available
    #[cfg(feature = "python")]
    fn value_of(&self, param: &StateParameter) -> Result<f64, NyxError> {
        self.value(param)
    }
}
