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
    md::prelude::GuidanceMode,
    Orbit, Spacecraft, State,
};
use crate::{md::StateParameter, NyxError};
use std::collections::HashMap;

use hifitime::Epoch;
use pyo3::class::basic::CompareOp;
use pyo3::prelude::*;
use pyo3::types::PyType;
use pythonize::{depythonize, pythonize};

#[pymethods]
impl Spacecraft {
    /// Initialize a new Spacecraft with optional thruster, mode, SRP, and Drag parameters.
    #[new]
    #[pyo3(
        text_signature = "(orbit, dry_mass_kg, fuel_mass_kg, srp_area_m2, drag_area_m2, cr, cd, thruster, mode)"
    )]
    pub fn py_new(
        orbit: Option<Orbit>,
        dry_mass_kg: Option<f64>,
        fuel_mass_kg: Option<f64>,
        srp: Option<SrpConfig>,
        drag: Option<DragConfig>,
        thruster: Option<Thruster>,
        mode: Option<GuidanceMode>,
    ) -> Result<Self, NyxError> {
        if orbit.is_none() && dry_mass_kg.is_none() && fuel_mass_kg.is_none() {
            Ok(Self::default())
        } else if orbit.is_none() || dry_mass_kg.is_none() || fuel_mass_kg.is_none() {
            Err(NyxError::CustomError(format!(
                "orbit and dry_mass_kg must be specified"
            )))
        } else {
            Ok(Self {
                orbit: orbit.unwrap(),
                dry_mass_kg: dry_mass_kg.unwrap(),
                fuel_mass_kg: fuel_mass_kg.unwrap_or(0.0),
                thruster,
                mode: mode.unwrap_or(GuidanceMode::Coast),
                stm: None,
                srp: srp.unwrap_or_else(|| SrpConfig::default()),
                drag: drag.unwrap_or_else(|| DragConfig::default()),
            })
        }
    }

    #[classmethod]
    fn load(_cls: &PyType, path: &str) -> Result<Self, ConfigError> {
        <Self as ConfigRepr>::load(path)
    }

    #[classmethod]
    fn load_many(_cls: &PyType, path: &str) -> Result<Vec<Self>, ConfigError> {
        <Self as ConfigRepr>::load_many(path)
    }

    #[classmethod]
    fn load_named(_cls: &PyType, path: &str) -> Result<BTreeMap<String, Self>, ConfigError> {
        <Self as ConfigRepr>::load_named(path)
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }

    fn __str__(&self) -> String {
        format!("{self}\n{self:x}")
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp) -> Result<bool, NyxError> {
        match op {
            CompareOp::Eq => Ok(self == other),
            CompareOp::Ne => Ok(self != other),
            _ => Err(NyxError::CustomError(format!("{op:?} not available"))),
        }
    }

    fn dumps(&self, py: Python) -> Result<PyObject, NyxError> {
        pythonize(py, &self).map_err(|e| NyxError::CustomError(e.to_string()))
    }

    fn __getstate__(&self, py: Python) -> Result<PyObject, NyxError> {
        self.dumps(py)
    }

    #[classmethod]
    /// Loads the SpacecraftDynamics from its YAML representation
    fn loads(_cls: &PyType, state: &PyAny) -> Result<Self, ConfigError> {
        depythonize(state).map_err(|e| ConfigError::InvalidConfig(e.to_string()))
    }

    fn __setstate__(&mut self, state: &PyAny) -> Result<(), ConfigError> {
        *self = depythonize(state).map_err(|e| ConfigError::InvalidConfig(e.to_string()))?;
        Ok(())
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
    fn value_of(&self, param: StateParameter) -> Result<f64, NyxError> {
        self.value(param)
    }

    fn srp(&self) -> SrpConfig {
        self.srp
    }

    fn drag(&self) -> DragConfig {
        self.drag
    }
}

#[pymethods]
impl SrpConfig {
    #[new]
    #[pyo3(text_signature = "(area_m2, cr=1.8)")]
    pub fn py_new(area_m2: Option<f64>, cr: Option<f64>) -> Self {
        Self {
            area_m2: area_m2.unwrap_or(0.0),
            cr: cr.unwrap_or(1.8),
        }
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp) -> Result<bool, NyxError> {
        match op {
            CompareOp::Eq => Ok(self == other),
            CompareOp::Ne => Ok(self != other),
            _ => Err(NyxError::CustomError(format!("{op:?} not available"))),
        }
    }

    fn dumps(&self, py: Python) -> Result<PyObject, NyxError> {
        pythonize(py, &self).map_err(|e| NyxError::CustomError(e.to_string()))
    }

    fn __getstate__(&self, py: Python) -> Result<PyObject, NyxError> {
        self.dumps(py)
    }

    fn __setstate__(&mut self, state: &PyAny) -> Result<(), ConfigError> {
        *self = depythonize(state).map_err(|e| ConfigError::InvalidConfig(e.to_string()))?;
        Ok(())
    }

    #[getter]
    fn get_area_m2(&self) -> f64 {
        self.area_m2
    }

    #[getter]
    fn get_cr(&self) -> f64 {
        self.cr
    }

    #[setter]
    fn set_area_m2(&mut self, area_m2: f64) -> PyResult<()> {
        self.area_m2 = area_m2;
        Ok(())
    }

    #[setter]
    fn set_cr(&mut self, cr: f64) -> PyResult<()> {
        self.cr = cr;
        Ok(())
    }
}

#[pymethods]
impl DragConfig {
    #[new]
    #[pyo3(text_signature = "(area_m2, cd=2.2)")]
    pub fn py_new(area_m2: Option<f64>, cd: Option<f64>) -> Self {
        Self {
            area_m2: area_m2.unwrap_or(0.0),
            cd: cd.unwrap_or(1.8),
        }
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp) -> Result<bool, NyxError> {
        match op {
            CompareOp::Eq => Ok(self == other),
            CompareOp::Ne => Ok(self != other),
            _ => Err(NyxError::CustomError(format!("{op:?} not available"))),
        }
    }

    fn dumps(&self, py: Python) -> Result<PyObject, NyxError> {
        pythonize(py, &self).map_err(|e| NyxError::CustomError(e.to_string()))
    }

    fn __getstate__(&self, py: Python) -> Result<PyObject, NyxError> {
        self.dumps(py)
    }

    fn __setstate__(&mut self, state: &PyAny) -> Result<(), ConfigError> {
        *self = depythonize(state).map_err(|e| ConfigError::InvalidConfig(e.to_string()))?;
        Ok(())
    }

    #[getter]
    fn get_area_m2(&self) -> f64 {
        self.area_m2
    }

    #[getter]
    fn get_cd(&self) -> f64 {
        self.cd
    }

    #[setter]
    fn set_area_m2(&mut self, area_m2: f64) -> PyResult<()> {
        self.area_m2 = area_m2;
        Ok(())
    }

    #[setter]
    fn set_cr(&mut self, cd: f64) -> PyResult<()> {
        self.cd = cd;
        Ok(())
    }
}
