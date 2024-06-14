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

pub use crate::io::ConfigError;
pub use crate::od::simulator::{Scheduler, Strand, TrkConfig};
use crate::python::PythonError;
use crate::{io::ConfigRepr, NyxError};
use hifitime::Duration;
use pyo3::basic::CompareOp;
use pyo3::prelude::*;
use pyo3::types::PyType;
use pythonize::{depythonize, pythonize};
use std::collections::BTreeMap;
use std::str::FromStr;

#[pymethods]
impl TrkConfig {
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

    #[new]
    #[pyo3(text_signature = "(sampling=None, strands=None, scheduler=None)")]
    fn py_new(
        sampling: Option<String>,
        strands: Option<Vec<Strand>>,
        scheduler: Option<Scheduler>,
    ) -> Result<Self, ConfigError> {
        let mut me = Self::default();

        if let Some(sampling) = sampling {
            me.sampling =
                Duration::from_str(&sampling).map_err(|e| ConfigError::InvalidConfig {
                    msg: format!("{e} invalid format for sampling"),
                })?;
        }

        me.strands = strands;

        if scheduler.is_some() {
            me.scheduler = scheduler;
        }

        Ok(me)
    }

    fn __repr__(&self) -> String {
        serde_yaml::to_string(&self).unwrap()
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp) -> Result<bool, PythonError> {
        match op {
            CompareOp::Eq => Ok(self == other),
            CompareOp::Ne => Ok(self != other),
            _ => Err(PythonError::OperationError { op }),
        }
    }

    fn dumps(&self, py: Python) -> Result<PyObject, NyxError> {
        pythonize(py, &self).map_err(|e| NyxError::CustomError { msg: e.to_string() })
    }

    fn __getstate__(&self, py: Python) -> Result<PyObject, NyxError> {
        self.dumps(py)
    }

    fn __setstate__(&mut self, state: &PyAny) -> Result<(), ConfigError> {
        *self =
            depythonize(state).map_err(|e| ConfigError::InvalidConfig { msg: e.to_string() })?;
        Ok(())
    }

    #[classmethod]
    /// Loads the TrkConfig from its YAML representation
    fn loads(_cls: &PyType, state: &PyAny) -> Result<Self, ConfigError> {
        depythonize(state).map_err(|e| ConfigError::InvalidConfig { msg: e.to_string() })
    }

    #[getter]
    fn get_sampling(&self) -> PyResult<Duration> {
        Ok(self.sampling)
    }

    #[setter]
    fn set_sampling(&mut self, sampling: Duration) -> PyResult<()> {
        self.sampling = sampling;
        Ok(())
    }
}
