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
pub use crate::io::ConfigError;
pub use crate::od::simulator::{Schedule, TrkConfig};
use crate::{io::ConfigRepr, od::simulator::Availability, NyxError};
use hifitime::{Duration, Epoch};
use pyo3::basic::CompareOp;
use pyo3::prelude::*;
use pyo3::types::PyType;
use pythonize::{depythonize, pythonize};
use std::{collections::HashMap, str::FromStr};

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
    fn load_named(_cls: &PyType, path: &str) -> Result<HashMap<String, Self>, ConfigError> {
        <Self as ConfigRepr>::load_named(path)
    }

    #[new]
    fn py_new(
        start: Option<String>,
        end: Option<String>,
        schedule_on: Option<String>,
        schedule_off: Option<String>,
        sampling: Option<String>,
    ) -> Result<Self, ConfigError> {
        let mut me = Self::default();

        if schedule_on.is_some() || schedule_off.is_some() {
            me.schedule = Schedule::Intermittent {
                on: Duration::from_str(schedule_on.unwrap().as_str()).map_err(|e| {
                    ConfigError::InvalidConfig(format!(
                        "{e} invalid format for schedule on (must be specified if schedule off is)"
                    ))
                })?,
                off: Duration::from_str(schedule_off.unwrap().as_str()).map_err(|e| {
                    ConfigError::InvalidConfig(format!(
                        "{e} invalid format for schedule off (must be specified if schedule on is)"
                    ))
                })?,
            };
        }

        if let Some(start) = start {
            if start.to_ascii_lowercase() == "visible" {
                me.start = Availability::Visible
            } else {
                me.start = Availability::Epoch(Epoch::from_str(&start).map_err(|e| {
                    ConfigError::InvalidConfig(format!("{e} invalid format for start availability"))
                })?)
            }
        }

        if let Some(end) = end {
            if end.to_ascii_lowercase() == "visible" {
                me.end = Availability::Visible
            } else {
                me.end = Availability::Epoch(Epoch::from_str(&end).map_err(|e| {
                    ConfigError::InvalidConfig(format!("{e} invalid format for end availability"))
                })?)
            }
        }

        if let Some(sampling) = sampling {
            me.sampling = Duration::from_str(&sampling).map_err(|e| {
                ConfigError::InvalidConfig(format!("{e} invalid format for sampling"))
            })?;
        }

        Ok(me)
    }

    fn __repr__(&self) -> String {
        serde_yaml::to_string(&self).unwrap()
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
    fn get_sampling(&self) -> PyResult<Duration> {
        Ok(self.sampling)
    }

    #[setter]
    fn set_sampling(&mut self, sampling: Duration) -> PyResult<()> {
        self.sampling = sampling;
        Ok(())
    }

    /// Allows setting the start and end availabilities and the sampling.
    /// Availabilities must be either `Visible` or an Epoch as a string.
    /// The sampling must be a Duration object.
    /// Example usage: `cfg.set(start='Visible', end='2020-01-01 15.26.30 UTC')`
    fn set(
        &mut self,
        start: Option<String>,
        end: Option<String>,
        sampling: Option<Duration>,
    ) -> Result<(), NyxError> {
        if let Some(start) = start {
            if start.to_ascii_lowercase() == "visible" {
                self.start = Availability::Visible
            } else {
                self.start = Availability::Epoch(Epoch::from_str(&start).map_err(|e| {
                    NyxError::CustomError(format!("{e} invalid format for start availability"))
                })?)
            }
        }
        if let Some(end) = end {
            if end.to_ascii_lowercase() == "visible" {
                self.end = Availability::Visible
            } else {
                self.end = Availability::Epoch(Epoch::from_str(&end).map_err(|e| {
                    NyxError::CustomError(format!("{e} invalid format for end availability"))
                })?)
            }
        }
        if let Some(sampling) = sampling {
            self.sampling = sampling;
        }
        Ok(())
    }
}
