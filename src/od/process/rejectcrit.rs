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

use crate::cosmic::Cosm;
use crate::io::{ConfigError, ConfigRepr, Configurable};
#[cfg(feature = "python")]
use crate::NyxError;
#[cfg(feature = "python")]
use pyo3::class::basic::CompareOp;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pythonize::{depythonize, pythonize};
use serde_derive::{Deserialize, Serialize};
use std::sync::Arc;

/// Reject measurements with a residual ratio greater than the provided sigmas values. Will only be turned used if at least min_accepted measurements have been processed so far.
/// If unsure, use the default: `FltResid::default()` in Rust, and `FltResid()` in Python (i.e. construct without arguments).
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(module = "nyx_space.orbit_determination"))]
pub struct FltResid {
    /// Minimum number of accepted measurements before applying the rejection criteria.
    pub min_accepted: usize,
    /// Number of sigmas for a measurement to be considered an outlier.
    pub num_sigmas: f64,
}

#[cfg(feature = "python")]
#[pymethods]
impl FltResid {
    #[new]
    #[pyo3(text_signature = "(min_accepted=None, num_sigmas=None)")]
    fn py_new(min_accepted: Option<usize>, num_sigmas: Option<f64>) -> Self {
        let mut me = Self::default();
        if let Some(min_accepted) = min_accepted {
            me.min_accepted = min_accepted;
        }
        if let Some(num_sigmas) = num_sigmas {
            me.num_sigmas = num_sigmas;
        }
        me
    }

    #[getter]
    fn get_min_accepted(&self) -> usize {
        self.min_accepted
    }

    #[setter(orbit)]
    fn py_set_min_accepted(&mut self, min_accepted: usize) -> PyResult<()> {
        self.min_accepted = min_accepted;
        Ok(())
    }

    #[getter]
    fn get_num_sigmas(&self) -> f64 {
        self.num_sigmas
    }

    #[setter(orbit)]
    fn py_set_num_sigmas(&mut self, num_sigmas: f64) -> PyResult<()> {
        self.num_sigmas = num_sigmas;
        Ok(())
    }

    fn __repr__(&self) -> String {
        format!("{self:?}")
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn dumps(&self, py: Python) -> Result<PyObject, NyxError> {
        pythonize(py, &self).map_err(|e| NyxError::CustomError { msg: e.to_string() })
    }

    fn __getstate__(&self, py: Python) -> Result<PyObject, NyxError> {
        self.dumps(py)
    }

    fn __setstate__(&mut self, state: &PyAny) -> Result<(), ConfigError> {
        *self = depythonize(state).map_err(|e| ConfigError::InvalidConfig(e.to_string()))?;
        Ok(())
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp) -> Result<bool, NyxError> {
        match op {
            CompareOp::Eq => Ok(self == other),
            CompareOp::Ne => Ok(self != other),
            _ => Err(NyxError::CustomError(format!("{op:?} not available"))),
        }
    }
}

impl Default for FltResid {
    fn default() -> Self {
        Self {
            min_accepted: 10,
            num_sigmas: 3.0,
        }
    }
}

impl ConfigRepr for FltResid {}

impl Configurable for FltResid {
    type IntermediateRepr = Self;

    fn from_config(cfg: Self, _cosm: Arc<Cosm>) -> Result<Self, ConfigError>
    where
        Self: Sized,
    {
        Ok(cfg)
    }

    fn to_config(&self) -> Result<Self::IntermediateRepr, ConfigError> {
        Ok(*self)
    }
}
