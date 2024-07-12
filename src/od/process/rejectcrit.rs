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

use crate::io::ConfigRepr;
#[cfg(feature = "python")]
use crate::python::PythonError;
#[cfg(feature = "python")]
use crate::NyxError;
#[cfg(feature = "python")]
use pyo3::class::basic::CompareOp;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pythonize::{depythonize, pythonize};
use serde_derive::{Deserialize, Serialize};

/// Reject measurements if the prefit is greater than the provided sigmas deviation from the measurement noise.
///
/// # Important
/// Some software, like ODTK, processes each measurement as a scalar. Nyx processes the measurements together.
/// As such, if the prefit on range is bad, then the Doppler measurement with the same time stamp will also be rejected.
/// This leads to better convergence of the filter, and more appropriate results.
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(module = "nyx_space.orbit_determination"))]
pub struct ResidRejectCrit {
    /// Number of sigmas for a measurement to be considered an outlier.
    pub num_sigmas: f64,
}

#[cfg(feature = "python")]
#[pymethods]
impl ResidRejectCrit {
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
        *self =
            depythonize(state).map_err(|e| ConfigError::InvalidConfig { msg: e.to_string() })?;
        Ok(())
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp) -> Result<bool, PythonError> {
        match op {
            CompareOp::Eq => Ok(self == other),
            CompareOp::Ne => Ok(self != other),
            _ => Err(PythonError::OperationError { op }),
        }
    }
}

impl Default for ResidRejectCrit {
    /// By default, a measurement is rejected if its prefit residual is greater the 4-sigma value of the measurement noise at that time step.
    /// This corresponds to [1 chance in in 15,787](https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule).
    fn default() -> Self {
        Self { num_sigmas: 4.0 }
    }
}

impl ConfigRepr for ResidRejectCrit {}
