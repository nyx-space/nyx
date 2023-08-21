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

use std::{collections::HashMap, sync::Arc};

use crate::{
    cosmic::Cosm,
    io::{estimate::OrbitEstimateSerde, ConfigRepr, Configurable},
    od::estimate::KfEstimate,
    NyxError, Orbit,
};
use nalgebra::Matrix6;
use numpy::PyReadonlyArrayDyn;
use pyo3::class::basic::CompareOp;
use pyo3::prelude::*;
use pyo3::types::PyType;
use pythonize::depythonize;

use super::ConfigError;

/// An estimate of an orbit with its covariance, the latter should be a numpy array of size 36.
#[derive(Debug, Clone, PartialEq)]
#[pyclass]
#[pyo3(text_signature = "(nominal_orbit, covariance)")]
pub(crate) struct OrbitEstimate(pub(crate) KfEstimate<Orbit>);

impl Configurable for OrbitEstimate {
    type IntermediateRepr = OrbitEstimateSerde;

    fn from_config(cfg: Self::IntermediateRepr, _cosm: Arc<Cosm>) -> Result<Self, ConfigError>
    where
        Self: Sized,
    {
        Ok(Self(KfEstimate::from_covar(
            Orbit::from(cfg.nominal),
            cfg.covar.to_matrix(),
        )))
    }

    fn to_config(&self) -> Result<Self::IntermediateRepr, ConfigError> {
        todo!()
    }
}

#[pymethods]
impl OrbitEstimate {
    #[new]
    fn new(nominal: Orbit, covar: PyReadonlyArrayDyn<f64>) -> Result<Self, NyxError> {
        // Check the shape of the input
        let mat6 = match covar.shape() {
            &[36] | &[36, 1] | &[6, 6] => {
                let data = covar
                    .as_slice()
                    .map_err(|e| NyxError::CustomError(format!("{e}")))?;
                let mut mat = Matrix6::zeros();
                for i in 0..6 {
                    for j in 0..6 {
                        mat[(i, j)] = data[6 * i + j];
                    }
                }
                mat
            }
            _ => {
                return Err(NyxError::CustomError(format!(
                    "covar must be 6x6 or 36x1 but is {:?}",
                    covar.shape()
                )))
            }
        };
        Ok(Self(KfEstimate::from_covar(nominal, mat6)))
    }

    #[classmethod]
    fn load(_cls: &PyType, path: &str) -> Result<Self, ConfigError> {
        let serde = OrbitEstimateSerde::load(path)?;

        let cosm = Cosm::de438();

        Self::from_config(serde, cosm)
    }

    #[classmethod]
    fn load_many(_cls: &PyType, path: &str) -> Result<Vec<Self>, ConfigError> {
        let stations = OrbitEstimateSerde::load_many(path)?;

        let cosm = Cosm::de438();

        let mut selves = Vec::with_capacity(stations.len());

        for serde in stations {
            selves.push(Self::from_config(serde, cosm.clone())?);
        }

        Ok(selves)
    }

    #[classmethod]
    fn load_named(_cls: &PyType, path: &str) -> Result<HashMap<String, Self>, ConfigError> {
        let orbits = OrbitEstimateSerde::load_named(path)?;

        let cosm = Cosm::de438();

        let mut selves = HashMap::with_capacity(orbits.len());

        for (k, v) in orbits {
            selves.insert(k, Self::from_config(v, cosm.clone())?);
        }

        Ok(selves)
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp) -> Result<bool, NyxError> {
        match op {
            CompareOp::Eq => Ok(self.__repr__() == other.__repr__()),
            CompareOp::Ne => Ok(self.__repr__() != other.__repr__()),
            _ => Err(NyxError::CustomError(format!("{op:?} not available"))),
        }
    }

    #[classmethod]
    /// Loads the SpacecraftDynamics from its YAML representation
    fn loads(_cls: &PyType, state: &PyAny) -> Result<Self, ConfigError> {
        <Self as Configurable>::from_config(
            depythonize(state).map_err(|e| ConfigError::InvalidConfig(e.to_string()))?,
            Cosm::de438(),
        )
    }

    // Manual getter/setters -- waiting on https://github.com/PyO3/pyo3/pull/2786

    #[getter]
    fn get_orbit(&self) -> PyResult<Orbit> {
        Ok(self.0.nominal_state)
    }

    #[getter]
    fn get_is_predicted(&self) -> PyResult<bool> {
        Ok(self.0.predicted)
    }

    fn __str__(&self) -> String {
        format!("{}", self.0)
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.0)
    }
}
