use super::{AccelModels, ForceModels, PropagatorConfig, SpacecraftSequence, Thruster};
use crate::md::prelude::FrameUid;
use crate::propagators::{IntegratorMethod, IntegratorOptions};
use crate::{dynamics::PointMasses, io::gravity::GravityFieldConfig};
use pyo3::exceptions::PyException;
use {
    crate::Spacecraft,
    anise::almanac::Almanac,
    pyo3::prelude::*,
    std::{collections::HashMap, sync::Arc},
};

#[cfg(feature = "python")]
#[pymethods]
impl SpacecraftSequence {
    #[new]
    fn py_new() -> Self {
        SpacecraftSequence::default()
    }

    #[classmethod]
    #[pyo3(name = "from_dhall")]
    fn py_from_dhall(_cls: &Bound<'_, pyo3::types::PyType>, dhall_str: &str) -> PyResult<Self> {
        serde_dhall::from_str(dhall_str)
            .parse()
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[classmethod]
    #[pyo3(name = "from_yaml")]
    fn py_from_yaml(_cls: &Bound<'_, pyo3::types::PyType>, yaml_str: &str) -> PyResult<Self> {
        serde_yml::from_str(yaml_str)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[pyo3(name = "setup")]
    fn py_setup(&mut self, py: Python<'_>, almanac: Py<Almanac>) -> PyResult<()> {
        let almanac_ref = almanac.borrow(py);
        self.setup(Arc::new(almanac_ref.clone()))
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[pyo3(name = "propagate")]
    fn py_propagate(
        &self,
        py: Python<'_>,
        state: Spacecraft,
        until_phase: Option<String>,
        almanac: Py<Almanac>,
    ) -> PyResult<Vec<(Option<String>, Vec<Spacecraft>)>> {
        let almanac_ref = almanac.borrow(py);
        let trajs = self
            .propagate(state, until_phase, Arc::new(almanac_ref.clone()))
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        let result = trajs
            .into_iter()
            .map(|traj| (traj.name, traj.states))
            .collect();
        Ok(result)
    }

    #[getter]
    fn get_thruster_sets(&self) -> HashMap<String, Thruster> {
        self.thruster_sets.clone()
    }

    fn thruster_set_insert(&mut self, name: String, thruster: Thruster) {
        self.thruster_sets.insert(name, thruster);
    }
    fn thruster_set_remove(&mut self, name: String) -> PyResult<()> {
        if self.thruster_sets.remove(&name).is_none() {
            Err(PyException::new_err(format!("{name} not in thruster set")))
        } else {
            Ok(())
        }
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl AccelModels {
    #[new]
    fn py_new(
        point_masses: Option<PointMasses>,
        gravity_field: Option<(GravityFieldConfig, FrameUid)>,
    ) -> Self {
        Self {
            point_masses,
            gravity_field,
        }
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl ForceModels {
    #[new]
    fn py_new() -> Self {
        // TODO Support building ForceModels
        Self {
            solar_pressure: None,
            drag: None,
        }
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl PropagatorConfig {
    #[new]
    fn py_new(
        accel_models: AccelModels,
        force_models: ForceModels,
        method: IntegratorMethod,
        options: IntegratorOptions,
    ) -> Self {
        Self {
            accel_models,
            force_models,
            method,
            options,
        }
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}
