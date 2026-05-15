#[cfg(feature = "python")]
use {
    super::SpacecraftSequence,
    crate::Spacecraft,
    anise::almanac::Almanac,
    pyo3::prelude::*,
    std::sync::Arc,
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
        serde_dhall::from_str(dhall_str).parse().map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[classmethod]
    #[pyo3(name = "from_yaml")]
    fn py_from_yaml(_cls: &Bound<'_, pyo3::types::PyType>, yaml_str: &str) -> PyResult<Self> {
        serde_yml::from_str(yaml_str).map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[pyo3(name = "setup")]
    fn py_setup(&mut self, py: Python<'_>, almanac: Py<Almanac>) -> PyResult<()> {
        let almanac_ref = almanac.borrow(py);
        self.setup(Arc::new(almanac_ref.clone())).map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
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
        let trajs = self.propagate(state, until_phase, Arc::new(almanac_ref.clone()))
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        let result = trajs.into_iter().map(|traj| (traj.name, traj.states)).collect();
        Ok(result)
    }
}
