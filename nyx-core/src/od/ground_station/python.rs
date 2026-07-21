use super::super::msr::MeasurementType;
use super::super::noise::StochasticNoise;
use super::GroundStation;
use anise::astro::Location;
use hifitime::Duration;
use indexmap::{IndexMap, IndexSet};
use pyo3::prelude::*;
use std::collections::HashMap;

#[cfg(feature = "python")]
#[pymethods]
impl GroundStation {
    #[new]
    #[pyo3(signature = (name, location, stochastic_noises, integration_time=None, light_time_correction=false, timestamp_noise_s=None))]
    fn py_new(
        name: String,
        location: Location,
        stochastic_noises: HashMap<MeasurementType, StochasticNoise>,
        integration_time: Option<Duration>,
        light_time_correction: Option<bool>,
        timestamp_noise_s: Option<StochasticNoise>,
    ) -> Self {
        Self {
            name,
            location,
            measurement_types: IndexSet::from_iter(
                stochastic_noises
                    .keys()
                    .copied()
                    .collect::<Vec<MeasurementType>>(),
            ),
            integration_time,
            light_time_correction: light_time_correction.unwrap_or(false),
            timestamp_noise_s,
            stochastic_noises: Some(stochastic_noises.into_iter().collect()),
        }
    }

    #[classmethod]
    #[pyo3(name = "from_yaml")]
    fn py_from_yaml(_cls: &Bound<'_, pyo3::types::PyType>, yaml_str: &str) -> PyResult<Self> {
        serde_yml::from_str(yaml_str)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[pyo3(name = "to_yaml")]
    fn py_to_yaml(&self) -> PyResult<String> {
        serde_yml::to_string(self)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[classmethod]
    #[pyo3(name = "load_many_yaml")]
    fn py_load_many_yaml(_cls: &Bound<'_, pyo3::types::PyType>, path: &str) -> PyResult<Vec<Self>> {
        use crate::io::ConfigRepr;
        Self::load_many(path).map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[classmethod]
    #[pyo3(name = "loads_many_yaml")]
    fn py_loads_many_yaml(
        _cls: &Bound<'_, pyo3::types::PyType>,
        yaml_str: &str,
    ) -> PyResult<Vec<Self>> {
        use crate::io::ConfigRepr;
        Self::loads_many(yaml_str)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[classmethod]
    #[pyo3(name = "dump_many_yaml")]
    fn py_dump_many_yaml(
        _cls: &Bound<'_, pyo3::types::PyType>,
        stations: Vec<Self>,
        path: &str,
    ) -> PyResult<()> {
        let s = serde_yml::to_string(&stations)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        std::fs::write(path, s).map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    #[classmethod]
    #[pyo3(name = "dumps_many_yaml")]
    fn py_dumps_many_yaml(
        _cls: &Bound<'_, pyo3::types::PyType>,
        stations: Vec<Self>,
    ) -> PyResult<String> {
        serde_yml::to_string(&stations)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[getter]
    pub fn get_name(&self) -> String {
        self.name.clone()
    }

    #[setter]
    pub fn set_name(&mut self, name: String) {
        self.name = name;
    }

    #[getter]
    pub fn get_location(&self) -> Location {
        self.location.clone()
    }

    #[setter]
    pub fn set_location(&mut self, location: Location) {
        self.location = location;
    }

    #[getter]
    pub fn get_integration_time(&self) -> Option<Duration> {
        self.integration_time
    }

    #[setter]
    pub fn set_integration_time(&mut self, integration_time: Option<Duration>) {
        self.integration_time = integration_time;
    }

    #[getter]
    pub fn get_light_time_correction(&self) -> bool {
        self.light_time_correction
    }

    #[setter]
    pub fn set_light_time_correction(&mut self, light_time_correction: bool) {
        self.light_time_correction = light_time_correction;
    }

    #[getter]
    pub fn get_timestamp_noise_s(&self) -> Option<StochasticNoise> {
        self.timestamp_noise_s
    }

    #[setter]
    pub fn set_timestamp_noise_s(&mut self, noise: Option<StochasticNoise>) {
        self.timestamp_noise_s = noise;
    }

    #[getter]
    pub fn get_measurement_types(&self) -> Vec<MeasurementType> {
        self.measurement_types.iter().cloned().collect()
    }

    #[setter]
    pub fn set_measurement_types(&mut self, types: Vec<MeasurementType>) {
        self.measurement_types = types.into_iter().collect();
    }

    pub fn add_measurement_type(&mut self, msr_type: MeasurementType, noise: StochasticNoise) {
        self.measurement_types.insert(msr_type);
        self.stochastic_noises
            .get_or_insert_with(IndexMap::new)
            .insert(msr_type, noise);
    }

    pub fn remove_measurement_type(&mut self, msr_type: &MeasurementType) -> bool {
        // (Note: Requires IndexSet to be used with the `shift_remove` method to maintain order,
        // fallback to `.remove()` if order preservation upon deletion is not strictly required)
        self.measurement_types.shift_remove(msr_type)
    }

    pub fn clear_measurement_types(&mut self) {
        self.measurement_types.clear();
    }

    #[getter]
    pub fn get_stochastic_noises(&self) -> Option<Vec<(MeasurementType, StochasticNoise)>> {
        self.stochastic_noises
            .as_ref()
            .map(|map| map.iter().map(|(k, v)| (*k, *v)).collect())
    }

    pub fn get_stochastic_noise(&self, m_type: &MeasurementType) -> Option<StochasticNoise> {
        self.stochastic_noises
            .as_ref()
            .and_then(|map| map.get(m_type).cloned())
    }

    pub fn set_stochastic_noise(&mut self, m_type: MeasurementType, noise: StochasticNoise) {
        self.stochastic_noises
            .get_or_insert_with(indexmap::IndexMap::new)
            .insert(m_type, noise);
    }

    pub fn remove_stochastic_noise(&mut self, m_type: &MeasurementType) -> Option<StochasticNoise> {
        self.stochastic_noises
            .as_mut()
            .and_then(|map| map.shift_remove(m_type))
    }

    pub fn clear_stochastic_noises(&mut self) {
        self.stochastic_noises = None;
    }

    fn __str__(&self) -> String {
        format!("{self}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }

    fn __eq__(&self, other: &Self) -> bool {
        self == other
    }
}
