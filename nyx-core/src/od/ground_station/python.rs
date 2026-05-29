use super::super::msr::MeasurementType;
use super::super::noise::StochasticNoise;
use super::GroundStation;
use anise::astro::Location;
use hifitime::Duration;
use indexmap::IndexSet;
use pyo3::prelude::*;
use std::collections::HashMap;

#[cfg(feature = "python")]
#[pymethods]
impl GroundStation {
    #[new]
    #[pyo3(signature = (name, location, measurement_types, integration_time=None, light_time_correction=false, timestamp_noise_s=None, stochastic_noises=None))]
    fn py_new(
        name: String,
        location: Location,
        measurement_types: Vec<MeasurementType>,
        integration_time: Option<Duration>,
        light_time_correction: Option<bool>,
        timestamp_noise_s: Option<StochasticNoise>,
        stochastic_noises: Option<HashMap<MeasurementType, StochasticNoise>>,
    ) -> Self {
        Self {
            name,
            location,
            measurement_types: IndexSet::from_iter(measurement_types),
            integration_time,
            light_time_correction: light_time_correction.unwrap_or(false),
            timestamp_noise_s,
            stochastic_noises: stochastic_noises.map(|sn| sn.into_iter().collect()),
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
