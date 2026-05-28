#[cfg(feature = "python")]
use super::GroundStation;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use anise::astro::Location;
#[cfg(feature = "python")]
use indexmap::{IndexSet, IndexMap};
#[cfg(feature = "python")]
use super::super::msr::MeasurementType;
#[cfg(feature = "python")]
use super::super::noise::StochasticNoise;
#[cfg(feature = "python")]
use hifitime::Duration;

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
        stochastic_noises: Option<std::collections::HashMap<MeasurementType, StochasticNoise>>,
    ) -> Self {
        let mut stochastics = None;
        if let Some(sn) = stochastic_noises {
            let mut im = IndexMap::new();
            for (k, v) in sn {
                im.insert(k, v);
            }
            stochastics = Some(im);
        }

        Self {
            name,
            location,
            measurement_types: IndexSet::from_iter(measurement_types),
            integration_time,
            light_time_correction: light_time_correction.unwrap_or(false),
            timestamp_noise_s,
            stochastic_noises: stochastics,
        }
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

    #[pyo3(name = "to_dhall")]
    fn py_to_dhall(&self) -> PyResult<String> {
        serde_dhall::serialize(self)
            .to_string()
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
        Self::load_many(path)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[classmethod]
    #[pyo3(name = "loads_many_yaml")]
    fn py_loads_many_yaml(_cls: &Bound<'_, pyo3::types::PyType>, yaml_str: &str) -> PyResult<Vec<Self>> {
        use crate::io::ConfigRepr;
        Self::loads_many(yaml_str)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[classmethod]
    #[pyo3(name = "dump_many_yaml")]
    fn py_dump_many_yaml(_cls: &Bound<'_, pyo3::types::PyType>, stations: Vec<Self>, path: &str) -> PyResult<()> {
        let s = serde_yml::to_string(&stations)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;
        std::fs::write(path, s)
            .map_err(|e| pyo3::exceptions::PyIOError::new_err(e.to_string()))
    }

    #[classmethod]
    #[pyo3(name = "dumps_many_yaml")]
    fn py_dumps_many_yaml(_cls: &Bound<'_, pyo3::types::PyType>, stations: Vec<Self>) -> PyResult<String> {
        serde_yml::to_string(&stations)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}
