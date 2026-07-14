use super::super::msr::MeasurementType;
use super::super::noise::StochasticNoise;
use super::PositionDevice;
use indexmap::IndexSet;
use pyo3::prelude::*;
use std::collections::HashMap;

#[cfg(feature = "python")]
#[pymethods]
impl PositionDevice {
    #[new]
    #[pyo3(signature = (name, stochastic_noises))]
    fn py_new(
        name: String,
        stochastic_noises: HashMap<MeasurementType, StochasticNoise>,
    ) -> Self {
        Self {
            name,
            measurement_types: IndexSet::from_iter(
                stochastic_noises
                    .keys()
                    .copied()
                    .collect::<Vec<MeasurementType>>(),
            ),
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

    fn __str__(&self) -> String {
        format!("{self}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}
