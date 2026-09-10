use super::super::msr::MeasurementType;
use super::super::noise::StochasticNoise;
use super::{DopplerConfig, GroundStation};
use anise::astro::Location;
use anise::frames::FrameUid;
use indexmap::{IndexMap, IndexSet};
use pyo3::prelude::*;
use std::collections::HashMap;

#[cfg(feature = "python")]
#[pymethods]
impl GroundStation {
    /// Create a new Ground Station.
    ///
    /// :type name: str
    /// :type location: Location
    /// :type stochastic_noises: dict[MeasurementType, StochasticNoise]
    /// :type doppler_config: DopplerConfig | None
    /// :type light_time_correction: bool | None
    /// :type timestamp_noise_s: StochasticNoise | None
    /// :type obstruction_body: FrameUid | None
    #[new]
    #[pyo3(signature = (name, location, stochastic_noises, doppler_config=None, light_time_correction=false, timestamp_noise_s=None, obstruction_body=None))]
    fn py_new(
        name: String,
        location: Location,
        stochastic_noises: HashMap<MeasurementType, StochasticNoise>,
        doppler_config: Option<DopplerConfig>,
        light_time_correction: Option<bool>,
        timestamp_noise_s: Option<StochasticNoise>,
        obstruction_body: Option<FrameUid>,
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
            doppler_config: Some(doppler_config.unwrap_or_default()),
            light_time_correction: light_time_correction.unwrap_or(false),
            timestamp_noise_s,
            stochastic_noises: Some(stochastic_noises.into_iter().collect()),
            obstruction_body,
        }
    }

    /// Load GroundStation from a YAML string.
    ///
    /// :type yaml_str: str
    /// :rtype: GroundStation
    #[classmethod]
    #[pyo3(name = "from_yaml")]
    fn py_from_yaml(_cls: &Bound<'_, pyo3::types::PyType>, yaml_str: &str) -> PyResult<Self> {
        serde_yml::from_str(yaml_str)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// :rtype: str
    #[pyo3(name = "to_yaml")]
    fn py_to_yaml(&self) -> PyResult<String> {
        serde_yml::to_string(self)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Load multiple GroundStations from a YAML file.
    ///
    /// :type path: str
    /// :rtype: list[GroundStation]
    #[classmethod]
    #[pyo3(name = "load_many_yaml")]
    fn py_load_many_yaml(_cls: &Bound<'_, pyo3::types::PyType>, path: &str) -> PyResult<Vec<Self>> {
        use crate::io::ConfigRepr;
        Self::load_many(path).map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    /// Load multiple GroundStations from a YAML string.
    ///
    /// :type yaml_str: str
    /// :rtype: list[GroundStation]
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

    /// Dump multiple GroundStations to a YAML file.
    ///
    /// :type stations: list[GroundStation]
    /// :type path: str
    /// :rtype: None
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

    /// Dump multiple GroundStations to a YAML string.
    ///
    /// :type stations: list[GroundStation]
    /// :rtype: str
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
    pub fn get_doppler_config(&self) -> Option<DopplerConfig> {
        self.doppler_config
    }

    #[setter]
    pub fn set_doppler_config(&mut self, doppler_config: Option<DopplerConfig>) {
        self.doppler_config = doppler_config;
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
    pub fn get_obstruction_body(&self) -> Option<FrameUid> {
        self.obstruction_body
    }

    #[setter]
    pub fn set_obstruction_body(&mut self, obstruction_body: Option<FrameUid>) {
        self.obstruction_body = obstruction_body;
    }

    #[getter]
    pub fn get_measurement_types(&self) -> Vec<MeasurementType> {
        self.measurement_types.iter().cloned().collect()
    }

    #[setter]
    pub fn set_measurement_types(&mut self, types: Vec<MeasurementType>) {
        self.measurement_types = types.into_iter().collect();
    }

    /// Add a measurement type with stochastic noise.
    ///
    /// :type msr_type: MeasurementType
    /// :type noise: StochasticNoise
    /// :rtype: None
    pub fn add_measurement_type(&mut self, msr_type: MeasurementType, noise: StochasticNoise) {
        self.measurement_types.insert(msr_type);
        self.stochastic_noises
            .get_or_insert_with(IndexMap::new)
            .insert(msr_type, noise);
    }

    /// Remove a measurement type.
    ///
    /// :type msr_type: MeasurementType
    /// :rtype: bool
    pub fn remove_measurement_type(&mut self, msr_type: &MeasurementType) -> bool {
        // (Note: Requires IndexSet to be used with the `shift_remove` method to maintain order,
        // fallback to `.remove()` if order preservation upon deletion is not strictly required)
        self.measurement_types.shift_remove(msr_type)
    }

    /// Clear all measurement types
    ///
    /// :rtype: None
    pub fn clear_measurement_types(&mut self) {
        self.measurement_types.clear();
    }

    #[getter]
    pub fn get_stochastic_noises(&self) -> Option<Vec<(MeasurementType, StochasticNoise)>> {
        self.stochastic_noises
            .as_ref()
            .map(|map| map.iter().map(|(k, v)| (*k, *v)).collect())
    }

    /// Get stochastic noise for a measurement type.
    ///
    /// :type m_type: MeasurementType
    /// :rtype: StochasticNoise | None
    pub fn get_stochastic_noise(&self, m_type: &MeasurementType) -> Option<StochasticNoise> {
        self.stochastic_noises
            .as_ref()
            .and_then(|map| map.get(m_type).cloned())
    }

    /// Set stochastic noise for a measurement type.
    ///
    /// :type m_type: MeasurementType
    /// :type noise: StochasticNoise
    /// :rtype: None
    pub fn set_stochastic_noise(&mut self, m_type: MeasurementType, noise: StochasticNoise) {
        self.stochastic_noises
            .get_or_insert_with(indexmap::IndexMap::new)
            .insert(m_type, noise);
    }

    /// Remove stochastic noise for a measurement type.
    ///
    /// :type m_type: MeasurementType
    /// :rtype: StochasticNoise | None
    pub fn remove_stochastic_noise(&mut self, m_type: &MeasurementType) -> Option<StochasticNoise> {
        self.stochastic_noises
            .as_mut()
            .and_then(|map| map.shift_remove(m_type))
    }

    /// Clear stochastic noises
    ///
    /// :rtype: None
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
