use crate::od::msr::MeasurementType;
use crate::od::noise::StochasticNoise;
use anise::prelude::Frame;
use indexmap::IndexMap;
use pyo3::prelude::*;

use super::PositionDevice;

#[cfg(feature = "python")]
#[pymethods]
impl PositionDevice {
    /// Create a new Position Tracking Device.
    ///
    /// :type name: str
    /// :type frame: Frame | None
    #[new]
    #[pyo3(signature = (name, frame=None))]
    fn py_new(name: String, frame: Option<Frame>) -> Self {
        Self::new(name, frame)
    }

    /// Add a measurement type with stochastic noise.
    ///
    /// :type msr_type: MeasurementType
    /// :type noise: StochasticNoise
    /// :rtype: PositionDevice
    #[pyo3(name = "with_noise")]
    pub fn py_with_noise(
        mut slf: PyRefMut<'_, Self>,
        msr_type: MeasurementType,
        noise: StochasticNoise,
    ) -> PyResult<PyRefMut<'_, Self>> {
        if slf.stochastic_noises.is_none() {
            slf.stochastic_noises = Some(IndexMap::new());
        }
        slf.stochastic_noises
            .as_mut()
            .unwrap()
            .insert(msr_type, noise);
        slf.measurement_types.insert(msr_type);
        Ok(slf)
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
    pub fn get_frame(&self) -> Option<Frame> {
        self.frame
    }

    #[setter]
    pub fn set_frame(&mut self, frame: Option<Frame>) {
        self.frame = frame;
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
