/// nyx_space.time
use pyo3::prelude::*;
use pyo3::types::PyType;

use crate::time::Epoch as EpochRs;

#[pyclass]
pub struct Epoch {
    pub inner: EpochRs
}

#[pymethods]
impl Epoch {
    #[classmethod]
    pub fn from_tai_seconds(_cls: &PyType, seconds: f64) -> PyResult<Epoch> {
        Ok(Epoch { inner: EpochRs::from_tai_seconds(seconds) })
    }

    pub fn as_tai_seconds(&self) -> f64
    {
        self.inner.as_tai_seconds()
    }

    #[classmethod]
    pub fn from_gregorian_utc(_cls: &PyType, year: i32,
        month: u8,
        day: u8,
        hour: u8,
        minute: u8,
        second: u8,
        nanos: u32
    ) -> PyResult<Epoch>
    {
        Ok(Epoch { inner: EpochRs::from_gregorian_utc(year, month, day, hour, minute, second, nanos) })
    }
}

#[pymodule]
pub fn time(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<Epoch>()?;
    Ok(())
}
