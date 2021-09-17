use pyo3::{prelude::*, types::PyType, wrap_pymodule};

use crate::io::gravity::HarmonicsMem as HarmonicsMemRs;
use crate::io::scenario::ScenarioSerde as ScenarioSerdeRs;

/// nyx_space.io
#[pymodule]
fn io(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(gravity))?;
    m.add_class::<ScenarioSerde>()?;
    Ok(())
}

/// nyx_space.io.gravity
#[pymodule]
fn gravity(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<HarmonicsMem>()?;
    Ok(())
}

// nyx_space.io.gravity.HarmonicsMem
#[pyclass]
pub struct HarmonicsMem {
    pub inner: HarmonicsMemRs,
}

#[pymethods]
impl HarmonicsMem {
    #[classmethod]
    pub fn from_cof(
        _cls: &PyType,
        filepath: &str,
        degree: usize,
        order: usize,
        gunzipped: bool,
    ) -> PyResult<Self> {
        Ok(Self {
            inner: HarmonicsMemRs::from_cof(filepath, degree, order, gunzipped)?,
        })
    }
}

#[pyclass]
pub struct ScenarioSerde {
    pub inner: ScenarioSerdeRs,
}

#[pymethods]
impl ScenarioSerde {
    #[classmethod]
    pub fn from_toml_str(_cls: &PyType, toml_str: &str) -> PyResult<Self> {
        Ok(Self {
            inner: toml::from_str(toml_str).unwrap(),
        })
    }
}
