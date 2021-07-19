use pyo3::{prelude::*, wrap_pymodule, types::PyType};

use crate::io::gravity::HarmonicsMem as HarmonicsMemRs;

/// nyx_space.io
#[pymodule]
fn io(_py: Python, m: &PyModule) -> PyResult<()>
{
    m.add_wrapped(wrap_pymodule!(gravity))?;
    Ok(())
}

/// nyx_space.io.gravity
#[pymodule]
fn gravity(_py: Python, m: &PyModule) -> PyResult<()>
{
    m.add_class::<HarmonicsMem>()?;
    Ok(())
}

// nyx_space.io.gravity.HarmonicsMem
#[pyclass]
pub struct HarmonicsMem
{
    pub inner: HarmonicsMemRs
}

#[pymethods]
impl HarmonicsMem {
    #[classmethod]
    pub fn from_cof(
        _cls: &PyType, 
        filepath: &str,
        degree: usize,
        order: usize,
        gunzipped: bool
    ) -> PyResult<Self>
    {
        Ok(Self { inner: HarmonicsMemRs::from_cof(filepath, degree, order, gunzipped)? })
    }
}
