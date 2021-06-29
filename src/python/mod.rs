/// This package provides python bindings for the rust crate
/// [nyx-space](http://gitlab.com/nyx-space/nyx) by building
/// a native Python extension using [PyO3](https://github.com/pyO3/pyO3)

use pyo3::prelude::*;
use pyo3::{wrap_pyfunction, wrap_pymodule};
/// Formats the sum of two numbers as string.
#[pyfunction]
fn sum_as_string(a: usize, b: usize) -> PyResult<String> {
    Ok((a + b).to_string())
}

/// A Python module implemented in Rust.

pub mod cosmic;

/// nyx_space.dynamics
#[pymodule]
fn dynamics(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sum_as_string, m)?)?;
    Ok(())
}

#[pymodule]
fn nyx_space(_py: Python, m: &PyModule) -> PyResult<()>
{
    m.add_wrapped(wrap_pymodule!(dynamics))?;
    Ok(())
}
