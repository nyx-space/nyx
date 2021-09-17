use crate::cosmic::Spacecraft as SpacecraftRs;
use crate::errors::NyxError;
use pyo3::exceptions::PyException;
/// This package provides python bindings for the rust crate
/// [nyx-space](http://gitlab.com/nyx-space/nyx) by building
/// a native Python extension using [PyO3](https://github.com/pyO3/pyO3)
use pyo3::prelude::*;
use pyo3::wrap_pymodule;

/// A Python module implemented in Rust.
mod io;
use io::PyInit_io;
mod time;
use time::PyInit_time;
mod cosmic;
use cosmic::PyInit_cosmic;
mod md;
use md::PyInit_md;

// mod dynamics;
// use dynamics::PyInit_dynamics;
// mod propagators;
// use propagators::PyInit_propagators;
// mod od;
// use od::PyInit_od;

impl std::convert::From<NyxError> for PyErr {
    fn from(err: NyxError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}

#[pymodule]
fn nyx_space(_py: Python, m: &PyModule) -> PyResult<()> {
    // m.add_wrapped(wrap_pymodule!(propagators))?;
    // m.add_wrapped(wrap_pymodule!(dynamics))?;
    m.add_wrapped(wrap_pymodule!(cosmic))?;
    m.add_wrapped(wrap_pymodule!(time))?;
    m.add_wrapped(wrap_pymodule!(io))?;
    m.add_wrapped(wrap_pymodule!(md))?;
    m.add_class::<Spacecraft>()?;
    Ok(())
}

#[pyclass]
pub struct Spacecraft {
    pub inner: SpacecraftRs,
}

#[pymethods]
impl Spacecraft {
    #[new]
    pub fn new(
        orbit: PyRef<crate::python::cosmic::Orbit>,
        dry_mass_kg: f64,
        fuel_mass_kg: f64,
        srp_area_m2: f64,
        drag_area_m2: f64,
        cr: f64,
        cd: f64,
    ) -> Self {
        Self {
            inner: SpacecraftRs::new(
                orbit.inner,
                dry_mass_kg,
                fuel_mass_kg,
                srp_area_m2,
                drag_area_m2,
                cr,
                cd,
            ),
        }
    }
}
