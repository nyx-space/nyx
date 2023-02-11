use pyo3::{prelude::*, types::PyType, types::PyDict};
use std::time::Duration;
use crate::propagators::Propagator as PropagatorRs;
use crate::propagators::PropInstance as PropInstanceRs;

use crate::python::od::State;
use crate::python::dynamics::Dynamics;

#[pymodule]
fn propagators(py: Python, m: &PyModule) -> PyResult<()>
{
    m.add_class::<Propagator>();
    m.add_class::<PropInstance>();
    Ok(())   
}

#[pyclass]
pub struct Propagator {
    pub inner: PropagatorRs<'static>
}
#[pymethods]
impl Propagator {
    #[classmethod]
    pub fn default(_cls: &PyType, dynamics: PyRef<Dynamics>) -> Self {
        Self { inner: PropagatorRs::default(dynamics.inner.clone()) }
    }
    pub fn with(&self, state: State) -> PropInstance
    {
        PropInstance { inner: self.inner.with(state.inner) }
    }
}

#[pyclass]
pub struct PropInstance {
    pub inner: PropInstanceRs<'static>
}
#[pymethods]
impl PropInstance {
    pub fn for_duration_sec(&self, duration_sec: f64) -> PyResult<State>
    {
        let result = self.inner
                                      .for_duration(
                                          Duration::from_millis(duration_sec * 1000.0)
                                       )?;
        Ok(State { inner: result })
    }
}
