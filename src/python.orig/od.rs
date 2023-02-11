use pyo3::{prelude::*, types::PyDict};
use std::sync::Arc;
use crate::od::State as StateRs;

#[pymodule]
fn od(py: Python, m: &PyModule) -> PyResult<()>
{
    m.add_class::<State>();
    Ok(())   
}

#[pyclass]
pub struct State {
    pub inner: Arc<dyn StateRs>
}
#[pymethods]
impl State {}
