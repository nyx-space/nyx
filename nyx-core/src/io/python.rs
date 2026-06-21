use super::ExportCfg;
use pyo3::prelude::*;

#[pymethods]
impl ExportCfg {
    #[pyo3(signature=(timestamped = false))]
    #[new]
    fn py_new(timestamped: bool) -> Self {
        if timestamped {
            Self::timestamped()
        } else {
            Self::default()
        }
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}
