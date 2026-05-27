use super::ExportCfg;
use pyo3::prelude::*;

#[pymethods]
impl ExportCfg {
    #[new]
    fn py_new(timestamped: Option<bool>) -> Self {
        if timestamped.unwrap_or_default() {
            Self::default()
        } else {
            Self::timestamped()
        }
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}
