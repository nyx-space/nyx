use super::IntegratorOptions;
use hifitime::Duration;
use pyo3::prelude::*;

#[pymethods]
impl IntegratorOptions {
    #[pyo3(signature=(min_step=None, max_step=None, tolerance=None))]
    #[new]
    fn py_new(
        min_step: Option<Duration>,
        max_step: Option<Duration>,
        tolerance: Option<f64>,
    ) -> Self {
        let mut opts = Self::default();
        if let Some(min_step) = min_step {
            opts.min_step = min_step;
        }
        if let Some(max_step) = max_step {
            opts.max_step = max_step;
        }
        if let Some(tolerance) = tolerance {
            opts.tolerance = tolerance;
        }
        opts.fixed_step = opts.min_step == opts.max_step;
        opts
    }

    #[getter]
    fn get_min_step(&self) -> Duration {
        self.min_step
    }
    #[getter]
    fn get_max_step(&self) -> Duration {
        self.max_step
    }
    #[getter]
    fn get_tolerance(&self) -> f64 {
        self.tolerance
    }
    fn __str__(&self) -> String {
        format!("{self}")
    }
    fn __repr__(&self) -> String {
        format!("{self} @ {self:p}")
    }
}
