use super::IntegratorOptions;
use hifitime::{Duration, Unit};
use pyo3::prelude::*;

#[pymethods]
impl IntegratorOptions {
    #[pyo3(signature=(min_step_s=None, max_step_s=None,tolerance=None))]
    #[new]
    fn py_new(min_step_s: Option<f64>, max_step_s: Option<f64>, tolerance: Option<f64>) -> Self {
        let mut opts = Self::default();
        if let Some(min_step_s) = min_step_s {
            opts.min_step = Unit::Second * min_step_s;
        }
        if let Some(max_step_s) = max_step_s {
            opts.max_step = Unit::Second * max_step_s;
        }
        if let Some(tolerance) = tolerance {
            opts.tolerance = tolerance;
        }
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
}
