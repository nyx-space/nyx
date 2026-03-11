use crate::io::InputOutputError;

use pyo3::{exceptions::PyException, prelude::*};

impl From<InputOutputError> for PyErr {
    fn from(err: InputOutputError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}
