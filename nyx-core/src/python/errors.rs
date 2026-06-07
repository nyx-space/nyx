use crate::{
    NyxError,
    io::{ConfigError, InputOutputError},
    md::trajectory::TrajError,
};
use pyo3::{exceptions::PyException, prelude::*};

impl From<InputOutputError> for PyErr {
    fn from(err: InputOutputError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}
impl From<ConfigError> for PyErr {
    fn from(err: ConfigError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}
impl From<TrajError> for PyErr {
    fn from(err: TrajError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}
impl From<NyxError> for PyErr {
    fn from(err: NyxError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}
