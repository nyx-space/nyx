use crate::{
    NyxError,
    dynamics::DynamicsError,
    io::{ConfigError, InputOutputError},
    md::trajectory::TrajError,
    od::ODError,
    propagators::PropagationError,
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
impl From<PropagationError> for PyErr {
    fn from(err: PropagationError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}
impl From<DynamicsError> for PyErr {
    fn from(err: DynamicsError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}
impl From<ODError> for PyErr {
    fn from(err: ODError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}
