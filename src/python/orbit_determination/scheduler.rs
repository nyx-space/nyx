/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
pub use crate::io::ConfigError;
pub use crate::od::simulator::{Cadence, Handoff, Scheduler, Strand};
use crate::python::PythonError;

use hifitime::Duration;
use pyo3::basic::CompareOp;
use pyo3::prelude::*;
use std::str::FromStr;

#[pymethods]
impl Scheduler {
    #[new]
    #[pyo3(text_signature = "(handoff, cadence_on=None, cadence_off=None)")]
    fn py_new(
        handoff: Handoff,
        cadence_on: Option<String>,
        cadence_off: Option<String>,
        min_samples: Option<u32>,
    ) -> Result<Self, ConfigError> {
        let mut me = Self::builder().handoff(handoff).build();

        if let Some(min_samples) = min_samples {
            me.min_samples = min_samples;
        }

        if cadence_on.is_some() || cadence_off.is_some() {
            me.cadence = Cadence::Intermittent {
                on: Duration::from_str(cadence_on.unwrap().as_str()).map_err(|e| {
                    ConfigError::InvalidConfig {
                        msg: format!(
                        "{e} invalid format for schedule on (must be specified if schedule off is)"
                    ),
                    }
                })?,
                off: Duration::from_str(cadence_off.unwrap().as_str()).map_err(|e| {
                    ConfigError::InvalidConfig {
                        msg: format!(
                        "{e} invalid format for schedule off (must be specified if schedule on is)"
                    ),
                    }
                })?,
            };
        }

        Ok(me)
    }

    fn __repr__(&self) -> String {
        format!("{self:?}")
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __richcmp__(&self, other: &Self, op: CompareOp) -> Result<bool, PythonError> {
        match op {
            CompareOp::Eq => Ok(self == other),
            CompareOp::Ne => Ok(self != other),
            _ => Err(PythonError::OperationError { op }),
        }
    }
}
