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

use pyo3::prelude::*;

use crate::md::{Event, StateParameter};
use hifitime::Unit;

#[pymethods]
impl Event {
    /// Initializes a new event. Arguments are "parameter: StateParameter" and "desired_value: float".
    #[new]
    fn py_new(
        parameter: StateParameter,
        desired_value: f64,
        epoch_precision: Option<Unit>,
        value_precision: Option<f64>,
    ) -> Self {
        if let Some(value_precision) = value_precision {
            if let Some(epoch_precision) = epoch_precision {
                Self::specific(parameter, desired_value, value_precision, epoch_precision)
            } else {
                Self::within_tolerance(parameter, desired_value, value_precision)
            }
        } else if let Some(epoch_precision) = epoch_precision {
            Self::specific(
                parameter,
                desired_value,
                parameter.default_event_precision(),
                epoch_precision,
            )
        } else {
            Self::new(parameter, desired_value)
        }
    }

    #[cfg(feature = "python")]
    fn __str__(&self) -> String {
        format!("{self}")
    }
}
