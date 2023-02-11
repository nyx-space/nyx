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

use crate::od::ui::*;
use crate::NyxError;
use hifitime::leap_seconds::{LatestLeapSeconds, LeapSecondsFile};
use hifitime::prelude::*;
use hifitime::ut1::Ut1Provider;
use pyo3::{exceptions::PyException, prelude::*};

impl std::convert::From<NyxError> for PyErr {
    fn from(err: NyxError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}

#[pymodule]
fn nyx_space(py: Python, m: &PyModule) -> PyResult<()> {
    register_time_module(py, m)?;
    register_od(py, m)?;

    Ok(())
}

/// Reexport hifitime as nyx_space.time
fn register_time_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()> {
    let sm = PyModule::new(py, "time")?;

    sm.add_class::<Epoch>()?;
    sm.add_class::<TimeScale>()?;
    sm.add_class::<TimeSeries>()?;
    sm.add_class::<Duration>()?;
    sm.add_class::<Unit>()?;
    sm.add_class::<LatestLeapSeconds>()?;
    sm.add_class::<LeapSecondsFile>()?;
    sm.add_class::<Ut1Provider>()?;

    parent_module.add_submodule(sm)?;
    Ok(())
}

/// Orbit determination
fn register_od(py: Python<'_>, parent_module: &PyModule) -> PyResult<()> {
    let sm = PyModule::new(py, "od")?;

    sm.add_class::<GroundStation>()?;

    parent_module.add_submodule(sm)?;
    Ok(())
}
