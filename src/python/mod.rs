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

use crate::io::ConfigError;
use crate::NyxError;
use hifitime::leap_seconds::{LatestLeapSeconds, LeapSecondsFile};
use hifitime::prelude::*;
use hifitime::ut1::Ut1Provider;
use pyo3::py_run;
use pyo3::{exceptions::PyException, prelude::*};

pub(crate) mod cosmic;
pub(crate) mod mission_design;
mod monte_carlo;
mod orbit_determination;

impl From<NyxError> for PyErr {
    fn from(err: NyxError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}

impl From<ConfigError> for PyErr {
    fn from(err: ConfigError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}

#[pymodule]
fn _nyx_space(py: Python, m: &PyModule) -> PyResult<()> {
    pyo3_log::init();

    register_time_module(py, m)?;
    orbit_determination::register_od(py, m)?;
    mission_design::register_md(py, m)?;
    cosmic::register_cosmic(py, m)?;
    monte_carlo::register_mc(py, m)?;

    Ok(())
}

/// Reexport hifitime as nyx_space.time
fn register_time_module(py: Python<'_>, parent_module: &PyModule) -> PyResult<()> {
    let sm = PyModule::new(py, "_nyx_space.time")?;

    sm.add_class::<Epoch>()?;
    sm.add_class::<TimeScale>()?;
    sm.add_class::<TimeSeries>()?;
    sm.add_class::<Duration>()?;
    sm.add_class::<Unit>()?;
    sm.add_class::<LatestLeapSeconds>()?;
    sm.add_class::<LeapSecondsFile>()?;
    sm.add_class::<Ut1Provider>()?;

    py_run!(py, sm, "import sys; sys.modules['nyx_space.time'] = sm");
    parent_module.add_submodule(sm)?;
    Ok(())
}
