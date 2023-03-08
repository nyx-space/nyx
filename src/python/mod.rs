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

use crate::io::tracking_data::DynamicTrackingArc;
use crate::io::trajectory_data::DynamicTrajectory;
use crate::io::ConfigError;
use crate::NyxError;
use hifitime::leap_seconds::{LatestLeapSeconds, LeapSecondsFile};
use hifitime::prelude::*;
use hifitime::ut1::Ut1Provider;
use pyo3::{exceptions::PyException, prelude::*};

pub(crate) mod cosmic;
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
fn nyx_space(py: Python, m: &PyModule) -> PyResult<()> {
    pyo3_log::init();

    register_time_module(py, m)?;
    register_od(py, m)?;
    register_md(py, m)?;
    register_cosmic(py, m)?;

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
    let sm = PyModule::new(py, "orbit_determination")?;

    sm.add_class::<orbit_determination::GroundStation>()?;
    sm.add_class::<orbit_determination::GroundTrackingArcSim>()?;
    sm.add_class::<DynamicTrackingArc>()?;
    sm.add_class::<orbit_determination::TrkConfig>()?;

    parent_module.add_submodule(sm)?;
    Ok(())
}

/// Mission design
fn register_md(py: Python<'_>, parent_module: &PyModule) -> PyResult<()> {
    let sm = PyModule::new(py, "mission_design")?;

    sm.add_class::<DynamicTrajectory>()?;

    parent_module.add_submodule(sm)?;
    Ok(())
}

/// nyx_space.cosmic
fn register_cosmic(py: Python<'_>, parent_module: &PyModule) -> PyResult<()> {
    let sm = PyModule::new(py, "cosmic")?;
    sm.add_class::<cosmic::Cosm>()?;
    sm.add_class::<cosmic::Bodies>()?;
    sm.add_class::<cosmic::Frame>()?;
    sm.add_class::<cosmic::Orbit>()?;

    parent_module.add_submodule(sm)?;
    Ok(())
}
