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
use crate::io::ExportCfg;
use crate::od::noise::GaussMarkov;
use crate::od::process::{FltResid, IterationConf};
pub use crate::od::simulator::TrkConfig;
pub use crate::{io::ConfigError, od::prelude::GroundStation};
use pyo3::{prelude::*, py_run};
mod arc;
pub(crate) mod estimate;
mod ground_station;
mod process;
mod trkconfig;

use estimate::OrbitEstimate;
use process::{predictor, process_tracking_arc};

pub(crate) fn register_od(py: Python<'_>, parent_module: &PyModule) -> PyResult<()> {
    let sm = PyModule::new(py, "_nyx_space.orbit_determination")?;

    sm.add_class::<GroundStation>()?;
    sm.add_class::<arc::GroundTrackingArcSim>()?;
    sm.add_class::<DynamicTrackingArc>()?;
    sm.add_class::<TrkConfig>()?;
    sm.add_class::<OrbitEstimate>()?;
    sm.add_class::<GaussMarkov>()?;
    sm.add_class::<FltResid>()?;
    sm.add_class::<IterationConf>()?;
    sm.add_class::<ExportCfg>()?;
    sm.add_function(wrap_pyfunction!(process_tracking_arc, sm)?)?;
    sm.add_function(wrap_pyfunction!(predictor, sm)?)?;

    py_run!(
        py,
        sm,
        "import sys; sys.modules['nyx_space.orbit_determination'] = sm"
    );
    parent_module.add_submodule(sm)?;
    Ok(())
}
