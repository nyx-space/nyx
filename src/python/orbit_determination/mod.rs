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

use std::collections::HashMap;

use crate::cosmic::Cosm;
use crate::io::tracking_data::DynamicTrackingArc;
use crate::io::trajectory_data::DynamicTrajectory;
use crate::od::msr::StdMeasurement;
use crate::od::simulator::arc::TrackingArcSim;
pub use crate::od::simulator::TrkConfig;
pub use crate::{io::ConfigError, od::prelude::GroundStation};
use crate::{NyxError, Spacecraft};
use pyo3::{prelude::*, py_run};
pub(crate) mod estimate;
mod ground_station;
mod process;
mod trkconfig;

use estimate::OrbitEstimate;
use process::process_tracking_arc;

pub(crate) fn register_od(py: Python<'_>, parent_module: &PyModule) -> PyResult<()> {
    let sm = PyModule::new(py, "nyx_space.orbit_determination")?;

    sm.add_class::<GroundStation>()?;
    sm.add_class::<GroundTrackingArcSim>()?;
    sm.add_class::<DynamicTrackingArc>()?;
    sm.add_class::<TrkConfig>()?;
    sm.add_class::<OrbitEstimate>()?;
    sm.add_function(wrap_pyfunction!(process_tracking_arc, sm)?)?;

    py_run!(
        py,
        sm,
        "import sys; sys.modules['nyx_space.orbit_determination'] = sm"
    );
    parent_module.add_submodule(sm)?;
    Ok(())
}

#[derive(Clone)]
#[pyclass]
pub struct GroundTrackingArcSim {
    inner: TrackingArcSim<Spacecraft, StdMeasurement, GroundStation>,
}

#[pymethods]
impl GroundTrackingArcSim {
    /// Initializes a new tracking arc simulation from the provided devices, trajectory, and the random number generator seed.
    #[new]
    pub fn with_seed(
        devices: Vec<GroundStation>,
        trajectory: DynamicTrajectory,
        configs: HashMap<String, TrkConfig>,
        seed: u64,
    ) -> Result<Self, NyxError> {
        // Try to convert the dynamic trajectory into an Orbit trajectory
        let traj = trajectory
            .to_traj()
            .map_err(|e| NyxError::CustomError(e.to_string()))?;

        let inner = TrackingArcSim::with_seed(devices, traj, configs, seed)
            .map_err(NyxError::ConfigError)?;

        Ok(Self { inner })
    }

    /// Generates simulated measurements and returns the path where the parquet file containing said measurements is stored.
    /// You may specify a metadata dictionary to be stored in the parquet file and whether the filename should include the timestamp.
    #[pyo3(text_signature = "(path, metadata=None, timestamp=False)")]
    pub fn generate_measurements(
        &mut self,
        path: String,
        metadata: Option<HashMap<String, String>>,
        timestamp: bool,
    ) -> Result<String, NyxError> {
        let cosm = Cosm::de438();
        let arc = self.inner.generate_measurements(cosm)?;

        // Save the tracking arc
        let maybe = arc.to_parquet(path, metadata, timestamp);

        match maybe {
            Ok(path) => Ok(format!("{}", path.to_str().unwrap())),
            Err(e) => Err(NyxError::CustomError(e.to_string())),
        }
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }
}
