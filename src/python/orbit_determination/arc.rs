/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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

use crate::cosmic::Cosm;
use crate::io::trajectory_data::TrajectoryLoader;
use crate::io::ExportCfg;
use crate::od::msr::RangeDoppler;
use crate::od::simulator::TrackingArcSim;
pub use crate::od::simulator::TrkConfig;
pub use crate::{io::ConfigError, od::prelude::GroundStation};
use crate::{NyxError, Orbit, Spacecraft};
use either::Either;
use pyo3::prelude::*;
use std::collections::BTreeMap;

#[derive(Clone)]
#[pyclass]
pub struct GroundTrackingArcSim {
    inner: Either<
        TrackingArcSim<Spacecraft, RangeDoppler, GroundStation>,
        TrackingArcSim<Orbit, RangeDoppler, GroundStation>,
    >,
}

#[pymethods]
impl GroundTrackingArcSim {
    /// Initializes a new tracking arc simulation from the provided devices, trajectory, and the random number generator seed.
    #[new]
    pub fn with_seed(
        devices: Vec<GroundStation>,
        trajectory: TrajectoryLoader,
        configs: BTreeMap<String, TrkConfig>,
        seed: u64,
    ) -> Result<Self, ConfigError> {
        // Try to convert the dynamic trajectory into a trajectory
        let inner = if let Ok(sc_traj) = trajectory.to_traj::<Spacecraft>() {
            let inner = TrackingArcSim::with_seed(devices, sc_traj, configs, seed)?;

            Either::Left(inner)
        } else if let Ok(traj) = trajectory.to_traj::<Orbit>() {
            let inner = TrackingArcSim::with_seed(devices, traj, configs, seed)?;

            Either::Right(inner)
        } else {
            return Err(ConfigError::InvalidConfig {
                msg: "Trajectory could neither be parsed as an orbit nor spacecraft trajectory"
                    .to_string(),
            });
        };

        Ok(Self { inner })
    }

    /// Generates simulated measurements and returns the path where the parquet file containing said measurements is stored.
    /// You may specify a metadata dictionary to be stored in the parquet file and whether the filename should include the timestamp.
    pub fn generate_measurements(
        &mut self,
        path: String,
        export_cfg: ExportCfg,
    ) -> Result<String, NyxError> {
        let cosm = Cosm::de438();
        let arc = match &mut self.inner {
            Either::Left(arc_sim) => arc_sim.generate_measurements(cosm)?,
            Either::Right(arc_sim) => arc_sim.generate_measurements(cosm)?,
        };

        // Save the tracking arc
        let maybe = arc.to_parquet(path, export_cfg);

        match maybe {
            Ok(path) => Ok(format!("{}", path.to_str().unwrap())),
            Err(e) => Err(NyxError::CustomError { msg: e.to_string() }),
        }
    }

    /// Generates a tracking schedule
    pub fn generate_schedule(&self) -> Result<BTreeMap<String, TrkConfig>, NyxError> {
        let cosm = Cosm::de438();
        match &self.inner {
            Either::Left(arc_sim) => arc_sim.generate_schedule(cosm),
            Either::Right(arc_sim) => arc_sim.generate_schedule(cosm),
        }
    }

    /// Builds a tracking schedule by generating it and storing it in this object.
    pub fn build_schedule(&mut self) -> Result<(), NyxError> {
        let cosm = Cosm::de438();
        match &mut self.inner {
            Either::Left(arc_sim) => arc_sim.build_schedule(cosm),
            Either::Right(arc_sim) => arc_sim.build_schedule(cosm),
        }
    }

    pub fn __repr__(&self) -> String {
        format!("{}", self.inner)
    }
}
