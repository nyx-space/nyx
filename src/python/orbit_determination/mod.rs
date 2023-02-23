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

use crate::io::odp::Cosm;
use crate::io::stations::StationSerde;
use crate::io::{ConfigRepr, Configurable};
use crate::od::msr::StdMeasurement;
use crate::od::simulator::arc::TrackingArcSim;
use crate::Orbit;
pub use crate::{io::ConfigError, od::ui::GroundStation};
use pyo3::prelude::*;

#[pymethods]
impl GroundStation {
    #[staticmethod]
    fn load_yaml(path: &str) -> Result<Self, ConfigError> {
        let serde = StationSerde::load_yaml(path)?;

        // Create a new Cosm until ANISE switch
        let cosm = Cosm::de438();

        GroundStation::from_config(&serde, cosm)
    }

    #[staticmethod]
    fn load_many_yaml(path: &str) -> Result<Vec<Self>, ConfigError> {
        let stations = StationSerde::load_many_yaml(path)?;

        // Create a new Cosm until ANISE switch
        let cosm = Cosm::de438();

        let mut selves = Vec::with_capacity(stations.len());

        for serde in stations {
            selves.push(GroundStation::from_config(&serde, cosm.clone())?);
        }

        Ok(selves)
    }

    fn __repr__(&self) -> String {
        format!("{self:?}")
    }

    fn __str__(&self) -> String {
        format!("{self}")
    }
}

#[pyclass]
pub struct GroundTrackingArcSim {
    inner: TrackingArcSim<Orbit, StdMeasurement, GroundStation>,
}

// #[pymethod]
// impl GroundTrackingArcSim {
//     pub fn with_seed(devices: Vec<D>, trajectory: Traj<Msr::State>, seed: u64) -> Self {
//         let rng = StdRng::seed_from_u64(seed);

//         Self::with_rng(devices, trajectory, rng)
//     }

//     pub fn new(devices: HashMap<String, D>, trajectory: Traj<Msr::State>) -> Self {
//         let rng = StdRng::from_entropy();

//         Self::with_rng(devices, trajectory, rng)
//     }
// }
