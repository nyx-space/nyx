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

use std::collections::BTreeMap;

use super::py_md::PyTrajectory;
use nyx_space::{
    Spacecraft,
    io::ConfigError,
    od::{
        GroundStation,
        prelude::{TrackingArcSim, TrkConfig},
    },
};

use pyo3::prelude::*;

#[derive(Clone)]
#[pyclass(from_py_object)]
pub struct GroundSpacecratTrackingArcSim {
    inner: TrackingArcSim<Spacecraft, GroundStation>,
}

#[pymethods]
impl GroundSpacecratTrackingArcSim {
    #[new]
    fn py_new(
        devices: BTreeMap<String, GroundStation>,
        trajectory: PyTrajectory,
        configs: BTreeMap<String, TrkConfig>,
        seed: Option<u64>,
    ) -> Result<Self, ConfigError> {
        let inner = match seed {
            Some(seed) => TrackingArcSim::with_seed(devices, trajectory.inner, configs, seed)?,
            None => TrackingArcSim::new(devices, trajectory.inner, configs)?,
        };
        Ok(Self { inner })
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }
}
