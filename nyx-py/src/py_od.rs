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

use super::py_md::PyTrajectory;
use anise::analysis::AnalysisError;
use anise::prelude::Almanac;
use nyx_space::{
    Spacecraft,
    io::ConfigError,
    od::{
        GroundStation,
        msr::TrackingDataArc,
        prelude::{TrackingArcSim, TrkConfig},
    },
};
use std::collections::BTreeMap;

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

    /// Generates measurements for the tracking arc using the defined strands
    ///
    /// # Warning
    /// This function will return an error if any of the devices defines as a scheduler.
    /// You must create the schedule first using `build_schedule` first.
    ///
    /// # Notes
    /// Although mutable, this function may be called several times to generate different measurements.
    ///
    /// # Algorithm
    /// For each tracking device, and for each strand within that device, sample the trajectory at the sample
    /// rate of the tracking device, adding a measurement whenever the spacecraft is visible.
    /// Build the measurements as a vector, ordered chronologically.
    ///
    fn generate_measurements(&mut self, almanac: Almanac) -> Result<TrackingDataArc, ConfigError> {
        self.inner.generate_measurements(almanac.into())
    }

    /// Builds the schedule provided the config. Requires the tracker to be a ground station.
    ///
    /// # Algorithm
    ///
    /// 1. For each tracking device:
    /// 2. Find when the vehicle trajectory has an elevation greater or equal to zero, and use that as the first start of the first tracking arc for this station
    /// 3. Find when the vehicle trajectory has an elevation less than zero (i.e. disappears below the horizon), after that initial epoch
    /// 4. Repeat 2, 3 until the end of the trajectory
    /// 5. Build each of these as "tracking strands" for this tracking device.
    /// 6. Organize all of the built tracking strands chronologically.
    /// 7. Iterate through all of the strands:
    ///    7.a. if that tracker is marked as `Greedy` and it ends after the start of the next strand, change the start date of the next strand.
    ///    7.b. if that tracker is marked as `Eager` and it ends after the start of the next strand, change the end date of the current strand.
    pub fn generate_schedule(
        &self,
        almanac: Almanac,
    ) -> Result<BTreeMap<String, TrkConfig>, AnalysisError> {
        self.inner.generate_schedule(almanac.into())
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }
}
