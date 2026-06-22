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
use nyx_space::od::estimate::KfEstimate;
use nyx_space::od::SpacecraftKalmanOD;
use nyx_space::od::prelude::ODError;

#[derive(Clone)]
#[pyclass(from_py_object, name = "SpacecraftODProcess")]
pub struct PySpacecraftODProcess {
    pub(crate) inner: SpacecraftKalmanOD,
}

#[pymethods]
impl PySpacecraftODProcess {
    #[pyo3(name = "process_arc")]
    fn py_process_arc(&self, initial_estimate: PyKfEstimate, arc: &TrackingDataArc) -> Result<PySpacecraftODSolution, PyErr> {
        let inner_res = self.inner.process_arc(initial_estimate.inner, arc).map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("{}", e)))?;
        Ok(PySpacecraftODSolution {
            inner: inner_res
        })
    }

    #[pyo3(name = "predict_until")]
    fn py_predict_until(&self, initial_estimate: PyKfEstimate, end_epoch: hifitime::Epoch) -> Result<PySpacecraftODSolution, PyErr> {
        let inner_res = self.inner.predict_until(initial_estimate.inner, end_epoch).map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("{}", e)))?;
        Ok(PySpacecraftODSolution {
            inner: inner_res
        })
    }

    #[pyo3(name = "predict_for")]
    fn py_predict_for(&self, initial_estimate: PyKfEstimate, duration: hifitime::Duration) -> Result<PySpacecraftODSolution, PyErr> {
        let inner_res = self.inner.predict_for(initial_estimate.inner, duration).map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("{}", e)))?;
        Ok(PySpacecraftODSolution {
            inner: inner_res
        })
    }
}

#[derive(Clone)]
#[pyclass(from_py_object, name = "KfEstimate")]
pub struct PyKfEstimate {
    pub(crate) inner: KfEstimate<Spacecraft>,
}

#[derive(Clone)]
#[pyclass(from_py_object, name = "Residual")]
pub struct PyResidual {
    pub(crate) inner: nyx_space::od::estimate::residual::Residual<nyx_space::linalg::Const<2>>,
}

#[pymethods]
impl PyResidual {
    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    fn __repr__(&self) -> String {
        format!("{} @ {self:p}", self.inner)
    }
}

#[derive(Clone)]
#[pyclass(from_py_object, name = "SpacecraftODSolution")]
pub struct PySpacecraftODSolution {
    pub(crate) inner: nyx_space::od::prelude::ODSolution<Spacecraft, KfEstimate<Spacecraft>, nyx_space::linalg::Const<2>, GroundStation>,
}

#[pymethods]
impl PySpacecraftODSolution {
    #[pyo3(name = "is_filter_run")]
    fn py_is_filter_run(&self) -> bool {
        self.inner.is_filter_run()
    }

    #[pyo3(name = "is_smoother_run")]
    fn py_is_smoother_run(&self) -> bool {
        self.inner.is_smoother_run()
    }

    #[pyo3(name = "to_traj")]
    fn py_to_traj(&self) -> Result<PyTrajectory, PyErr> {
        let traj = self.inner.to_traj().map_err(|e| pyo3::exceptions::PyValueError::new_err(format!("{}", e)))?;
        Ok(PyTrajectory { inner: traj })
    }

    #[pyo3(name = "accepted_residuals")]
    fn py_accepted_residuals(&self) -> Vec<PyResidual> {
        self.inner.accepted_residuals().into_iter().map(|r| PyResidual { inner: r }).collect()
    }

    #[pyo3(name = "rejected_residuals")]
    fn py_rejected_residuals(&self) -> Vec<PyResidual> {
        self.inner.rejected_residuals().into_iter().map(|r| PyResidual { inner: r }).collect()
    }
}

#[pymethods]
impl PyKfEstimate {
    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    fn __repr__(&self) -> String {
        format!("{} @ {self:p}", self.inner)
    }
}

#[derive(Clone)]
#[pyclass(from_py_object)]
pub struct GroundTrackingArcSim {
    inner: TrackingArcSim<Spacecraft, GroundStation>,
}

#[pymethods]
impl GroundTrackingArcSim {
    /// Initializes a simulated tracking architecture for a spacecraft.
    ///
    /// Flight dynamics operations demand strict determinism for Monte Carlo covariance analysis.
    /// When a seed is provided, the underlying pseudo-random number generator guarantees repeatable
    /// measurement noise characteristics. Omission of the seed defaults to system entropy, mimicking
    /// stochastic operational data streams.
    ///
    /// :param devices: Mapping of ground station identifiers to their respective physical definitions.
    /// :type devices: dict[str, GroundStation]
    /// :param trajectory: The deterministic ephemeris of the target spacecraft.
    /// :type trajectory: PyTrajectory
    /// :param configs: Tasking and measurement constraints per device (e.g., sample rates, cadences).
    /// :type configs: dict[str, TrkConfig]
    /// :param seed: Initialization seed for the underlying PRNG.
    /// :type seed: int | None
    /// :rtype: SpacecraftTrackingArcSim
    /// :raises ConfigError: If initialization fails due to malformed configurations.
    #[pyo3(signature=(devices, trajectory, configs, seed=None))]
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

    /// Simulates operational tracking data across predefined tracking strands.
    ///
    /// This function strictly demands that a schedule already exists (stored in the `config` field).
    /// If a device is configured as a scheduler but lacks pre-computed
    /// strands, this function will raise an error rather than implicitly hallucinating a tracking
    /// pass. Call `generate_schedule` to build a schedule first.
    /// For each tracking device, the trajectory is sampled at the specific hardware rate,
    /// synthesizing measurements only when the spacecraft is visible.
    ///
    /// :type almanac: Almanac
    /// :rtype: TrackingDataArc
    /// :raises ConfigError: If a scheduling configuration is present but the schedule was not built prior to execution.
    fn generate_measurements(&mut self, almanac: Almanac) -> Result<TrackingDataArc, ConfigError> {
        self.inner.generate_measurements(almanac.into())
    }

    /// Builds the schedule provided the config.
    ///
    /// # Algorithm
    ///
    /// 1. For each tracking device:
    /// 2. Find when the vehicle elevation above ground station mask is greater or equal to zero, and use that as the first start of the first tracking arc for this station
    /// 3. Find when the vehicle drops below the mask, after that initial epoch
    /// 4. Repeat 2, 3 until the end of the trajectory
    /// 5. Build each of these as "tracking strands" for this tracking device.
    /// 6. Organize all of the built tracking strands chronologically.
    /// 7. Iterate through all of the strands to adjust for tracker Greedy/Eager configuration.
    /// `Greedy` trackers will delay the start of subsequent station contacts, whereas `Eager` trackers will terminate
    /// current tracking strands prematurely to allow the next station to acquire.
    /// :type almanac: Almanac
    /// :rtype: dict[str, TrkConfig]
    /// :raises AnalysisError: If underlying location dataset injection or visibility computation fails.
    pub fn generate_schedule(
        &self,
        almanac: Almanac,
    ) -> Result<BTreeMap<String, TrkConfig>, AnalysisError> {
        self.inner.generate_schedule(almanac.into())
    }

    /// Builds a schedule using the generate_schedule function, and set that schedule in this instance's configuration.
    /// :type almanac: Almanac
    pub fn build_schedule(&mut self, almanac: Almanac) -> Result<(), AnalysisError> {
        self.inner.build_schedule(almanac.into())
    }

    #[getter]
    fn configs(&self) -> BTreeMap<String, TrkConfig> {
        self.inner.configs.clone()
    }

    #[getter]
    fn devices(&self) -> BTreeMap<String, GroundStation> {
        self.inner.devices.clone()
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    fn __repr__(&self) -> String {
        format!("{} @ {self:p}", self.inner)
    }
}
