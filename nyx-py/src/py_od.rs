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
use hifitime::{Duration, Epoch};
use nyx_space::linalg::Const;
use nyx_space::od::estimate::KfEstimate;
use nyx_space::od::prelude::{
    ODError, ODSolution, Residual, SpacecraftKalmanScalarOD, TrackingArcSim, TrkConfig,
};
use nyx_space::{
    Spacecraft,
    io::ConfigError,
    od::{
        GroundStation,
        msr::{MeasurementType, TrackingDataArc},
    },
};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use std::collections::BTreeMap;

#[derive(Clone)]
#[pyclass(from_py_object, name = "SpacecraftODProcess")]
pub struct PySpacecraftODProcess {
    pub(crate) inner: SpacecraftKalmanScalarOD,
}

#[pymethods]
impl PySpacecraftODProcess {
    #[pyo3(name = "process_arc")]
    fn py_process_arc(
        &self,
        initial_estimate: PyKfEstimate,
        arc: &TrackingDataArc,
    ) -> Result<PySpacecraftODSolution, PyErr> {
        let inner_res = self
            .inner
            .process_arc(initial_estimate.inner, arc)
            .map_err(|e| PyValueError::new_err(format!("{}", e)))?;
        Ok(PySpacecraftODSolution { inner: inner_res })
    }

    #[pyo3(name = "predict_until")]
    fn py_predict_until(
        &self,
        initial_estimate: PyKfEstimate,
        end_epoch: Epoch,
    ) -> Result<PySpacecraftODSolution, PyErr> {
        let inner_res = self
            .inner
            .predict_until(initial_estimate.inner, end_epoch)
            .map_err(|e| PyValueError::new_err(format!("{}", e)))?;
        Ok(PySpacecraftODSolution { inner: inner_res })
    }

    #[pyo3(name = "predict_for")]
    fn py_predict_for(
        &self,
        initial_estimate: PyKfEstimate,
        duration: Duration,
    ) -> Result<PySpacecraftODSolution, PyErr> {
        let inner_res = self
            .inner
            .predict_for(initial_estimate.inner, duration)
            .map_err(|e| PyValueError::new_err(format!("{}", e)))?;
        Ok(PySpacecraftODSolution { inner: inner_res })
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
    pub(crate) inner: Residual<Const<1>>,
}

#[pymethods]
impl PyResidual {
    #[getter]
    fn epoch(&self) -> Epoch {
        self.inner.epoch
    }

    #[getter]
    fn prefit(&self) -> Vec<f64> {
        self.inner.prefit.as_slice().to_vec()
    }

    #[getter]
    fn postfit(&self) -> Vec<f64> {
        self.inner.postfit.as_slice().to_vec()
    }

    #[getter]
    fn ratio(&self) -> f64 {
        self.inner.ratio
    }

    #[getter]
    fn rejected(&self) -> bool {
        self.inner.rejected
    }

    #[getter]
    fn tracker(&self) -> Option<String> {
        self.inner.tracker.clone()
    }

    #[getter]
    fn nis(&self) -> f64 {
        self.inner.nis()
    }

    fn whitened_residual(&self, msr_type: MeasurementType) -> Option<f64> {
        self.inner.whitened_resid(msr_type)
    }

    fn real_obs(&self, msr_type: MeasurementType) -> Option<f64> {
        self.inner.real_obs(msr_type)
    }

    fn computed_obs(&self, msr_type: MeasurementType) -> Option<f64> {
        self.inner.computed_obs(msr_type)
    }

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
    pub(crate) inner: ODSolution<Spacecraft, KfEstimate<Spacecraft>, Const<1>, GroundStation>,
}

#[pymethods]
impl PySpacecraftODSolution {
    fn is_filter_run(&self) -> bool {
        self.inner.is_filter_run()
    }

    fn is_smoother_run(&self) -> bool {
        self.inner.is_smoother_run()
    }

    fn to_traj(&self) -> Result<PyTrajectory, PyErr> {
        let traj = self
            .inner
            .to_traj()
            .map_err(|e| PyValueError::new_err(format!("{}", e)))?;
        Ok(PyTrajectory { inner: traj })
    }

    fn accepted_residuals(&self) -> Vec<PyResidual> {
        self.inner
            .accepted_residuals()
            .into_iter()
            .map(|r| PyResidual { inner: r })
            .collect()
    }

    fn rejected_residuals(&self) -> Vec<PyResidual> {
        self.inner
            .rejected_residuals()
            .into_iter()
            .map(|r| PyResidual { inner: r })
            .collect()
    }
    /// Returns the root mean square of the prefit residuals
    pub fn rms_prefit_residuals(&self) -> f64 {
        self.inner.rms_prefit_residuals()
    }

    /// Returns the root mean square of the postfit residuals
    pub fn rms_postfit_residuals(&self) -> f64 {
        self.inner.rms_postfit_residuals()
    }

    /// Returns the root mean square of the prefit residual ratios
    pub fn rms_residual_ratios(&self) -> f64 {
        self.inner.rms_residual_ratios()
    }

    /// Computes the fraction of residual ratios that lie within ±threshold.
    pub fn residual_ratio_within_threshold(&self, threshold: f64) -> Result<f64, ODError> {
        self.inner.residual_ratio_within_threshold(threshold)
    }

    /// Computes the Kolmogorov–Smirnov statistic for the aggregated residual ratios of the accepted residuals.
    ///
    /// Returns Ok(ks_statistic) if residuals are available.
    pub fn ks_test_normality(&self) -> Result<f64, ODError> {
        self.inner.ks_test_normality()
    }

    /// Checks whether the whitened residuals of the accepted residuals pass a normality test at a given significance level `alpha`, default to 0.05.
    ///
    /// This uses a simplified KS-test threshold: D_alpha = c(α) / √n.
    /// For example, for α = 0.05, c(α) is approximately 1.36.
    /// α=0.05 means a 5% probability of Type I error (incorrectly rejecting the null hypothesis when it is true).
    /// This threshold is standard in many statistical tests to balance sensitivity and false positives
    ///
    /// The computation of the c(alpha) is from https://en.wikipedia.org/w/index.php?title=Kolmogorov%E2%80%93Smirnov_test&oldid=1280260701#Two-sample_Kolmogorov%E2%80%93Smirnov_test
    ///
    /// Returns Ok(true) if the residuals are consistent with a normal distribution,
    /// Ok(false) if not, or None if no residuals are available.
    #[pyo3(signature=(alpha=None))]
    pub fn is_normal(&self, alpha: Option<f64>) -> Result<bool, ODError> {
        self.inner.is_normal(alpha)
    }

    /// Checks whether the filter innovations are statistically consistent
    /// by performing a Chi-squared test on the Normalized Innovation Squared (NIS).
    ///
    /// For each accepted residual, NIS is computed as:
    /// ```text
    ///     prefit^T * S_k^-1 * prefit
    /// ```
    ///
    /// The sum of NIS values should fall within the confidence interval of a
    /// Chi-squared distribution with degrees of freedom `k = n * m`, where `n`
    /// is the number of residuals and `m` is the measurement dimension.
    ///
    /// Returns Ok(true) if the filter is consistent, Ok(false) if the filter
    /// is over-confident or under-confident, or an error if no residuals are available.
    #[pyo3(signature=(alpha=None))]
    pub fn is_nis_consistent(&self, alpha: Option<f64>) -> Result<bool, ODError> {
        self.inner.is_nis_consistent(alpha)
    }

    /// Checks whether the filter estimates are statistically consistent
    /// by performing a Chi-squared test on the Normalized Estimation Error Squared (NEES).
    ///
    /// For each estimate, NEES is computed as:
    /// ```text
    ///     error^T * P^-1 * error
    /// ```
    /// where `error` is the difference between the estimated state and the true state,
    /// and `P` is the estimated state covariance matrix.
    ///
    /// The sum of NEES values should fall within the confidence interval of a
    /// Chi-squared distribution with degrees of freedom `k = n * dim`, where `n`
    /// is the number of estimates and `dim` is the state dimension.
    ///
    /// Returns Ok(true) if the filter is consistent, Ok(false) if the filter
    /// is over-confident or under-confident, or an error if no estimates are available.
    #[pyo3(signature=(truth_traj, alpha=None))]
    pub fn is_nees_consistent(
        &self,
        truth_traj: &PyTrajectory,
        alpha: Option<f64>,
    ) -> Result<bool, ODError> {
        self.inner.is_nees_consistent(&truth_traj.inner, alpha)
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
