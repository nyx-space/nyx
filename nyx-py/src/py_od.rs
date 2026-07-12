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

use super::py_md::{Propagator, PyTrajectory};
use anise::analysis::AnalysisError;
use anise::ephemerides::ephemeris::{Ephemeris, LocalFrame};
use anise::prelude::Almanac;
use hifitime::{Duration, Epoch};
use ndarray::Array2;
use numpy::{PyArray2, PyReadonlyArray1};
use nyx_space::NyxError;
use nyx_space::io::{ExportCfg, InputOutputError};
use nyx_space::linalg::{Const, OVector};
use nyx_space::mc::{MvnSpacecraft, StateDispersion};
use nyx_space::od::blse::Estimate;
use nyx_space::od::estimate::KfEstimate;
use nyx_space::od::prelude::{
    KalmanVariant, ODError, ODSolution, Residual, SigmaRejection, SpacecraftKalmanScalarOD,
    TrackingArcSim, TrkConfig,
};
use nyx_space::od::snc::ProcessNoise3D;
use nyx_space::propagators::PropagationError;
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
use pyo3::types::PyType;
use std::collections::BTreeMap;

#[derive(Clone)]
#[pyclass(from_py_object, name = "ProcessNoise")]
pub struct PyProcessNoise {
    pub(crate) inner: ProcessNoise3D,
}

#[pymethods]
impl PyProcessNoise {
    #[classmethod]
    fn from_velocity_m_s(
        _cls: &Bound<'_, PyType>,
        vx_m_s: f64,
        vy_m_s: f64,
        vz_m_s: f64,
        noise_duration: Duration,
        disable_time: Duration,
        local_frame: Option<LocalFrame>,
    ) -> Self {
        let inner = ProcessNoise3D::from_velocity_km_s(
            &[vx_m_s * 1e-3, vy_m_s * 1e-3, vz_m_s * 1e-3],
            noise_duration,
            disable_time,
            local_frame,
        );
        Self { inner }
    }

    #[classmethod]
    fn from_accel_m_s2(
        _cls: &Bound<'_, PyType>,
        ax_m_s2: f64,
        ay_m_s2: f64,
        az_m_s2: f64,
        disable_time: Duration,
        local_frame: Option<LocalFrame>,
        x_decay_s: Option<f64>,
        y_decay_s: Option<f64>,
        z_decay_s: Option<f64>,
    ) -> Self {
        let mut inner = ProcessNoise3D::from_diagonal(
            &[ax_m_s2 * 1e-3, ay_m_s2 * 1e-3, az_m_s2 * 1e-3],
            disable_time,
            local_frame,
        );
        if x_decay_s.is_some() || y_decay_s.is_some() || z_decay_s.is_some() {
            inner.decay_diag = Some(vec![
                x_decay_s.unwrap_or_default(),
                y_decay_s.unwrap_or_default(),
                z_decay_s.unwrap_or_default(),
            ]);
        }

        Self { inner }
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    fn __repr__(&self) -> String {
        format!("{} @ {self:p}", self.inner)
    }

    fn __eq__(&self, other: Self) -> bool {
        self.inner == other.inner
    }
}

#[derive(Clone)]
#[pyclass(from_py_object, name = "SpacecraftODProcess")]
pub struct PySpacecraftODProcess {
    pub(crate) inner: SpacecraftKalmanScalarOD,
}

#[pymethods]
impl PySpacecraftODProcess {
    #[new]
    #[pyo3(signature=(prop, kf_variant, devices, sigma_reject=Some(SigmaRejection::default()), process_noise=None))]
    fn py_new(
        prop: Propagator,
        kf_variant: KalmanVariant,
        devices: BTreeMap<String, GroundStation>,
        sigma_reject: Option<SigmaRejection>,
        process_noise: Option<PyProcessNoise>,
    ) -> Result<Self, PropagationError> {
        let almanac = prop.almanac.clone();
        let mut inner = SpacecraftKalmanScalarOD::new(
            prop.build()?,
            kf_variant,
            sigma_reject,
            devices,
            almanac,
        );

        if let Some(pn) = process_noise {
            inner = inner.with_process_noise(pn.inner);
        }

        Ok(Self { inner })
    }
    /// Process the provided tracking arc for this orbit determination process.
    fn process_arc(
        &self,
        initial_estimate: PySpacecraftEstimate,
        arc: &TrackingDataArc,
    ) -> Result<PySpacecraftODSolution, PyErr> {
        let inner_res = self
            .inner
            .process_arc(initial_estimate.inner, arc)
            .map_err(|e| PyValueError::new_err(format!("{e}")))?;
        Ok(PySpacecraftODSolution { inner: inner_res })
    }

    /// Perform a time update. Continuously predicts the trajectory until the provided end epoch, with covariance mapping at each step.
    fn predict_until(
        &self,
        initial_estimate: PySpacecraftEstimate,
        end_epoch: Epoch,
    ) -> Result<PySpacecraftODSolution, PyErr> {
        let inner_res = self
            .inner
            .predict_until(initial_estimate.inner, end_epoch)
            .map_err(|e| PyValueError::new_err(format!("{e}")))?;
        Ok(PySpacecraftODSolution { inner: inner_res })
    }

    /// Perform a time update. Continuously predicts the trajectory for the provided duration, with covariance mapping at each step.
    fn predict_for(
        &self,
        initial_estimate: PySpacecraftEstimate,
        duration: Duration,
    ) -> Result<PySpacecraftODSolution, PyErr> {
        let inner_res = self
            .inner
            .predict_for(initial_estimate.inner, duration)
            .map_err(|e| PyValueError::new_err(format!("{e}")))?;
        Ok(PySpacecraftODSolution { inner: inner_res })
    }

    #[getter]
    fn sigma_rejection(&self) -> Option<SigmaRejection> {
        self.inner.sigma_reject
    }

    #[setter]
    fn set_sigma_rejection(&mut self, sigma_reject: Option<SigmaRejection>) {
        self.inner.sigma_reject = sigma_reject;
    }

    #[getter]
    fn variant(&self) -> KalmanVariant {
        self.inner.kf_variant
    }

    #[setter]
    fn set_variant(&mut self, kf_variant: KalmanVariant) {
        self.inner.kf_variant = kf_variant
    }
}

#[derive(Clone)]
#[pyclass(from_py_object, name = "SpacecraftEstimate")]
pub struct PySpacecraftEstimate {
    pub(crate) inner: KfEstimate<Spacecraft>,
}

#[pymethods]
impl PySpacecraftEstimate {
    /// Initializes a new filter estimate from the nominal state (not dispersed) and the diagonal of the covariance
    #[classmethod]
    fn from_diag(
        _cls: &Bound<'_, PyType>,
        nominal: Spacecraft,
        diag: PyReadonlyArray1<f64>,
    ) -> PyResult<Self> {
        if diag.len()? != 9 {
            return Err(PyValueError::new_err(format!(
                "Expected a numpy array of size 9, got {}",
                diag.len()?
            )));
        }

        let diag_vec = OVector::<f64, Const<9>>::from_iterator(diag.as_slice()?.iter().copied());

        let inner = KfEstimate::from_diag(nominal, diag_vec);

        Ok(Self { inner })
    }

    /// Generates an initial Kalman filter state estimate dispersed from the nominal state using the provided standard deviation parameters.
    ///
    /// The resulting estimate will have a diagonal covariance matrix constructed from the variances of each parameter.
    #[pyo3(signature=(nominal_state, dispersions, seed=None))]
    #[classmethod]
    fn from_dispersions(
        _cls: &Bound<'_, PyType>,
        nominal_state: Spacecraft,
        dispersions: Vec<StateDispersion>,
        seed: Option<u128>,
    ) -> Result<Self, NyxError> {
        let inner = KfEstimate::from_dispersions(nominal_state, dispersions, seed)
            .map_err(|e| NyxError::GenericError { msg: e.to_string() })?;

        Ok(Self { inner })
    }

    /// Builds a multivariate random variable spacecraft from this estimate's nominal state and covariance, zero mean.
    fn to_random_variable(&self) -> Result<MvnSpacecraft, NyxError> {
        self.inner.to_random_variable().map_err(|e| *e)
    }

    #[getter]
    fn predicted(&self) -> bool {
        self.inner.predicted
    }

    #[getter]
    fn nominal_state(&self) -> Spacecraft {
        self.inner.nominal_state
    }

    #[getter]
    fn state(&self) -> Spacecraft {
        self.inner.state()
    }

    #[getter]
    fn state_deviations(&self) -> Vec<f64> {
        self.inner
            .state_deviation
            .iter()
            .copied()
            .collect::<Vec<f64>>()
    }

    /// Returns the 9x9 covariance of this estimate.
    /// :rtype: numpy.ndarray
    #[getter]
    fn covariance<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        // Extract data from SMatrix (column-major order, hence the transpose)
        let data: Vec<f64> = self.inner.covar.transpose().iter().copied().collect();

        // Create an ndarray Array2 (row-major order)
        let state_dcm =
            Array2::from_shape_vec((9, 9), data).expect("9x9 matrix always has 81 elements");

        let pt_state_dcm = PyArray2::<f64>::from_owned_array(py, state_dcm);

        Ok(pt_state_dcm)
    }

    /// Returns whether this estimate is within some bound
    /// The 68-95-99.7 rule is a good way to assess whether the filter is operating normally
    fn within_sigma(&self, sigma: f64) -> bool {
        self.inner.within_sigma(sigma)
    }

    /// Returns whether this estimate is within three sigmas
    fn within_3sigma(&self) -> bool {
        self.inner.within_3sigma()
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }

    fn __repr__(&self) -> String {
        format!("{} @ {self:p}", self.inner)
    }

    fn __eq__(&self, other: &Self) -> bool {
        self.inner == other.inner
    }
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

    /// Returns the normalized innovation squared (NIS) as the norm squares of the whitened residual
    fn nis(&self) -> f64 {
        self.inner.nis()
    }

    /// Returns the whitened residual for this measurement type, if available
    fn whitened_residual(&self, msr_type: MeasurementType) -> Option<f64> {
        self.inner.whitened_resid(msr_type)
    }

    /// Returns the real observation for this measurement type, if available
    fn real_obs(&self, msr_type: MeasurementType) -> Option<f64> {
        self.inner.real_obs(msr_type)
    }

    /// Returns the computed/expected observation for this measurement type, if available
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
            .map_err(|e| PyValueError::new_err(format!("{e}")))?;
        Ok(PyTrajectory { inner: traj })
    }

    /// Export OD solutions, gains, ratios, residuals, sigmas, etc. to parquet
    fn to_parquet(&self, path: &str, cfg: ExportCfg) -> Result<String, ODError> {
        self.inner
            .to_parquet(path, cfg)
            .map(|path| path.to_string_lossy().to_string())
    }

    /// Export to an ANISE ephemeris, which can be converted to a CCSDS OEM
    fn to_ephemeris(&self, object_id: String) -> Ephemeris {
        self.inner.to_ephemeris(object_id)
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

    /// Smoothes this OD solution, returning a new OD solution and the filter-smoother consistency ratios, with updated **postfit** residuals, and where the ratio now represents the filter-smoother consistency ratio.
    ///
    /// Notes:
    ///  1. Gains will be scrubbed because the smoother process does not recompute the gain.
    ///  2. Prefit residuals, ratios, and measurement covariances are not updated, as these depend on the filtering process.
    ///  3. Note: this function consumes the current OD solution to prevent reusing the wrong one.
    ///
    ///
    /// To assess whether the smoothing process improved the solution, compare the RMS of the postfit residuals from the filter and the smoother process.
    ///
    /// # Filter-Smoother consistency ratio
    ///
    /// The **filter-smoother consistency ratio** is used to evaluate the consistency between the state estimates produced by a filter (e.g., Kalman filter) and a smoother.
    /// This ratio is called "filter smoother consistency test" in the ODTK MathSpec.
    ///
    /// It is computed as follows:
    ///
    /// #### 1. Define the State Estimates
    /// **Filter state estimate**:
    /// $ \hat{X}_{f,k} $
    /// This is the state estimate at time step $ k $ from the filter.
    ///
    /// **Smoother state estimate**:
    /// $ \hat{X}_{s,k} $
    /// This is the state estimate at time step $ k $ from the smoother.
    ///
    /// #### 2. Define the Covariances
    ///
    /// **Filter covariance**:
    /// $ P_{f,k} $
    /// This is the covariance of the state estimate at time step $ k $ from the filter.
    ///
    /// **Smoother covariance**:
    /// $ P_{s,k} $
    /// This is the covariance of the state estimate at time step $ k $ from the smoother.
    ///
    /// #### 3. Compute the Differences
    ///
    /// **State difference**:
    /// $ \Delta X_k = \hat{X}_{s,k} - \hat{X}_{f,k} $
    ///
    /// **Covariance difference**:
    /// $ \Delta P_k = P_{s,k} - P_{f,k} $
    ///
    /// #### 4. Calculate the Consistency Ratio
    /// For each element $ i $ of the state vector, compute the ratio:
    ///
    /// $$
    /// R_{i,k} = \frac{\Delta X_{i,k}}{\sqrt{\Delta P_{i,k}}}
    /// $$
    ///
    /// #### 5. Evaluate Consistency
    /// - If $ |R_{i,k}| \leq 3 $ for all $ i $ and $ k $, the filter-smoother consistency test is satisfied, indicating good consistency.
    /// - If $ |R_{i,k}| > 3 $ for any $ i $ or $ k $, the test fails, suggesting potential modeling inconsistencies or issues with the estimation process.
    ///
    fn smooth(&self, almanac: &Almanac) -> Result<Self, ODError> {
        let inner = self.clone().inner.smooth(almanac)?;

        Ok(Self { inner })
    }

    #[classmethod]
    fn from_parquet(
        _cls: &Bound<'_, PyType>,
        path: &str,
        devices: BTreeMap<String, GroundStation>,
    ) -> Result<Self, InputOutputError> {
        let inner = ODSolution::from_parquet(path, devices)?;

        Ok(Self { inner })
    }

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
    fn generate_measurements(&mut self, almanac: &Almanac) -> Result<TrackingDataArc, ConfigError> {
        self.inner.generate_measurements(almanac)
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
        almanac: &Almanac,
    ) -> Result<BTreeMap<String, TrkConfig>, AnalysisError> {
        self.inner.generate_schedule(almanac)
    }

    /// Builds a schedule using the generate_schedule function, and set that schedule in this instance's configuration.
    /// :type almanac: Almanac
    pub fn build_schedule(&mut self, almanac: &Almanac) -> Result<(), AnalysisError> {
        self.inner.build_schedule(almanac)
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
