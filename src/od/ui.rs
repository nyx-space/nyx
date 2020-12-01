use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, DimName};

pub use super::estimate::*;
pub use super::kalman::*;
pub use super::ranging::*;
pub use super::residual::*;
pub use super::snc::*;
pub use super::srif::*;
pub use super::*;

use crate::propagators::error_ctrl::ErrorCtrl;
use crate::propagators::PropInstance;
use crate::time::Duration;
use crate::State;

use std::fmt;
use std::marker::PhantomData;
use std::sync::mpsc::channel;

/// Defines the stopping condition for the smoother
#[derive(Clone, Copy, Debug)]
pub enum SmoothingArc {
    /// Stop smoothing when the gap between estimate is the provided floating point number in seconds
    TimeGap(Duration),
    /// Stop smoothing at the provided Epoch.
    Epoch(Epoch),
    /// Stop smoothing at the first prediction
    Prediction,
    /// Only stop once all estimates have been processed
    All,
}

impl fmt::Display for SmoothingArc {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            SmoothingArc::All => write!(f, "all estimates smoothed"),
            SmoothingArc::Epoch(e) => write!(f, "{}", e.as_gregorian_tai_str()),
            SmoothingArc::TimeGap(g) => write!(f, "time gap of {}", g),
            SmoothingArc::Prediction => write!(f, "first prediction"),
        }
    }
}

/// An orbit determination process. Note that everything passed to this structure is moved.
pub struct ODProcess<
    'a,
    // D: Dynamics<StateSize = Msr::StateSize>,
    D: Dynamics,
    E: ErrorCtrl,
    Msr: Measurement<StateSize = D::StateSize>,
    N: MeasurementDevice<D::StateType, Msr>,
    T: EkfTrigger,
    A: DimName,
    // S: D::StateType,
    K: Filter<D::StateSize, A, Msr::MeasurementSize, D::StateType>,
    // MsrIn,
> where
    DefaultAllocator: Allocator<f64, D::StateSize>
        + Allocator<f64, D::PropVecSize>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::StateSize>
        + Allocator<f64, Msr::StateSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, D::StateSize>
        + Allocator<f64, D::StateSize, Msr::MeasurementSize>
        + Allocator<f64, D::StateSize, D::StateSize>
        + Allocator<f64, A>
        + Allocator<f64, A, A>
        + Allocator<f64, D::StateSize, A>
        + Allocator<f64, A, D::StateSize>,
{
    /// PropInstance used for the estimation
    pub prop: PropInstance<'a, D, E>,
    /// Kalman filter itself
    pub kf: K,
    /// List of measurement devices used
    pub devices: Vec<N>,
    /// Whether or not these devices can make simultaneous measurements of the spacecraft
    pub simultaneous_msr: bool,
    /// Vector of estimates available after a pass
    pub estimates: Vec<K::Estimate>,
    /// Vector of residuals available after a pass
    pub residuals: Vec<Residual<Msr::MeasurementSize>>,
    pub ekf_trigger: T,
    init_state: D::StateType,
    _marker: PhantomData<A>,
}

impl<
        'a,
        D: Dynamics,
        E: ErrorCtrl,
        Msr: Measurement<StateSize = D::StateSize>,
        N: MeasurementDevice<D::StateType, Msr>,
        T: EkfTrigger,
        A: DimName,
        // S: State<D::StateSize>,
        K: Filter<D::StateSize, A, Msr::MeasurementSize, D::StateType>,
        // MsrIn,
    > ODProcess<'a, D, E, Msr, N, T, A, K>
where
    DefaultAllocator: Allocator<f64, D::StateSize>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::StateSize>
        + Allocator<f64, Msr::StateSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, D::StateSize>
        + Allocator<f64, D::StateSize, Msr::MeasurementSize>
        + Allocator<f64, D::StateSize, D::StateSize>
        + Allocator<f64, D::PropVecSize>
        + Allocator<f64, A>
        + Allocator<f64, A, A>
        + Allocator<f64, D::StateSize, A>
        + Allocator<f64, A, D::StateSize>,
{
    pub fn ekf(
        prop: PropInstance<'a, D, E>,
        kf: K,
        devices: Vec<N>,
        simultaneous_msr: bool,
        num_expected_msr: usize,
        trigger: T,
    ) -> Self {
        let init_state = prop.state();
        let mut estimates = Vec::with_capacity(num_expected_msr + 1);
        estimates.push(kf.previous_estimate().clone());
        Self {
            prop,
            kf,
            devices,
            simultaneous_msr,
            estimates,
            residuals: Vec::with_capacity(num_expected_msr),
            ekf_trigger: trigger,
            init_state,
            _marker: PhantomData::<A>,
        }
    }

    pub fn default_ekf(prop: PropInstance<'a, D, E>, kf: K, devices: Vec<N>, trigger: T) -> Self {
        let init_state = prop.state();
        let mut estimates = Vec::with_capacity(10_001);
        estimates.push(kf.previous_estimate().clone());
        Self {
            prop,
            kf,
            devices,
            simultaneous_msr: false,
            estimates,
            residuals: Vec::with_capacity(10_000),
            ekf_trigger: trigger,
            init_state,
            _marker: PhantomData::<A>,
        }
    }

    /// Allows to smooth the provided estimates. Returns the smoothed estimates or an error.
    ///
    /// Estimates must be ordered in chronological order. This function will smooth the
    /// estimates from the last in the list to the first one.
    pub fn smooth(&mut self, condition: SmoothingArc) -> Result<Vec<K::Estimate>, NyxError> {
        let l = self.estimates.len() - 1;
        let mut k = l - 1;

        info!("Smoothing {} estimates until {}", l + 1, condition);
        let mut smoothed = Vec::with_capacity(l + 1);
        // Set the first item of the smoothed estimates to the last estimate (we cannot smooth the very last estimate)
        smoothed.push(self.estimates[k + 1].clone());

        loop {
            // Borrow the previously smoothed estimate of the k+1 estimate
            let sm_est_kp1 = &smoothed[l - k - 1];
            let x_kp1_l = sm_est_kp1.state_deviation();
            let p_kp1_l = sm_est_kp1.covar();
            // Borrow the k-th estimate, which we're smoothing with the next estimate
            let est_k = &self.estimates[k];
            let x_k_k = &est_k.state_deviation();
            let p_k_k = &est_k.covar();
            // Borrow the k+1-th estimate, which we're smoothing with the next estimate
            let est_kp1 = &self.estimates[k + 1];

            let phi_kp1_k = est_kp1.stm();
            // let p_kp1_k = phi_kp1_k * p_k_k * phi_kp1_k.transpose(); // TODO: Add SNC here, which is effectively covar_bar!
            let p_kp1_k = est_kp1.predicted_covar();
            let p_kp1_k_inv = &p_kp1_k
                .try_inverse()
                .ok_or_else(|| NyxError::SingularCovarianceMatrix)?;
            // Compute Sk
            let sk = p_k_k * phi_kp1_k.transpose() * p_kp1_k_inv;
            // Compute smoothed estimate
            let x_k_l = x_k_k + &sk * (x_kp1_l - phi_kp1_k * x_k_k);
            // Compute smoothed covariance
            let p_k_l = p_k_k + &sk * (p_kp1_l - &est_kp1.covar()) * &sk.transpose();
            // Store into vector
            let mut smoothed_est_k = est_k.clone();
            // Compute the smoothed state deviation
            smoothed_est_k.set_state_deviation(x_k_l);
            // Compute the smoothed covariance
            smoothed_est_k.set_covar(p_k_l);
            // Move on
            smoothed.push(smoothed_est_k);
            if k == 0 {
                break;
            }
            k -= 1;
            // Check the smoother stopping condition
            let next_est = &self.estimates[k];
            match condition {
                SmoothingArc::Epoch(e) => {
                    // If the epoch of the next estimate is _before_ the stopping time, stop smoothing
                    if next_est.epoch() < e {
                        break;
                    }
                }
                SmoothingArc::TimeGap(gap_s) => {
                    if est_k.epoch() - next_est.epoch() > gap_s {
                        break;
                    }
                }
                SmoothingArc::Prediction => {
                    if next_est.predicted() {
                        break;
                    }
                }
                SmoothingArc::All => {}
            }
        }

        info!(
            "Condition reached after smoothing {} estimates ",
            smoothed.len()
        );

        // Now, let's add all of the other estimates so that the same indexing can be done
        // between all the estimates and the smoothed estimates
        if k > 0 {
            // Add the estimate that might have been skipped.
            // smoothed.push(self.estimates[k + 1].clone());
            loop {
                smoothed.push(self.estimates[k].clone());
                if k == 0 {
                    break;
                }
                k -= 1;
            }
        }

        // And reverse to maintain the order of estimates
        smoothed.reverse();
        Ok(smoothed)
    }

    /// Allows iterating on the filter solution. Requires specifying a smoothing condition to know where to stop the smoothing.
    pub fn iterate(
        &mut self,
        measurements: &[Msr],
        condition: SmoothingArc,
    ) -> Result<(), NyxError> {
        // First, smooth the estimates
        let smoothed = match self.smooth(condition) {
            Ok(smoothed) => smoothed,
            Err(e) => return Err(e),
        };
        // Reset the propagator
        self.prop.state = self.init_state;
        // Empty the estimates and add the first smoothed estimate as the initial estimate
        self.estimates = Vec::new();
        self.estimates.push(smoothed[0].clone());
        self.kf.set_previous_estimate(&smoothed[0]);
        // And re-run the filter
        self.process_measurements(measurements)
    }

    /// Allows processing all measurements with covariance mapping.
    ///
    /// Important notes:
    /// + the measurements have be to mapped to a fixed time corresponding to the step of the propagator
    pub fn process_measurements(&mut self, measurements: &[Msr]) -> Result<(), NyxError> {
        let (tx, rx) = channel();
        self.prop.tx_chan = Some(tx);
        assert!(
            !measurements.is_empty(),
            "must have at least one measurement"
        );
        // Start by propagating the estimator (on the same thread).
        let num_msrs = measurements.len();

        let prop_time = measurements[num_msrs - 1].epoch() - self.kf.previous_estimate().epoch();
        info!(
            "Navigation propagating for a total of {} seconds (~ {:.3} days)",
            prop_time,
            prop_time / 86_400.0
        );

        // Push the initial estimate
        let prev = self.kf.previous_estimate().clone();
        let mut prev_dt = prev.epoch();

        let mut reported = vec![false; 11];
        let mut arc_warned = false;

        info!(
            "Processing {} measurements with covariance mapping",
            num_msrs
        );

        for (msr_cnt, msr) in measurements.iter().enumerate() {
            let next_msr_epoch = msr.epoch();

            let delta_t = next_msr_epoch - prev_dt;
            self.prop.until_time_elapsed(delta_t)?;

            while let Ok(nominal_state) = rx.try_recv() {
                // Get the datetime and info needed to compute the theoretical measurement according to the model
                // let meas_input = self.prop.dynamics.to_measurement(&nominal_state);
                let dt = nominal_state.epoch();

                // Update the STM of the KF (needed between each measurement or time update)
                let stm = nominal_state.stm().unwrap();
                self.kf.update_stm(stm);

                // Check if we should do a time update or a measurement update
                if next_msr_epoch > dt {
                    if msr_cnt == 0 && !arc_warned {
                        warn!("OD arc starts prior to first measurement");
                        arc_warned = true;
                    }
                    // No measurement can be used here, let's just do a time update
                    debug!("time update {}", dt.as_gregorian_tai_str());
                    match self.kf.time_update(nominal_state) {
                        Ok(est) => {
                            if self.kf.is_extended() {
                                // self.prop.state
                                // self.prop.dynamics.set_estimated_state(
                                //     self.prop.dynamics.estimated_state() + est.state_deviation(),
                                // );
                            }
                            self.estimates.push(est);
                        }
                        Err(e) => return Err(e),
                    }
                } else {
                    // The epochs match, so this is a valid measurement to use
                    // Get the computed observations
                    for device in self.devices.iter() {
                        if let Some(computed_meas) = device.measure(&nominal_state) {
                            if computed_meas.visible() {
                                self.kf.update_h_tilde(computed_meas.sensitivity());

                                // Switch back from extended if necessary
                                if self.kf.is_extended() && self.ekf_trigger.disable_ekf(dt) {
                                    self.kf.set_extended(false);
                                    info!("EKF disabled @ {}", dt.as_gregorian_tai_str());
                                }

                                match self.kf.measurement_update(
                                    nominal_state,
                                    msr.observation(),
                                    computed_meas.observation(),
                                ) {
                                    Ok((est, res)) => {
                                        debug!(
                                            "msr update msr #{} {}",
                                            msr_cnt,
                                            dt.as_gregorian_tai_str()
                                        );

                                        // Switch to EKF if necessary, and update the dynamics and such
                                        // Note: we call enable_ekf first to ensure that the trigger gets
                                        // called in case it needs to save some information (e.g. the
                                        // StdEkfTrigger needs to store the time of the previous measurement).
                                        if self.ekf_trigger.enable_ekf(&est)
                                            && !self.kf.is_extended()
                                        {
                                            self.kf.set_extended(true);
                                            if !est.within_3sigma() {
                                                warn!(
                                                    "EKF enabled @ {} but filter DIVERGING",
                                                    dt.as_gregorian_tai_str()
                                                );
                                            } else {
                                                info!(
                                                    "EKF enabled @ {}",
                                                    dt.as_gregorian_tai_str()
                                                );
                                            }
                                        }
                                        if self.kf.is_extended() {
                                            self.prop.state =
                                                self.prop.state + est.state_deviation();
                                            // self.prop.dynamics.set_estimated_state(
                                            //     self.prop
                                            //         .dynamics
                                            //         .extract_estimated_state(&nominal_state)
                                            //         + est.state_deviation(),
                                            // );
                                        }
                                        self.estimates.push(est);
                                        self.residuals.push(res);
                                    }
                                    Err(e) => return Err(e),
                                }

                                // If we do not have simultaneous measurements from different devices
                                // then we don't need to check the visibility from other devices
                                // if one is in visibility.
                                if !self.simultaneous_msr {
                                    break;
                                }
                            }
                        }
                    }

                    let msr_prct = (10.0 * (msr_cnt as f64) / (num_msrs as f64)) as usize;
                    if !reported[msr_prct] {
                        info!(
                            "{:>3}% done ({:.0} measurements processed)",
                            10 * msr_prct,
                            msr_cnt
                        );
                        reported[msr_prct] = true;
                    }
                }
            }

            // Update the prev_dt for the next pass
            prev_dt = msr.epoch();
        }

        // Always report the 100% mark
        if !reported[10] {
            info!("{:>3}% done ({:.0} measurements processed)", 100, num_msrs);
        }

        Ok(())
    }

    /// Allows for covariance mapping without processing measurements
    pub fn map_covar(&mut self, end_epoch: Epoch) -> Result<(), NyxError> {
        let (tx, rx) = channel();
        self.prop.tx_chan = Some(tx);
        // Start by propagating the estimator (on the same thread).
        let prop_time = end_epoch - self.kf.previous_estimate().epoch();
        info!("Propagating for {} seconds", prop_time);

        self.prop.until_time_elapsed(prop_time)?;

        info!("Mapping covariance");

        while let Ok(nominal_state) = rx.try_recv() {
            // Update the STM of the KF (needed between each measurement or time update)
            // let stm = self.prop.dynamics.extract_stm(&nominal_state);
            self.kf.update_stm(nominal_state.stm()?);
            info!("final time update {:?}", nominal_state.epoch());
            match self.kf.time_update(nominal_state) {
                Ok(est) => {
                    if self.kf.is_extended() {
                        // let est_state = est.state_deviation().clone();
                        self.prop.state = self.prop.state + est.state_deviation();
                        // self.prop.dynamics.set_estimated_state(
                        //     self.prop.dynamics.extract_estimated_state(&nominal_state) + est_state,
                        // );
                    }
                    self.estimates.push(est);
                }
                Err(e) => return Err(e),
            }
        }

        Ok(())
    }
}

impl<
        'a,
        // D: Dynamics<StateSize = Msr::StateSize>,
        D: Dynamics,
        E: ErrorCtrl,
        Msr: Measurement<StateSize = D::StateSize>,
        N: MeasurementDevice<D::StateType, Msr>,
        A: DimName,
        K: Filter<D::StateSize, A, Msr::MeasurementSize, D::StateType>,
        // MsrIn,
    > ODProcess<'a, D, E, Msr, N, CkfTrigger, A, K>
where
    DefaultAllocator: Allocator<f64, D::StateSize>
        + Allocator<f64, D::PropVecSize>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::StateSize>
        + Allocator<f64, Msr::StateSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, D::StateSize>
        + Allocator<f64, D::StateSize, Msr::MeasurementSize>
        + Allocator<f64, D::StateSize, D::StateSize>
        + Allocator<f64, A>
        + Allocator<f64, A, A>
        + Allocator<f64, D::StateSize, A>
        + Allocator<f64, A, D::StateSize>,
{
    pub fn ckf(
        prop: PropInstance<'a, D, E>,
        kf: K,
        devices: Vec<N>,
        simultaneous_msr: bool,
        num_expected_msr: usize,
    ) -> Self {
        let init_state = prop.state();
        let mut estimates = Vec::with_capacity(num_expected_msr + 1);
        estimates.push(kf.previous_estimate().clone());
        Self {
            prop,
            kf,
            devices,
            simultaneous_msr,
            estimates,
            residuals: Vec::with_capacity(num_expected_msr),
            ekf_trigger: CkfTrigger {},
            init_state,
            _marker: PhantomData::<A>,
        }
    }

    pub fn default_ckf(prop: PropInstance<'a, D, E>, kf: K, devices: Vec<N>) -> Self {
        let init_state = prop.state();
        let mut estimates = Vec::with_capacity(10_001);
        estimates.push(kf.previous_estimate().clone());
        Self {
            prop,
            kf,
            devices,
            simultaneous_msr: false,
            estimates,
            residuals: Vec::with_capacity(10_000),
            ekf_trigger: CkfTrigger {},
            init_state,
            _marker: PhantomData::<A>,
        }
    }
}
/// A trait detailing when to switch to from a CKF to an EKF
pub trait EkfTrigger {
    fn enable_ekf<S, E, T: State<S>>(&mut self, est: &E) -> bool
    where
        S: DimName,
        E: Estimate<S, T>,
        DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>;

    /// Return true if the filter should not longer be as extended.
    /// By default, this returns false, i.e. when a filter has been switched to an EKF, it will
    /// remain as such.
    fn disable_ekf(&mut self, _epoch: Epoch) -> bool {
        false
    }
}

/// CkfTrigger will never switch a KF to an EKF
pub struct CkfTrigger;

impl EkfTrigger for CkfTrigger {
    fn enable_ekf<S, E, T: State<S>>(&mut self, _est: &E) -> bool
    where
        S: DimName,
        E: Estimate<S, T>,
        DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
    {
        false
    }
}

/// An EkfTrigger on the number of measurements processed and a time between measurements.
pub struct StdEkfTrigger {
    pub num_msrs: usize,
    /// In seconds!
    pub disable_time: Duration,
    /// Set to the sigma number needed to switch to the EKF (cf. 68–95–99.7 rule). If number is negative, this is ignored.
    pub within_sigma: f64,
    prev_msr_dt: Option<Epoch>,
    cur_msrs: usize,
}

impl StdEkfTrigger {
    pub fn new(num_msrs: usize, disable_time: Duration) -> Self {
        Self {
            num_msrs,
            disable_time,
            within_sigma: -1.0,
            prev_msr_dt: None,
            cur_msrs: 0,
        }
    }
}

impl EkfTrigger for StdEkfTrigger {
    fn enable_ekf<S, E, T: State<S>>(&mut self, est: &E) -> bool
    where
        S: DimName,
        E: Estimate<S, T>,
        DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
    {
        if !est.predicted() {
            // If this isn't a prediction, let's update the previous measurement time
            self.prev_msr_dt = Some(est.epoch());
        }
        self.cur_msrs += 1;
        self.cur_msrs >= self.num_msrs
            && ((self.within_sigma > 0.0 && est.within_sigma(self.within_sigma))
                || self.within_sigma <= 0.0)
    }

    fn disable_ekf(&mut self, epoch: Epoch) -> bool {
        // Return true if there is a prev msr dt, and the next measurement time is more than the disable time seconds away
        match self.prev_msr_dt {
            Some(prev_dt) => {
                if (epoch - prev_dt).abs() > self.disable_time {
                    self.cur_msrs = 0;
                    true
                } else {
                    false
                }
            }
            None => false,
        }
    }
}
