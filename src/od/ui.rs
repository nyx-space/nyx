use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, DimName};

pub use super::estimate::*;
pub use super::kalman::*;
pub use super::ranging::*;
pub use super::residual::*;
pub use super::srif::*;
pub use super::*;

use crate::propagators::error_ctrl::ErrorCtrl;
use crate::propagators::Propagator;

use std::marker::PhantomData;
use std::sync::mpsc::channel;

/// An orbit determination process. Note that everything passed to this structure is moved.
pub struct ODProcess<
    'a,
    D: Estimable<MsrIn, LinStateSize = Msr::StateSize>,
    E: ErrorCtrl,
    Msr: Measurement,
    N: MeasurementDevice<MsrIn, Msr>,
    T: EkfTrigger,
    A: DimName,
    K: Filter<D::LinStateSize, A, Msr::MeasurementSize, D::StateType>,
    MsrIn,
> where
    D::StateType: EstimableState<Msr::StateSize>,
    DefaultAllocator: Allocator<f64, D::StateSize>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::StateSize>
        + Allocator<f64, Msr::StateSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, D::LinStateSize>
        + Allocator<f64, D::LinStateSize, Msr::MeasurementSize>
        + Allocator<f64, D::LinStateSize, D::LinStateSize>
        + Allocator<f64, A, A>
        + Allocator<f64, D::LinStateSize, A>
        + Allocator<f64, A, D::LinStateSize>,
{
    /// Propagator used for the estimation
    pub prop: Propagator<'a, D, E>,
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
    _marker: PhantomData<A>,
}

impl<
        'a,
        D: Estimable<MsrIn, LinStateSize = Msr::StateSize>,
        E: ErrorCtrl,
        Msr: Measurement,
        N: MeasurementDevice<MsrIn, Msr>,
        T: EkfTrigger,
        A: DimName,
        K: Filter<D::LinStateSize, A, Msr::MeasurementSize, D::StateType>,
        MsrIn,
    > ODProcess<'a, D, E, Msr, N, T, A, K, MsrIn>
where
    D::StateType: EstimableState<Msr::StateSize>,
    DefaultAllocator: Allocator<f64, D::StateSize>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::StateSize>
        + Allocator<f64, Msr::StateSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, D::LinStateSize>
        + Allocator<f64, D::LinStateSize, Msr::MeasurementSize>
        + Allocator<f64, D::LinStateSize, D::LinStateSize>
        + Allocator<f64, A, A>
        + Allocator<f64, D::LinStateSize, A>
        + Allocator<f64, A, D::LinStateSize>,
{
    pub fn ekf(
        prop: Propagator<'a, D, E>,
        kf: K,
        devices: Vec<N>,
        simultaneous_msr: bool,
        num_expected_msr: usize,
        trigger: T,
    ) -> Self {
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
            _marker: PhantomData::<A>,
        }
    }

    pub fn default_ekf(prop: Propagator<'a, D, E>, kf: K, devices: Vec<N>, trigger: T) -> Self {
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
            _marker: PhantomData::<A>,
        }
    }

    /// Allows to smooth the provided estimates. If there is a result, it'll be a filter error.
    ///
    /// Estimates must be ordered in chronological order. This function will smooth the
    /// estimates from the last in the list to the first one.
    pub fn smooth(&mut self) -> Option<FilterError> {
        info!("Smoothing {} estimates", self.estimates.len());
        let mut smoothed = Vec::with_capacity(self.estimates.len());

        // Is this smoothing correct ? I should try using the previous STM and start are estimates N-2
        for estimate in self.estimates.iter().rev() {
            let mut sm_est = estimate.clone();
            // TODO: Ensure that SNC was _not_ enabled
            let mut stm_inv = estimate.stm().clone();
            if !stm_inv.try_inverse_mut() {
                return Some(FilterError::StateTransitionMatrixSingular);
            }
            sm_est.set_covar(&stm_inv * estimate.covar() * &stm_inv.transpose());
            sm_est.set_state_deviation(&stm_inv * estimate.state_deviation());
            smoothed.push(sm_est);
        }

        // And reverse to maintain order
        smoothed.reverse();
        // And store
        self.estimates = smoothed;

        None
    }

    /// Allows iterating on the filter solution
    pub fn iterate(&mut self, measurements: &[Msr], map_covar: bool) -> Option<FilterError> {
        // First, smooth the estimates
        self.smooth();
        // Get the first estimate post-smoothing
        let mut init_smoothed = self.estimates[0].clone();
        println!("{}", init_smoothed.epoch().as_gregorian_tai_str());
        // Reset the propagator
        self.prop.reset();
        let mut iterated_state = self.prop.dynamics.state_vector();
        for (i, x) in init_smoothed.state_deviation().iter().enumerate() {
            iterated_state[i] += x;
        }
        self.prop
            .dynamics
            .set_state(self.prop.dynamics.time(), &iterated_state);
        // Set the filter's initial state to this smoothed estimate
        init_smoothed.set_state_deviation(VectorN::<f64, Msr::StateSize>::zeros());
        self.kf.set_previous_estimate(&init_smoothed);
        // And re-run the filter
        if map_covar {
            self.process_measurements_covar(measurements)?;
        } else {
            self.process_measurements(measurements)?;
        }
        None
    }

    /// Allows processing all measurements without covariance mapping.
    pub fn process_measurements(&mut self, measurements: &[Msr]) -> Option<FilterError> {
        info!("Processing {} measurements", measurements.len());

        let mut prev_dt = self.kf.previous_estimate().epoch();
        let mut reported = vec![false; 11];
        let num_msrs = measurements.len();

        for (msr_cnt, real_meas) in measurements.iter().enumerate() {
            let next_epoch = real_meas.epoch();
            // Propagate the dynamics to the measurement, and then start the filter.
            let delta_time = next_epoch - prev_dt;
            prev_dt = next_epoch; // Update the epoch for the next computation
            self.prop.until_time_elapsed(delta_time);
            // Update the STM of the KF
            self.kf.update_stm(self.prop.dynamics.stm());
            let nominal_state = self.prop.state();
            let meas_input = self.prop.dynamics.to_measurement(&nominal_state);
            // Get the computed observations
            for device in self.devices.iter() {
                if let Some(computed_meas) = device.measure(&meas_input) {
                    if computed_meas.visible() {
                        self.kf.update_h_tilde(computed_meas.sensitivity());

                        // Switch back from extended if necessary
                        if self.kf.is_extended() && self.ekf_trigger.disable_ekf(real_meas.epoch())
                        {
                            self.kf.set_extended(false);
                            info!(
                                "EKF disabled @ {}",
                                real_meas.epoch().as_gregorian_tai_str()
                            );
                        }

                        match self.kf.measurement_update(
                            nominal_state,
                            real_meas.observation(),
                            computed_meas.observation(),
                        ) {
                            Ok((est, res)) => {
                                // Switch to extended if necessary, and update the dynamics and such
                                // Note: we call enable_ekf first to ensure that the trigger gets
                                // called in case it needs to save some information (e.g. the
                                // StdEkfTrigger needs to store the time of the previous measurement).
                                if self.ekf_trigger.enable_ekf(&est) && !self.kf.is_extended() {
                                    self.kf.set_extended(true);
                                    if !est.within_3sigma() {
                                        warn!(
                                            "EKF enabled @ {} but filter DIVERGING",
                                            real_meas.epoch().as_gregorian_tai_str()
                                        );
                                    } else {
                                        info!(
                                            "EKF enabled @ {}",
                                            real_meas.epoch().as_gregorian_tai_str()
                                        );
                                    }
                                }
                                if self.kf.is_extended() {
                                    self.prop.dynamics.set_estimated_state(
                                        self.prop.dynamics.estimated_state()
                                            + est.state_deviation(),
                                    );
                                }
                                self.estimates.push(est);
                                self.residuals.push(res);
                            }
                            Err(e) => return Some(e),
                        }
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
        // Always report the 100% mark
        if !reported[10] {
            info!("{:>3}% done ({:.0} measurements processed)", 100, num_msrs);
        }

        None
    }

    /// Allows processing all measurements with covariance mapping.
    ///
    /// Important notes:
    /// + the measurements have be to mapped to a fixed time corresponding to the step of the propagator
    pub fn process_measurements_covar(&mut self, measurements: &[Msr]) -> Option<FilterError> {
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
        self.estimates.push(prev);
        for msr in measurements {
            let delta_t = msr.epoch() - prev_dt;
            self.prop.until_time_elapsed(delta_t);
            prev_dt = msr.epoch();
        }
        info!(
            "Processing {} measurements with covariance mapping",
            num_msrs
        );

        let mut msr_cnt = 0_usize;
        let mut reported = vec![false; 11];

        while let Ok(nominal_state) = rx.try_recv() {
            // Get the datetime and info needed to compute the theoretical measurement according to the model
            let meas_input = self.prop.dynamics.to_measurement(&nominal_state);

            let dt = nominal_state.epoch();

            let mut num_msr_processed = 0_u32;
            loop {
                // Update the STM of the KF (needed between each measurement or time update)
                let stm = self.prop.dynamics.extract_stm(&nominal_state);
                self.kf.update_stm(stm);
                if msr_cnt < num_msrs {
                    // Get the next measurement
                    let real_meas = &measurements[msr_cnt];
                    let next_epoch = real_meas.epoch();
                    if next_epoch < dt {
                        // We missed a measurement! Let's try to catch up.
                        error!(
                            "Skipping msr #{}: nav and msr not in sync: msr.dt = {} \t nav.dt = {}",
                            msr_cnt,
                            next_epoch.as_gregorian_tai_str(),
                            dt.as_gregorian_tai_str()
                        );
                        msr_cnt += 1;
                        continue;
                    } else if next_epoch > dt {
                        // No measurement can be used here, let's just do a time update (unless we have already done a time update)
                        if num_msr_processed == 0 {
                            debug!("time update {}", dt.as_gregorian_tai_str());
                            match self.kf.time_update(nominal_state) {
                                Ok(est) => {
                                    if self.kf.is_extended() {
                                        self.prop.dynamics.set_estimated_state(
                                            self.prop.dynamics.estimated_state()
                                                + est.state_deviation(),
                                        );
                                    }
                                    self.estimates.push(est);
                                }
                                Err(e) => return Some(e),
                            }
                        }
                        break; // Move on to the next propagator output
                    } else {
                        // The epochs match, so this is a valid measurement to use
                        // Get the computed observations
                        for device in self.devices.iter() {
                            if let Some(computed_meas) = device.measure(&meas_input) {
                                if computed_meas.visible() {
                                    self.kf.update_h_tilde(computed_meas.sensitivity());

                                    // Switch back from extended if necessary
                                    if self.kf.is_extended()
                                        && self.ekf_trigger.disable_ekf(real_meas.epoch())
                                    {
                                        self.kf.set_extended(false);
                                        info!(
                                            "EKF disabled @ {}",
                                            real_meas.epoch().as_gregorian_tai_str()
                                        );
                                    }

                                    match self.kf.measurement_update(
                                        nominal_state,
                                        real_meas.observation(),
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
                                                        real_meas.epoch().as_gregorian_tai_str()
                                                    );
                                                } else {
                                                    info!(
                                                        "EKF enabled @ {}",
                                                        real_meas.epoch().as_gregorian_tai_str()
                                                    );
                                                }
                                            }
                                            if self.kf.is_extended() {
                                                self.prop.dynamics.set_estimated_state(
                                                    self.prop
                                                        .dynamics
                                                        .extract_estimated_state(&nominal_state)
                                                        + est.state_deviation(),
                                                );
                                            }
                                            self.estimates.push(est);
                                            self.residuals.push(res);
                                        }
                                        Err(e) => return Some(e),
                                    }
                                    if !self.simultaneous_msr {
                                        break;
                                    }
                                }
                            }
                        }
                        // And increment the measurement counter
                        msr_cnt += 1;
                        num_msr_processed += 1;
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
                } else {
                    // No more measurements, we can only do a time update
                    debug!("final time update {:?}", dt.as_gregorian_tai_str());
                    match self.kf.time_update(nominal_state) {
                        Ok(est) => {
                            if self.kf.is_extended() {
                                self.prop.dynamics.set_estimated_state(
                                    self.prop.dynamics.extract_estimated_state(&nominal_state)
                                        + est.state_deviation(),
                                );
                            }
                            self.estimates.push(est);
                        }
                        Err(e) => return Some(e),
                    }
                    // Leave the loop {...} which processes several measurements at once
                    break;
                }
            }
        }

        // Always report the 100% mark
        if !reported[10] {
            info!("{:>3}% done ({:.0} measurements processed)", 100, num_msrs);
        }

        None
    }

    /// Allows for covariance mapping without processing measurements
    pub fn map_covar(&mut self, end_epoch: Epoch) -> Option<FilterError> {
        let (tx, rx) = channel();
        self.prop.tx_chan = Some(tx);
        // Start by propagating the estimator (on the same thread).
        let prop_time = end_epoch - self.kf.previous_estimate().epoch();
        info!("Propagating for {} seconds", prop_time);

        self.prop.until_time_elapsed(prop_time);
        info!("Mapping covariance");

        while let Ok(nominal_state) = rx.try_recv() {
            // Update the STM of the KF (needed between each measurement or time update)
            let stm = self.prop.dynamics.extract_stm(&nominal_state);
            self.kf.update_stm(stm);
            info!("final time update {:?}", nominal_state.epoch());
            match self.kf.time_update(nominal_state) {
                Ok(est) => {
                    if self.kf.is_extended() {
                        let est_state = est.state_deviation().clone();
                        self.prop.dynamics.set_estimated_state(
                            self.prop.dynamics.extract_estimated_state(&nominal_state) + est_state,
                        );
                    }
                    self.estimates.push(est);
                }
                Err(e) => return Some(e),
            }
        }

        None
    }
}

impl<
        'a,
        D: Estimable<MsrIn, LinStateSize = Msr::StateSize>,
        E: ErrorCtrl,
        Msr: Measurement,
        N: MeasurementDevice<MsrIn, Msr>,
        A: DimName,
        K: Filter<D::LinStateSize, A, Msr::MeasurementSize, D::StateType>,
        MsrIn,
    > ODProcess<'a, D, E, Msr, N, CkfTrigger, A, K, MsrIn>
where
    D::StateType: EstimableState<Msr::StateSize>,
    DefaultAllocator: Allocator<f64, D::StateSize>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::StateSize>
        + Allocator<f64, Msr::StateSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, D::LinStateSize>
        + Allocator<f64, D::LinStateSize, Msr::MeasurementSize>
        + Allocator<f64, D::LinStateSize, D::LinStateSize>
        + Allocator<f64, A, A>
        + Allocator<f64, D::LinStateSize, A>
        + Allocator<f64, A, D::LinStateSize>,
{
    pub fn ckf(
        prop: Propagator<'a, D, E>,
        kf: K,
        devices: Vec<N>,
        simultaneous_msr: bool,
        num_expected_msr: usize,
    ) -> Self {
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
            _marker: PhantomData::<A>,
        }
    }

    pub fn default_ckf(prop: Propagator<'a, D, E>, kf: K, devices: Vec<N>) -> Self {
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
            _marker: PhantomData::<A>,
        }
    }
}
/// A trait detailing when to switch to from a CKF to an EKF
pub trait EkfTrigger {
    fn enable_ekf<S, E, T: EstimableState<S>>(&mut self, est: &E) -> bool
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
    fn enable_ekf<S, E, T: EstimableState<S>>(&mut self, _est: &E) -> bool
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
    pub disable_time: f64,
    prev_msr_dt: Option<Epoch>,
    cur_msrs: usize,
}

impl StdEkfTrigger {
    pub fn new(num_msrs: usize, disable_time: f64) -> Self {
        Self {
            num_msrs,
            disable_time,
            prev_msr_dt: None,
            cur_msrs: 0,
        }
    }
}

impl EkfTrigger for StdEkfTrigger {
    fn enable_ekf<S, E, T: EstimableState<S>>(&mut self, est: &E) -> bool
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
