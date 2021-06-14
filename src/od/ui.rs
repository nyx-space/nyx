/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, DimName};
use crate::md::trajectory::Traj;

pub use super::estimate::*;
pub use super::kalman::*;
pub use super::ranging::*;
pub use super::residual::*;
pub use super::snc::*;
pub use super::*;

use crate::propagators::error_ctrl::ErrorCtrl;
use crate::propagators::PropInstance;
pub use crate::time::{Duration, TimeUnit};
use crate::State;

use std::convert::TryFrom;
use std::default::Default;
use std::fmt;
use std::marker::PhantomData;
use std::ops::Add;
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
            SmoothingArc::All => write!(f, "all estimates"),
            SmoothingArc::Epoch(e) => write!(f, "{}", e),
            SmoothingArc::TimeGap(g) => write!(f, "time gap of {}", g),
            SmoothingArc::Prediction => write!(f, "first prediction"),
        }
    }
}

/// Defines a filter iteration configuration. Allows iterating on an OD solution until convergence criteria is met.
/// The root mean squared of the postfit residuals is used to assess convergence between iterations.
#[derive(Clone, Copy, Debug)]
pub struct IterationConf {
    /// The number of measurements to account for in the iteration
    pub smoother: SmoothingArc,
    /// The absolute tolerance of the RMS postfit residual
    pub absolute_tol: f64,
    /// The relative tolerance between the latest RMS postfit residual and the best RMS postfit residual so far
    pub relative_tol: f64,
    /// The maximum number of iterations to allow (will raise an error if the filter has not converged after this many iterations)
    pub max_iterations: usize,
    /// The maximum number of subsequent divergences in RMS.
    pub max_divergences: usize,
    /// Set to true to force an ODP failure when the convergence criteria is not met
    pub force_failure: bool,
    /// Set to true to use the RMS prefit instead of postfit
    pub use_prefit: bool,
}

impl Default for IterationConf {
    /// The default absolute tolerance is 1e-2 (calibrated on an EKF with error).
    fn default() -> Self {
        Self {
            smoother: SmoothingArc::All,
            absolute_tol: 1e-2,
            relative_tol: 1e-3,
            max_iterations: 15,
            max_divergences: 3,
            force_failure: false,
            use_prefit: false,
        }
    }
}

impl fmt::Display for IterationConf {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let kind = if self.use_prefit { "prefit" } else { "postfit" };
        write!(f, "Iterate {} residuals until abs = {:.2e}, or rel = {:.2e}, or {} iterations, or {} subsequent divergences with smoothing condition of {}",
            kind,
            self.absolute_tol,
            self.relative_tol,
            self.max_iterations,
            self.max_divergences,
            self.smoother)
    }
}

impl TryFrom<SmoothingArc> for IterationConf {
    type Error = NyxError;

    /// Converts a smoother into an interation configuration to iterate just once without failing
    fn try_from(smoother: SmoothingArc) -> Result<Self, Self::Error> {
        Ok(Self {
            smoother,
            force_failure: false,
            ..Default::default()
        })
    }
}

/// An orbit determination process. Note that everything passed to this structure is moved.
#[allow(clippy::upper_case_acronyms)]
pub struct ODProcess<
    'a,
    D: Dynamics,
    E: ErrorCtrl,
    Msr: Measurement<StateSize = <S as State>::Size>,
    N: MeasurementDevice<S, Msr>,
    T: EkfTrigger,
    A: DimName,
    S: EstimateFrom<D::StateType>,
    K: Filter<S, A, Msr::MeasurementSize>,
> where
    D::StateType: Add<OVector<f64, <S as State>::Size>, Output = D::StateType>,
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <S as State>::Size>
        + Allocator<f64, <S as State>::VecLength>
        + Allocator<f64, <D::StateType as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::StateSize>
        + Allocator<f64, Msr::StateSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, Msr::MeasurementSize>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<f64, <S as State>::Size, <S as State>::Size>
        + Allocator<f64, A>
        + Allocator<f64, A, A>
        + Allocator<f64, <D::StateType as State>::Size, A>
        + Allocator<f64, A, <D::StateType as State>::Size>
        + Allocator<f64, <S as State>::Size, A>
        + Allocator<f64, A, <S as State>::Size>,
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
        Msr: Measurement<StateSize = <S as State>::Size>,
        N: MeasurementDevice<S, Msr>,
        T: EkfTrigger,
        A: DimName,
        S: EstimateFrom<D::StateType>,
        K: Filter<S, A, Msr::MeasurementSize>,
    > ODProcess<'a, D, E, Msr, N, T, A, S, K>
where
    D::StateType: Add<OVector<f64, <S as State>::Size>, Output = D::StateType>,
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::StateSize>
        + Allocator<f64, Msr::StateSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, <D::StateType as State>::Size>
        + Allocator<f64, Msr::MeasurementSize, <S as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, Msr::MeasurementSize>
        + Allocator<f64, <S as State>::Size, Msr::MeasurementSize>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>
        + Allocator<f64, A>
        + Allocator<f64, A, A>
        + Allocator<f64, <D::StateType as State>::Size, A>
        + Allocator<f64, A, <D::StateType as State>::Size>
        + Allocator<f64, <S as State>::Size>
        + Allocator<f64, <S as State>::VecLength>
        + Allocator<f64, <S as State>::Size, <S as State>::Size>
        + Allocator<f64, <S as State>::Size, A>
        + Allocator<f64, A, <S as State>::Size>,
{
    pub fn ekf(
        prop: PropInstance<'a, D, E>,
        kf: K,
        devices: Vec<N>,
        simultaneous_msr: bool,
        num_expected_msr: usize,
        trigger: T,
    ) -> Self {
        let init_state = prop.state;
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
        let init_state = prop.state;
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
    pub fn smooth(&self, condition: SmoothingArc) -> Result<Vec<K::Estimate>, NyxError> {
        let l = self.estimates.len() - 1;

        info!("Smoothing {} estimates until {}", l + 1, condition);
        let mut smoothed = Vec::with_capacity(self.estimates.len());
        // Set the first item of the smoothed estimates to the last estimate (we cannot smooth the very last estimate)
        smoothed.push(self.estimates.last().unwrap().clone());

        loop {
            // Borrow the previously smoothed estimate of the k+1 estimate
            let sm_est_kp1 = &self.estimates[l - smoothed.len() + 1].clone();
            let x_kp1_l = sm_est_kp1.state_deviation();
            let p_kp1_l = sm_est_kp1.covar();
            // Borrow the k-th estimate, which we're smoothing with the next estimate
            let est_k = &self.estimates[l - smoothed.len()];
            let x_k_k = &est_k.state_deviation();
            let p_k_k = &est_k.covar();
            // Borrow the k+1-th estimate, which we're smoothing with the next estimate
            let est_kp1 = &self.estimates[l - smoothed.len() + 1];

            // Check the smoother stopping condition
            match condition {
                SmoothingArc::Epoch(e) => {
                    // If the epoch of the next estimate is _before_ the stopping time, stop smoothing
                    if est_kp1.epoch() < e {
                        break;
                    }
                }
                SmoothingArc::TimeGap(gap_s) => {
                    if est_k.epoch() - est_kp1.epoch() > gap_s {
                        break;
                    }
                }
                SmoothingArc::Prediction => {
                    if est_kp1.predicted() {
                        break;
                    }
                }
                SmoothingArc::All => {}
            }

            let phi_kp1_k = est_kp1.stm();
            // let p_kp1_k = phi_kp1_k * p_k_k * phi_kp1_k.transpose(); // TODO: Add SNC here, which is effectively covar_bar!
            let p_kp1_k = est_kp1.predicted_covar();
            let p_kp1_k_inv = &p_kp1_k
                .try_inverse()
                .ok_or(NyxError::SingularCovarianceMatrix)?;
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
            if smoothed.len() == self.estimates.len() {
                break;
            }
        }

        info!(
            "Smoothing condition reached after {} estimates ",
            smoothed.len()
        );

        // Now, let's add all of the other estimates so that the same indexing can be done
        // between all the estimates and the smoothed estimates
        if smoothed.len() < self.estimates.len() {
            // Add the estimates that might have been skipped.
            let mut k = self.estimates.len() - smoothed.len();
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

    /// Returns the root mean square of the prefit residuals
    pub fn rms_prefit_residual(&self) -> f64 {
        let mut sum = 0.0;
        for residual in &self.residuals {
            // sum += residual.prefit.dot(&residual.prefit);
            let mut msr_noise_item_inv = self.kf.measurement_noise(residual.dt).diagonal().clone();
            msr_noise_item_inv.apply(|m| 1.0 / m);
            sum += residual.prefit.dot(&msr_noise_item_inv).powi(2);
        }
        (sum / (self.estimates.len() as f64)).sqrt()
    }

    /// Returns the root mean square of the postfit residuals
    pub fn rms_postfit_residual(&self) -> f64 {
        let mut sum = 0.0;
        for residual in &self.residuals {
            sum += residual.postfit.dot(&residual.postfit);
        }
        (sum / (self.estimates.len() as f64)).sqrt()
    }

    /// Allows iterating on the filter solution. Requires specifying a smoothing condition to know where to stop the smoothing.
    pub fn iterate(&mut self, measurements: &[Msr], config: IterationConf) -> Result<(), NyxError> {
        // Compute the initial RMS
        let mut best_rms = if config.use_prefit {
            self.rms_prefit_residual()
        } else {
            self.rms_postfit_residual()
        };
        let mut previous_rms = best_rms;
        let mut divergence_cnt = 0;
        let mut iter_cnt = 0;
        loop {
            if best_rms <= config.absolute_tol {
                info!("*****************");
                info!("*** CONVERGED ***");
                info!("*****************");

                info!(
                    "Filter converged to absolute tolerance ({:.2e} < {:.2e}) after {} iterations",
                    best_rms, config.absolute_tol, iter_cnt
                );
                return Ok(());
            }

            iter_cnt += 1;

            info!("***************************");
            info!("*** Iteration number {} ***", iter_cnt);
            info!("***************************");

            // First, smooth the estimates
            let smoothed = self.smooth(config.smoother)?;
            // Reset the propagator
            self.prop.state = self.init_state;
            // Empty the estimates and add the first smoothed estimate as the initial estimate
            self.estimates = Vec::with_capacity(measurements.len());
            self.estimates.push(smoothed[0].clone());
            self.kf.set_previous_estimate(&smoothed[0]);
            // And re-run the filter
            self.process_measurements(measurements)?;

            // Compute the new RMS
            let new_rms = if config.use_prefit {
                self.rms_prefit_residual()
            } else {
                self.rms_postfit_residual()
            };
            let cur_rel_rms = (new_rms - best_rms).abs() / best_rms;
            if cur_rel_rms < config.relative_tol {
                info!("*****************");
                info!("*** CONVERGED ***");
                info!("*****************");
                info!(
                    "New RMS: {:.5}\tPrevious RMS: {:.5}\tBest RMS: {:.5}",
                    new_rms, previous_rms, best_rms
                );
                info!(
                    "Filter converged to relative tolerance ({:.2e} < {:.2e}) after {} iterations",
                    cur_rel_rms, config.relative_tol, iter_cnt
                );
                return Ok(());
            }

            if new_rms > previous_rms {
                warn!(
                    "New RMS: {:.5}\tPrevious RMS: {:.5}\tBest RMS: {:.5}",
                    new_rms, previous_rms, best_rms
                );
                divergence_cnt += 1;
                previous_rms = new_rms;
                if divergence_cnt >= config.max_divergences {
                    let msg = format!(
                        "Filter iterations have continuously diverged {} times: {}",
                        config.max_divergences, config
                    );
                    if config.force_failure {
                        return Err(NyxError::MaxIterReached(msg));
                    } else {
                        error!("{}", msg);
                        return Ok(());
                    }
                } else {
                    warn!("Filter iteration caused divergence {} of {} acceptable subsequent divergences", divergence_cnt, config.max_divergences);
                }
            } else {
                info!(
                    "New RMS: {:.5}\tPrevious RMS: {:.5}\tBest RMS: {:.5}",
                    new_rms, previous_rms, best_rms
                );
                // Reset the counter
                divergence_cnt = 0;
                previous_rms = new_rms;
                if previous_rms < best_rms {
                    best_rms = previous_rms;
                }
            }

            if iter_cnt >= config.max_iterations {
                let msg = format!(
                    "Filter has iterated {} times but failed to reach filter convergence criteria: {}",
                    config.max_iterations, config
                );
                if config.force_failure {
                    return Err(NyxError::MaxIterReached(msg));
                } else {
                    error!("{}", msg);
                    return Ok(());
                }
            }
        }
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
        info!("Navigation propagating for a total of {}", prop_time);

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
            self.prop.for_duration(delta_t)?;

            while let Ok(prop_state) = rx.try_recv() {
                let nominal_state = S::extract(&prop_state);

                // Get the datetime and info needed to compute the theoretical measurement according to the model
                let dt = nominal_state.epoch();

                // Update the STM of the KF (needed between each measurement or time update)
                let stm = nominal_state.stm()?;
                self.kf.update_stm(stm);

                // Check if we should do a time update or a measurement update
                if next_msr_epoch > dt {
                    if msr_cnt == 0 && !arc_warned {
                        warn!("OD arc starts prior to first measurement");
                        arc_warned = true;
                    }
                    // No measurement can be used here, let's just do a time update
                    trace!("time update {}", dt);
                    match self.kf.time_update(nominal_state) {
                        Ok(est) => {
                            // State deviation is always zero for an EKF time update
                            // therefore we don't do anything different for an extended filter
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
                                    info!("EKF disabled @ {}", dt);
                                }

                                match self.kf.measurement_update(
                                    nominal_state,
                                    &msr.observation(),
                                    &computed_meas.observation(),
                                ) {
                                    Ok((est, res)) => {
                                        trace!("msr update #{} @ {}", msr_cnt, dt);
                                        // Switch to EKF if necessary, and update the dynamics and such
                                        // Note: we call enable_ekf first to ensure that the trigger gets
                                        // called in case it needs to save some information (e.g. the
                                        // StdEkfTrigger needs to store the time of the previous measurement).
                                        if self.ekf_trigger.enable_ekf(&est)
                                            && !self.kf.is_extended()
                                        {
                                            self.kf.set_extended(true);
                                            if !est.within_3sigma() {
                                                warn!("EKF enabled @ {} but filter DIVERGING", dt);
                                            } else {
                                                info!("EKF enabled @ {}", dt);
                                            }
                                        }
                                        if self.kf.is_extended() {
                                            self.prop.state =
                                                self.prop.state + est.state_deviation();
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

        self.prop.for_duration(prop_time)?;

        info!("Mapping covariance");

        while let Ok(prop_state) = rx.try_recv() {
            let nominal_state = S::extract(&prop_state);
            // Update the STM of the KF (needed between each measurement or time update)
            self.kf.update_stm(nominal_state.stm()?);
            info!("final time update {}", nominal_state.epoch());
            match self.kf.time_update(nominal_state) {
                Ok(est) => {
                    if self.kf.is_extended() {
                        self.prop.state = self.prop.state + est.state_deviation();
                    }
                    self.estimates.push(est);
                }
                Err(e) => return Err(e),
            }
        }

        Ok(())
    }

    /// Builds the navigation trajectory for the estimated state only (no covariance until https://gitlab.com/nyx-space/nyx/-/issues/199!)
    pub fn to_nav_traj(&self) -> Result<Traj<S>, NyxError>
    where
        DefaultAllocator: Allocator<f64, <S as State>::VecLength>,
    {
        if self.estimates.is_empty() {
            Err(NyxError::NoStateData(
                "No navigation trajectory to generate: run the OD process first".to_string(),
            ))
        } else {
            let (tx, rx) = channel();
            let start_state = self.estimates[0].state();
            for estimate in &self.estimates {
                tx.send(estimate.state()).unwrap();
            }
            Traj::new(start_state, rx)
        }
    }
}

impl<
        'a,
        D: Dynamics,
        E: ErrorCtrl,
        Msr: Measurement<StateSize = <S as State>::Size>,
        N: MeasurementDevice<S, Msr>,
        A: DimName,
        S: EstimateFrom<D::StateType>,
        K: Filter<S, A, Msr::MeasurementSize>,
    > ODProcess<'a, D, E, Msr, N, CkfTrigger, A, S, K>
where
    D::StateType: Add<OVector<f64, <S as State>::Size>, Output = D::StateType>,
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::StateSize>
        + Allocator<f64, Msr::StateSize>
        + Allocator<f64, Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, Msr::MeasurementSize>
        + Allocator<f64, <S as State>::Size, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, <S as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<f64, <S as State>::Size>
        + Allocator<f64, <S as State>::VecLength>
        + Allocator<f64, <S as State>::Size, <S as State>::Size>
        + Allocator<f64, A>
        + Allocator<f64, A, A>
        + Allocator<f64, <D::StateType as State>::Size, A>
        + Allocator<f64, A, <D::StateType as State>::Size>
        + Allocator<f64, <S as State>::Size, A>
        + Allocator<f64, A, <S as State>::Size>,
{
    pub fn ckf(
        prop: PropInstance<'a, D, E>,
        kf: K,
        devices: Vec<N>,
        simultaneous_msr: bool,
        num_expected_msr: usize,
    ) -> Self {
        let init_state = prop.state;
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
        let init_state = prop.state;
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
    fn enable_ekf<T: State, E>(&mut self, est: &E) -> bool
    where
        E: Estimate<T>,
        DefaultAllocator: Allocator<f64, <T as State>::Size>
            + Allocator<f64, <T as State>::VecLength>
            + Allocator<f64, <T as State>::Size, <T as State>::Size>;

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
    fn enable_ekf<T: State, E>(&mut self, _est: &E) -> bool
    where
        E: Estimate<T>,
        DefaultAllocator: Allocator<f64, <T as State>::Size>
            + Allocator<f64, <T as State>::VecLength>
            + Allocator<f64, <T as State>::Size, <T as State>::Size>,
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
    fn enable_ekf<T: State, E>(&mut self, est: &E) -> bool
    where
        E: Estimate<T>,
        DefaultAllocator: Allocator<f64, <T as State>::Size>
            + Allocator<f64, <T as State>::VecLength>
            + Allocator<f64, <T as State>::Size, <T as State>::Size>,
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
