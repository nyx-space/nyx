/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName};
use crate::md::trajectory::{Interpolatable, Traj};

pub use crate::od::estimate::*;
pub use crate::od::kalman::*;
pub use crate::od::measurement::*;
pub use crate::od::residual::*;
pub use crate::od::snc::*;
pub use crate::od::*;

use crate::propagators::error_ctrl::ErrorCtrl;
use crate::propagators::PropInstance;
pub use crate::time::{Duration, Unit};
use crate::State;
mod conf;
pub use conf::{IterationConf, SmoothingArc};
mod trigger;
pub use trigger::{CkfTrigger, EkfTrigger, KfTrigger};

use std::marker::PhantomData;
use std::ops::Add;

use self::msr::arc::TrackingArc;

/// An orbit determination process. Note that everything passed to this structure is moved.
#[allow(clippy::upper_case_acronyms)]
pub struct ODProcess<
    'a,
    D: Dynamics,
    E: ErrorCtrl,
    Msr: Measurement,
    T: KfTrigger,
    A: DimName,
    S: EstimateFrom<D::StateType, Msr> + Interpolatable,
    K: Filter<S, A, Msr::MeasurementSize>,
> where
    D::StateType: Interpolatable + Add<OVector<f64, <S as State>::Size>, Output = D::StateType>,
    <DefaultAllocator as Allocator<f64, <D::StateType as State>::VecLength>>::Buffer: Send,
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <S as State>::Size>
        + Allocator<f64, <S as State>::VecLength>
        + Allocator<f64, <D::StateType as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, S::Size>
        + Allocator<f64, S::Size>
        + Allocator<f64, Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, Msr::MeasurementSize>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<usize, <D::StateType as State>::Size, <D::StateType as State>::Size>
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
    /// Vector of estimates available after a pass
    pub estimates: Vec<K::Estimate>,
    /// Vector of residuals available after a pass
    pub residuals: Vec<Residual<Msr::MeasurementSize>>,
    pub ekf_trigger: T,
    pub cosm: Arc<Cosm>,
    init_state: D::StateType,
    _marker: PhantomData<A>,
}

impl<
        'a,
        D: Dynamics,
        E: ErrorCtrl,
        Msr: Measurement,
        T: KfTrigger,
        A: DimName,
        S: EstimateFrom<D::StateType, Msr> + Interpolatable,
        K: Filter<S, A, Msr::MeasurementSize>,
    > ODProcess<'a, D, E, Msr, T, A, S, K>
where
    D::StateType: Interpolatable + Add<OVector<f64, <S as State>::Size>, Output = D::StateType>,
    <DefaultAllocator as Allocator<f64, <D::StateType as State>::VecLength>>::Buffer: Send,
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, S::Size>
        + Allocator<f64, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, <D::StateType as State>::Size>
        + Allocator<f64, Msr::MeasurementSize, <S as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, Msr::MeasurementSize>
        + Allocator<f64, <S as State>::Size, Msr::MeasurementSize>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<usize, <D::StateType as State>::Size, <D::StateType as State>::Size>
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
    pub fn ekf(prop: PropInstance<'a, D, E>, kf: K, trigger: T, cosm: Arc<Cosm>) -> Self {
        let init_state = prop.state;
        let mut estimates = Vec::with_capacity(10_001);
        estimates.push(kf.previous_estimate().clone());
        Self {
            prop,
            kf,
            estimates,
            residuals: Vec::with_capacity(10_000),
            ekf_trigger: trigger,
            cosm,
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
            let sm_est_kp1 = &self.estimates[l - smoothed.len() + 1];
            let x_kp1_l = sm_est_kp1.state_deviation();
            let p_kp1_l = sm_est_kp1.covar();
            // Borrow the k-th estimate, which we're smoothing with the next estimate
            let est_k = &self.estimates[l - smoothed.len()];
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
                    if est_kp1.epoch() - est_k.epoch() > gap_s {
                        break;
                    }
                }
                SmoothingArc::Prediction => {
                    if est_kp1.predicted() {
                        break;
                    }
                }
                SmoothingArc::All => {
                    // if est_k.epoch() == est_kp1.epoch() {
                    //     break;
                    // }
                }
            }

            // Compute the STM between both steps taken by the filter
            // The filter will reset the STM between each estimate it computes, time update or measurement update.
            // Therefore, the STM is simply the inverse of the one we used previously.
            // est_kp1 is the estimate that used the STM from time k to time k+1. So the STM stored there
            // is \Phi_{k \to k+1}. Let's invert that.
            let phi_kp1_k = &est_kp1
                .stm()
                .clone()
                .try_inverse()
                .ok_or(NyxError::SingularStateTransitionMatrix)?;

            // Compute smoothed state deviation
            let x_k_l = phi_kp1_k * x_kp1_l;
            // Compute smoothed covariance
            let p_k_l = phi_kp1_k * p_kp1_l * phi_kp1_k.transpose();
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
            let mut msr_noise_item_inv: OVector<f64, Msr::MeasurementSize> =
                self.kf.measurement_noise(residual.dt).diagonal().clone();
            for i in 0..msr_noise_item_inv.len() {
                msr_noise_item_inv[i] = 1.0 / msr_noise_item_inv[i];
            }
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
    pub fn iterate_arc<Dev>(
        &mut self,
        arc: &TrackingArc<Msr>,
        config: IterationConf,
    ) -> Result<(), NyxError>
    where
        Dev: TrackingDeviceSim<S, Msr>,
    {
        let measurements = &arc.measurements;

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
            self.process_arc::<Dev>(arc)?;

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

    /// Process the provided tracking arc for this orbit determination process.
    #[allow(clippy::erasing_op)]
    pub fn process_arc<Dev>(&mut self, arc: &TrackingArc<Msr>) -> Result<(), NyxError>
    where
        Dev: TrackingDeviceSim<S, Msr>,
    {
        let mut devices = arc.rebuild_devices::<S, Dev>(self.cosm.clone()).unwrap();

        let measurements = &arc.measurements;
        assert!(
            measurements.len() >= 2,
            "must have at least two measurements"
        );
        // Start by propagating the estimator (on the same thread).
        let num_msrs = measurements.len();
        let step_size = arc.min_duration_sep().unwrap();
        // Update the step size of the navigation propagator if it isn't already fixed step
        if !self.prop.fixed_step {
            self.prop.set_step(step_size, false);
        }
        let prop_time = measurements[num_msrs - 1].1.epoch() - self.kf.previous_estimate().epoch();
        info!("Navigation propagating for a total of {prop_time} with step size {step_size}");

        let mut epoch = self.prop.state.epoch();

        let mut reported = [false; 11];
        info!("Processing {num_msrs} measurements with covariance mapping");

        // In the following, we'll continuously tell the propagator to advance by a single step until the next measurement.
        // This is effectively a clone of the "for_duration" function of the propagator.
        // However, after every step of the propagator, we will perform a time update and reset the STM.
        // This ensures that we have a time update every time the propagator error function says that the integration would
        // be wrong if the step would be larger. Moreover, we will reset the STM to ensure that the state transition matrix
        // stays linear (the error control function may verify this).

        for (msr_cnt, (device_name, msr)) in measurements.iter().enumerate() {
            let next_msr_epoch = msr.epoch();

            // Advance the propagator
            loop {
                let delta_t = next_msr_epoch - epoch;

                // Propagator for the minimum time between the step size and the duration to the next measurement.
                // Ensure that we don't go backward if the previous step we took was indeed backward.
                let next_step_size = delta_t.min(if self.prop.details.step.is_negative() {
                    step_size
                } else {
                    self.prop.details.step
                });
                debug!("advancing propagator by {next_step_size} (Δt to next msr: {delta_t})");
                let (_, traj_covar) = self.prop.for_duration_with_traj(next_step_size)?;

                // Now that we've advanced the propagator, let's see whether we're at the time of the next measurement.

                // Extract the state and update the STM in the filter.
                let nominal_state = S::extract(self.prop.state);
                // Get the datetime and info needed to compute the theoretical measurement according to the model
                epoch = nominal_state.epoch();

                if nominal_state.epoch() == next_msr_epoch {
                    // Perform a measurement update

                    // First, let's build a trajectory of the states without the covariances.
                    let mut traj = Traj::new();
                    for state in traj_covar.states {
                        traj.states.push(S::extract(state));
                    }

                    // Get the computed observations
                    match devices.get_mut(device_name) {
                        Some(device) => {
                            if let Some((computed_meas, device_loc)) =
                                device.measure_as_seen(epoch, &traj, None, self.cosm.clone())
                            {
                                // Switch back from extended if necessary
                                if self.kf.is_extended() && self.ekf_trigger.disable_ekf(epoch) {
                                    self.kf.set_extended(false);
                                    info!("EKF disabled @ {epoch}");
                                }

                                let h_tilde = S::sensitivity(msr, nominal_state, device_loc);

                                self.kf.update_h_tilde(h_tilde);

                                match self.kf.measurement_update(
                                    nominal_state,
                                    &msr.observation(),
                                    &computed_meas.observation(),
                                ) {
                                    Ok((estimate, residual)) => {
                                        debug!("msr update #{msr_cnt} @ {epoch}");
                                        // Switch to EKF if necessary, and update the dynamics and such
                                        // Note: we call enable_ekf first to ensure that the trigger gets
                                        // called in case it needs to save some information (e.g. the
                                        // StdEkfTrigger needs to store the time of the previous measurement).
                                        if self.ekf_trigger.enable_ekf(&estimate)
                                            && !self.kf.is_extended()
                                        {
                                            self.kf.set_extended(true);
                                            if !estimate.within_3sigma() {
                                                warn!("EKF enabled @ {epoch} but filter DIVERGING");
                                            } else {
                                                info!("EKF enabled @ {epoch}");
                                            }
                                        }
                                        if self.kf.is_extended() {
                                            self.prop.state =
                                                self.prop.state + estimate.state_deviation();
                                        }
                                        self.prop.state.reset_stm();
                                        self.estimates.push(estimate);
                                        self.residuals.push(residual);
                                    }
                                    Err(e) => return Err(e),
                                }
                            } else {
                                warn!("Real observation exists @ {epoch} but simulated {device_name} does not see it -- ignoring measurement");
                            }
                        }
                        None => {
                            error!("Tracking arc references {device_name} which is not in the list of configured devices")
                        }
                    }

                    let msr_prct = (10.0 * (msr_cnt as f64) / (num_msrs as f64)) as usize;
                    if !reported[msr_prct] {
                        info!(
                            "{:>3}% done ({msr_cnt:.0} measurements processed)",
                            10 * msr_prct,
                        );
                        reported[msr_prct] = true;
                    }

                    break;
                } else {
                    // No measurement can be used here, let's just do a time update and continue advancing the propagator.
                    debug!("time update {epoch}");
                    match self.kf.time_update(nominal_state) {
                        Ok(est) => {
                            // State deviation is always zero for an EKF time update
                            // therefore we don't do anything different for an extended filter
                            self.estimates.push(est);
                        }
                        Err(e) => return Err(e),
                    }
                    self.prop.state.reset_stm();
                }
            }
        }

        // Always report the 100% mark
        if !reported[10] {
            info!("100% done ({num_msrs:.0} measurements processed)");
        }

        Ok(())
    }

    /// Allows for covariance mapping without processing measurements
    pub fn map_covar(&mut self, end_epoch: Epoch) -> Result<(), NyxError> {
        let prop_time = end_epoch - self.kf.previous_estimate().epoch();
        info!("Propagating for {prop_time} seconds and mapping covariance",);

        loop {
            let mut epoch = self.prop.state.epoch();
            if epoch + self.prop.details.step > end_epoch {
                self.prop.until_epoch(end_epoch)?;
            } else {
                self.prop.single_step()?;
            }
            // Perform time update

            // Extract the state and update the STM in the filter.
            let prop_state = self.prop.state;
            let nominal_state = S::extract(prop_state);
            // Get the datetime and info needed to compute the theoretical measurement according to the model
            epoch = nominal_state.epoch();
            // No measurement can be used here, let's just do a time update
            debug!("time update {epoch}");
            match self.kf.time_update(nominal_state) {
                Ok(est) => {
                    // State deviation is always zero for an EKF time update
                    // therefore we don't do anything different for an extended filter
                    self.estimates.push(est);
                }
                Err(e) => return Err(e),
            }
            self.prop.state.reset_stm();
            if epoch == end_epoch {
                break;
            }
        }

        Ok(())
    }

    /// Builds the navigation trajectory for the estimated state only (no covariance until https://gitlab.com/nyx-space/nyx/-/issues/199!)
    pub fn to_traj(&self) -> Result<Traj<S>, NyxError>
    where
        DefaultAllocator: Allocator<f64, <S as State>::VecLength>,
        S: Interpolatable,
    {
        if self.estimates.is_empty() {
            Err(NyxError::NoStateData(
                "No navigation trajectory to generate: run the OD process first".to_string(),
            ))
        } else {
            Ok(Traj {
                states: self
                    .estimates
                    .iter()
                    .map(|est| est.nominal_state())
                    .collect(),
                name: None,
            })
        }
    }
}

impl<
        'a,
        D: Dynamics,
        E: ErrorCtrl,
        Msr: Measurement,
        A: DimName,
        S: EstimateFrom<D::StateType, Msr> + Interpolatable,
        K: Filter<S, A, Msr::MeasurementSize>,
    > ODProcess<'a, D, E, Msr, CkfTrigger, A, S, K>
where
    D::StateType: Interpolatable + Add<OVector<f64, <S as State>::Size>, Output = D::StateType>,
    <DefaultAllocator as Allocator<f64, <D::StateType as State>::VecLength>>::Buffer: Send,
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, S::Size>
        + Allocator<f64, S::Size>
        + Allocator<f64, Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, Msr::MeasurementSize>
        + Allocator<f64, <S as State>::Size, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, <S as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<usize, <D::StateType as State>::Size, <D::StateType as State>::Size>
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
    pub fn ckf(prop: PropInstance<'a, D, E>, kf: K, cosm: Arc<Cosm>) -> Self {
        let init_state = prop.state;
        let mut estimates = Vec::with_capacity(10_001);
        estimates.push(kf.previous_estimate().clone());
        Self {
            prop,
            kf,
            estimates,
            residuals: Vec::with_capacity(10_000),
            cosm,
            ekf_trigger: CkfTrigger {},
            init_state,
            _marker: PhantomData::<A>,
        }
    }
}
