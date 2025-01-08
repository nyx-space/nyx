/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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
pub use crate::od::ground_station::*;
pub use crate::od::snc::*;
pub use crate::od::*;
use crate::propagators::PropInstance;
pub use crate::time::{Duration, Unit};
use anise::prelude::Almanac;
use indexmap::IndexSet;
use msr::sensitivity::TrackerSensitivity;
use snafu::prelude::*;
mod conf;
pub use conf::{IterationConf, SmoothingArc};
mod trigger;
pub use trigger::EkfTrigger;
mod rejectcrit;
use self::msr::TrackingDataArc;
pub use self::rejectcrit::ResidRejectCrit;
use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::ops::Add;
mod export;

/// An orbit determination process. Note that everything passed to this structure is moved.
#[allow(clippy::upper_case_acronyms)]
pub struct ODProcess<
    'a,
    D: Dynamics,
    MsrSize: DimName,
    Accel: DimName,
    K: Filter<D::StateType, Accel, MsrSize>,
    Trk: TrackerSensitivity<D::StateType, D::StateType>,
> where
    D::StateType:
        Interpolatable + Add<OVector<f64, <D::StateType as State>::Size>, Output = D::StateType>,
    <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
    DefaultAllocator: Allocator<<D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::VecLength>
        + Allocator<MsrSize>
        + Allocator<MsrSize, <D::StateType as State>::Size>
        + Allocator<MsrSize, MsrSize>
        + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<Accel>
        + Allocator<Accel, Accel>
        + Allocator<<D::StateType as State>::Size, Accel>
        + Allocator<Accel, <D::StateType as State>::Size>,
{
    /// PropInstance used for the estimation
    pub prop: PropInstance<'a, D>,
    /// Kalman filter itself
    pub kf: K,
    /// Tracking devices
    pub devices: BTreeMap<String, Trk>,
    /// Vector of estimates available after a pass
    pub estimates: Vec<K::Estimate>,
    /// Vector of residuals available after a pass
    pub residuals: Vec<Option<Residual<MsrSize>>>,
    pub ekf_trigger: Option<EkfTrigger>,
    /// Residual rejection criteria allows preventing bad measurements from affecting the estimation.
    pub resid_crit: Option<ResidRejectCrit>,
    pub almanac: Arc<Almanac>,
    init_state: D::StateType,
    _marker: PhantomData<Accel>,
}

impl<
        'a,
        D: Dynamics,
        MsrSize: DimName,
        Accel: DimName,
        K: Filter<D::StateType, Accel, MsrSize>,
        Trk: TrackerSensitivity<D::StateType, D::StateType>,
    > ODProcess<'a, D, MsrSize, Accel, K, Trk>
where
    D::StateType:
        Interpolatable + Add<OVector<f64, <D::StateType as State>::Size>, Output = D::StateType>,
    <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
    DefaultAllocator: Allocator<<D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::VecLength>
        + Allocator<MsrSize>
        + Allocator<MsrSize, <D::StateType as State>::Size>
        + Allocator<MsrSize, MsrSize>
        + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<Accel>
        + Allocator<Accel, Accel>
        + Allocator<<D::StateType as State>::Size, Accel>
        + Allocator<Accel, <D::StateType as State>::Size>,
{
    /// Initialize a new orbit determination process with an optional trigger to switch from a CKF to an EKF.
    pub fn new(
        prop: PropInstance<'a, D>,
        kf: K,
        devices: BTreeMap<String, Trk>,
        ekf_trigger: Option<EkfTrigger>,
        resid_crit: Option<ResidRejectCrit>,
        almanac: Arc<Almanac>,
    ) -> Self {
        let init_state = prop.state;
        Self {
            prop: prop.quiet(),
            kf,
            devices,
            estimates: Vec::with_capacity(10_000),
            residuals: Vec::with_capacity(10_000),
            ekf_trigger,
            resid_crit,
            almanac,
            init_state,
            _marker: PhantomData::<Accel>,
        }
    }

    /// Initialize a new orbit determination process with an Extended Kalman filter. The switch from a classical KF to an EKF is based on the provided trigger.
    pub fn ekf(
        prop: PropInstance<'a, D>,
        kf: K,
        devices: BTreeMap<String, Trk>,
        trigger: EkfTrigger,
        resid_crit: Option<ResidRejectCrit>,
        almanac: Arc<Almanac>,
    ) -> Self {
        let init_state = prop.state;
        Self {
            prop: prop.quiet(),
            kf,
            devices,
            estimates: Vec::with_capacity(10_000),
            residuals: Vec::with_capacity(10_000),
            ekf_trigger: Some(trigger),
            resid_crit,
            almanac,
            init_state,
            _marker: PhantomData::<Accel>,
        }
    }

    /// Allows to smooth the provided estimates. Returns the smoothed estimates or an error.
    ///
    /// Estimates must be ordered in chronological order. This function will smooth the
    /// estimates from the last in the list to the first one.
    pub fn smooth(&self, condition: SmoothingArc) -> Result<Vec<K::Estimate>, ODError> {
        let l = self.estimates.len() - 1;

        info!("Smoothing {} estimates until {}", l + 1, condition);
        let mut smoothed = Vec::with_capacity(self.estimates.len());
        // Set the first item of the smoothed estimates to the last estimate (we cannot smooth the very last estimate)
        smoothed.push(self.estimates.last().unwrap().clone());

        loop {
            let k = l - smoothed.len();
            // Borrow the previously smoothed estimate of the k+1 estimate
            let sm_est_kp1 = &self.estimates[k + 1];
            let x_kp1_l = sm_est_kp1.state_deviation();
            let p_kp1_l = sm_est_kp1.covar();
            // Borrow the k-th estimate, which we're smoothing with the next estimate
            let est_k = &self.estimates[k];
            // Borrow the k+1-th estimate, which we're smoothing with the next estimate
            let est_kp1 = &self.estimates[k + 1];

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
                SmoothingArc::All => {}
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
                .ok_or(ODError::SingularStateTransitionMatrix)?;

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

        // Note that we have yet to reverse the list, so we print them backward
        info!(
            "Smoothed {} estimates (from {} to {})",
            smoothed.len(),
            smoothed.last().unwrap().epoch(),
            smoothed[0].epoch(),
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

    /// Returns the root mean square of the prefit residual ratios
    pub fn rms_residual_ratios(&self) -> f64 {
        let mut sum = 0.0;
        for residual in self.residuals.iter().flatten() {
            sum += residual.ratio.powi(2);
        }
        (sum / (self.residuals.len() as f64)).sqrt()
    }

    /// Allows iterating on the filter solution. Requires specifying a smoothing condition to know where to stop the smoothing.
    pub fn iterate_arc(
        &mut self,
        arc: &TrackingDataArc,
        config: IterationConf,
    ) -> Result<(), ODError> {
        let mut best_rms = self.rms_residual_ratios();
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
                break;
            }

            iter_cnt += 1;

            // Prevent infinite loop when iterating prior to turning on the EKF.
            if let Some(trigger) = &mut self.ekf_trigger {
                trigger.reset();
            }

            info!("***************************");
            info!("*** Iteration number {iter_cnt:02} ***");
            info!("***************************");

            // First, smooth the estimates
            let smoothed = self.smooth(config.smoother)?;
            // Reset the propagator
            self.prop.state = self.init_state;
            // Empty the estimates and add the first smoothed estimate as the initial estimate
            self.estimates = Vec::with_capacity(arc.measurements.len().max(self.estimates.len()));
            self.residuals = Vec::with_capacity(arc.measurements.len().max(self.estimates.len()));

            self.kf.set_previous_estimate(&smoothed[0]);
            // And re-run the filter
            self.process_arc(arc)?;

            // Compute the new RMS
            let new_rms = self.rms_residual_ratios();
            let cur_rms_num = (new_rms - previous_rms).abs();
            let cur_rel_rms = cur_rms_num / previous_rms;
            if cur_rel_rms < config.relative_tol || cur_rms_num < config.absolute_tol * best_rms {
                if previous_rms < best_rms {
                    best_rms = previous_rms;
                }
                info!("*****************");
                info!("*** CONVERGED ***");
                info!("*****************");
                info!(
                    "New residual RMS: {:.5}\tPrevious RMS: {:.5}\tBest RMS: {:.5}",
                    new_rms, previous_rms, best_rms
                );
                if cur_rel_rms < config.relative_tol {
                    info!(
                        "Filter converged on relative tolerance ({:.2e} < {:.2e}) after {} iterations",
                        cur_rel_rms, config.relative_tol, iter_cnt
                    );
                } else {
                    info!(
                        "Filter converged on relative change ({:.2e} < {:.2e} * {:.2e}) after {} iterations",
                        cur_rms_num, config.absolute_tol, best_rms, iter_cnt
                    );
                }
                break;
            } else if new_rms > previous_rms {
                warn!(
                    "New residual RMS: {:.5}\tPrevious RMS: {:.5}\tBest RMS: {:.5} ({cur_rel_rms:.2e} > {:.2e})",
                    new_rms, previous_rms, best_rms, config.relative_tol
                );
                divergence_cnt += 1;
                previous_rms = new_rms;
                if divergence_cnt >= config.max_divergences {
                    let msg = format!(
                        "Filter iterations have continuously diverged {} times: {}",
                        config.max_divergences, config
                    );
                    if config.force_failure {
                        return Err(ODError::Diverged {
                            loops: config.max_divergences,
                        });
                    } else {
                        error!("{}", msg);
                        break;
                    }
                } else {
                    warn!("Filter iteration caused divergence {} of {} acceptable subsequent divergences", divergence_cnt, config.max_divergences);
                }
            } else {
                info!(
                    "New residual RMS: {:.5}\tPrevious RMS: {:.5}\tBest RMS: {:.5} ({cur_rel_rms:.2e} > {:.2e})",
                    new_rms, previous_rms, best_rms, config.relative_tol
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
                    return Err(ODError::Diverged {
                        loops: config.max_divergences,
                    });
                } else {
                    error!("{}", msg);
                    break;
                }
            }
        }

        Ok(())
    }

    /// Process the provided measurements for this orbit determination process given the associated devices.
    ///
    /// # Argument details
    /// + The measurements must be a list mapping the name of the measurement device to the measurement itself.
    /// + The name of all measurement devices must be present in the provided devices, i.e. the key set of `devices` must be a superset of the measurement device names present in the list.
    /// + The maximum step size to ensure we don't skip any measurements.
    #[allow(clippy::erasing_op)]
    pub fn process_arc(&mut self, arc: &TrackingDataArc) -> Result<(), ODError> {
        let measurements = &arc.measurements;
        ensure!(
            measurements.len() >= 2,
            TooFewMeasurementsSnafu {
                need: 2_usize,
                action: "running a Kalman filter"
            }
        );

        let max_step = match arc.min_duration_sep() {
            Some(step_size) => step_size,
            None => {
                return Err(ODError::TooFewMeasurements {
                    action: "determining the minimum step size",
                    need: 2,
                })
            }
        };

        ensure!(
            !max_step.is_negative() && max_step != Duration::ZERO,
            StepSizeSnafu { step: max_step }
        );

        // Check proper configuration.
        if MsrSize::USIZE > arc.unique_types().len() {
            error!("Filter misconfigured: expect high rejection count!");
            error!(
                "Arc only contains {} measurement types, but filter configured for {}.",
                arc.unique_types().len(),
                MsrSize::USIZE
            );
            error!("Filter should be configured for these numbers to match.");
            error!("Consider running subsequent arcs if ground stations provide different measurements.")
        }

        // Start by propagating the estimator.
        let num_msrs = measurements.len();

        // Update the step size of the navigation propagator if it isn't already fixed step
        if !self.prop.fixed_step {
            self.prop.set_step(max_step, false);
        }

        let prop_time = arc.end_epoch().unwrap() - self.kf.previous_estimate().epoch();
        info!("Navigation propagating for a total of {prop_time} with step size {max_step}");

        let mut epoch = self.prop.state.epoch();

        let mut reported = [false; 11];
        reported[0] = true; // Prevent showing "0% done"
        info!("Processing {num_msrs} measurements with covariance mapping");

        // We'll build a trajectory of the estimated states. This will be used to compute the measurements.
        let mut traj: Traj<D::StateType> = Traj::new();

        let mut msr_accepted_cnt: usize = 0;
        let tick = Epoch::now().unwrap();

        for (msr_cnt, (epoch_ref, msr)) in measurements.iter().enumerate() {
            let next_msr_epoch = *epoch_ref;

            // Advance the propagator
            loop {
                let delta_t = next_msr_epoch - epoch;

                // Propagator for the minimum time between the maximum step size, the next step size, and the duration to the next measurement.
                let next_step_size = delta_t.min(self.prop.step_size).min(max_step);

                // Remove old states from the trajectory
                // This is a manual implementation of `retaint` because we know it's a sorted vec, so no need to resort every time
                let mut index = traj.states.len();
                while index > 0 {
                    index -= 1;
                    if traj.states[index].epoch() >= epoch {
                        break;
                    }
                }
                traj.states.truncate(index);

                debug!("propagate for {next_step_size} (Δt to next msr: {delta_t})");
                let (_, traj_covar) = self
                    .prop
                    .for_duration_with_traj(next_step_size)
                    .context(ODPropSnafu)?;

                for state in traj_covar.states {
                    // NOTE: At the time being, only spacecraft estimation is possible, and the trajectory will always be the exact state
                    // that was propagated. Even once ground station biases are estimated, these won't go through the propagator.
                    traj.states.push(state);
                }

                // Now that we've advanced the propagator, let's see whether we're at the time of the next measurement.

                // Extract the state and update the STM in the filter.
                let nominal_state = self.prop.state;
                // Get the datetime and info needed to compute the theoretical measurement according to the model
                epoch = nominal_state.epoch();

                // Perform a measurement update
                if nominal_state.epoch() == next_msr_epoch {
                    // Get the computed observations
                    match self.devices.get_mut(&msr.tracker) {
                        Some(device) => {
                            if let Some(computed_meas) =
                                device.measure(epoch, &traj, None, self.almanac.clone())?
                            {
                                let msr_types = device.measurement_types();

                                // Switch back from extended if necessary
                                if let Some(trigger) = &mut self.ekf_trigger {
                                    if self.kf.is_extended() && trigger.disable_ekf(epoch) {
                                        self.kf.set_extended(false);
                                        info!("EKF disabled @ {epoch}");
                                    }
                                }

                                // Perform several measurement updates to ensure the desired dimensionality.
                                let windows = msr_types.len() / MsrSize::USIZE;
                                let mut msr_rejected = false;
                                for wno in 0..=windows {
                                    let mut cur_msr_types = IndexSet::new();
                                    for msr_type in msr_types
                                        .iter()
                                        .copied()
                                        .skip(wno * MsrSize::USIZE)
                                        .take(MsrSize::USIZE)
                                    {
                                        cur_msr_types.insert(msr_type);
                                    }

                                    if cur_msr_types.is_empty() {
                                        // We've processed all measurements.
                                        break;
                                    }

                                    // Check that the observation is valid.
                                    for val in
                                        msr.observation::<MsrSize>(&cur_msr_types).iter().copied()
                                    {
                                        ensure!(
                                            val.is_finite(),
                                            InvalidMeasurementSnafu {
                                                epoch: next_msr_epoch,
                                                val
                                            }
                                        );
                                    }

                                    let h_tilde = device
                                        .h_tilde::<MsrSize>(
                                            msr,
                                            &cur_msr_types,
                                            &nominal_state,
                                            self.almanac.clone(),
                                        )
                                        .unwrap();

                                    self.kf.update_h_tilde(h_tilde);

                                    let real_obs = msr.observation(&cur_msr_types);

                                    let computed_obs = computed_meas
                                        .observation::<MsrSize>(&cur_msr_types)
                                        - device.measurement_bias_vector::<MsrSize>(
                                            &cur_msr_types,
                                            epoch,
                                        )?;

                                    match self.kf.measurement_update(
                                        nominal_state,
                                        &real_obs,
                                        &computed_obs,
                                        device.measurement_covar_matrix(&cur_msr_types, epoch)?,
                                        self.resid_crit,
                                    ) {
                                        Ok((estimate, mut residual)) => {
                                            debug!("processed measurement #{msr_cnt} for {cur_msr_types:?} @ {epoch} from {}", device.name());

                                            residual.tracker = Some(device.name());
                                            residual.msr_types = cur_msr_types;

                                            if residual.rejected {
                                                msr_rejected = true;
                                            }

                                            // Switch to EKF if necessary, and update the dynamics and such
                                            // Note: we call enable_ekf first to ensure that the trigger gets
                                            // called in case it needs to save some information (e.g. the
                                            // StdEkfTrigger needs to store the time of the previous measurement).

                                            if let Some(trigger) = &mut self.ekf_trigger {
                                                if trigger.enable_ekf(&estimate)
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
                                                    self.prop.state = self.prop.state
                                                        + estimate.state_deviation();
                                                }
                                            }

                                            self.prop.state.reset_stm();

                                            self.estimates.push(estimate);
                                            self.residuals.push(Some(residual));
                                        }
                                        Err(e) => return Err(e),
                                    }
                                }
                                if !msr_rejected {
                                    msr_accepted_cnt += 1;
                                }
                            } else {
                                warn!("Ignoring observation @ {epoch} because simulated {} does not expect it", msr.tracker);
                            }
                        }
                        None => {
                            error!(
                                "Tracker {} is not in the list of configured devices",
                                msr.tracker
                            )
                        }
                    }

                    let msr_prct = (10.0 * (msr_cnt as f64) / (num_msrs as f64)) as usize;
                    if !reported[msr_prct] {
                        let num_rejected = msr_cnt - msr_accepted_cnt.saturating_sub(1);
                        let msg = format!(
                            "{:>3}% done - {msr_accepted_cnt:.0} measurements accepted, {:.0} rejected",
                            10 * msr_prct, num_rejected
                        );
                        if msr_accepted_cnt < num_rejected {
                            warn!("{msg}");
                        } else {
                            info!("{msg}");
                        }
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
                            // We push None so that the residuals and estimates are aligned
                            self.residuals.push(None);
                        }
                        Err(e) => return Err(e),
                    }
                    self.prop.state.reset_stm();
                }
            }
        }

        // Always report the 100% mark
        if !reported[10] {
            let tock_time = Epoch::now().unwrap() - tick;
            info!(
                "100% done - {msr_accepted_cnt:.0} measurements accepted, {:.0} rejected (done in {tock_time})",
                num_msrs - msr_accepted_cnt
            );
        }

        Ok(())
    }

    /// Continuously predicts the trajectory until the provided end epoch, with covariance mapping at each step. In other words, this performs a time update.
    pub fn predict_until(&mut self, step: Duration, end_epoch: Epoch) -> Result<(), ODError> {
        let prop_time = end_epoch - self.kf.previous_estimate().epoch();
        info!("Mapping covariance for {prop_time} with {step} step");

        loop {
            let mut epoch = self.prop.state.epoch();
            if epoch + self.prop.details.step > end_epoch {
                self.prop.until_epoch(end_epoch).context(ODPropSnafu)?;
            } else {
                self.prop.for_duration(step).context(ODPropSnafu)?;
            }

            // Perform time update

            // Extract the state and update the STM in the filter.
            let nominal_state = self.prop.state;
            // Get the datetime and info needed to compute the theoretical measurement according to the model
            epoch = nominal_state.epoch();
            // No measurement can be used here, let's just do a time update
            debug!("time update {epoch}");
            match self.kf.time_update(nominal_state) {
                Ok(est) => {
                    // State deviation is always zero for an EKF time update
                    // therefore we don't do anything different for an extended filter
                    self.estimates.push(est);
                    self.residuals.push(None);
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

    /// Continuously predicts the trajectory for the provided duration, with covariance mapping at each step. In other words, this performs a time update.
    pub fn predict_for(&mut self, step: Duration, duration: Duration) -> Result<(), ODError> {
        let end_epoch = self.kf.previous_estimate().epoch() + duration;
        self.predict_until(step, end_epoch)
    }

    /// Builds the navigation trajectory for the estimated state only
    pub fn to_traj(&self) -> Result<Traj<D::StateType>, NyxError>
    where
        DefaultAllocator: Allocator<<D::StateType as State>::VecLength>,
    {
        if self.estimates.is_empty() {
            Err(NyxError::NoStateData {
                msg: "No navigation trajectory to generate: run the OD process first".to_string(),
            })
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
        MsrSize: DimName,
        Accel: DimName,
        K: Filter<D::StateType, Accel, MsrSize>,
        Trk: TrackerSensitivity<D::StateType, D::StateType>,
    > ODProcess<'a, D, MsrSize, Accel, K, Trk>
where
    D::StateType:
        Interpolatable + Add<OVector<f64, <D::StateType as State>::Size>, Output = D::StateType>,
    <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
    DefaultAllocator: Allocator<<D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::VecLength>
        + Allocator<MsrSize>
        + Allocator<MsrSize, <D::StateType as State>::Size>
        + Allocator<MsrSize, MsrSize>
        + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<Accel>
        + Allocator<Accel, Accel>
        + Allocator<<D::StateType as State>::Size, Accel>
        + Allocator<Accel, <D::StateType as State>::Size>,
{
    pub fn ckf(
        prop: PropInstance<'a, D>,
        kf: K,
        devices: BTreeMap<String, Trk>,
        resid_crit: Option<ResidRejectCrit>,
        almanac: Arc<Almanac>,
    ) -> Self {
        let init_state = prop.state;
        Self {
            prop: prop.quiet(),
            kf,
            devices,
            estimates: Vec::with_capacity(10_000),
            residuals: Vec::with_capacity(10_000),
            resid_crit,
            ekf_trigger: None,
            init_state,
            almanac,
            _marker: PhantomData::<Accel>,
        }
    }
}
