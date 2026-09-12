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
use crate::od::msr::IntegrationRef;
pub use crate::od::snc::*;
pub use crate::od::*;
use crate::propagators::Propagator;
pub use crate::time::{Duration, Unit};
use anise::prelude::Almanac;
use indexmap::IndexSet;
use log::{debug, error, info, warn};
use msr::sensitivity::TrackerSensitivity;
use snafu::prelude::*;
use solution::kalman::KalmanVariant;
use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::ops::Add;
use typed_builder::TypedBuilder;

mod rejectcrit;
use self::kalman::KalmanFilter;
use self::msr::TrackingDataArc;
pub use self::rejectcrit::SigmaRejection;
mod solution;
pub use solution::{NormalizedConsistency, ODSolution};
mod initializers;

/// An orbit determination process (ODP) which filters OD measurements through a Kalman filter.
#[derive(Clone, TypedBuilder)]
#[builder(doc)]
#[allow(clippy::upper_case_acronyms)]
pub struct KalmanODProcess<
    D: Dynamics,
    MsrSize: DimName,
    Accel: DimName,
    Trk: TrackerSensitivity<D::StateType, D::StateType>,
> where
    D::StateType:
        Interpolatable + Add<OVector<f64, <D::StateType as State>::Size>, Output = D::StateType>,
    <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
    <DefaultAllocator as Allocator<<D::StateType as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>>::Buffer<f64>: Copy,
    DefaultAllocator: Allocator<<D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::VecLength>
        + Allocator<MsrSize>
        + Allocator<MsrSize, <D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::Size, MsrSize>
        + Allocator<MsrSize, MsrSize>
        + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<Accel>
        + Allocator<Accel, Accel>
        + Allocator<<D::StateType as State>::Size, Accel>
        + Allocator<Accel, <D::StateType as State>::Size>,
{
    /// Propagator used for the estimation
    pub prop: Propagator<D>,
    /// Kalman filter variant
    #[builder(default)]
    pub kf_variant: KalmanVariant,
    /// Residual rejection criteria allows preventing bad measurements from affecting the estimation.
    #[builder(default, setter(strip_option))]
    pub sigma_reject: Option<SigmaRejection>,
    /// Tracking devices
    #[builder(default_code = "BTreeMap::new()")]
    pub devices: BTreeMap<String, Trk>,
    /// A sets of process noise (usually noted Q), must be ordered chronologically
    #[builder(default_code = "vec![]")]
    pub process_noise: Vec<ProcessNoise<Accel>>,
    /// Maximum step size where the STM linearization is assumed correct (1 minute is usually fine)
    #[builder(default_code = "1 * Unit::Minute")]
    pub max_step: Duration,
    /// Precision of the measurement epoch when processing measurements.
    #[builder(default_code = "1 * Unit::Microsecond")]
    pub epoch_precision: Duration,
    pub almanac: Arc<Almanac>,
    #[builder(default_code = "PhantomData::<MsrSize>")]
    _msr_size: PhantomData<MsrSize>,
}

impl<
        D: Dynamics,
        MsrSize: DimName,
        Accel: DimName,
        Trk: TrackerSensitivity<D::StateType, D::StateType>,
    > KalmanODProcess<D, MsrSize, Accel, Trk>
where
    D::StateType:
        Interpolatable + Add<OVector<f64, <D::StateType as State>::Size>, Output = D::StateType>,
    <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
    <DefaultAllocator as Allocator<<D::StateType as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>>::Buffer<f64>: Copy,
    DefaultAllocator: Allocator<<D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::VecLength>
        + Allocator<MsrSize>
        + Allocator<MsrSize, <D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::Size, MsrSize>
        + Allocator<MsrSize, MsrSize>
        + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<Accel>
        + Allocator<Accel, Accel>
        + Allocator<<D::StateType as State>::Size, Accel>
        + Allocator<Accel, <D::StateType as State>::Size>
        + Allocator<nalgebra::Const<1>, MsrSize>,
{
    /// Process the provided tracking arc for this orbit determination process.
    #[allow(clippy::erasing_op)]
    pub fn process_arc(
        &self,
        initial_estimate: KfEstimate<D::StateType>,
        arc: &TrackingDataArc,
    ) -> Result<ODSolution<D::StateType, KfEstimate<D::StateType>, MsrSize, Trk>, ODError> {
        // Initialize the solution.
        let mut od_sol = ODSolution::new(self.devices.clone(), arc.unique_types());

        let measurements = &arc.measurements;
        ensure!(
            measurements.len() >= 2,
            TooFewMeasurementsSnafu {
                need: 2_usize,
                action: "running a Kalman filter"
            }
        );

        ensure!(
            !self.max_step.is_negative() && self.max_step != Duration::ZERO,
            StepSizeSnafu { step: self.max_step }
        );

        // Check proper configuration.
        if MsrSize::DIM > arc.unique_types().len() {
            error!("Filter misconfigured: expect high rejection count!");
            error!(
                "Arc only contains {} measurement types, but filter configured for {}.",
                arc.unique_types().len(),
                MsrSize::DIM
            );
            error!("Filter should be configured for these numbers to match.");
            error!("Consider running subsequent arcs if ground stations provide different measurements or switch to Scalar processing.");
        }

        let mut cfg_errors = vec![];
        for tracker in &arc.unique_aliases() {
            if let Some(device) = self.devices.get(tracker) {
                if let Err(e) =
                device.is_compatible(tracker, &arc.clone().filter_by_tracker(tracker.clone())) {
                    cfg_errors.push(e);
                }
            } else {
                error!("Tracker `{tracker}` from TrackingDataArc is not configured in the OD Process.");
                error!("Measurements from `{tracker}` will be ignored!");
            }
        }
        if !cfg_errors.is_empty() {
            // Return all the errors at once.
            let msg = cfg_errors.iter().map(|e| e.to_string()).collect::<Vec<String>>().join("\n---\n");
            return Err(ODError::ODConfigError { source: ConfigError::InvalidConfig { msg } })
        }

        // Start by propagating the estimator.
        let num_msrs = measurements.len();

        // Set up the propagator instance.
        let prop = self.prop.clone();
        let mut prop_instance = prop.with(initial_estimate.nominal_state().with_stm(), self.almanac.clone()).quiet();

        // Update the step size of the navigation propagator if it isn't already fixed step
        if !prop_instance.fixed_step {
            prop_instance.set_step(self.max_step, false);
        }

        let prop_time = arc.end_epoch().unwrap() - initial_estimate.epoch();
        info!("Navigation propagating for a total of {prop_time} with step size {}", self.max_step);

        let resid_crit = if arc.force_reject {
            warn!("Rejecting all measurements from {arc} (force_reject is True in TrackingDataArc)");
            Some(SigmaRejection { num_sigmas: 0.0 })
        } else {
            self.sigma_reject
        };

        let mut epoch = prop_instance.state.epoch();

        let mut reported = [false; 11];
        reported[0] = true; // Prevent showing "0% done"
        info!(
            "Processing {num_msrs} measurement epochs from {:?}",
            arc.unique_aliases()
        );

        // Set up the Kalman filter.
        let mut kf = KalmanFilter::<D::StateType, Accel> {
            prev_estimate: initial_estimate,
            process_noise: self.process_noise.clone(),
            variant: self.kf_variant,
            prev_used_snc: 0,
        };

        kf.initialize_process_noises();

        let mut devices = self.devices.clone();

        // We'll build a trajectory of the estimated states. This will be used to compute the measurements.
        let mut traj: Traj<D::StateType> = Traj::new();

        let mut msr_accepted_cnt: usize = 0;
        let mut msr_rejected_cnt: usize = 0;
        let mut unknown_trackers = IndexSet::new();
        let tick = Epoch::now().unwrap();

        for (msr_cnt, msr) in measurements.iter().enumerate() {
            let next_msr_epoch = msr.epoch;

            // Advance the propagator
            loop {
                let delta_t = next_msr_epoch - epoch;

                // Propagate for the minimum time between the maximum step size, the next step size, and the duration to the next measurement.
                let next_step_size = delta_t.min(prop_instance.step_size).min(self.max_step);

                // Remove any states at or after the current propagation epoch before appending new steps
                let keep_idx = traj.states.partition_point(|s| s.epoch() < epoch);
                traj.states.truncate(keep_idx);

                debug!("propagate for {next_step_size} (Δt to next msr: {delta_t})");
                let (latest_state, traj_covar) = prop_instance
                    .for_duration_with_traj(next_step_size)
                    .context(ODPropSnafu)?;

                for state in traj_covar.states {
                    // NOTE: At the time being, only spacecraft estimation is possible, and the trajectory will always be the exact state
                    // that was propagated. Even once ground station biases are estimated, these won't go through the propagator.
                    traj.states.push(state);
                }

                // Now that we've advanced the propagator, let's see whether we're at the time of the next measurement.

                // Extract the state and update the STM in the filter.
                let nominal_state = prop_instance.state;
                // Get the datetime and info needed to compute the theoretical measurement according to the model
                epoch = nominal_state.epoch();

                // Perform a measurement update, accounting for possible errors in measurement timestamps
                if (nominal_state.epoch() - next_msr_epoch).abs() < self.epoch_precision {
                    // Force the state epoch to match the measurement epoch exactly.
                    // This prevents infinite loops where the propagator (especially if fixed step)
                    // fails to step a tiny amount (drift) to reach the exact measurement time.
                    prop_instance.state.set_epoch(next_msr_epoch);

                    if msr.rejected {
                        debug!("Skipping manually rejected measurement at {epoch}");
                        let est = kf.time_update(nominal_state)?;
                        od_sol.push_time_update(est);
                        prop_instance.state.reset_stm();
                        msr_rejected_cnt += 1;
                    } else {
                        // Get the computed observations
                        match devices.get_mut(&msr.tracker) {
                            Some(device) => {
                                let msr_types = device.measurement_types().clone();

                                // Inspect the measurement to see if look-ahead is required for light-time computation
                                let pre_lookahead_len = traj.states.len();
                                if let Some(dop_cfg) = msr.doppler_config && dop_cfg.integration_ref != IntegrationRef::End {
                                    let lookahead_by = match dop_cfg.integration_ref {
                                        IntegrationRef::Start => dop_cfg.integration_time,
                                        IntegrationRef::Middle => dop_cfg.integration_time * 0.5,
                                        _ => unreachable!()
                                    };
                                    let mut lookahead_prop = prop.with(latest_state, self.almanac.clone()).quiet();
                                    let (_, lookahead_traj) = lookahead_prop
                                        .for_duration_with_traj(lookahead_by)
                                        .context(ODPropSnafu)?;

                                    for state in lookahead_traj.states.into_iter().filter(|s| s.epoch() > latest_state.epoch()) {
                                        traj.states.push(state);
                                    }
                                }

                                // Current nominal prior (needed for separate processing of simultaneous measurements)
                                let prior_nominal_state = prop_instance.state;
                                let mut current_state_estimate = prior_nominal_state;
                                let mut any_measurement_accepted = false;

                                // Perform several measurement updates to ensure the desired dimensionality.
                                let windows = msr_types.len() / MsrSize::DIM;
                                for wno in 0..=windows {
                                    // Update the nominal state in case we're ingesting several measurements
                                    // sequentially for the same epoch.
                                    let cur_msr_types = msr_types
                                        .iter()
                                        .copied()
                                        .skip(wno * MsrSize::DIM)
                                        .take(MsrSize::DIM)
                                        .collect::<IndexSet<_>>();

                                    if cur_msr_types.is_empty() {
                                        // We've processed all measurements.
                                        break;
                                    }

                                    // If this measurement type is unavailable, continue to the next one.
                                    if !msr.availability(&cur_msr_types)
                                        .iter()
                                        .any(|avail| *avail)
                                    {
                                        continue;
                                    }

                                    // Grab the un-modulo'd real observation
                                    let mut real_obs: OVector<f64, MsrSize> =
                                        msr.observation(&cur_msr_types);

                                    // Check that the observation is valid.
                                    for val in real_obs.iter().copied() {
                                        ensure!(
                                            val.is_finite(),
                                            InvalidMeasurementSnafu {
                                                epoch: next_msr_epoch,
                                                val
                                            }
                                        );
                                    }

                                    // Compute device specific matrices
                                    // Sensitivity (H tilde) is computed on the _pristine_ nominal estimate
                                    // i.e. it is not polluted by a partial measurement update if there are
                                    // multiple concurrent measurement processed sequentially.
                                    let h_tilde = device.h_tilde::<MsrSize>(
                                        msr,
                                        &cur_msr_types,
                                        &prior_nominal_state,
                                        &self.almanac,
                                    )?;

                                    let measurement_covar = device
                                        .measurement_covar_matrix(&cur_msr_types, epoch)?;

                                    // Evaluate the observation from the trajectory with the look-ahead states
                                    // but it does not include any of the states from the measurement update.
                                    let computed_meas_res = device.measure(epoch, &traj, None, &self.almanac);

                                    if let Some(computed_meas) = computed_meas_res?
                                    {
                                        // Apply any biases on the computed observation
                                        let obs_bias = device.measurement_bias_vector::<MsrSize>(
                                            &cur_msr_types,
                                            epoch,
                                        )?;

                                        let mut computed_obs = computed_meas
                                            .observation::<MsrSize>(&cur_msr_types)
                                            - obs_bias;

                                        // Apply the modulo to the real obs
                                        if let Some(moduli) = &arc.moduli {
                                            let mut obs_ambiguity =
                                                OVector::<f64, MsrSize>::zeros();

                                            for (i, msr_type) in cur_msr_types.iter().enumerate() {
                                                if let Some(modulus) = moduli.get(msr_type) {
                                                    let k = computed_obs[i].div_euclid(*modulus);
                                                    // real_obs = measured_obs + k * modulus
                                                    obs_ambiguity[i] = k * *modulus;
                                                }
                                            }
                                            real_obs += obs_ambiguity;
                                        }

                                        // Map prior shifts to account for partial state update: h_eff = h(x0) + H * (x_curr - x0)
                                        let delta_state = current_state_estimate.to_state_vector() - prior_nominal_state.to_state_vector();
                                        let obs_shift = &h_tilde * delta_state;
                                        computed_obs += obs_shift;

                                        // Kalman measurement update on the filter covariance and state
                                        let (estimate, mut residual, gain) = kf.measurement_update(
                                            current_state_estimate,
                                            real_obs,
                                            computed_obs,
                                            measurement_covar,
                                            h_tilde,
                                            resid_crit,
                                        )?;

                                        debug!(
                                            "processed measurement #{msr_cnt} for {cur_msr_types:?} @ {epoch} from {}",
                                            device.name()
                                        );

                                        if residual.rejected {
                                            msr_rejected_cnt += 1;
                                        } else {
                                            msr_accepted_cnt += 1;
                                            current_state_estimate = estimate.state();
                                            any_measurement_accepted = true;
                                        }


                                        residual.tracker = Some(device.name());
                                        residual.msr_types = cur_msr_types;
                                        od_sol.push_measurement_update(estimate, residual, gain);
                                    } else {
                                        debug!(
                                            "Device {} does not expect measurement at {epoch}, skipping",
                                            msr.tracker
                                        );
                                        msr_rejected_cnt += 1;
                                    }
                                }

                                // Strip the temporary states to maintain trajectory causality
                                traj.states.truncate(pre_lookahead_len);

                                if any_measurement_accepted && kf.replace_state() {
                                    // Only update the state of the EKF if at least one residual was not rejected.
                                    prop_instance.state = current_state_estimate;
                                    traj.states.pop();
                                    traj.states.push(prop_instance.state);
                                }

                                // Reset the STM strictly once per epoch, after all updates have been absorbed
                                prop_instance.state.reset_stm();
                            }
                            None => {
                                if !unknown_trackers.contains(&msr.tracker) {
                                    error!(
                                        "Tracker {} is not in the list of configured devices",
                                        msr.tracker
                                    );
                                    unknown_trackers.insert(msr.tracker.clone());
                                }
                            }
                        }
                    }

                    let msr_prct = (10.0 * (msr_cnt as f64) / (num_msrs as f64)) as usize;
                    if !reported[msr_prct] {
                        let msg = format!(
                            "{:>3}% done - {epoch} - {msr_accepted_cnt:.0} measurements accepted, {:.0} rejected",
                            10 * msr_prct, msr_rejected_cnt
                        );
                        if msr_accepted_cnt < msr_rejected_cnt {
                            warn!("{msg}");
                        } else {
                            info!("{msg}");
                        }
                        reported[msr_prct] = true;
                    }

                    break;
                } else {
                    // No measurement can be used here, let's just do a time update and continue advancing the propagator.
                    // State deviation is always zero for an EKF time update so we don't do anything different than for a CKF.
                    let est = kf.time_update(nominal_state)?;
                    od_sol.push_time_update(est);
                    prop_instance.state.reset_stm();
                }
            }
        }

        // Always report the 100% mark
        if !reported[10] {
            let tock_time = Epoch::now().unwrap() - tick;
            info!(
                "100% done - {epoch} - {msr_accepted_cnt} measurements accepted, {msr_rejected_cnt} rejected (done in {tock_time})",
            );
        }

        Ok(od_sol)
    }

    /// Perform a time update. Continuously predicts the trajectory until the provided end epoch, with covariance mapping at each step.
    pub fn predict_until(
        &self,
        initial_estimate: KfEstimate<D::StateType>,
        end_epoch: Epoch,
    ) -> Result<ODSolution<D::StateType, KfEstimate<D::StateType>, MsrSize, Trk>, ODError> {
        // Initialize the solution with no measurement types.
        let mut od_sol = ODSolution::new(self.devices.clone(), IndexSet::new());

        od_sol.push_time_update(initial_estimate);

        // Set up the propagator instance.
        let prop = self.prop.clone();
        let mut prop_instance = prop.with(initial_estimate.nominal_state().with_stm(), self.almanac.clone()).quiet();


        // Set up the Kalman filter.
        let mut kf = KalmanFilter::<D::StateType, Accel> {
            prev_estimate: initial_estimate,
            process_noise: self.process_noise.clone(),
            variant: self.kf_variant,
            prev_used_snc: 0,
        };

        let prop_time = end_epoch - kf.previous_estimate().epoch();
        info!("Mapping covariance for {prop_time} every {} until {end_epoch}", self.max_step);

        loop {
            let nominal_state = prop_instance.for_duration(self.max_step).context(ODPropSnafu)?;
            // Extract the state and update the STM in the filter.
            // Get the datetime and info needed to compute the theoretical measurement according to the model
            let epoch = nominal_state.epoch();
            // No measurement can be used here, let's just do a time update
            debug!("time update {epoch}");
            let est = kf.time_update(nominal_state)?;
            od_sol.push_time_update(est);
            prop_instance.state.reset_stm();
            if epoch >= end_epoch {
                break;
            }
        }

        Ok(od_sol)
    }

    /// Perform a time update. Continuously predicts the trajectory for the provided duration, with covariance mapping at each step.
    pub fn predict_for(
        &self,
        initial_estimate: KfEstimate<D::StateType>,
        duration: Duration,
    ) -> Result<ODSolution<D::StateType, KfEstimate<D::StateType>, MsrSize, Trk>, ODError> {
        let end_epoch = initial_estimate.nominal_state().epoch() + duration;
        self.predict_until(initial_estimate, end_epoch)
    }
}
