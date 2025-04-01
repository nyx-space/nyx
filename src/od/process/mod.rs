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
use crate::propagators::Propagator;
pub use crate::time::{Duration, Unit};
use anise::prelude::Almanac;
use indexmap::IndexSet;
use msr::sensitivity::TrackerSensitivity;
use snafu::prelude::*;
use solution::kalman::KalmanVariant;
use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::ops::Add;

mod conf;
pub use conf::{IterationConf, SmoothingArc};
mod rejectcrit;
use self::kalman::KalmanFilter;
use self::msr::TrackingDataArc;
pub use self::rejectcrit::ResidRejectCrit;
mod solution;
pub use solution::ODSolution;

/// An orbit determination process (ODP) which filters OD measurements through a Kalman filter.
#[derive(Clone)]
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
    pub kf_variant: KalmanVariant,
    /// Residual rejection criteria allows preventing bad measurements from affecting the estimation.
    pub resid_crit: Option<ResidRejectCrit>,
    /// Tracking devices
    pub devices: BTreeMap<String, Trk>,
    /// A sets of process noise (usually noted Q), must be ordered chronologically
    pub process_noise: Vec<ProcessNoise<Accel>>,
    pub almanac: Arc<Almanac>,
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
    /// Initialize a new orbit determination process with an optional trigger to switch from a CKF to an EKF.
    pub fn new(
        prop: Propagator<D>,
        kf_variant: KalmanVariant,
        resid_crit: Option<ResidRejectCrit>,
        devices: BTreeMap<String, Trk>,
        almanac: Arc<Almanac>,
    ) -> Self {
        Self {
            prop,
            kf_variant,
            devices,
            resid_crit,
            process_noise: vec![],
            almanac,
            _msr_size: PhantomData::<MsrSize>,
        }
    }

    /// Set (or replaces) the existing process noise configuration.
    pub fn from_process_noise(
        prop: Propagator<D>,
        kf_variant: KalmanVariant,
        devices: BTreeMap<String, Trk>,
        resid_crit: Option<ResidRejectCrit>,
        process_noise: ProcessNoise<Accel>,
        almanac: Arc<Almanac>,
    ) -> Self {
        let mut me = Self::new(
            prop,
            kf_variant,
            resid_crit,
            devices,
            almanac
        );
        me.process_noise.push(process_noise);
        
        me
    }

    /// Set (or replaces) the existing process noise configuration.
    pub fn with_process_noise(mut self, process_noise: ProcessNoise<Accel>) -> Self {
        self.process_noise.clear();
        self.process_noise.push(process_noise);
        self
    }

    /// Pushes the provided process noise to the list the existing process noise configurations.
    pub fn and_with_process_noise(mut self, process_noise: ProcessNoise<Accel>) -> Self {
        self.process_noise.push(process_noise);
        self
    }

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

        // Set up the propagator instance.
        let prop = self.prop.clone();
        let mut prop_instance = prop.with(initial_estimate.nominal_state().with_stm(), self.almanac.clone()).quiet();

        // Update the step size of the navigation propagator if it isn't already fixed step
        if !prop_instance.fixed_step {
            prop_instance.set_step(max_step, false);
        }

        
        let prop_time = arc.end_epoch().unwrap() - initial_estimate.epoch();
        info!("Navigation propagating for a total of {prop_time} with step size {max_step}");

        let resid_crit = if arc.force_reject {
            warn!("Rejecting all measurements from {arc} as requested");
            Some(ResidRejectCrit { num_sigmas: 0.0 })
        } else {
            self.resid_crit
        };

        let mut epoch = prop_instance.state.epoch();

        let mut reported = [false; 11];
        reported[0] = true; // Prevent showing "0% done"
        info!(
            "Processing {num_msrs} measurements from {:?}",
            arc.unique_aliases()
        );

        // Set up the Kalman filter.
        let mut kf = KalmanFilter::<D::StateType, Accel> {
            prev_estimate: initial_estimate,
            process_noise: self.process_noise.clone(),
            variant: self.kf_variant,
            prev_used_snc: 0,
        };

        let mut devices = self.devices.clone();

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
                let next_step_size = delta_t.min(prop_instance.step_size).min(max_step);

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
                let (_, traj_covar) = prop_instance
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
                // TODO: Move epoch precision to process configuration.
                if (nominal_state.epoch() - next_msr_epoch).abs() < Unit::Microsecond * 1 {
                    // Get the computed observations
                    match devices.get_mut(&msr.tracker) {
                        Some(device) => {
                            if let Some(computed_meas) =
                                device.measure(epoch, &traj, None, self.almanac.clone())?
                            {
                                let msr_types = device.measurement_types();

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

                                    // If this measurement type is unavailable, continue to the next one.
                                    if !msr.availability(&cur_msr_types).iter().any(|avail| *avail)
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
                                                epoch: *epoch_ref,
                                                val
                                            }
                                        );
                                    }

                                    // Compute device specific matrices
                                    let h_tilde = device.h_tilde::<MsrSize>(
                                        msr,
                                        &cur_msr_types,
                                        &nominal_state,
                                        self.almanac.clone(),
                                    )?;

                                    let measurement_covar =
                                        device.measurement_covar_matrix(&cur_msr_types, epoch)?;

                                    // Apply any biases on the computed observation
                                    let computed_obs = computed_meas
                                        .observation::<MsrSize>(&cur_msr_types)
                                        - device.measurement_bias_vector::<MsrSize>(
                                            &cur_msr_types,
                                            epoch,
                                        )?;

                                    // Apply the modulo to the real obs
                                    if let Some(moduli) = &arc.moduli {
                                        let mut obs_ambiguity = OVector::<f64, MsrSize>::zeros();

                                        for (i, msr_type) in cur_msr_types.iter().enumerate() {
                                            if let Some(modulus) = moduli.get(msr_type) {
                                                let k = computed_obs[i].div_euclid(*modulus);
                                                // real_obs = measured_obs + k * modulus
                                                obs_ambiguity[i] = k * *modulus;
                                            }
                                        }
                                        real_obs += obs_ambiguity;
                                    }

                                    match kf.measurement_update(
                                        nominal_state,
                                        real_obs,
                                        computed_obs,
                                        measurement_covar,
                                        h_tilde,
                                        resid_crit,
                                    ) {
                                        Ok((estimate, mut residual, gain)) => {
                                            debug!("processed measurement #{msr_cnt} for {cur_msr_types:?} @ {epoch} from {}", device.name());

                                            residual.tracker = Some(device.name());
                                            residual.msr_types = cur_msr_types;

                                            if residual.rejected {
                                                msr_rejected = true;
                                            }

                                            if kf.replace_state() {
                                                prop_instance.state = estimate.state();
                                            }

                                            prop_instance.state.reset_stm();

                                            od_sol
                                                .push_measurement_update(estimate, residual, gain);
                                        }
                                        Err(e) => return Err(e),
                                    }
                                }
                                if !msr_rejected {
                                    msr_accepted_cnt += 1;
                                }
                            } else {
                                debug!("ignoring observation @ {epoch} because simulated {} does not expect it", msr.tracker);
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
                    debug!("time update {epoch:?}, next msr {next_msr_epoch:?}");
                    match kf.time_update(nominal_state) {
                        Ok(est) => {
                            // State deviation is always zero for an EKF time update so we don't do anything different than for a CKF.
                            od_sol.push_time_update(est);
                        }
                        Err(e) => return Err(e),
                    }
                    prop_instance.state.reset_stm();
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

        Ok(od_sol)
    }

    /// Perform a time update. Continuously predicts the trajectory until the provided end epoch, with covariance mapping at each step.
    pub fn predict_until(
        &self,
        initial_estimate: KfEstimate<D::StateType>,
        step: Duration,
        end_epoch: Epoch,
    ) -> Result<ODSolution<D::StateType, KfEstimate<D::StateType>, MsrSize, Trk>, ODError> {
        // Initialize the solution with no measurement types.
        let mut od_sol = ODSolution::new(self.devices.clone(), IndexSet::new());

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
        info!("Mapping covariance for {prop_time} until {end_epoch} with {step} step");

        loop {
            let nominal_state = prop_instance.for_duration(step).context(ODPropSnafu)?;
            // Extract the state and update the STM in the filter.
            // Get the datetime and info needed to compute the theoretical measurement according to the model
            let epoch = nominal_state.epoch();
            // No measurement can be used here, let's just do a time update
            debug!("time update {epoch}");
            match kf.time_update(nominal_state) {
                Ok(est) => {
                    od_sol.push_time_update(est);
                }
                Err(e) => return Err(e),
            }
            prop_instance.state.reset_stm();
            if epoch >= end_epoch {
                break;
            }
        }

        Ok(od_sol)
    }

    /// Perform a time update. Continuously predicts the trajectory for the provided duration, with covariance mapping at each step. In other words, this performs a time update.
    pub fn predict_for(
        &self,
        initial_estimate: KfEstimate<D::StateType>,
        step: Duration,
        duration: Duration,
    ) -> Result<ODSolution<D::StateType, KfEstimate<D::StateType>, MsrSize, Trk>, ODError> {
        let end_epoch = initial_estimate.nominal_state().epoch() + duration;
        self.predict_until(initial_estimate, step, end_epoch)
    }
}
