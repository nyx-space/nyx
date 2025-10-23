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

#![allow(clippy::type_complexity)] // Allow complex types for generics
#![allow(unused_imports)] // Keep imports for context even if slightly unused in snippet
use crate::linalg::allocator::Allocator;
use crate::linalg::{Const, DefaultAllocator, DimName, OMatrix, OVector, U1}; // Use U1 for MsrSize
use crate::md::trajectory::{Interpolatable, Traj}; // May not need Traj if we propagate point-to-point
pub use crate::od::estimate::*;
pub use crate::od::ground_station::*;
pub use crate::od::snc::*; // SNC not typically used in BLS, but keep context
pub use crate::od::*;
use crate::propagators::Propagator;
pub use crate::time::{Duration, Epoch, Unit};
use anise::prelude::Almanac;
use indexmap::IndexSet;
use log::error;
use log::{debug, info, trace, warn};
use msr::sensitivity::TrackerSensitivity; // Assuming this is the correct path
use nalgebra::{Cholesky, Dyn, Matrix, VecStorage};
use snafu::prelude::*;
use solution::msr::MeasurementType;
use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::ops::Add;
use std::sync::Arc;
use typed_builder::TypedBuilder;

mod solution;

pub use solution::BLSSolution;

use self::msr::TrackingDataArc;

/// Solver choice for the Batch Least Squares estimator
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BLSSolver {
    /// Standard Normal Equations: (H^T W H) dx = H^T W dy
    NormalEquations,
    /// Levenberg-Marquardt: (H^T W H + lambda * D^T D) dx = H^T W dy
    LevenbergMarquardt,
}

/// Configuration for the Batch Least Squares estimator
#[derive(Clone, TypedBuilder)]
#[builder(doc)]
pub struct BatchLeastSquares<
    D: Dynamics,
    Trk: TrackerSensitivity<D::StateType, D::StateType>, // Use the same TrackerSensitivity
> where
    D::StateType:
        Interpolatable + Add<OVector<f64, <D::StateType as State>::Size>, Output = D::StateType>,
    <D::StateType as State>::Size: DimName, // Add DimName bound for state size
    // Add Allocator constraints similar to KalmanODProcess, but using U1 for MsrSize
    <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
    <DefaultAllocator as Allocator<<D::StateType as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>>::Buffer<f64>: Copy,
    DefaultAllocator: Allocator<<D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::VecLength>
        + Allocator<U1> // MsrSize is U1
        + Allocator<U1, <D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::Size, U1>
        + Allocator<U1, U1>
        + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>,
{
    /// Propagator used for the estimation reference trajectory
    pub prop: Propagator<D>,
    /// Tracking devices
    pub devices: BTreeMap<String, Trk>,
    /// Solver method
    #[builder(default = BLSSolver::NormalEquations)]
    pub solver: BLSSolver,
    /// Convergence tolerance on the norm of the correction on position, in kilometers
    #[builder(default = 1e-4)]
    pub tolerance_pos_km: f64,
    /// Maximum number of iterations
    #[builder(default = 10)]
    pub max_iterations: usize,
    /// Maximum step size where the STM linearization is assumed correct
    /// (30 seconds is usually fine, but too large and info matrix could be singular)
    #[builder(default_code = "30 * Unit::Second")]
    pub max_step: Duration,
    /// Precision of the measurement epoch when processing measurements.
    #[builder(default_code = "1 * Unit::Microsecond")]
    pub epoch_precision: Duration,
    /// Initial damping factor for Levenberg-Marquardt
    #[builder(default = 10.0)]
    pub lm_lambda_init: f64,
    /// Factor to decrease lambda by in LM
    #[builder(default = 10.0)] // Decrease aggressively if step is good
    pub lm_lambda_decrease: f64,
    /// Factor to increase lambda by in LM
    #[builder(default = 10.0)] // Increase aggressively if step is bad
    pub lm_lambda_increase: f64,
    /// Minimum value for LM lambda
    #[builder(default = 1e-12)]
    pub lm_lambda_min: f64,
    /// Maximum value for LM lambda
    #[builder(default = 1e12)]
    pub lm_lambda_max: f64,
    /// Use diagonal scaling (D = sqrt(diag(H^T W H))) in LM
    #[builder(default = true)]
    pub lm_use_diag_scaling: bool,
    pub almanac: Arc<Almanac>,
}

#[allow(type_alias_bounds)]
type StateMatrix<D: Dynamics> =
    OMatrix<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>;

impl<D, Trk> BatchLeastSquares<D, Trk>
where
    D: Dynamics,
    Trk: TrackerSensitivity<D::StateType, D::StateType> + Clone, // Add Clone requirement for Trk
    D::StateType: Interpolatable
        + Add<OVector<f64, <D::StateType as State>::Size>, Output = D::StateType>
        + std::fmt::Debug, // Add Debug for logging
    <D::StateType as State>::Size: DimName,
    <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
    <DefaultAllocator as Allocator<<D::StateType as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>>::Buffer<f64>: Copy,
    DefaultAllocator: Allocator<<D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::VecLength>
        + Allocator<U1>
        + Allocator<U1, <D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::Size, U1>
        + Allocator<U1, U1>
        + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>,
{
    /// Processes a tracking data arc to estimate the state using Batch Least Squares.
    pub fn estimate(
        &self,
        initial_guess: D::StateType,
        arc: &TrackingDataArc,
    ) -> Result<BLSSolution<D::StateType>, ODError> {
        let measurements = &arc.measurements;
        let num_measurements = measurements.len();
        let mut devices = self.devices.clone();

        ensure!(
            num_measurements >= 2,
            TooFewMeasurementsSnafu {
                need: 2_usize,
                action: "BLSE"
            }
        );

        info!(
            "Using {:?} in the Batch Least Squares estimation with {num_measurements} measurements",
            self.solver
        );
        info!("Initial guess: {}", initial_guess.orbit());

        let mut current_estimate = initial_guess;
        let mut current_covariance = StateMatrix::<D>::zeros();
        let mut converged = false;
        let mut corr_pos_km = f64::MAX;
        let mut lambda = self.lm_lambda_init;
        let mut current_rms = f64::MAX;
        let mut iter: usize = 0;

        let mut unknown_trackers = IndexSet::new();

        // --- Iteration Loop ---
        while iter < self.max_iterations {
            iter += 1;
            info!("[{iter}/{}] Current estimate: {}", self.max_iterations, current_estimate.orbit());

            // Re-initialize matrices for this iteration
            // Information Matrix: Lambda = H^T * W * H
            let mut info_matrix = StateMatrix::<D>::identity();
            // Normal Matrix: N = H^T * W * dy
            let mut normal_matrix = OVector::<f64, <D::StateType as State>::Size>::zeros();
            // Sum of squares of weighted residuals for RMS calculation and LM cost
            let mut sum_sq_weighted_residuals = 0.0;

            // Set up a single propagator for the whole iteration.
            let mut prop_inst = self.prop.with(current_estimate.with_stm(), self.almanac.clone()).quiet();
            let mut epoch = current_estimate.epoch();

            // Store the STM to the start of the batch.
            let mut stm = StateMatrix::<D>::identity();

            for (epoch_ref, msr) in measurements.iter() {
                let msr_epoch = *epoch_ref;

                loop {
                    let delta_t = msr_epoch - epoch;
                    if delta_t <= Duration::ZERO {
                        // Move onto the next measurement.
                        break;
                    }

                    // Propagate for the minimum time between the maximum step size, the next step size, and the duration to the next measurement.
                    let next_step = delta_t.min(prop_inst.step_size).min(self.max_step);

                    // Propagate reference state from the previous state to msr_epoch
                    let this_state = prop_inst.for_duration(next_step).context(ODPropSnafu)?;
                    epoch = this_state.epoch();

                    // Grab the STM Phi(t_{i+1}, t_i) from the propagated state's STM.
                    let step_stm = this_state.stm().expect("STM unavailable");
                    // Compute the STM Phi(t_{i+1}, t_0) = Phi(t_{i+1}, t_i) * Phi(t_i, t_0)
                    stm = step_stm * stm;

                    if (epoch - msr_epoch).abs() < self.epoch_precision {
                        // Get the correct tracking device
                        let device = match devices.get_mut(&msr.tracker) {
                            Some(d) => d,
                            None => {
                                if !unknown_trackers.contains(&msr.tracker) {
                                    error!(
                                        "Tracker {} is not in the list of configured devices",
                                        msr.tracker
                                    );
                                }
                                unknown_trackers.insert(msr.tracker.clone());
                                continue;
                            }
                        };

                        for msr_type in msr.data.keys().copied() {
                            let mut msr_types = IndexSet::new();
                            msr_types.insert(msr_type);

                            let h_tilde = device
                            .h_tilde::<U1>(msr, &msr_types, &this_state, self.almanac.clone())?;

                            // Compute expected measurement H(X(t_i))
                            let computed_meas_opt = device
                                .measure_instantaneous(this_state, None, self.almanac.clone())?;

                            let computed_meas = match computed_meas_opt {
                                Some(cm) => cm,
                                None => {
                                    debug!("Device {} does not expect measurement at epoch {msr_epoch}, skipping", msr.tracker);
                                    continue;
                                }
                            };

                            // Get the computed observation value
                            let computed_obs = computed_meas.observation::<U1>(&msr_types)[0];

                            // Get real observation y_i
                            let real_obs = msr.observation::<U1>(&msr_types)[0];

                            // Sanity check measurement value
                            ensure!(
                                real_obs.is_finite(),
                                InvalidMeasurementSnafu {
                                    epoch: msr_epoch,
                                    val: real_obs
                                }
                            );

                            // Compute residual dy = y_i - H(X(t_i))
                            let residual = real_obs - computed_obs;

                            // Get measurement variance R (assuming 1x1 matrix) and weight W = 1/R
                            let r_matrix = device
                                .measurement_covar_matrix::<U1>(&msr_types, msr_epoch)?;
                            let r_variance = r_matrix[(0, 0)];

                            ensure!(r_variance > 0.0, SingularNoiseRkSnafu);
                            let weight = 1.0 / r_variance;

                            // Compute H_matrix = H_tilde * Phi(t_i, t_0) (sensitivity wrt initial state X_0)
                            let h_matrix = h_tilde * stm;

                            // Accumulate Information Matrix: info_matrix += H^T * W * H
                            // Recall that the weight is a scalar, so we can move it to the end of the operation.
                            info_matrix += h_matrix.transpose() * &h_matrix * weight;

                            // Accumulate Normal Matrix: normal_matrix += H^T * W * y
                            normal_matrix += h_matrix.transpose() * residual * weight;

                            // Accumulate sum of squares of weighted residuals
                            sum_sq_weighted_residuals += weight * residual * residual;
                        }
                    }
                }
            }

            // --- Solve for State Correction dx ---
            let state_correction: OVector<f64, <D::StateType as State>::Size>;
            let iteration_cost_decreased; // For LM logic

            // Use num_measurements for consistency
            let current_iter_rms = (sum_sq_weighted_residuals / num_measurements as f64).sqrt();

            match self.solver {
                BLSSolver::NormalEquations => {
                    // Solve Lambda * dx = N => dx = Lambda^-1 * N
                    let info_matrix_chol = match info_matrix.cholesky() {
                         Some(chol) => chol,
                         None => return Err(ODError::SingularInformationMatrix)
                    };
                    state_correction = info_matrix_chol.solve(&normal_matrix);
                    // Assume NE always decreases cost locally
                    iteration_cost_decreased = true;
                    current_rms = current_iter_rms;
                }
                BLSSolver::LevenbergMarquardt => {
                     // Solve (Lambda + lambda * D^T D) * dx = N
                     // D^T D is a diagonal scaling matrix.
                     // Common choices: D^T D = I or D^T D = diag(Lambda)
                    let mut d_sq = StateMatrix::<D>::identity();
                    if self.lm_use_diag_scaling {
                        // Use D^T D = diag(Lambda)
                        for i in 0..6 {
                            d_sq[(i, i)] = info_matrix.diagonal()[i];
                        }
                        // Ensure diagonal elements are positive for stability
                        for i in 0..6 {
                            if d_sq[(i, i)] <= 0.0 {
                                d_sq[(i, i)] = 1e-6; // Set a small positive floor
                                warn!("LM Scaling: Found non-positive diagonal element {} in H^TWH, using floor.", info_matrix[(i, i)]);
                            }
                        }
                    } // else d_sq remains Identity

                    // Inner LM loop to find suitable lambda
                    let augmented_matrix = info_matrix + d_sq * lambda;

                    if let Some(aug_chol) = augmented_matrix.cholesky() {
                        state_correction = aug_chol.solve(&normal_matrix);

                        // --- LM Lambda Update Logic ---
                        // Simple strategy: Check if RMS decreased. More robust methods exist.
                        // For a simple check, we compare current_iter_rms with the previous iteration's RMS.
                        if current_iter_rms < current_rms || iter == 0 {
                            // Cost (approximated by RMS) decreased or first iteration
                            iteration_cost_decreased = true;
                            // Decrease damping
                            lambda /= self.lm_lambda_decrease;
                            // Clamp to min
                            lambda = lambda.max(self.lm_lambda_min);
                            debug!("LM: Cost decreased (RMS {current_rms} -> {current_iter_rms}). Decreasing lambda to {lambda}");
                            current_rms = current_iter_rms;
                        } else {
                             // Cost increased or stalled
                             iteration_cost_decreased = false;
                             // Increase damping
                             lambda *= self.lm_lambda_increase;
                             // Clamp to max
                             lambda = lambda.min(self.lm_lambda_max);
                             debug!("LM: Cost increased/stalled (RMS {current_rms} -> {current_iter_rms}). Increasing lambda to {lambda}");
                             // Don't update current_rms baseline if cost increased
                        }

                    } else {
                        // Augmented matrix is singular, increase lambda significantly and retry
                        warn!("LM: Augmented matrix (H^TWH + lambda*D^2) singular with lambda={lambda}. Increasing lambda.");
                        lambda *= self.lm_lambda_increase * 10.0; // Increase more aggressively
                        lambda = lambda.min(self.lm_lambda_max);
                        // Skip update in this iteration, force retry with larger lambda next time if possible
                        // Skip the rest of the loop and go to next iteration
                        continue;
                    }
                }
            }

            // --- Update State Estimate ---
            // Only update if the step is considered successful (esp. for LM)
            // Also hit if using normal equations because iteration_cost_decreased is forced to true
            if iteration_cost_decreased {
                current_estimate = current_estimate + state_correction;
                corr_pos_km = state_correction.fixed_rows::<3>(0).norm();

                let corr_vel_km_s = state_correction.fixed_rows::<3>(3).norm();
                info!(
                    "[{iter}/{}] RMS: {current_iter_rms:.3}; corrections: {:.3} m\t{:.3} m/s",
                    self.max_iterations,
                    corr_pos_km * 1e3,
                    corr_vel_km_s * 1e3
                );

                // Update the covariance
                current_covariance = match info_matrix.try_inverse() {
                    Some(cov) => cov,
                    None => {
                        warn!("Information matrix H^TWH is singular.");
                        StateMatrix::<D>::identity()
                    }
               };

                // --- Check Convergence ---
                if corr_pos_km < self.tolerance_pos_km {
                    info!("Converged in {iter} iterations.");
                    converged = true;
                    break;
                }
            } else if self.solver == BLSSolver::LevenbergMarquardt {
                 // LM step was rejected (cost increased)
                 info!("[{iter}/{}] LM: Step rejected, increasing lambda.", self.max_iterations);
                 // Reset correction norm as step was bad
                 corr_pos_km = f64::MAX;
                 // The loop will continue with the increased lambda
            }
        }

        if !converged {
            warn!("Not converged after {} iterations.", self.max_iterations);
        }

        info!("Batch Least Squares estimation completed.");
        Ok(BLSSolution {
            estimated_state: current_estimate,
            covariance: current_covariance,
            num_iterations: iter,
            final_rms: current_rms,
            final_corr_pos_km: corr_pos_km,
            converged,
        })
    }
}
