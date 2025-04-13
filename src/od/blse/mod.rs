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

// Define potential errors
#[derive(Debug, Snafu)]
#[snafu(visibility(pub(crate)))]
pub enum BLSError {
    #[snafu(display("Propagation failed: {source}"))]
    Propagation { source: ODError },

    #[snafu(display("Device '{device_name}' not found in configuration"))]
    DeviceNotFound { device_name: String },

    #[snafu(display("Measurement error from device: {source}"))]
    MeasurementError { source: ODError },

    #[snafu(display("Sensitivity (H_tilde) computation error: {source}"))]
    SensitivityError { source: ODError },

    #[snafu(display("Could not get measurement covariance: {source}"))]
    MeasurementCovarError { source: ODError },

    #[snafu(display("State does not contain STM at epoch {epoch:?}"))]
    MissingSTM { epoch: Epoch },

    #[snafu(display("Measurement at epoch {epoch:?} could not be computed"))]
    MeasurementComputationFailed { epoch: Epoch },

    #[snafu(display("Singular matrix encountered: {details}"))]
    SingularMatrix { details: String },

    #[snafu(display("Too few measurements ({count}) to estimate state"))]
    TooFewMeasurements { count: usize },

    #[snafu(display("Maximum iterations ({max_iter}) reached without convergence"))]
    MaxIterationsReached { max_iter: usize },

    #[snafu(display("Levenberg-Marquardt failed: {details}"))]
    LMFailure { details: String },

    #[snafu(display("Invalid measurement value {val} at epoch {epoch:?}"))]
    InvalidMeasurementValue { epoch: Epoch, val: f64 },
}

/// Solver choice for the Batch Least Squares estimator
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BLSSolver {
    /// Standard Normal Equations: (H^T W H) dx = H^T W dy
    NormalEquations,
    /// Levenberg-Marquardt: (H^T W H + lambda * D^T D) dx = H^T W dy
    LevenbergMarquardt,
    // Could add more later, e.g., QR decomposition based
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
    /// Convergence tolerance on the norm of the state correction vector (dx)
    #[builder(default = 1e-3)]
    pub tolerance: f64,
    /// Maximum number of iterations
    #[builder(default = 10)]
    pub max_iterations: usize,
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
    /// Precision for comparing epochs during trajectory generation (optional)
    #[builder(default = 1 * Unit::Microsecond)]
    pub epoch_precision: Duration,
    /// Almanac for environmental models
    pub almanac: Arc<Almanac>,
    // Keep PhantomData for Dynamics type D if needed, although Trk might imply it
    #[builder(default, setter(skip))]
    _dynamics: PhantomData<D>,
}

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
    ) -> Result<BLSSolution<D::StateType>, BLSError> {
        let measurements = &arc.measurements;
        let num_measurements = measurements.len();
        let mut devices = self.devices.clone();

        ensure!(
            num_measurements >= 2,
            TooFewMeasurementsSnafu {
                count: num_measurements
            }
        );

        info!(
            "Starting Batch Least Squares estimation with {num_measurements} measurements",
        );
        info!("Initial guess: {}", initial_guess.orbit());

        let mut current_estimate = initial_guess;
        let mut converged = false;
        let mut last_correction_norm = f64::MAX;
        let mut lambda = self.lm_lambda_init;
        let mut current_rms = f64::MAX;
        let mut iter: usize = 0;

        // --- Iteration Loop ---
        while iter < self.max_iterations {
            iter += 1;
            info!("[{iter}/{}] Current estimate: {}", self.max_iterations, current_estimate.orbit());

            // Re-initialize matrices for this iteration
            // Information Matrix: Lambda = H^T * W * H
            let mut info_matrix = OMatrix::<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>::identity();
            // Normal Matrix: N = H^T * W * dy
            let mut normal_matrix = OVector::<f64, <D::StateType as State>::Size>::zeros();
            // Sum of squares of weighted residuals for RMS calculation and LM cost
            let mut sum_sq_weighted_residuals = 0.0;

            for (msr_idx, (epoch_ref, msr)) in measurements.iter().enumerate() {
                let msr_epoch = *epoch_ref;

                // Propagate reference state from initial_epoch to msr_epoch
                // NOTE: This re-propagates from t0 each time, which is conceptually clean for BLS
                // but potentially slow. A full trajectory propagation per iteration is faster
                // but needs careful state storage/retrieval. Let's stick to this simpler way first.
                // Reset instance state to the start of the iteration's estimate
                let mut state_to_prop = current_estimate.with_stm();
                state_to_prop.reset_stm();

                let mut prop_to_msr = self.prop.with(state_to_prop, self.almanac.clone()).quiet();
                let state_at_msr = prop_to_msr.until_epoch(msr_epoch).expect("TODO:");

                // Get the STM Phi(t_i, t_0) from the propagated state
                let stm = state_at_msr.stm().expect("TODO:");

                // Get the correct tracking device
                let device = match devices.get_mut(&msr.tracker) {
                    Some(d) => d,
                    None => return Err(BLSError::DeviceNotFound { device_name: msr.tracker.clone() }),
                };

                // Compute H_tilde = dH/dX(t_i) at measurement time t_i
                // Create a temporary Traj containing only the state at measurement time,
                // as device.measure might expect it.
                // TODO: Check if device.measure can work with just a single StateType.
                // If not, we need this temporary Traj.

                for msr_type in msr.data.keys().copied() {
                    let mut msr_types = IndexSet::new();
                    msr_types.insert(msr_type);

                    let h_tilde = device
                    .h_tilde::<U1>(msr, &msr_types, &state_at_msr, self.almanac.clone())
                    .map_err(|e| BLSError::SensitivityError { source: e })?;

                    // Compute expected measurement H(X(t_i))
                    let computed_meas_opt = device
                        .measure_instantaneous(state_at_msr, None, self.almanac.clone())
                        .map_err(|e| BLSError::MeasurementError { source: e })?;

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
                        InvalidMeasurementValueSnafu { epoch: msr_epoch, val: real_obs }
                    );

                    // Compute residual dy = y_i - H(X(t_i))
                    let residual = real_obs - computed_obs;

                    // Get measurement variance R (assuming 1x1 matrix) and weight W = 1/R
                    let r_matrix = device
                        .measurement_covar_matrix(&msr_types, msr_epoch)
                        .map_err(|e| BLSError::MeasurementCovarError { source: e })?;
                    let r_variance = r_matrix[(0, 0)];
                    ensure!(r_variance > 0.0, SingularMatrixSnafu {
                        details: format!("Zero measurement variance R for msr {} @ {}", msr_idx, msr_epoch)
                    });
                    let weight = 1.0; // / r_variance; // TODO: R^-1 to be fixed.

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

            // --- Solve for State Correction dx ---
            let orbit_correction: OVector<f64, <D::StateType as State>::Size>;
            let iteration_cost_decreased; // For LM logic

            // Use num_measurements for consistency
            let current_iter_rms = (sum_sq_weighted_residuals / num_measurements as f64).sqrt();

            match self.solver {
                BLSSolver::NormalEquations => {
                    // Solve Lambda * dx = N => dx = Lambda^-1 * N
                    let info_matrix_chol = match info_matrix.cholesky() {
                         Some(chol) => chol,
                         None => return Err(BLSError::SingularMatrix{ details: format!("Information matrix H^TWH is singular in iteration {}", iter+1)})
                    };
                    orbit_correction = info_matrix_chol.solve(&normal_matrix);
                    // Assume NE always decreases cost locally
                    iteration_cost_decreased = true;
                    current_rms = current_iter_rms;
                }
                BLSSolver::LevenbergMarquardt => {
                     // Solve (Lambda + lambda * D^T D) * dx = N
                     // D^T D is a diagonal scaling matrix.
                     // Common choices: D^T D = I or D^T D = diag(Lambda)
                    let mut d_sq = OMatrix::<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>::identity();
                    if self.lm_use_diag_scaling {
                        // Use D^T D = diag(Lambda)
                        for i in 0..6 {
                            d_sq[(i, i)] = info_matrix.diagonal()[i];
                        }
                        // Ensure diagonal elements are positive for stability
                        for i in 0..6 {
                            if d_sq[(i, i)] <= 0.0 {
                                d_sq[(i, i)] = 1e-6; // Set a small positive floor
                                warn!("LM Scaling: Found non-positive diagonal element {} in H^TWH, using floor.", info_matrix[(i,i)]);
                            }
                        }
                    } // else d_sq remains Identity

                    // Inner LM loop to find suitable lambda
                    let augmented_matrix = info_matrix + d_sq * lambda;

                    if let Some(aug_chol) = augmented_matrix.cholesky() {
                        orbit_correction = aug_chol.solve(&normal_matrix);

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
                            debug!("LM: Cost decreased (RMS {} -> {}). Decreasing lambda to {}", current_rms, current_iter_rms, lambda);
                            current_rms = current_iter_rms; // Update RMS baseline
                        } else {
                             // Cost increased or stalled
                             iteration_cost_decreased = false;
                             lambda *= self.lm_lambda_increase; // Increase damping
                             lambda = lambda.min(self.lm_lambda_max); // Clamp to max
                             debug!("LM: Cost increased/stalled (RMS {} -> {}). Increasing lambda to {}", current_rms, current_iter_rms, lambda);
                             // Don't update current_rms baseline if cost increased
                        }

                    } else {
                        // Augmented matrix is singular, increase lambda significantly and retry
                        warn!("LM: Augmented matrix (H^TWH + lambda*D^2) singular with lambda={}. Increasing lambda.", lambda);
                        lambda *= self.lm_lambda_increase * 10.0; // Increase more aggressively
                        lambda = lambda.min(self.lm_lambda_max);
                        // Skip update in this iteration, force retry with larger lambda next time if possible
                        continue; // Skip the rest of the loop and go to next iteration
                    }
                }
            }


            // --- Update State Estimate ---
            // Only update if the step is considered successful (esp. for LM)
            // Also hit if using normal equations because iteration_cost_decreased is forced to true
            if iteration_cost_decreased {
                let mut state_correction = OVector::<f64, <D::StateType as State>::Size>::zeros();
                for i in 0..6 {
                    state_correction[i] = orbit_correction[i];
                }
                current_estimate = current_estimate + state_correction;
                info!(
                    "[{iter}/{}] RMS: {current_iter_rms:.3}; corrections: {:.3} m\t{:.3} m/s",
                    self.max_iterations,
                    orbit_correction.fixed_rows::<3>(0).norm() * 1e3,
                    orbit_correction.fixed_rows::<3>(3).norm() * 1e3
                );
                debug!("Updated estimate: {}", current_estimate.orbit());

                 last_correction_norm = orbit_correction.norm();

                // --- Check Convergence ---
                if last_correction_norm < self.tolerance {
                    info!(
                        "Converged in {iter} iterations. Final correction norm: {:.3e}",
                        last_correction_norm
                    );
                    converged = true;
                    break;
                }
            } else if self.solver == BLSSolver::LevenbergMarquardt {
                 // LM step was rejected (cost increased)
                 info!("LM: Step rejected in iteration {iter}. Retrying with increased lambda.");
                 // Reset correction norm as step was bad
                 last_correction_norm = f64::MAX;
                 // The loop will continue with the increased lambda
            }

        }

        if !converged {
            warn!(
                "Maximum iterations ({}) reached without convergence. Last correction norm: {:.3e}",
                self.max_iterations, last_correction_norm
            );
            // Optionally return error or the last state
            // return Err(BLSError::MaxIterationsReached { max_iter: self.max_iterations });
        }

        // Compute final covariance P = (H^T W H)^-1
        // Need to recompute H^T W H using the *final* estimate
        info!("Computing final covariance matrix");
        let mut final_info_matrix = OMatrix::<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>::identity();

         // Propagate final reference trajectory
         let final_prop = self.prop.clone();
         let mut final_prop_instance = final_prop.with(current_estimate.with_stm(), self.almanac.clone()).quiet();
         final_prop_instance.state.reset_stm(); // Ensure STM starts as identity

        for (epoch_ref, msr) in measurements.iter() {
            let msr_epoch = *epoch_ref;

            // Similar propagation as in the iteration loop, but using the final converged estimate
            let mut state_to_prop = current_estimate.with_stm();
            state_to_prop.reset_stm();

            let mut prop_to_msr = self.prop.with(state_to_prop, self.almanac.clone()).quiet();
            let state_at_msr = prop_to_msr.until_epoch(msr_epoch).expect("TODO:");

            for msr_type in msr.data.keys().copied() {
                let mut msr_types = IndexSet::new();
                msr_types.insert(msr_type);

                let stm = state_at_msr.stm().expect("TODO:");
                let device = devices.get(&msr.tracker).ok_or_else(|| BLSError::DeviceNotFound { device_name: msr.tracker.clone() })?;


                let h_tilde = device
                    .h_tilde::<U1>(msr, &msr_types, &state_at_msr, self.almanac.clone())
                    .map_err(|e| BLSError::SensitivityError { source: e })?;

                let r_matrix = device
                    .measurement_covar_matrix(&msr_types, msr_epoch)
                    .map_err(|e| BLSError::MeasurementCovarError { source: e })?;
                let weight = 1.0 / r_matrix[(0, 0)];
                let h_matrix = h_tilde * stm;
                final_info_matrix.gemm(weight, &h_matrix.transpose(), &h_matrix, 1.0);
            }
        }

        let final_covariance = match final_info_matrix.try_inverse() {
             Some(cov) => cov,
             None => {
                 warn!("Final information matrix H^TWH is singular. Returning identity covariance.");
                 OMatrix::<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>::identity()
             }
        };

        let mut covariance = OMatrix::<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>::zeros();
        for i in 0..6 {
            for j in 0..6 {
                covariance[(i, i)] = final_covariance[(i, j)];
            }
        }

        info!("Batch Least Squares estimation completed.");
        Ok(BLSSolution {
            estimated_state: current_estimate,
            covariance,
            num_iterations: iter,
            final_rms: current_rms,
            final_correction_norm: last_correction_norm,
            converged,
        })
    }
}
