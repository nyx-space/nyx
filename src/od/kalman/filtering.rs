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

pub use crate::errors::NyxError;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector};
pub use crate::od::estimate::{Estimate, KfEstimate, Residual};
use crate::od::prelude::KalmanVariant;
use crate::od::process::ResidRejectCrit;
pub use crate::od::snc::ProcessNoise;
use crate::od::{ODDynamicsSnafu, ODError, State};
pub use crate::time::{Epoch, Unit};
use log::info;
use snafu::prelude::*;

use super::KalmanFilter;

impl<T, A> KalmanFilter<T, A>
where
    A: DimName,
    T: State,
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<<T as State>::VecLength>
        + Allocator<A>
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<A, A>
        + Allocator<<T as State>::Size, A>
        + Allocator<A, <T as State>::Size>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    /// Returns the previous estimate
    pub fn previous_estimate(&self) -> &KfEstimate<T> {
        &self.prev_estimate
    }

    pub fn set_previous_estimate(&mut self, est: &KfEstimate<T>) {
        self.prev_estimate = *est;
    }

    /// Computes a time update/prediction (i.e. advances the filter estimate with the updated STM).
    ///
    /// May return a FilterError if the STM was not updated.
    pub fn time_update(&mut self, nominal_state: T) -> Result<KfEstimate<T>, ODError> {
        let stm = nominal_state.stm().context(ODDynamicsSnafu)?;
        let mut covar_bar = stm * self.prev_estimate.covar * stm.transpose();

        // Apply any process noise as in a normal time update, if applicable
        for (i, snc) in self.process_noise.iter().enumerate().rev() {
            if let Some(snc_contrib) = snc.propagate::<<T as State>::Size>(
                nominal_state.orbit(),
                nominal_state.epoch() - self.prev_estimate.epoch(),
            )? {
                if self.prev_used_snc != i {
                    info!("Switched to {i}-th {snc}");
                    self.prev_used_snc = i;
                }
                // Let's add the process noise
                covar_bar += snc_contrib;
                // And break so we don't add any more process noise
                break;
            }
        }

        let state_bar = if matches!(self.variant, KalmanVariant::DeviationTracking) {
            stm * self.prev_estimate.state_deviation
        } else {
            OVector::<f64, <T as State>::Size>::zeros()
        };
        let estimate = KfEstimate {
            nominal_state,
            state_deviation: state_bar,
            covar: covar_bar,
            covar_bar,
            stm,
            predicted: true,
        };
        self.prev_estimate = estimate;
        // Update the prev epoch for all SNCs
        for snc in &mut self.process_noise {
            snc.prev_epoch = Some(self.prev_estimate.epoch());
        }
        Ok(estimate)
    }

    /// Computes the measurement update with a provided real observation and computed observation.
    ///
    /// May return a FilterError if the STM or sensitivity matrices were not updated.
    pub fn measurement_update<M: DimName>(
        &mut self,
        nominal_state: T,
        real_obs: OVector<f64, M>,
        computed_obs: OVector<f64, M>,
        r_k: OMatrix<f64, M, M>,
        h_tilde: OMatrix<f64, M, <T as State>::Size>,
        resid_rejection: Option<ResidRejectCrit>,
    ) -> Result<
        (
            KfEstimate<T>,
            Residual<M>,
            Option<OMatrix<f64, <T as State>::Size, M>>,
        ),
        ODError,
    >
    where
        DefaultAllocator: Allocator<M>
            + Allocator<M, M>
            + Allocator<M, <T as State>::Size>
            + Allocator<<T as State>::Size, M>
            + Allocator<nalgebra::Const<1>, M>,
    {
        let epoch = nominal_state.epoch();

        // Grab the state transition matrix.
        let stm = nominal_state.stm().context(ODDynamicsSnafu)?;

        // Propagate the covariance.
        let mut covar_bar = stm * self.prev_estimate.covar * stm.transpose();

        // Apply any process noise as in a normal time update, if applicable
        for (i, snc) in self.process_noise.iter().enumerate().rev() {
            if let Some(snc_contrib) = snc.propagate::<<T as State>::Size>(
                nominal_state.orbit(),
                nominal_state.epoch() - self.prev_estimate.epoch(),
            )? {
                if self.prev_used_snc != i {
                    info!("Switched to {i}-th {snc}");
                    self.prev_used_snc = i;
                }
                // Let's add the process noise
                covar_bar += snc_contrib;
                // And break so we don't add any more process noise
                break;
            }
        }

        // Project the propagated covariance into the measurement space.
        let p_ht = covar_bar * h_tilde.transpose();
        let h_p_ht = &h_tilde * &p_ht;

        // Compute the innovation matrix (S_k).
        let s_k = &h_p_ht + &r_k;

        // Compute observation deviation/error (usually marked as y_i)
        let prefit = real_obs.clone() - computed_obs.clone();

        // Compute the prefit ratio for the automatic rejection.
        // The measurement covariance is the square of the measurement itself.
        // So we compute its Cholesky decomposition to return to the non squared values.
        let s_k_chol = match s_k.clone().cholesky() {
            Some(r_k_clone) => r_k_clone,
            None => {
                // In very rare case, when there isn't enough noise in the measurements,
                // the inverting of S_k fails. If so, we revert back to the nominal Kalman derivation.
                r_k.clone().cholesky().ok_or(ODError::SingularNoiseRk)?
            }
        };

        // Get the L factor from the Cholesky decomposition
        let l_matrix = s_k_chol.l();

        // Solve L * v = prefit for the whitened residual vector v. This is an O(n^2) triangular solve, faster than a full Cholesky solve.
        let whitened_resid = l_matrix.solve_lower_triangular(&prefit).unwrap();

        // Compute the RMS ratio using the norm of the whitened vector This is the true Mahalanobis-based N-sigma ratio.
        let ratio = (whitened_resid.norm_squared() / (M::DIM as f64)).sqrt();

        // Compute the physical 1-sigma envelop. Using the diagonal of S_k (not L) is correct for physical innovation plots.
        // let innovation_trend = l_matrix.diagonal();
        let innovation_trend = s_k.diagonal().map(|x| x.sqrt());

        if let Some(resid_reject) = resid_rejection {
            if ratio > resid_reject.num_sigmas {
                // Reject this whole measurement and perform only a time update
                let pred_est = self.time_update(nominal_state)?;
                let resid = Residual::rejected(
                    epoch,
                    prefit,
                    ratio,
                    innovation_trend,
                    real_obs,
                    computed_obs,
                );

                return Ok((pred_est, resid, None));
            }
        }

        // Instead of inverting the innovation matrix S_k, we will use the (super short) arXiv paper 1111.4144
        // which shows how to use the Cholesky decomposition to invert a matrix, core tenets repeated here for my reference.
        // \forall A \ in \mathbb{R}^{n\times n}, X=A^{-1} <=> A*X=I
        // Cholesky: A = L*L^T
        // Therefore, L*L^T*X = I
        // 1. Solve L * Y = I  => Y = L^{-1} (via forward sub)
        // 2. Solve L^T * X = Y => X = (L^T)^{-1} * L^{-1} = A^{-1}
        //
        // _However_, we can be more clever still!
        // Instead of explicitly inverting the innovation matrix S_k, we solve the linear system
        // S_k * K^T = H * P using Cholesky decomposition.
        // This avoids the numerical instability of computing S_k^-1 directly.
        // Math context:
        // We want to solve A * X = B for X.
        // 1. Decompose A into L * L^T (Cholesky).
        // 2. Solve L * Y = B for Y (Forward substitution).
        // 3. Solve L^T * X = Y for X (Backward substitution).

        // Prepare the RHS of the linear system: (P * H^T)^T = H * P
        // We want to solve: S_k * K^T = H * P
        // So K = (S_k \ (H * P))^T
        let rhs = p_ht.transpose();

        // Solve for Gain using Cholesky
        // We try standard Cholesky first.
        let gain = match s_k.clone().cholesky() {
            Some(chol) => {
                // SOLVE, don't invert.
                // chol.solve(B) computes S_k^{-1} * B more stably than inv(S_k) * B
                let k_t = chol.solve(&rhs);
                k_t.transpose()
            }
            None => {
                // If this fails, revert the LU decomposition of nalgebra
                // Invert the innovation covariance.
                match s_k.try_inverse() {
                    Some(s_k_inv) => covar_bar * &h_tilde.transpose() * &s_k_inv,
                    None => {
                        eprintln!(
                            "SINGULAR GAIN\nr = {r_k}\nh = {h_tilde:.3e}\ncovar = {covar_bar:.3e}"
                        );
                        return Err(ODError::SingularKalmanGain);
                    }
                }
            }
        };

        // Compute the state estimate, depends on the variant.
        let (state_hat, res) = match self.variant {
            KalmanVariant::ReferenceUpdate => {
                // In EKF, the state hat is actually the state deviation. We trust the gain to be correct,
                // so we just apply it directly to the prefit residual.
                let state_hat = &gain * &prefit;
                let postfit = &prefit - (&h_tilde * state_hat);
                let resid = Residual::accepted(
                    epoch,
                    prefit,
                    postfit,
                    ratio,
                    innovation_trend,
                    real_obs,
                    computed_obs,
                );
                (state_hat, resid)
            }
            KalmanVariant::DeviationTracking => {
                // Time update
                let state_bar = stm * self.prev_estimate.state_deviation;
                let postfit = &prefit - (&h_tilde * state_bar);
                (
                    state_bar + &gain * &postfit,
                    Residual::accepted(
                        epoch,
                        prefit,
                        postfit,
                        ratio,
                        innovation_trend,
                        real_obs,
                        computed_obs,
                    ),
                )
            }
        };

        // Compute covariance (Joseph update)
        let first_term =
            OMatrix::<f64, <T as State>::Size, <T as State>::Size>::identity() - &gain * &h_tilde;
        let covar =
            first_term * covar_bar * first_term.transpose() + &gain * &r_k * &gain.transpose();

        // Force symmetry on the covariance
        let covar = 0.5 * (covar + covar.transpose());

        // And wrap up
        let estimate = KfEstimate {
            nominal_state,
            state_deviation: state_hat,
            covar,
            covar_bar,
            stm,
            predicted: false,
        };

        self.prev_estimate = estimate;
        // Update the prev epoch for all SNCs
        for snc in &mut self.process_noise {
            snc.prev_epoch = Some(self.prev_estimate.epoch());
        }

        Ok((estimate, res, Some(gain)))
    }

    pub fn replace_state(&self) -> bool {
        matches!(self.variant, KalmanVariant::ReferenceUpdate)
    }

    /// Overwrites all of the process noises to the one provided
    pub fn set_process_noise(&mut self, snc: ProcessNoise<A>) {
        self.process_noise = vec![snc];
    }
}
