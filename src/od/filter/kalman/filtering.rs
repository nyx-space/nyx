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
use crate::od::process::ResidRejectCrit;
pub use crate::od::snc::SNC;
use crate::od::{Filter, ODDynamicsSnafu, ODError, State};
pub use crate::time::{Epoch, Unit};
use snafu::prelude::*;

use super::KF;

impl<T, A, M> Filter<T, A, M> for KF<T, A, M>
where
    A: DimName,
    M: DimName,
    T: State,
    DefaultAllocator: Allocator<M>
        + Allocator<<T as State>::Size>
        + Allocator<<T as State>::VecLength>
        + Allocator<A>
        + Allocator<M, M>
        + Allocator<M, <T as State>::Size>
        + Allocator<<T as State>::Size, M>
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<A, A>
        + Allocator<<T as State>::Size, A>
        + Allocator<A, <T as State>::Size>
        + Allocator<nalgebra::Const<1>, M>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    type Estimate = KfEstimate<T>;

    /// Returns the previous estimate
    fn previous_estimate(&self) -> &Self::Estimate {
        &self.prev_estimate
    }

    fn set_previous_estimate(&mut self, est: &Self::Estimate) {
        self.prev_estimate = *est;
    }

    /// Update the sensitivity matrix (or "H tilde"). This function **must** be called prior to each
    /// call to `measurement_update`.
    fn update_h_tilde(&mut self, h_tilde: OMatrix<f64, M, <T as State>::Size>) {
        self.h_tilde = h_tilde;
        self.h_tilde_updated = true;
    }

    /// Computes a time update/prediction (i.e. advances the filter estimate with the updated STM).
    ///
    /// May return a FilterError if the STM was not updated.
    fn time_update(&mut self, nominal_state: T) -> Result<Self::Estimate, ODError> {
        let stm = nominal_state.stm().context(ODDynamicsSnafu)?;
        let mut covar_bar = stm * self.prev_estimate.covar * stm.transpose();

        // Try to apply an SNC, if applicable
        for (i, snc) in self.process_noise.iter().enumerate().rev() {
            if let Some(snc_matrix) = snc.to_matrix(nominal_state.epoch()) {
                // Check if we're using another SNC than the one before
                if self.prev_used_snc != i {
                    info!("Switched to {}-th {}", i, snc);
                    self.prev_used_snc = i;
                }

                // Let's compute the Gamma matrix, an approximation of the time integral
                // which assumes that the acceleration is constant between these two measurements.
                let mut gamma = OMatrix::<f64, <T as State>::Size, A>::zeros();
                let delta_t = (nominal_state.epoch() - self.prev_estimate.epoch()).to_seconds();
                for blk in 0..A::dim() / 3 {
                    for i in 0..3 {
                        let idx_i = i + A::dim() * blk;
                        let idx_j = i + 3 * blk;
                        let idx_k = i + 3 + A::dim() * blk;
                        // For first block
                        // (0, 0) (1, 1) (2, 2) <=> \Delta t^2/2
                        // (3, 0) (4, 1) (5, 2) <=> \Delta t
                        // Second block
                        // (6, 3) (7, 4) (8, 5) <=> \Delta t^2/2
                        // (9, 3) (10, 4) (11, 5) <=> \Delta t
                        // * \Delta t^2/2
                        // (i, i) when blk = 0
                        // (i + A::dim() * blk, i + 3) when blk = 1
                        // (i + A::dim() * blk, i + 3 * blk)
                        // * \Delta t
                        // (i + 3, i) when blk = 0
                        // (i + 3, i + 9) when blk = 1 (and I think i + 12 + 3)
                        // (i + 3 + A::dim() * blk, i + 3 * blk)
                        gamma[(idx_i, idx_j)] = delta_t.powi(2) / 2.0;
                        gamma[(idx_k, idx_j)] = delta_t;
                    }
                }
                // Let's add the process noise
                covar_bar += &gamma * snc_matrix * &gamma.transpose();
                // And break so we don't add any more process noise
                break;
            }
        }

        let state_bar = if self.ekf {
            OVector::<f64, <T as State>::Size>::zeros()
        } else {
            stm * self.prev_estimate.state_deviation
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
    fn measurement_update(
        &mut self,
        nominal_state: T,
        real_obs: OVector<f64, M>,
        computed_obs: OVector<f64, M>,
        r_k: OMatrix<f64, M, M>,
        resid_rejection: Option<ResidRejectCrit>,
    ) -> Result<(Self::Estimate, Residual<M>), ODError> {
        if !self.h_tilde_updated {
            return Err(ODError::SensitivityNotUpdated);
        }

        let epoch = nominal_state.epoch();

        // Grab the state transition matrix.
        let stm = nominal_state.stm().context(ODDynamicsSnafu)?;

        // Propagate the covariance.
        let covar_bar = stm * self.prev_estimate.covar * stm.transpose();

        // Project the propagated covariance into the measurement space.
        let h_p_ht = &self.h_tilde * covar_bar * &self.h_tilde.transpose();

        // Compute the innovation matrix (S_k).
        let s_k = &h_p_ht + &r_k;

        // Compute observation deviation/error (usually marked as y_i)
        let prefit = real_obs.clone() - computed_obs.clone();

        // Compute the prefit ratio for the automatic rejection.
        // The measurement covariance is the square of the measurement itself.
        // So we compute its Cholesky decomposition to return to the non squared values.
        let r_k_chol = match s_k.clone().cholesky() {
            Some(r_k_clone) => r_k_clone.l(),
            None => {
                // In very rare case, when there isn't enough noise in the measurements,
                // the inverting of S_k fails. If so, we revert back to the nominal Kalman derivation.
                r_k.clone().cholesky().ok_or(ODError::SingularNoiseRk)?.l()
            }
        };

        // Compute the ratio as the average of each component of the prefit over the square root of the measurement
        // matrix r_k. Refer to ODTK MathSpec equation 4.10.
        let ratio = s_k
            .diagonal()
            .iter()
            .copied()
            .enumerate()
            .map(|(idx, r)| prefit[idx] / r.sqrt())
            .sum::<f64>()
            / (M::USIZE as f64);

        if let Some(resid_reject) = resid_rejection {
            if ratio.abs() > resid_reject.num_sigmas {
                // Reject this whole measurement and perform only a time update
                let pred_est = self.time_update(nominal_state)?;
                return Ok((
                    pred_est,
                    Residual::rejected(
                        epoch,
                        prefit,
                        ratio,
                        r_k_chol.diagonal(),
                        real_obs,
                        computed_obs,
                    ),
                ));
            }
        }

        // Invert the innovation covariance.
        let s_k_inv = match s_k.try_inverse() {
            Some(s_k_inv) => s_k_inv,
            None => return Err(ODError::SingularKalmanGain),
        };

        let gain = covar_bar * &self.h_tilde.transpose() * &s_k_inv;

        // Compute the state estimate
        let (state_hat, res) = if self.ekf {
            // In EKF, the state hat is actually the state deviation. We trust the gain to be correct,
            // so we just apply it directly to the prefit residual.
            let state_hat = &gain * &prefit;
            let postfit = &prefit - (&self.h_tilde * state_hat);
            (
                state_hat,
                Residual::accepted(
                    epoch,
                    prefit,
                    postfit,
                    ratio,
                    r_k_chol.diagonal(),
                    real_obs,
                    computed_obs,
                ),
            )
        } else {
            // Time update
            let state_bar = stm * self.prev_estimate.state_deviation;
            let postfit = &prefit - (&self.h_tilde * state_bar);
            (
                state_bar + &gain * &postfit,
                Residual::accepted(
                    epoch,
                    prefit,
                    postfit,
                    ratio,
                    r_k_chol.diagonal(),
                    real_obs,
                    computed_obs,
                ),
            )
        };

        // Compute covariance (Joseph update)
        let first_term = OMatrix::<f64, <T as State>::Size, <T as State>::Size>::identity()
            - &gain * &self.h_tilde;
        let covar =
            first_term * covar_bar * first_term.transpose() + &gain * &r_k * &gain.transpose();

        // And wrap up
        let estimate = KfEstimate {
            nominal_state,
            state_deviation: state_hat,
            covar,
            covar_bar,
            stm,
            predicted: false,
        };

        self.h_tilde_updated = false;
        self.prev_estimate = estimate;
        // Update the prev epoch for all SNCs
        for snc in &mut self.process_noise {
            snc.prev_epoch = Some(self.prev_estimate.epoch());
        }
        Ok((estimate, res))
    }

    fn is_extended(&self) -> bool {
        self.ekf
    }

    fn set_extended(&mut self, status: bool) {
        self.ekf = status;
    }

    /// Overwrites all of the process noises to the one provided
    fn set_process_noise(&mut self, snc: SNC<A>) {
        self.process_noise = vec![snc];
    }
}
