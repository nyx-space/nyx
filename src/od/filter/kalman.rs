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
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector, U3};
pub use crate::od::estimate::{Estimate, KfEstimate, Residual};
use crate::od::process::ResidRejectCrit;
pub use crate::od::snc::SNC;
use crate::od::{Filter, ODDynamicsSnafu, ODError, State};
pub use crate::time::{Epoch, Unit};
use snafu::prelude::*;

/// Defines both a Classical and an Extended Kalman filter (CKF and EKF)
/// T: Type of state
/// A: Acceleration size (for SNC)
/// M: Measurement size (used for the sensitivity matrix)
#[derive(Debug, Clone)]
#[allow(clippy::upper_case_acronyms)]
pub struct KF<T, A, M>
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
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<A, A>
        + Allocator<<T as State>::Size, A>
        + Allocator<A, <T as State>::Size>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    /// The previous estimate used in the KF computations.
    pub prev_estimate: KfEstimate<T>,
    /// A sets of process noise (usually noted Q), must be ordered chronologically
    pub process_noise: Vec<SNC<A>>,
    /// Determines whether this KF should operate as a Conventional/Classical Kalman filter or an Extended Kalman Filter.
    /// Recall that one should switch to an Extended KF only once the estimate is good (i.e. after a few good measurement updates on a CKF).
    pub ekf: bool,
    h_tilde: OMatrix<f64, M, <T as State>::Size>,
    h_tilde_updated: bool,
    prev_used_snc: usize,
}

impl<T, A, M> KF<T, A, M>
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
        + Allocator<A, <T as State>::Size>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    /// Initializes this KF with an initial estimate, measurement noise, and one process noise
    pub fn new(initial_estimate: KfEstimate<T>, process_noise: SNC<A>) -> Self {
        assert_eq!(
            A::dim() % 3,
            0,
            "SNC can only be applied to accelerations multiple of 3"
        );

        // Set the initial epoch of the SNC
        let mut process_noise = process_noise;
        process_noise.init_epoch = Some(initial_estimate.epoch());

        Self {
            prev_estimate: initial_estimate,
            process_noise: vec![process_noise],
            ekf: false,
            h_tilde: OMatrix::<f64, M, <T as State>::Size>::zeros(),
            h_tilde_updated: false,
            prev_used_snc: 0,
        }
    }

    /// Initializes this KF with an initial estimate, measurement noise, and several process noise
    /// WARNING: SNCs MUST be ordered chronologically! They will be selected automatically by walking
    /// the list of SNCs backward until one can be applied!
    pub fn with_sncs(initial_estimate: KfEstimate<T>, process_noises: Vec<SNC<A>>) -> Self {
        assert_eq!(
            A::dim() % 3,
            0,
            "SNC can only be applied to accelerations multiple of 3"
        );
        let mut process_noises = process_noises;
        // Set the initial epoch of the SNC
        for snc in &mut process_noises {
            snc.init_epoch = Some(initial_estimate.epoch());
        }

        Self {
            prev_estimate: initial_estimate,
            process_noise: process_noises,
            ekf: false,
            h_tilde: OMatrix::<f64, M, <T as State>::Size>::zeros(),
            h_tilde_updated: false,
            prev_used_snc: 0,
        }
    }
}

impl<T, M> KF<T, U3, M>
where
    M: DimName,
    T: State,
    DefaultAllocator: Allocator<M>
        + Allocator<<T as State>::Size>
        + Allocator<<T as State>::VecLength>
        + Allocator<M, M>
        + Allocator<M, <T as State>::Size>
        + Allocator<<T as State>::Size, M>
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<U3, U3>
        + Allocator<<T as State>::Size, U3>
        + Allocator<U3, <T as State>::Size>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    /// Initializes this KF without SNC
    pub fn no_snc(initial_estimate: KfEstimate<T>) -> Self {
        Self {
            prev_estimate: initial_estimate,
            process_noise: Vec::new(),
            ekf: false,
            h_tilde: OMatrix::<f64, M, <T as State>::Size>::zeros(),
            h_tilde_updated: false,
            prev_used_snc: 0,
        }
    }
}

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
        real_obs: &OVector<f64, M>,
        computed_obs: &OVector<f64, M>,
        measurement_covar: OMatrix<f64, M, M>,
        resid_rejection: Option<ResidRejectCrit>,
    ) -> Result<(Self::Estimate, Residual<M>), ODError> {
        if !self.h_tilde_updated {
            return Err(ODError::SensitivityNotUpdated);
        }

        let stm = nominal_state.stm().context(ODDynamicsSnafu)?;

        let epoch = nominal_state.epoch();

        let covar_bar = stm * self.prev_estimate.covar * stm.transpose();

        let h_tilde_t = &self.h_tilde.transpose();
        let h_p_ht = &self.h_tilde * covar_bar * h_tilde_t;
        // Account for state uncertainty in the measurement noise. Equation 4.10 of ODTK MathSpec.
        let r_k = &h_p_ht + measurement_covar;

        // Compute observation deviation (usually marked as y_i)
        let prefit = real_obs - computed_obs;

        // Compute the prefit ratio for the automatic rejection.
        // The measurement covariance is the square of the measurement itself.
        // So we compute its Cholesky decomposition to return to the non squared values.
        let r_k_chol_inv = r_k
            .clone()
            .cholesky()
            .ok_or(ODError::SingularNoiseRk)?
            .l()
            .try_inverse()
            .ok_or(ODError::SingularNoiseRk)?;

        // Then we multiply the inverted r_k matrix with prefit, giving us the vector of
        // ratios of the prefit over the r_k matrix diagonal.
        let ratio_vec = r_k_chol_inv * &prefit;
        // Compute the ratio as the dot product.
        let ratio = ratio_vec.dot(&ratio_vec);

        if let Some(resid_reject) = resid_rejection {
            if ratio > resid_reject.num_sigmas {
                // Reject this whole measurement and perform only a time update
                let pred_est = self.time_update(nominal_state)?;
                return Ok((
                    pred_est,
                    Residual::rejected(epoch, prefit, ratio, r_k.diagonal()),
                ));
            }
        }

        // Compute the Kalman gain but first adding the measurement noise to H⋅P⋅H^T
        let mut innovation_covar = h_p_ht + &r_k;
        if !innovation_covar.try_inverse_mut() {
            return Err(ODError::SingularKalmanGain);
        }

        let gain = covar_bar * h_tilde_t * &innovation_covar;

        // Compute the state estimate
        let (state_hat, res) = if self.ekf {
            let state_hat = &gain * &prefit;
            let postfit = &prefit - (&self.h_tilde * state_hat);
            (
                state_hat,
                Residual::accepted(epoch, prefit, postfit, ratio, r_k.diagonal()),
            )
        } else {
            // Must do a time update first
            let state_bar = stm * self.prev_estimate.state_deviation;
            let postfit = &prefit - (&self.h_tilde * state_bar);
            (
                state_bar + &gain * &postfit,
                Residual::accepted(epoch, prefit, postfit, ratio, r_k.diagonal()),
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
