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
use crate::md::trajectory::Interpolatable;
pub use crate::od::estimate::*;
pub use crate::od::*;
use anise::prelude::Almanac;
use log::info;
use msr::sensitivity::TrackerSensitivity;
use std::ops::Add;

use super::ODSolution;

impl<StateType, EstType, MsrSize, Trk> ODSolution<StateType, EstType, MsrSize, Trk>
where
    StateType: Interpolatable + Add<OVector<f64, <StateType as State>::Size>, Output = StateType>,
    EstType: Estimate<StateType>,
    MsrSize: DimName,
    Trk: TrackerSensitivity<StateType, StateType>,
    <DefaultAllocator as Allocator<<StateType as State>::VecLength>>::Buffer<f64>: Send,
    DefaultAllocator: Allocator<<StateType as State>::Size>
        + Allocator<<StateType as State>::VecLength>
        + Allocator<MsrSize>
        + Allocator<MsrSize, <StateType as State>::Size>
        + Allocator<MsrSize, MsrSize>
        + Allocator<<StateType as State>::Size, <StateType as State>::Size>
        + Allocator<<StateType as State>::Size, MsrSize>,
{
    /// Smoothes this OD solution, returning a new OD solution and the filter-smoother consistency ratios, with updated **postfit** residuals, and where the ratio now represents the filter-smoother consistency ratio.
    ///
    /// Notes:
    ///  1. Gains will be scrubbed because the smoother process does not recompute the gain.
    ///  2. Prefit residuals, ratios, and measurement covariances are not updated, as these depend on the filtering process.
    ///  3. Note: this function consumes the current OD solution to prevent reusing the wrong one.
    ///
    ///
    /// To assess whether the smoothing process improved the solution, compare the RMS of the postfit residuals from the filter and the smoother process.
    ///
    /// # Filter-Smoother consistency ratio
    ///
    /// The **filter-smoother consistency ratio** is used to evaluate the consistency between the state estimates produced by a filter (e.g., Kalman filter) and a smoother.
    /// This ratio is called "filter smoother consistency test" in the ODTK MathSpec.
    ///
    /// It is computed as follows:
    ///
    /// #### 1. Define the State Estimates
    /// **Filter state estimate**:
    /// $ \hat{X}_{f,k} $
    /// This is the state estimate at time step $ k $ from the filter.
    ///
    /// **Smoother state estimate**:
    /// $ \hat{X}_{s,k} $
    /// This is the state estimate at time step $ k $ from the smoother.
    ///
    /// #### 2. Define the Covariances
    ///
    /// **Filter covariance**:
    /// $ P_{f,k} $
    /// This is the covariance of the state estimate at time step $ k $ from the filter.
    ///
    /// **Smoother covariance**:
    /// $ P_{s,k} $
    /// This is the covariance of the state estimate at time step $ k $ from the smoother.
    ///
    /// #### 3. Compute the Differences
    ///
    /// **State difference**:
    /// $ \Delta X_k = \hat{X}_{s,k} - \hat{X}_{f,k} $
    ///
    /// **Covariance difference**:
    /// $ \Delta P_k = P_{s,k} - P_{f,k} $
    ///
    /// #### 4. Calculate the Consistency Ratio
    /// For each element $ i $ of the state vector, compute the ratio:
    ///
    /// $$
    /// R_{i,k} = \frac{\Delta X_{i,k}}{\sqrt{\Delta P_{i,k}}}
    /// $$
    ///
    /// #### 5. Evaluate Consistency
    /// - If $ |R_{i,k}| \leq 3 $ for all $ i $ and $ k $, the filter-smoother consistency test is satisfied, indicating good consistency.
    /// - If $ |R_{i,k}| > 3 $ for any $ i $ or $ k $, the test fails, suggesting potential modeling inconsistencies or issues with the estimation process.
    ///
    pub fn smooth(self, almanac: Arc<Almanac>) -> Result<Self, ODError> {
        let l = self.estimates.len() - 1;

        let mut smoothed = Self {
            estimates: Vec::with_capacity(self.estimates.len()),
            residuals: Vec::with_capacity(self.residuals.len()),
            gains: Vec::with_capacity(self.estimates.len()),
            filter_smoother_ratios: Vec::with_capacity(self.estimates.len()),
            devices: self.devices.clone(),
            measurement_types: self.measurement_types.clone(),
        };

        // Set the first item of the smoothed estimates to the last estimate (we cannot smooth the very last estimate)
        smoothed
            .estimates
            .push(self.estimates.last().unwrap().clone());
        smoothed
            .residuals
            .push(self.residuals.last().unwrap().clone());
        smoothed.gains.push(None);
        smoothed.filter_smoother_ratios.push(None);

        loop {
            let k = l - smoothed.estimates.len();
            // Borrow the previously smoothed estimate of the k+1 estimate
            let sm_est_kp1 = &self.estimates[k + 1];
            let x_kp1_l = sm_est_kp1.state_deviation();
            let p_kp1_l = sm_est_kp1.covar();
            // Borrow the k-th estimate, which we're smoothing with the next estimate
            let est_k = &self.estimates[k];
            // Borrow the k+1-th estimate, which we're smoothing with the next estimate
            let est_kp1 = &self.estimates[k + 1];

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
            // Recompute the residual if available.
            if let Some(mut residual) = self.residuals[k + 1].clone() {
                let tracker = residual
                    .tracker
                    .as_ref()
                    .expect("tracker unset in smoother process");

                let device = smoothed
                    .devices
                    .get_mut(tracker)
                    .expect("unknown tracker in smoother process");

                let new_state_est = smoothed_est_k.state();
                let epoch = new_state_est.epoch();

                if let Some(computed_meas) =
                    device.measure_instantaneous(new_state_est, None, almanac.clone())?
                {
                    // Only recompute the computed observation from the update state estimate.
                    residual.computed_obs = computed_meas
                        .observation::<MsrSize>(&residual.msr_types)
                        - device.measurement_bias_vector::<MsrSize>(&residual.msr_types, epoch)?;

                    // Update the postfit residual.
                    residual.postfit = &residual.real_obs - &residual.computed_obs;

                    // Store the updated data.
                    smoothed.residuals.push(Some(residual));
                } else {
                    smoothed.residuals.push(None);
                }
            } else {
                smoothed.residuals.push(None);
            }

            // Compute the filter-smoother consistency ratio.
            let delta_covar = est_k.covar() - smoothed_est_k.covar();
            let delta_state =
                est_k.state().to_state_vector() - smoothed_est_k.state().to_state_vector();

            let fs_ratios = OVector::<f64, <StateType as State>::Size>::from_iterator(
                delta_state
                    .iter()
                    .enumerate()
                    .map(|(i, dx)| dx / delta_covar[(i, i)].sqrt()),
            );

            smoothed.estimates.push(smoothed_est_k);
            smoothed.filter_smoother_ratios.push(Some(fs_ratios));
            // Set all gains to None.
            smoothed.gains.push(None);

            if smoothed.estimates.len() == self.estimates.len() {
                break;
            }
        }

        // Note that we have yet to reverse the list, so we print them backward
        info!(
            "Smoothed {} estimates (from {} to {})",
            smoothed.estimates.len(),
            smoothed.estimates.last().unwrap().epoch(),
            smoothed.estimates[0].epoch(),
        );

        // Now, let's add all of the other estimates so that the same indexing can be done
        // between all the estimates and the smoothed estimates
        if smoothed.estimates.len() < self.estimates.len() {
            // Add the estimates that might have been skipped.
            let mut k = self.estimates.len() - smoothed.estimates.len();
            loop {
                smoothed.estimates.push(self.estimates[k].clone());
                if k == 0 {
                    break;
                }
                k -= 1;
            }
        }

        // And reverse to maintain the order of estimates
        smoothed.estimates.reverse();
        smoothed.residuals.reverse();
        smoothed.filter_smoother_ratios.reverse();

        Ok(smoothed)
    }
}
