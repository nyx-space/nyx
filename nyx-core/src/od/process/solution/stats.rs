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
use crate::linalg::{Const, DMatrix, DefaultAllocator, DimMin, DimName, DimSub, OVector};
use crate::md::trajectory::{Interpolatable, Traj};
pub use crate::od::estimate::*;
pub use crate::od::*;
use log::{info, warn};
use msr::sensitivity::TrackerSensitivity;
use snafu::ResultExt;
use statrs::distribution::{ContinuousCDF, Normal};
use std::fmt;
use std::ops::Add;

use super::ODSolution;

#[cfg(feature = "python")]
use pyo3::prelude::*;

#[cfg_attr(feature = "python", pyclass(from_py_object, get_all))]
#[derive(Copy, Clone, Debug)]
pub struct NormalizedConsistency {
    pub normalized_sum: f64,
    pub k: f64,
    pub lower_bound: f64,
    pub upper_bound: f64,
    pub is_nees: bool,
}

#[cfg_attr(feature = "python", pymethods)]
impl NormalizedConsistency {
    pub fn name(&self) -> &'static str {
        if self.is_nees {
            "NEES"
        } else {
            "NIS"
        }
    }

    /// Returns true if there are more than 35 degrees of freedom
    pub fn has_statistical_power(&self) -> bool {
        self.k > 35.0
    }

    /// Returns true if the sum metric is within the lower and upper bounds
    pub fn is_consistent(&self) -> bool {
        self.normalized_sum < self.upper_bound && self.normalized_sum > self.lower_bound
    }

    /// Returns true if the normalized sum is less than the lower bound
    pub fn is_underconfident(&self) -> bool {
        self.normalized_sum < self.lower_bound
    }

    /// Returns true if the normalized sum is greater than the upper bound
    pub fn is_overconfident(&self) -> bool {
        self.normalized_sum > self.upper_bound
    }

    /// Log the status to the logger
    pub fn log(&self) {
        let name = self.name();
        let sum = self.normalized_sum;
        let upper_bound = self.upper_bound;
        let lower_bound = self.lower_bound;
        if sum > upper_bound {
            warn!(
                "{name} consistency test failed high: {name} sum {sum:.3} > upper bound {upper_bound:.3}."
            );
            if self.is_nees {
                warn!("Estimation errors are larger than the covariance suggests.");
            } else {
                warn!("Innovations are larger than expected");
            }
            warn!(
                "Filter is overconfident: P, R, or Q may be too small, or the dynamics/measurement model may be biased."
            );
        } else if sum < lower_bound {
            warn!(
                "{name} consistency test failed low: {name} sum {sum:.3} < lower bound {lower_bound:.3}."
            );
            if self.is_nees {
                warn!("Estimation errors are smaller than expected.");
            } else {
                warn!("Innovations are smaller than expected.")
            }
            warn!("Filter is underconfident: P, R, or Q may be too large.");
        } else {
            info!("{name} passed")
        }
    }

    #[cfg(feature = "python")]
    fn __str__(&self) -> String {
        format!("{self}")
    }

    #[cfg(feature = "python")]
    fn __repr__(&self) -> String {
        format!("{self} @ {self:p}")
    }
}

impl fmt::Display for NormalizedConsistency {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let is_ok = self.is_consistent();
        write!(
            f,
            "{} consistency {} (k={}; bounds: {:.3} < {:.3} < {:.3})",
            self.name(),
            if is_ok { "PASSED" } else { "FAILED" },
            self.k,
            self.lower_bound,
            self.normalized_sum,
            self.upper_bound
        )
    }
}

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
    /// Returns the root mean square of the prefit residuals
    pub fn rms_prefit_residuals(&self) -> f64 {
        let mut sum = 0.0;
        for residual in self.residuals.iter().flatten() {
            sum += residual.prefit.dot(&residual.prefit);
        }
        (sum / (self.residuals.len() as f64)).sqrt()
    }

    /// Returns the root mean square of the postfit residuals
    pub fn rms_postfit_residuals(&self) -> f64 {
        let mut sum = 0.0;
        for residual in self.residuals.iter().flatten() {
            sum += residual.postfit.dot(&residual.postfit);
        }
        (sum / (self.residuals.len() as f64)).sqrt()
    }

    /// Returns the root mean square of the prefit residual ratios
    pub fn rms_residual_ratios(&self) -> f64 {
        let mut sum = 0.0;
        for residual in self.residuals.iter().flatten() {
            sum += residual.ratio.powi(2);
        }
        (sum / (self.residuals.len() as f64)).sqrt()
    }

    /// Computes the fraction of residual ratios that lie within ±threshold.
    pub fn residual_ratio_within_threshold(&self, threshold: f64) -> Result<f64, ODError> {
        let mut total = 0;
        let mut count_within = 0;
        for residual in self.residuals.iter().flatten() {
            total += 1;
            if residual.ratio.abs() <= threshold {
                count_within += 1;
            }
        }
        if total > 0 {
            Ok(count_within as f64 / total as f64)
        } else {
            Err(ODError::ODNoResiduals {
                action: "percentage of residuals within threshold",
            })
        }
    }

    /// Computes the Kolmogorov–Smirnov statistic for the aggregated residual ratios of the accepted residuals.
    ///
    /// Returns Ok(ks_statistic) if residuals are available.
    pub fn ks_test_normality(&self) -> Result<f64, ODError> {
        let mut residual_ratios = self
            .accepted_residuals()
            .iter()
            .flat_map(|resid| resid.whitened_resid.into_iter().copied())
            .collect::<Vec<f64>>();

        if residual_ratios.is_empty() {
            return Err(ODError::ODNoResiduals {
                action: "percentage of residuals within threshold",
            });
        }
        residual_ratios.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let n = residual_ratios.len() as f64;
        let mean = residual_ratios.iter().sum::<f64>() / n;
        let variance = residual_ratios
            .iter()
            .map(|v| (v - mean).powi(2))
            .sum::<f64>()
            / n;
        let std_dev = variance.sqrt();

        // Create a normal distribution based on the empirical mean and std deviation.
        let normal_dist = Normal::new(mean, std_dev).unwrap();

        // Compute the maximum difference between the empirical CDF and the normal CDF.
        let mut d_stat = 0.0;
        for (i, &value) in residual_ratios.iter().enumerate() {
            let empirical_cdf = (i + 1) as f64 / n;
            let model_cdf = normal_dist.cdf(value);
            let diff = (empirical_cdf - model_cdf).abs();
            if diff > d_stat {
                d_stat = diff;
            }
        }
        Ok(d_stat)
    }

    /// Checks whether the whitened residuals of the accepted residuals pass a normality test at a given significance level `alpha`, default to 0.05.
    ///
    /// This uses a simplified KS-test threshold: D_alpha = c(α) / √n.
    /// For example, for α = 0.05, c(α) is approximately 1.36.
    /// α=0.05 means a 5% probability of Type I error (incorrectly rejecting the null hypothesis when it is true).
    /// This threshold is standard in many statistical tests to balance sensitivity and false positives
    ///
    /// The computation of the c(alpha) is from https://en.wikipedia.org/w/index.php?title=Kolmogorov%E2%80%93Smirnov_test&oldid=1280260701#Two-sample_Kolmogorov%E2%80%93Smirnov_test
    ///
    /// Returns Ok(true) if the residuals are consistent with a normal distribution,
    /// Ok(false) if not, or None if no residuals are available.
    pub fn is_normal(&self, alpha: Option<f64>) -> Result<bool, ODError> {
        let n = self.accepted_residuals().len();
        if n == 0 {
            return Err(ODError::ODNoResiduals {
                action: "evaluating residual normality",
            });
        } else if n < 35 {
            warn!("KS normality test unreliable for n={n} < 35; skipping");
        }

        let ks_stat = self.ks_test_normality()?;

        // Default to 5% probability
        let alpha = alpha.unwrap_or(0.05);

        // Compute the c_alpha
        let c_alpha = (-(alpha * 0.5).ln() * 0.5).sqrt();

        let d_critical = c_alpha / (n as f64).sqrt();

        Ok(ks_stat <= d_critical)
    }

    /// Checks whether the filter innovations are statistically consistent
    /// by performing a Chi-squared test on the Normalized Innovation Squared (NIS).
    ///
    /// For each accepted residual, NIS is computed as:
    /// ```text
    ///     prefit^T * S_k^-1 * prefit
    /// ```
    ///
    /// The sum of NIS values should fall within the confidence interval of a
    /// Chi-squared distribution with degrees of freedom `k = n * m`, where `n`
    /// is the number of residuals and `m` is the measurement dimension.
    ///
    /// Returns Ok(true) if the filter is consistent, Ok(false) if the filter
    /// is over-confident or under-confident, or an error if no residuals are available.
    pub fn nis_consistency(&self, alpha: Option<f64>) -> Result<NormalizedConsistency, ODError> {
        let n = self.accepted_residuals().len();

        if n == 0 {
            return Err(ODError::ODNoResiduals {
                action: "evaluating NIS consistency",
            });
        }

        // Sum the NIS across all residuals.
        // Assuming your Residual struct has an `nis` field or a method to compute it
        // from the ratio: `ratio.powi(2) * measurement_dim`
        let nis_sum: f64 = self.accepted_residuals().iter().map(|r| r.nis()).sum();

        // Total degrees of freedom: number of residuals * measurement dimension.
        // Adjust `M::DIM` based on how you access the dimension in this scope.
        let k = (n * MsrSize::DIM) as f64;
        if k < 35.0 {
            warn!("NIS consistency test lacks statistical power for n*MsrSize={k} < 35");
        }
        // Default to a 5% probability of Type I error (95% confidence interval)
        let alpha = alpha.unwrap_or(0.05);

        // For a two-sided test, we need the standard normal quantile for 1 - (alpha / 2).
        // If alpha = 0.05, the critical z-score is approximately 1.95996.
        let z_critical = Normal::new(0.0, 1.0)
            .unwrap()
            .inverse_cdf(1.0 - alpha / 2.0);

        // Use the Wilson-Hilferty transformation to approximate the Chi-squared
        // lower and upper critical bounds.
        let factor = 2.0 / (9.0 * k);
        let lower_bound = k * (1.0 - factor - z_critical * factor.sqrt()).powi(3);
        let upper_bound = k * (1.0 - factor + z_critical * factor.sqrt()).powi(3);

        Ok(NormalizedConsistency {
            normalized_sum: nis_sum,
            lower_bound,
            upper_bound,
            k,
            is_nees: false,
        })
    }

    /// Checks whether the filter estimates are statistically consistent
    /// by performing a Chi-squared test on the Normalized Estimation Error Squared (NEES).
    ///
    /// For each estimate, NEES is computed as:
    /// ```text
    ///     error^T * P^-1 * error
    /// ```
    /// where `error` is the difference between the estimated state and the true state,
    /// and `P` is the estimated state covariance matrix.
    ///
    /// The sum of NEES values should fall within the confidence interval of a
    /// Chi-squared distribution with degrees of freedom `k = n * dim`, where `n`
    /// is the number of estimates and `dim` is the actively estimated state dimension.
    pub fn nees_consistency(
        &self,
        truth_traj: &Traj<StateType>,
        alpha: Option<f64>,
    ) -> Result<NormalizedConsistency, ODError>
    where
        StateType::Size: DimMin<StateType::Size>,
        <StateType::Size as DimMin<StateType::Size>>::Output: DimSub<Const<1>>,
        <<StateType as State>::Size as DimMin<<StateType as State>::Size>>::Output:
            DimSub<Const<1>>,
        DefaultAllocator: Allocator<StateType::Size, Const<1>>
            + Allocator<Const<1>, <StateType as State>::Size>
            + Allocator<<StateType::Size as DimMin<StateType::Size>>::Output, StateType::Size>
            + Allocator<StateType::Size, <StateType::Size as DimMin<StateType::Size>>::Output>
            + Allocator<<StateType::Size as DimMin<StateType::Size>>::Output>
            + Allocator<
                <<StateType::Size as DimMin<StateType::Size>>::Output as DimSub<Const<1>>>::Output,
            >,
    {
        let n_total = self.estimates.len();

        // We will skip the apriori estimate.
        if n_total <= 1 {
            return Err(ODError::ODNoResiduals {
                action: "evaluating NEES consistency",
            });
        }

        // The exact number of epochs evaluated in the NEES sum
        let n = n_total - 1;

        let mut nees_sum = 0.0;
        let mut est_size = None;

        for est in self.estimates.iter().skip(1) {
            let epoch = est.epoch();

            let truth_state = truth_traj.at(epoch).context(ODTrajSnafu {
                details: "interpolating truth during NEES test",
            })?;

            // Extract the state vectors
            let x_est = est.state().to_state_vector();
            let x_true = truth_state.to_state_vector();
            let error = x_est - x_true;
            let cov = est.covar();
            let cov = 0.5 * (&cov + cov.transpose());

            if est_size.is_none() {
                // Deduce est_size by physically inspecting the diagonal variances.
                // Unestimated deterministic parameters will have zero variance.
                let mut active_count = 0;
                for i in 0..StateType::Size::dim() {
                    // Check if variance is strictly positive
                    if cov[(i, i)] > 0.0 {
                        active_count += 1;
                    } else {
                        // Estimated parameters are contiguous from index 0
                        break;
                    }
                }
                if active_count < 6 {
                    warn!("Detected est_size of {active_count} < 6; forcing bound to 6");
                    active_count = 6;
                }
                est_size = Some(active_count);
            }

            // At this point, est_size is set and is valid
            let true_est_size = est_size.unwrap();

            // Dynamically slice the actively estimated block
            let cov_e = cov
                .view((0, 0), (true_est_size, true_est_size))
                .into_owned();
            let err_e = error.rows(0, true_est_size).into_owned();

            // Execute Symmetric Eigendecomposition.
            // This is mathematically safe because cov was forced to be symmetric via 0.5 * (P + P^T)
            let eigen = cov_e.symmetric_eigen();

            // Define the numerical noise floor relative to the largest valid physical variance
            let max_lambda = eigen.eigenvalues.max();
            let epsilon = max_lambda * (true_est_size as f64) * std::f64::EPSILON;

            // Construct the inverted eigenvalue spectrum
            let mut inv_eigenvalues = eigen.eigenvalues.clone();
            for i in 0..true_est_size {
                // By strictly checking > epsilon, we explicitly trap BOTH:
                // 1. Positive numerical noise (e.g., +1e-19)
                // 2. Severe negative numerical drift (e.g., -1e-5)
                if inv_eigenvalues[i] > epsilon {
                    inv_eigenvalues[i] = 1.0 / inv_eigenvalues[i];
                } else {
                    // Clamp corrupted non-PSD subspaces exactly to zero
                    inv_eigenvalues[i] = 0.0;
                }
            }

            // Reconstruct the symmetric positive semi-definite pseudo-inverse: P^+ = Q * Lambda^+ * Q^T
            let p_inv = &eigen.eigenvectors
                * DMatrix::from_diagonal(&inv_eigenvalues)
                * eigen.eigenvectors.transpose();

            // Compute quadratic form: e^T * P^+ * e
            let nees_matrix = err_e.transpose() * p_inv * err_e;
            nees_sum += nees_matrix[(0, 0)];
        }

        // Total degrees of freedom: N evaluated epochs * actively estimated dimension
        let k = (n * est_size.unwrap()) as f64;

        if k < 35.0 {
            warn!("NEES consistency test lacks statistical power for n*StateDim={k} < 35");
        }

        // Default to a 5% probability of Type I error (95% confidence interval)
        let alpha = alpha.unwrap_or(0.05);

        // For a two-sided test, we need the standard normal quantile for 1 - (alpha / 2).
        let z_critical = Normal::new(0.0, 1.0)
            .unwrap()
            .inverse_cdf(1.0 - alpha / 2.0);

        // Use the Wilson-Hilferty transformation to approximate the Chi-squared lower and upper critical bounds.
        let factor = 2.0 / (9.0 * k);
        let lower_bound = k * (1.0 - factor - z_critical * factor.sqrt()).powi(3);
        let upper_bound = k * (1.0 - factor + z_critical * factor.sqrt()).powi(3);

        Ok(NormalizedConsistency {
            normalized_sum: nees_sum,
            lower_bound,
            upper_bound,
            k,
            is_nees: true,
        })
    }
}
