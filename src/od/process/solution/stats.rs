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
use log::warn;
use msr::sensitivity::TrackerSensitivity;
use statrs::distribution::{ContinuousCDF, Normal};
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
        let n = self.residuals.iter().flatten().count();
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
    pub fn is_nis_consistent(&self, alpha: Option<f64>) -> Result<bool, ODError> {
        let n = self.accepted_residuals().iter().count();

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

        if nis_sum > upper_bound {
            warn!("NIS consistency test failed high: NIS sum {nis_sum:.6} > upper bound {upper_bound:.6}. Innovations are larger than expected.");
            warn!("Filter may be overconfident: P, R, or Q may be too small, or the dynamics/measurement model may be biased.");
            Ok(false)
        } else if nis_sum < lower_bound {
            warn!("NIS consistency test failed low: NIS sum {nis_sum:.6} < lower bound {lower_bound:.6}. Innovations are smaller than expected.");
            warn!("Filter may be underconfident: P, R, or Q may be too large.");
            Ok(false)
        } else {
            Ok(true)
        }
    }
}
