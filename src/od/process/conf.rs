/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use crate::time::{Duration, Epoch};
use crate::NyxError;
use std::convert::TryFrom;
use std::default::Default;
use std::fmt;
use typed_builder::TypedBuilder;

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Defines the stopping condition for the smoother
#[derive(Clone, Copy, Debug)]
pub enum SmoothingArc {
    /// Stop smoothing when the gap between estimate is the provided floating point number in seconds
    TimeGap(Duration),
    /// Stop smoothing at the provided Epoch.
    Epoch(Epoch),
    /// Stop smoothing at the first prediction
    Prediction,
    /// Only stop once all estimates have been processed
    All,
}

impl fmt::Display for SmoothingArc {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            SmoothingArc::All => write!(f, "all estimates"),
            SmoothingArc::Epoch(e) => write!(f, "{e}"),
            SmoothingArc::TimeGap(g) => write!(f, "time gap of {g}"),
            SmoothingArc::Prediction => write!(f, "first prediction"),
        }
    }
}

impl Default for SmoothingArc {
    fn default() -> Self {
        Self::All
    }
}

/// Defines a filter iteration configuration. Allows iterating on an OD solution until convergence criteria is met.
/// The root mean squared of the prefit residuals ratios is used to assess convergence between iterations.
#[derive(Clone, Copy, Debug, TypedBuilder)]
#[cfg_attr(feature = "python", pyclass)]
pub struct IterationConf {
    /// The number of measurements to account for in the iteration
    #[builder(default)]
    pub smoother: SmoothingArc,
    /// The absolute tolerance of the RMS prefit residual ratios
    #[builder(default = 1e-1)]
    pub absolute_tol: f64,
    /// The relative tolerance between the latest RMS prefit residual ratios and the best RMS prefit residual ratios so far
    #[builder(default = 1e-2)]
    pub relative_tol: f64,
    /// The maximum number of iterations to allow (will raise an error if the filter has not converged after this many iterations)
    #[builder(default = 15)]
    pub max_iterations: usize,
    /// The maximum number of subsequent divergences in RMS.
    #[builder(default = 3)]
    pub max_divergences: usize,
    /// Set to true to force an ODP failure when the convergence criteria is not met
    #[builder(default = false)]
    pub force_failure: bool,
}

impl IterationConf {
    /// Iterate and smooth only once
    pub fn once() -> Self {
        Self {
            max_iterations: 1,
            ..Default::default()
        }
    }
}

#[cfg(feature = "python")]
#[pymethods]
impl IterationConf {
    #[new]
    fn py_new(
        strategy: Option<String>,
        duration: Option<Duration>,
        epoch: Option<Epoch>,
        absolute_tol: Option<f64>,
        relative_tol: Option<f64>,
        max_iterations: Option<usize>,
        max_divergences: Option<usize>,
        force_failure: Option<bool>,
    ) -> Result<Self, NyxError> {
        let mut me = Self {
            smoother: if let Some(strategy) = strategy {
                match strategy.to_lowercase().trim() {
                    "all" => SmoothingArc::All,
                    "prediction" => SmoothingArc::Prediction,
                    _ => {
                        return Err(NyxError::CustomError {
                            msg: "strategy should be `all` or `prediction`".to_string(),
                        })
                    }
                }
            } else if let Some(duration) = duration {
                SmoothingArc::TimeGap(duration)
            } else if let Some(epoch) = epoch {
                SmoothingArc::Epoch(epoch)
            } else {
                return Err(NyxError::CustomError {
                    msg: "Smoothing arc strategy not specified".to_string(),
                });
            },
            ..Default::default()
        };

        if let Some(abs_tol) = absolute_tol {
            me.absolute_tol = abs_tol;
        }

        if let Some(rel_tol) = relative_tol {
            me.relative_tol = rel_tol;
        }

        if let Some(max_it) = max_iterations {
            me.max_iterations = max_it;
        }

        if let Some(max_div) = max_divergences {
            me.max_divergences = max_div;
        }

        if let Some(force_failure) = force_failure {
            me.force_failure = force_failure;
        }

        Ok(me)
    }

    fn __repr__(&self) -> String {
        format!("{self}")
    }

    fn __str__(&self) -> String {
        format!("Smoothing {self}")
    }
}

impl Default for IterationConf {
    /// The default absolute tolerance is 1e-2 (calibrated on an EKF with error).
    fn default() -> Self {
        Self {
            smoother: SmoothingArc::All,
            absolute_tol: 1e-1,
            relative_tol: 1e-2,
            max_iterations: 15,
            max_divergences: 3,
            force_failure: false,
        }
    }
}

impl fmt::Display for IterationConf {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Iterate until abs = {:.2e}, or rel = {:.2e}, or {} iterations, or {} subsequent divergences with smoothing condition of {}",
            self.absolute_tol,
            self.relative_tol,
            self.max_iterations,
            self.max_divergences,
            self.smoother)
    }
}

impl TryFrom<SmoothingArc> for IterationConf {
    type Error = NyxError;

    /// Converts a smoother into an interation configuration to iterate just once without failing
    fn try_from(smoother: SmoothingArc) -> Result<Self, Self::Error> {
        Ok(Self {
            smoother,
            force_failure: false,
            ..Default::default()
        })
    }
}
