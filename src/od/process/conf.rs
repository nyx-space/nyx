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

/// Defines a filter iteration configuration. Allows iterating on an OD solution until convergence criteria is met.
/// The root mean squared of the postfit residuals is used to assess convergence between iterations.
#[derive(Clone, Copy, Debug)]
pub struct IterationConf {
    /// The number of measurements to account for in the iteration
    pub smoother: SmoothingArc,
    /// The absolute tolerance of the RMS postfit residual
    pub absolute_tol: f64,
    /// The relative tolerance between the latest RMS postfit residual and the best RMS postfit residual so far
    pub relative_tol: f64,
    /// The maximum number of iterations to allow (will raise an error if the filter has not converged after this many iterations)
    pub max_iterations: usize,
    /// The maximum number of subsequent divergences in RMS.
    pub max_divergences: usize,
    /// Set to true to force an ODP failure when the convergence criteria is not met
    pub force_failure: bool,
    /// Set to true to use the RMS prefit instead of postfit
    pub use_prefit: bool,
}

impl Default for IterationConf {
    /// The default absolute tolerance is 1e-2 (calibrated on an EKF with error).
    fn default() -> Self {
        Self {
            smoother: SmoothingArc::All,
            absolute_tol: 1e-2,
            relative_tol: 1e-3,
            max_iterations: 15,
            max_divergences: 3,
            force_failure: false,
            use_prefit: false,
        }
    }
}

impl fmt::Display for IterationConf {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let kind = if self.use_prefit { "prefit" } else { "postfit" };
        write!(f, "Iterate {} residuals until abs = {:.2e}, or rel = {:.2e}, or {} iterations, or {} subsequent divergences with smoothing condition of {}",
            kind,
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
