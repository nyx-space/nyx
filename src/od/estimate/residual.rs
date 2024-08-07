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

use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OVector};
use hifitime::Epoch;
use std::fmt;

/// Stores an Estimate, as the result of a `time_update` or `measurement_update`.
#[derive(Debug, Clone, PartialEq)]
pub struct Residual<M>
where
    M: DimName,
    DefaultAllocator: Allocator<M>,
{
    /// Date time of this Residual
    pub epoch: Epoch,
    /// The prefit residual in the units of the measurement type
    pub prefit: OVector<f64, M>,
    /// The postfit residual in the units of the measurement type
    pub postfit: OVector<f64, M>,
    /// The prefit residual ratio, i.e. `r' * (H*P*H')^-1 * r`, where `r` is the prefit residual, `H` is the sensitivity matrix, and `P` is the covariance matrix.
    pub ratio: f64,
    /// The tracker measurement noise (variance)) for this tracker at this time.
    pub tracker_msr_noise: OVector<f64, M>,
    /// Whether or not this was rejected
    pub rejected: bool,
    /// Name of the tracker that caused this residual
    pub tracker: Option<String>,
}

impl<M> Residual<M>
where
    M: DimName,
    DefaultAllocator: Allocator<M>,
{
    /// An empty estimate. This is useful if wanting to store an estimate outside the scope of a filtering loop.
    pub fn zeros() -> Self {
        Self {
            epoch: Epoch::from_tai_seconds(0.0),
            prefit: OVector::<f64, M>::zeros(),
            postfit: OVector::<f64, M>::zeros(),
            tracker_msr_noise: OVector::<f64, M>::zeros(),
            ratio: 0.0,
            rejected: true,
            tracker: None,
        }
    }

    /// Flags a Residual as rejected.
    pub fn rejected(
        epoch: Epoch,
        prefit: OVector<f64, M>,
        ratio: f64,
        tracker_msr_noise: OVector<f64, M>,
    ) -> Self {
        Self {
            epoch,
            prefit,
            postfit: OVector::<f64, M>::zeros(),
            ratio,
            tracker_msr_noise,
            rejected: true,
            tracker: None,
        }
    }

    pub fn accepted(
        epoch: Epoch,
        prefit: OVector<f64, M>,
        postfit: OVector<f64, M>,
        ratio: f64,
        tracker_msr_noise: OVector<f64, M>,
    ) -> Self {
        Self {
            epoch,
            prefit,
            postfit,
            ratio,
            tracker_msr_noise,
            rejected: false,
            tracker: None,
        }
    }
}

impl<M> fmt::Display for Residual<M>
where
    M: DimName,
    DefaultAllocator: Allocator<M> + Allocator<M>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "Residual from {} at {}: ratio = {:.3}\nPrefit {} Postfit {}",
            self.tracker.as_ref().unwrap_or(&"Unknown".to_string()),
            self.epoch,
            self.ratio,
            &self.prefit,
            &self.postfit
        )
    }
}

impl<M> fmt::LowerExp for Residual<M>
where
    M: DimName,
    DefaultAllocator: Allocator<M> + Allocator<M>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Prefit {:e} Postfit {:e}", &self.prefit, &self.postfit)
    }
}
