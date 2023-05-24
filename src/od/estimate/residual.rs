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

use crate::hifitime::Epoch;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OVector};
use std::fmt;

/// Stores an Estimate, as the result of a `time_update` or `measurement_update`.
#[derive(Debug, Clone, PartialEq)]
pub struct Residual<M>
where
    M: DimName,
    DefaultAllocator: Allocator<f64, M>,
{
    /// Date time of this Residual
    pub epoch: Epoch,
    /// The prefit residual in the units of the measurement type
    pub prefit: OVector<f64, M>,
    /// The postfit residual in the units of the measurement type
    pub postfit: OVector<f64, M>,
    /// The prefit residual ratio, i.e. `r' * (H*P*H')^-1 * r`, where `r` is the prefit residual, `H` is the sensitivity matrix, and `P` is the covariance matrix.
    pub ratio: f64,
    /// Whether or not this was rejected
    pub rejected: bool,
}

impl<M> Residual<M>
where
    M: DimName,
    DefaultAllocator: Allocator<f64, M>,
{
    /// An empty estimate. This is useful if wanting to store an estimate outside the scope of a filtering loop.
    pub fn zeros() -> Self {
        Self {
            epoch: Epoch::from_tai_seconds(0.0),
            prefit: OVector::<f64, M>::zeros(),
            postfit: OVector::<f64, M>::zeros(),
            ratio: 0.0,
            rejected: true,
        }
    }

    /// Flags a Residual as rejected.
    pub fn rejected(epoch: Epoch, prefit: OVector<f64, M>, ratio: f64) -> Self {
        Self {
            epoch,
            prefit,
            postfit: OVector::<f64, M>::zeros(),
            ratio,
            rejected: true,
        }
    }

    pub fn new(
        epoch: Epoch,
        prefit: OVector<f64, M>,
        postfit: OVector<f64, M>,
        ratio: f64,
    ) -> Self {
        Self {
            epoch,
            prefit,
            postfit,
            ratio,
            rejected: false,
        }
    }
}

impl<M> fmt::Display for Residual<M>
where
    M: DimName,
    DefaultAllocator: Allocator<f64, M> + Allocator<usize, M>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Prefit {} Postfit {}", &self.prefit, &self.postfit)
    }
}

impl<M> fmt::LowerExp for Residual<M>
where
    M: DimName,
    DefaultAllocator: Allocator<f64, M> + Allocator<usize, M>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Prefit {:e} Postfit {:e}", &self.prefit, &self.postfit)
    }
}
