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
use crate::linalg::{DefaultAllocator, U3};
pub use crate::od::estimate::{Estimate, KfEstimate, Residual};
pub use crate::od::snc::ProcessNoise;
use crate::od::State;
pub use crate::time::{Epoch, Unit};

use super::{KalmanVariant, KF};

// Only process noises of size 3 are supported, so only these have initializers.

impl<T> KF<T, U3>
where
    T: State,
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<<T as State>::VecLength>
        + Allocator<U3>
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<U3, U3>
        + Allocator<<T as State>::Size, U3>
        + Allocator<U3, <T as State>::Size>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    /// Initializes this KF with an initial estimate, measurement noise, and one process noise
    pub fn new(initial_estimate: KfEstimate<T>, variant: KalmanVariant) -> Self {
        Self {
            prev_estimate: initial_estimate,
            process_noise: vec![],
            variant,
            prev_used_snc: 0,
        }
    }

    /// Set (or replaces) the existing process noise configuration.
    pub fn with_process_noise(mut self, mut process_noise: ProcessNoise<U3>) -> Self {
        process_noise.init_epoch = Some(self.prev_estimate.epoch());
        self.process_noise.clear();
        self.process_noise.push(process_noise);
        self
    }

    /// Pushes the provided process noise to the list the existing process noise configurations.
    pub fn and_with_process_noise(mut self, process_noise: ProcessNoise<U3>) -> Self {
        self.process_noise.push(process_noise);
        self
    }
}
