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
use crate::linalg::{DefaultAllocator, DimName, U3};
pub use crate::od::estimate::{Estimate, KfEstimate, Residual};
pub use crate::od::snc::SNC;
use crate::od::State;
pub use crate::time::{Epoch, Unit};

use super::KF;

impl<T, A> KF<T, A>
where
    A: DimName,
    T: State,
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<<T as State>::VecLength>
        + Allocator<A>
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
            prev_used_snc: 0,
        }
    }
}

impl<T> KF<T, U3>
where
    T: State,
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<<T as State>::VecLength>
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
            prev_used_snc: 0,
        }
    }
}
