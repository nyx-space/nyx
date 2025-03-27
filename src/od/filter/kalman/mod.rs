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
use crate::linalg::{DefaultAllocator, DimName};
pub use crate::od::estimate::{Estimate, KfEstimate, Residual};
pub use crate::od::snc::ProcessNoise;
use crate::od::State;
pub use crate::time::{Epoch, Unit};

pub mod filtering;
pub mod initializers;

/// Defines both a Classical and an Extended Kalman filter (CKF and EKF)
/// T: Type of state
/// A: Acceleration size (for SNC)
/// M: Measurement size (used for the sensitivity matrix)
#[derive(Debug, Clone)]
#[allow(clippy::upper_case_acronyms)]
pub struct KF<T, A>
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
    /// The previous estimate used in the KF computations.
    pub prev_estimate: KfEstimate<T>,
    /// A sets of process noise (usually noted Q), must be ordered chronologically
    pub process_noise: Vec<ProcessNoise<A>>,
    /// The variant of this Kalman filter.
    pub variant: KalmanVariant,
    pub prev_used_snc: usize,
}

#[derive(Copy, Clone, Debug, Default, PartialEq)]
pub enum KalmanVariant {
    /// Configures the filter as a standard Extended Kalman Filter (EKF) update,
    /// updating the full reference state in the process' propagator at each measurement update.
    #[default]
    ReferenceUpdate,
    /// Tracks the state deviation (formerly called Classical Kalman Filter (CKF)) and does not update the reference in the process' propagator.
    DeviationTracking,
    /// Configures the Kalman Filter as an Iterative Extended Kalman Filter (IEKF). This should only be used to trust the observations even more than a single EKF update would.
    ///
    /// The three Option fields control the iteration:
    /// - pos_km: if Some(val), then iterations stop when the RSS of the position part of the state error drops below `val` in km.
    /// - vel_km_s: if Some(val), then iterations stop when the RSS of the velocity part drops below `val` in km/s.
    /// - max_iter: if Some(n), force exit after n iterations.
    ///
    /// If all three are None, then the filter acts like the StandardEKF mode.
    IterativeUpdate {
        pos_km: Option<f64>,
        vel_km_s: Option<f64>,
        max_iter: Option<u8>,
    },
}
