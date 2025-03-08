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

use nalgebra::Const;

pub use crate::errors::NyxError;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector, U1, U3};
pub use crate::od::estimate::{Estimate, LsqEstimate, Residual};
use crate::od::process::ResidRejectCrit;
pub use crate::od::snc::SNC;
use crate::od::{Filter, ODError, State};
pub use crate::time::{Epoch, Unit};

/// [LSQ] Iterative Least-Square filter.
/// - T: type of [State]
#[derive(Clone)]
#[allow(clippy::upper_case_acronyms)]
pub struct LSQ<T>
where
    T: State,
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<T::Size, T::Size>
        + Allocator<<T as State>::VecLength>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    /// Internal state
    prev_estimate: LsqEstimate<T>,
    /// Measurement vector
    b: OVector<f64, T::Size>,
    /// Internal matrix
    h: OMatrix<f64, T::Size, T::Size>,
}

impl<T> LSQ<T>
where
    T: State,
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<T::Size, T::Size>
        + Allocator<<T as State>::VecLength>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    /// Creates a new [LSQ] filter with inital [State]
    pub fn new(initial_state: LsqEstimate<T>) -> Self {
        assert_eq!(
            T::Size::USIZE,
            T::VecLength::USIZE,
            "H matrix and Noise process should have identical dimensions in LSQ"
        );

        Self {
            prev_estimate: initial_state,
            b: OVector::<f64, T::Size>::zeros(),
            h: OMatrix::<f64, T::Size, T::Size>::zeros(),
        }
    }
}

impl<T> Filter<T, T::Size, T::Size> for LSQ<T>
where
    T: State,
    DefaultAllocator: Allocator<T::Size>
        + Allocator<<T as State>::VecLength, <T as State>::Size>
        + Allocator<<T as State>::Size, <T as State>::VecLength>
        + Allocator<<T as State>::VecLength>
        + Allocator<<T as State>::VecLength, <T as State>::VecLength>
        + Allocator<T::Size, T::Size>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    type Estimate = LsqEstimate<T>;

    fn previous_estimate(&self) -> &Self::Estimate {
        &self.prev_estimate
    }

    fn set_previous_estimate(&mut self, state: &Self::Estimate) {
        self.prev_estimate = *state;
    }

    /// Not applicable to LSQ
    fn is_extended(&self) -> bool {
        false
    }

    /// Not applicable to LSQ
    fn set_extended(&mut self, _: bool) {}

    fn set_process_noise(&mut self, snc: SNC<T::Size>) {
        let t = self.prev_estimate.epoch();

        // TODO: currently limited to one vector
        // <=> one noise process ?

        if let Some(matrix) = snc.to_matrix(t) {
            assert_eq!(matrix.nrows(), T::Size::USIZE, "dimension issue");

            // Measurement is always a vector in LSQ
            let col = matrix.column(0);

            // this cast should be 100% of the time feasible
            // reduce matrix[0] as column of b vector
            // self.b = col.into(); // this triggers a rather complex Matrix/Vector size issue
        }
    }

    fn update_h_tilde(&mut self, h_tilde: OMatrix<f64, T::Size, T::Size>) {
        self.h = h_tilde;
    }

    fn time_update(&mut self, nominal_state: T) -> Result<Self::Estimate, ODError> {
        Ok(self.prev_estimate)
    }

    fn measurement_update(
        &mut self,
        nominal_state: T,
        b: &OVector<f64, T::Size>,
        _: &OVector<f64, T::Size>,
        _: OMatrix<f64, T::Size, T::Size>,
        _: Option<ResidRejectCrit>,
    ) -> Result<(Self::Estimate, Residual<T::Size>), ODError> {
        // apply the LSQ equation
        let ht = self.h.transpose();
        let ht_h = ht.clone() * self.h.clone();
        let ht_h_inv = ht_h.try_inverse().ok_or(ODError::MatrixInversion)?;
        let ht_b = ht * self.b;

        let estimate = LsqEstimate {
            nominal_state,
            covar: ht_h_inv,
            stm: self.h,
            state_deviation: OMatrix::<f64, T::Size, T::Size>::zeros(),
        };

        self.prev_estimate = estimate;
        Ok(estimate)
    }
}
