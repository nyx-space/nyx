/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use crate::cosmic::Frame;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector, U3, U6};
use crate::time::{Duration, Epoch};

use std::fmt;

#[allow(clippy::upper_case_acronyms)]
pub type SNC3 = SNC<U3>;
#[allow(clippy::upper_case_acronyms)]
pub type SNC6 = SNC<U6>;

#[derive(Clone)]
#[allow(clippy::upper_case_acronyms)]
pub struct SNC<A: DimName>
where
    DefaultAllocator: Allocator<f64, A> + Allocator<f64, A, A>,
{
    /// Time at which this SNC starts to become applicable
    pub start_time: Option<Epoch>,
    /// Specify the frame of this SNC -- CURRENTLY UNIMPLEMENTED
    pub frame: Option<Frame>,
    /// Enables state noise compensation (process noise) only be applied if the time between measurements is less than the disable_time amount in seconds
    pub disable_time: Duration,
    // Stores the initial epoch when the SNC is requested, needed for decay. Kalman filter will edit this automatically.
    pub init_epoch: Option<Epoch>,
    diag: OVector<f64, A>,
    decay_diag: Option<Vec<f64>>,
    // Stores the previous epoch of the SNC request, needed for disable time
    pub prev_epoch: Option<Epoch>,
}

impl<A> fmt::Debug for SNC<A>
where
    A: DimName,
    DefaultAllocator: Allocator<f64, A> + Allocator<f64, A, A>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut fmt_cov = Vec::with_capacity(A::dim());
        if let Some(decay) = &self.decay_diag {
            for (i, dv) in decay.iter().enumerate() {
                fmt_cov.push(format!("{:.1e} × exp(- {:.1e} × t)", self.diag[i], dv));
            }
        } else {
            for i in 0..A::dim() {
                fmt_cov.push(format!("{:.1e}", self.diag[i]));
            }
        }
        write!(
            f,
            "SNC: diag({}) {}",
            fmt_cov.join(", "),
            if let Some(start) = self.start_time {
                format!("starting at {}", start.as_gregorian_utc_str())
            } else {
                "".to_string()
            }
        )
    }
}

impl<A> fmt::Display for SNC<A>
where
    A: DimName,
    DefaultAllocator: Allocator<f64, A> + Allocator<f64, A, A>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl<A: DimName> SNC<A>
where
    DefaultAllocator: Allocator<f64, A> + Allocator<f64, A, A>,
{
    /// Initialize a state noise compensation structure from the diagonal values
    pub fn from_diagonal(disable_time: Duration, values: &[f64]) -> Self {
        assert_eq!(
            values.len(),
            A::dim(),
            "Not enough values for the size of the SNC matrix"
        );

        let mut diag = OVector::zeros();
        for (i, v) in values.iter().enumerate() {
            diag[i] = *v;
        }

        Self {
            diag,
            disable_time,
            start_time: None,
            frame: None,
            decay_diag: None,
            init_epoch: None,
            prev_epoch: None,
        }
    }

    /// Initialize an SNC with a time at which it should start
    pub fn with_start_time(disable_time: Duration, values: &[f64], start_time: Epoch) -> Self {
        let mut me = Self::from_diagonal(disable_time, values);
        me.start_time = Some(start_time);
        me
    }

    /// Initialize an exponentially decaying SNC with initial SNC and decay constants.
    /// Decay constants in seconds since start of the tracking pass.
    pub fn with_decay(
        disable_time: Duration,
        initial_snc: &[f64],
        decay_constants_s: &[f64],
    ) -> Self {
        assert_eq!(
            decay_constants_s.len(),
            A::dim(),
            "Not enough decay constants for the size of the SNC matrix"
        );

        let mut me = Self::from_diagonal(disable_time, initial_snc);
        me.decay_diag = Some(decay_constants_s.to_vec());
        me
    }

    /// Returns the SNC matrix (_not_ incl. Gamma matrix approximation) at the provided Epoch.
    /// May be None if:
    ///  1. Start time of this matrix is _after_ epoch
    ///  2. Time between epoch and previous epoch (set in the Kalman filter!) is longer than disabling time
    pub fn to_matrix(&self, epoch: Epoch) -> Option<OMatrix<f64, A, A>> {
        if let Some(start_time) = self.start_time {
            if start_time > epoch {
                // This SNC applies only later
                debug!("@{} SNC starts at {}", epoch, start_time);
                return None;
            }
        }

        // Check the disable time, and return no SNC if the previous SNC was computed too long ago
        if let Some(prev_epoch) = self.prev_epoch {
            if epoch - prev_epoch > self.disable_time {
                debug!(
                    "@{} SNC disabled: prior use greater than {}",
                    epoch, self.disable_time
                );
                return None;
            }
        }
        // Build a static matrix
        let mut snc = OMatrix::<f64, A, A>::zeros();
        for i in 0..self.diag.nrows() {
            snc[(i, i)] = self.diag[i];
        }

        if let Some(decay) = &self.decay_diag {
            // Let's apply the decay to the diagonals
            let total_delta_t = (epoch - self.init_epoch.unwrap()).in_seconds();
            for i in 0..self.diag.nrows() {
                snc[(i, i)] *= (-decay[i] * total_delta_t).exp();
            }
        }

        trace!(
            "@{} SNC diag {:?}",
            epoch,
            snc.diagonal().iter().copied().collect::<Vec<f64>>()
        );

        Some(snc)
    }
}

#[test]
fn test_snc_init() {
    use crate::time::Unit;
    let snc_expo = SNC3::with_decay(
        2 * Unit::Minute,
        &[1e-6, 1e-6, 1e-6],
        &[3600.0, 3600.0, 3600.0],
    );
    println!("{}", snc_expo);

    let snc_std = SNC3::with_start_time(
        2 * Unit::Minute,
        &[1e-6, 1e-6, 1e-6],
        Epoch::from_et_seconds(3600.0),
    );
    println!("{}", snc_std);
}
