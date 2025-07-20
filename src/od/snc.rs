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

use crate::dynamics::guidance::LocalFrame;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector, U3, U6};
use crate::time::{Duration, Epoch};

use std::fmt;

#[allow(clippy::upper_case_acronyms)]
pub type ProcessNoise3D = ProcessNoise<U3>;

#[allow(clippy::upper_case_acronyms)]
pub type ProcessNoise6D = ProcessNoise<U6>;

#[deprecated = "SNC has been renamed to ProcessNoise"]
#[allow(clippy::upper_case_acronyms)]
pub type SNC<A> = ProcessNoise<A>;
#[deprecated = "SNC has been renamed to ProcessNoise"]
#[allow(clippy::upper_case_acronyms)]
pub type SNC3 = ProcessNoise<U3>;
#[deprecated = "SNC has been renamed to ProcessNoise"]
#[allow(clippy::upper_case_acronyms)]
pub type SNC6 = ProcessNoise<U6>;

#[derive(Clone)]
#[allow(clippy::upper_case_acronyms)]
pub struct ProcessNoise<A: DimName>
where
    DefaultAllocator: Allocator<A> + Allocator<A, A>,
{
    /// Time at which this SNC starts to become applicable
    pub start_time: Option<Epoch>,
    /// Specify the local frame of this SNC
    pub local_frame: Option<LocalFrame>,
    /// Enables state noise compensation (process noise) only be applied if the time between measurements is less than the disable_time
    pub disable_time: Duration,
    // Stores the initial epoch when the SNC is requested, needed for decay. Kalman filter will edit this automatically.
    pub init_epoch: Option<Epoch>,
    diag: OVector<f64, A>,
    decay_diag: Option<Vec<f64>>,
    // Stores the previous epoch of the SNC request, needed for disable time
    pub prev_epoch: Option<Epoch>,
}

impl<A> fmt::Debug for ProcessNoise<A>
where
    A: DimName,
    DefaultAllocator: Allocator<A> + Allocator<A, A>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut fmt_cov = Vec::with_capacity(A::dim());
        if let Some(decay) = &self.decay_diag {
            for (i, dv) in decay.iter().enumerate() {
                fmt_cov.push(format!(
                    "{:.1e} × exp(- {:.1e} × t)",
                    self.diag[i] * 1e6,
                    dv
                ));
            }
        } else {
            for i in 0..A::dim() {
                fmt_cov.push(format!("{:.1e}", self.diag[i] * 1e6));
            }
        }
        write!(
            f,
            "SNC: diag({}) mm/s^2 {} {}",
            fmt_cov.join(", "),
            if let Some(lf) = self.local_frame {
                format!("in {lf:?}")
            } else {
                "".to_string()
            },
            if let Some(start) = self.start_time {
                format!("starting at {start}")
            } else {
                "".to_string()
            }
        )
    }
}

impl<A> fmt::Display for ProcessNoise<A>
where
    A: DimName,
    DefaultAllocator: Allocator<A> + Allocator<A, A>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{self:?}")
    }
}

impl<A: DimName> ProcessNoise<A>
where
    DefaultAllocator: Allocator<A> + Allocator<A, A>,
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
            local_frame: None,
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
                debug!("@{epoch} SNC starts at {start_time}");
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
            let total_delta_t = (epoch - self.init_epoch.unwrap()).to_seconds();
            for i in 0..self.diag.nrows() {
                snc[(i, i)] *= (-decay[i] * total_delta_t).exp();
            }
        }

        debug!(
            "@{} SNC diag {:?}",
            epoch,
            snc.diagonal().iter().copied().collect::<Vec<f64>>()
        );

        Some(snc)
    }
}

impl ProcessNoise3D {
    /// Initialize the process noise from velocity errors over time
    pub fn from_velocity_km_s(
        velocity_noise: &[f64; 3],
        noise_duration: Duration,
        disable_time: Duration,
        local_frame: Option<LocalFrame>,
    ) -> Self {
        let mut diag = OVector::<f64, U3>::zeros();
        for (i, v) in velocity_noise.iter().enumerate() {
            diag[i] = *v / noise_duration.to_seconds();
        }

        Self {
            diag,
            disable_time,
            start_time: None,
            local_frame,
            decay_diag: None,
            init_epoch: None,
            prev_epoch: None,
        }
    }
}

#[test]
fn test_snc_init() {
    use crate::time::Unit;
    let snc_expo = ProcessNoise3D::with_decay(
        2 * Unit::Minute,
        &[1e-6, 1e-6, 1e-6],
        &[3600.0, 3600.0, 3600.0],
    );
    println!("{snc_expo}");

    let snc_std = ProcessNoise3D::with_start_time(
        2 * Unit::Minute,
        &[1e-6, 1e-6, 1e-6],
        Epoch::from_et_seconds(3600.0),
    );
    println!("{snc_std}");

    let snc_vel = ProcessNoise3D::from_velocity_km_s(
        &[1e-2, 1e-2, 1e-2],
        Unit::Hour * 2,
        Unit::Minute * 2,
        Some(LocalFrame::RIC),
    );

    let diag_val = 1e-2 / (Unit::Hour * 2).to_seconds();
    assert_eq!(
        snc_vel
            .to_matrix(Epoch::from_et_seconds(3600.0))
            .unwrap()
            .diagonal(),
        OVector::<f64, U3>::new(diag_val, diag_val, diag_val)
    );
    assert_eq!(snc_vel.local_frame, Some(LocalFrame::RIC));
}
