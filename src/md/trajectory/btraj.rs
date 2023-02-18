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

use super::spline::INTERPOLATION_SAMPLES;
use super::{InterpState, TrajError};
use crate::errors::NyxError;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::time::Epoch;

/// Store a trajectory of any State.
#[derive(Clone, Default)]
pub struct BTraj<S: InterpState>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    /// Optionally name this trajectory
    pub name: Option<String>,
    /// We use a vector because we know that the states are produced in a chronological manner (the direction does not matter).
    pub states: Vec<S>,
}

impl<S: InterpState> BTraj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    pub fn new() -> Self {
        Self {
            name: None,
            states: Vec::new(),
        }
    }
    /// Orders the states, can be used to store the states out of order
    pub fn finalize(&mut self) {
        self.states.sort_by(|a, b| a.epoch().cmp(&b.epoch()));
    }

    /// Evaluate the trajectory at this specific epoch.
    pub fn at(&self, epoch: Epoch) -> Result<S, NyxError> {
        if self.states.is_empty() {
            return Err(NyxError::Trajectory(TrajError::NoInterpolationData(epoch)));
        }
        match self
            .states
            .binary_search_by(|state| state.epoch().cmp(&epoch))
        {
            Ok(idx) => {
                // Oh wow, we actually had this exact state!
                Ok(self.states[idx])
            }
            Err(idx) => {
                // This is the closest index, so let's grab the items around it.
                // NOTE: This is essentially the same code as in ANISE for the Hermite SPK type 13

                // We didn't find it, so let's build an interpolation here.
                let num_left = INTERPOLATION_SAMPLES / 2;

                // Ensure that we aren't fetching out of the window
                let mut first_idx = idx.saturating_sub(num_left);
                let last_idx = self.states.len().min(first_idx + INTERPOLATION_SAMPLES);

                // Check that we have enough samples
                if last_idx == self.states.len() {
                    first_idx = last_idx - 2 * num_left;
                }

                let states = Vec::with_capacity(last_idx - first_idx);

                self.states[idx].interpolate(epoch, &states)
            }
        }
    }
}
