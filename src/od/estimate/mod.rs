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

use super::State;
use crate::cosmic::Orbit;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, OMatrix, OVector};
use crate::Spacecraft;
use hifitime::Epoch;
use std::cmp::PartialEq;
use std::fmt;

pub mod residual;
pub use residual::Residual;
pub mod kfestimate;
pub use kfestimate::KfEstimate;
mod sc_uncertainty;
pub use sc_uncertainty::SpacecraftUncertainty;

/// Stores an Estimate, as the result of a `time_update` or `measurement_update`.
pub trait Estimate<T: State>
where
    Self: Clone + PartialEq + Sized + fmt::Display,
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<<T as State>::VecLength>,
{
    /// An empty estimate. This is useful if wanting to store an estimate outside the scope of a filtering loop.
    fn zeros(state: T) -> Self;
    /// Epoch of this Estimate
    fn epoch(&self) -> Epoch {
        self.state().epoch()
    }
    // Sets the epoch
    fn set_epoch(&mut self, dt: Epoch) {
        self.state().set_epoch(dt);
    }
    /// The estimated state
    fn state(&self) -> T {
        self.nominal_state().add(self.state_deviation())
    }
    /// The state deviation as computed by the filter.
    fn state_deviation(&self) -> OVector<f64, <T as State>::Size>;
    /// The nominal state as reported by the filter dynamics
    fn nominal_state(&self) -> T;
    /// The Covariance of this estimate. Will return the predicted covariance if this is a time update/prediction.
    fn covar(&self) -> OMatrix<f64, <T as State>::Size, <T as State>::Size>;
    /// The predicted covariance of this estimate from the time update
    fn predicted_covar(&self) -> OMatrix<f64, <T as State>::Size, <T as State>::Size>;
    /// Sets the state deviation.
    fn set_state_deviation(&mut self, new_state: OVector<f64, <T as State>::Size>);
    /// Sets the Covariance of this estimate
    fn set_covar(&mut self, new_covar: OMatrix<f64, <T as State>::Size, <T as State>::Size>);
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    fn predicted(&self) -> bool;
    /// The STM used to compute this Estimate
    fn stm(&self) -> &OMatrix<f64, <T as State>::Size, <T as State>::Size>;
    /// Returns whether this estimate is within some bound
    /// The 68-95-99.7 rule is a good way to assess whether the filter is operating normally
    fn within_sigma(&self, sigma: f64) -> bool {
        let state = self.state_deviation();
        let covar = self.covar();
        for i in 0..state.len() {
            let bound = covar[(i, i)].sqrt() * sigma;
            if state[i] > bound || state[i] < -bound {
                return false;
            }
        }
        true
    }
    /// Returns whether this estimate is within 3 sigma, which represent 99.7% for a Normal distribution
    fn within_3sigma(&self) -> bool {
        self.within_sigma(3.0)
    }
}

/// A trait to store a navigation solution, can be used in conjunction with KfEstimate
pub trait NavSolution<T>: Estimate<Spacecraft>
where
    T: State,
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<<T as State>::VecLength>,
{
    fn orbital_state(&self) -> Orbit;
    /// Returns the nominal state as computed by the dynamics
    fn expected_state(&self) -> Orbit;
}

impl NavSolution<Spacecraft> for KfEstimate<Spacecraft> {
    fn orbital_state(&self) -> Orbit {
        self.state().orbit
    }
    fn expected_state(&self) -> Orbit {
        self.nominal_state().orbit
    }
}
