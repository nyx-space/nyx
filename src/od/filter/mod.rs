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

use self::kalman::Residual;

use super::estimate::Estimate;
use super::process::ResidRejectCrit;
use super::snc::ProcessNoise;
use super::ODError;
pub use crate::dynamics::Dynamics;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector};
pub use crate::{State, TimeTagged};
pub mod kalman;

/// Defines a Filter trait where S is the size of the estimated state, A the number of acceleration components of the EOMs (used for process noise matrix size), M the size of the measurements.
pub trait Filter<T, A, M>
where
    A: DimName,
    M: DimName,
    T: State,
    DefaultAllocator: Allocator<M>
        + Allocator<<T as State>::Size>
        + Allocator<<T as State>::VecLength>
        + Allocator<A>
        + Allocator<M, M>
        + Allocator<M, <T as State>::Size>
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<A, A>
        + Allocator<<T as State>::Size, A>
        + Allocator<A, <T as State>::Size>,
{
    type Estimate: Estimate<T>;

    /// Returns the previous estimate
    fn previous_estimate(&self) -> &Self::Estimate;

    /// Set the previous estimate
    fn set_previous_estimate(&mut self, est: &Self::Estimate);

    /// Computes a time update/prediction at the provided nominal state (i.e. advances the filter estimate with the updated STM).
    ///
    /// Returns an error if the STM was not updated.
    fn time_update(&mut self, nominal_state: T) -> Result<Self::Estimate, ODError>;

    /// Computes the measurement update with a provided real observation and computed observation.
    ///
    /// The nominal state is the state used for the computed observation.
    /// The real observation is the observation that was actually measured.
    /// The computed observation is the observation that was computed from the nominal state.
    ///
    /// Returns the updated estimate and the residual. The residual may be zero if the residual ratio check prevented the ingestion of this measurement.
    ///
    /// # Arguments
    ///
    /// * `nominal_state`: the nominal state at which the observation was computed.
    /// * `real_obs`: the real observation that was measured.
    /// * `computed_obs`: the computed observation from the nominal state.
    /// * `measurement_covar`: the measurement covariance associated with this time update (i./e. the square of the standard deviation)
    /// * `h_tilder`: the sensitivity matrix that maps the measurement(s) into the estimate state space
    /// * `resid_rejection`: the automatic residual rejection criteria, if enabled.
    fn measurement_update(
        &mut self,
        nominal_state: T,
        real_obs: OVector<f64, M>,
        computed_obs: OVector<f64, M>,
        measurement_covar: OMatrix<f64, M, M>,
        h_tilde: OMatrix<f64, M, <T as State>::Size>,
        resid_rejection: Option<ResidRejectCrit>,
    ) -> Result<(Self::Estimate, Residual<M>), ODError>;

    /// Returns whether the filter is an extended filter (e.g. EKF)
    fn is_extended(&self) -> bool;

    /// Sets the filter to be extended or not depending on the value of status
    fn set_extended(&mut self, status: bool);

    /// Sets the process noise matrix of the estimated state
    fn set_process_noise(&mut self, snc: ProcessNoise<A>);
}
