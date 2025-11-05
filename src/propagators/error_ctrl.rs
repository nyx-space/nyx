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

use serde::{Deserialize, Serialize};

use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OVector, U3};

// This determines when to take into consideration the magnitude of the state_delta and
// prevents dividing by too small of a number.
const REL_ERR_THRESH: f64 = 0.1;

/// Manages how a propagator computes the error in the current step.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub enum ErrorControl {
    /// An RSS state error control for Cartesian states (e.g., position and velocity).
    RSSCartesianState,
    /// An RSS step error control for Cartesian states (e.g., position and velocity).
    RSSCartesianStep,
    /// An RSS state error control. Recommended for high-accuracy simulations.
    RSSState,
    /// An RSS step error control. Should only be used with state slices that have the same units.
    RSSStep,
    /// Computes the largest error at each component. Not recommended for states with mixed units.
    LargestError,
    /// A largest state error control.
    LargestState,
    /// A largest step error control. Should only be used with state slices that have the same units.
    LargestStep,
}

impl ErrorControl {
    /// Computes the error of the current step.
    pub fn estimate<N: DimName>(
        self,
        error_est: &OVector<f64, N>,
        candidate: &OVector<f64, N>,
        cur_state: &OVector<f64, N>,
    ) -> f64
    where
        DefaultAllocator: Allocator<N>,
    {
        match self {
            ErrorControl::RSSCartesianState => {
                if N::dim() >= 6 {
                    let err_radius = RSSState::estimate::<U3>(
                        &error_est.fixed_rows::<3>(0).into_owned(),
                        &candidate.fixed_rows::<3>(0).into_owned(),
                        &cur_state.fixed_rows::<3>(0).into_owned(),
                    );
                    let err_velocity = RSSState::estimate::<U3>(
                        &error_est.fixed_rows::<3>(3).into_owned(),
                        &candidate.fixed_rows::<3>(3).into_owned(),
                        &cur_state.fixed_rows::<3>(3).into_owned(),
                    );
                    err_radius.max(err_velocity)
                } else {
                    RSSStep::estimate(error_est, candidate, cur_state)
                }
            }
            ErrorControl::RSSCartesianStep => {
                if N::dim() >= 6 {
                    let err_radius = RSSStep::estimate::<U3>(
                        &error_est.fixed_rows::<3>(0).into_owned(),
                        &candidate.fixed_rows::<3>(0).into_owned(),
                        &cur_state.fixed_rows::<3>(0).into_owned(),
                    );
                    let err_velocity = RSSStep::estimate::<U3>(
                        &error_est.fixed_rows::<3>(3).into_owned(),
                        &candidate.fixed_rows::<3>(3).into_owned(),
                        &cur_state.fixed_rows::<3>(3).into_owned(),
                    );
                    err_radius.max(err_velocity)
                } else {
                    RSSStep::estimate(error_est, candidate, cur_state)
                }
            }
            ErrorControl::RSSState => {
                let mag = 0.5 * (candidate + cur_state).norm();
                let err = error_est.norm();
                if mag > REL_ERR_THRESH {
                    err / mag
                } else {
                    err
                }
            }
            ErrorControl::RSSStep => {
                let mag = (candidate - cur_state).norm();
                let err = error_est.norm();
                if mag > REL_ERR_THRESH.sqrt() {
                    err / mag
                } else {
                    err
                }
            }
            ErrorControl::LargestError => {
                let state_delta = candidate - cur_state;
                let mut max_err = 0.0;
                for (i, prop_err_i) in error_est.iter().enumerate() {
                    let err = if state_delta[i] > REL_ERR_THRESH {
                        (prop_err_i / state_delta[i]).abs()
                    } else {
                        prop_err_i.abs()
                    };
                    if err > max_err {
                        max_err = err;
                    }
                }
                max_err
            }
            ErrorControl::LargestState => {
                let sum_state = candidate + cur_state;
                let mut mag = 0.0f64;
                let mut err = 0.0f64;
                for i in 0..N::dim() {
                    mag += 0.5 * sum_state[i].abs();
                    err += error_est[i].abs();
                }

                if mag > REL_ERR_THRESH {
                    err / mag
                } else {
                    err
                }
            }
            ErrorControl::LargestStep => {
                let state_delta = candidate - cur_state;
                let mut mag = 0.0f64;
                let mut err = 0.0f64;
                for i in 0..N::dim() {
                    mag += state_delta[i].abs();
                    err += error_est[i].abs();
                }

                if mag > REL_ERR_THRESH {
                    err / mag
                } else {
                    err
                }
            }
        }
    }
}

impl Default for ErrorControl {
    fn default() -> Self {
        Self::RSSCartesianStep
    }
}

#[derive(Clone, Copy)]
#[allow(clippy::upper_case_acronyms)]
struct RSSStep;
impl RSSStep {
    fn estimate<N: DimName>(
        error_est: &OVector<f64, N>,
        candidate: &OVector<f64, N>,
        cur_state: &OVector<f64, N>,
    ) -> f64
    where
        DefaultAllocator: Allocator<N>,
    {
        let mag = (candidate - cur_state).norm();
        let err = error_est.norm();
        if mag > REL_ERR_THRESH.sqrt() {
            err / mag
        } else {
            err
        }
    }
}

#[derive(Clone, Copy)]
#[allow(clippy::upper_case_acronyms)]
struct RSSState;
impl RSSState {
    fn estimate<N: DimName>(
        error_est: &OVector<f64, N>,
        candidate: &OVector<f64, N>,
        cur_state: &OVector<f64, N>,
    ) -> f64
    where
        DefaultAllocator: Allocator<N>,
    {
        let mag = 0.5 * (candidate + cur_state).norm();
        let err = error_est.norm();
        if mag > REL_ERR_THRESH {
            err / mag
        } else {
            err
        }
    }
}
