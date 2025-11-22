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

/// The Error Control manages how a propagator computes the error in the current step.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum ErrorControl {
    /// An RSS state error control which effectively for the provided vector composed of two vectors of the same unit, both of size 3 (e.g. position + velocity).
    RSSCartesianState,
    /// An RSS step error control which effectively for the provided vector composed of two vectors of the same unit, both of size 3 (e.g. position + velocity).
    #[default]
    RSSCartesianStep,
    /// An RSS state error control: when in doubt, use this error controller, especially for high accurracy.
    ///
    /// Here is the warning from GMAT R2016a on this error controller:
    /// > This is a more stringent error control method than [`rss_step`] that is often used as the default in other software such as STK.
    /// > If you set [the] accuracy to a very small number, 1e-13 for example, and set the error control to [`rss_step`], integrator
    /// > performance will be poor, for little if any improvement in the accuracy of the orbit integration.
    /// > For more best practices of these integrators (which clone those in GMAT), please refer to the
    /// > [GMAT reference](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/doc/help/src/Resource_NumericalIntegrators.xml#L1292).
    /// > (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/ODEModel.cpp#L3004]
    RSSState,
    /// An RSS step error control which effectively computes the L2 norm of the provided Vector of size 3
    ///
    /// Note that this error controller should be preferably be used only with slices of a state with the same units.
    /// For example, one should probably use this for position independently of using it for the velocity.
    /// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/ODEModel.cpp#L3045]
    RSSStep,
    /// A largest error control which effectively computes the largest error at each component
    ///
    /// This is a standard error computation algorithm, but it's arguably bad if the state's components have different units.
    /// It calculates the largest local estimate of the error from the integration (`error_est`)
    /// given the difference in the candidate state and the previous state (`state_delta`).
    /// This error estimator is from the physical model estimator of GMAT
    /// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/PhysicalModel.cpp#L987]
    LargestError,
    /// A largest state error control
    ///
    /// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/ODEModel.cpp#L3018]
    LargestState,

    /// A largest step error control which effectively computes the L1 norm of the provided Vector of size 3
    ///
    /// Note that this error controller should be preferably be used only with slices of a state with the same units.
    /// For example, one should probably use this for position independently of using it for the velocity.
    /// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/ODEModel.cpp#L3033]
    LargestStep,
}

impl ErrorControl {
    /// Computes the actual error of the current step.
    ///
    /// The `error_est` is the estimated error computed from the difference in the two stages of
    /// of the RK propagator. The `candidate` variable is the candidate state, and `cur_state` is
    /// the current state. This function must return the error.
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

/// An RSS step error control which effectively computes the L2 norm of the provided Vector of size 3
///
/// Note that this error controller should be preferably be used only with slices of a state with the same units.
/// For example, one should probably use this for position independently of using it for the velocity.
/// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/ODEModel.cpp#L3045]
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

/// An RSS state error control: when in doubt, use this error controller, especially for high accurracy.
///
/// Here is the warning from GMAT R2016a on this error controller:
/// > This is a more stringent error control method than [`rss_step`] that is often used as the default in other software such as STK.
/// > If you set [the] accuracy to a very small number, 1e-13 for example, and set the error control to [`rss_step`], integrator
/// > performance will be poor, for little if any improvement in the accuracy of the orbit integration.
/// > For more best practices of these integrators (which clone those in GMAT), please refer to the
/// > [GMAT reference](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/doc/help/src/Resource_NumericalIntegrators.xml#L1292).
/// > (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/ODEModel.cpp#L3004]
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
