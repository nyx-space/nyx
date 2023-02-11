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
extern crate rstats;

use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::md::trajectory::{InterpState, Traj};
use crate::md::StateParameter;
use crate::time::{Duration, Epoch};
use crate::NyxError;
pub use rstats::Stats;

use super::DispersedState;

/// A structure storing the result of a single Monte Carlo run
pub struct Run<S: InterpState, R>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
    <DefaultAllocator as Allocator<f64, S::VecLength>>::Buffer: Send,
{
    /// The index of this run
    pub index: usize,
    /// The original dispersed state
    pub dispersed_state: DispersedState<S>,
    /// The result from this run
    pub result: Result<R, NyxError>,
}

/// A structure of Monte Carlo results
pub struct Results<S: InterpState, R>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
    <DefaultAllocator as Allocator<f64, S::VecLength>>::Buffer: Send,
{
    /// Raw data from each run, sorted by run index for O(1) access to each run
    pub runs: Vec<Run<S, R>>,
    /// Name of this scenario
    pub scenario: String,
}

/// A structure that stores the result of a propagation segment of a Monte Carlo.
pub struct PropResult<S: InterpState>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
    <DefaultAllocator as Allocator<f64, S::VecLength>>::Buffer: Send,
{
    pub state: S,
    pub traj: Traj<S>,
}

impl<S: InterpState> Results<S, PropResult<S>>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
    <DefaultAllocator as Allocator<f64, S::VecLength>>::Buffer: Send,
{
    /// Returns the value of the requested state parameter for all trajectories from `start` to `end` every `step` and
    /// using the value of `value_if_run_failed` if set and skipping that run if the run failed
    pub fn every_value_of_between(
        &self,
        param: &StateParameter,
        step: Duration,
        start: Epoch,
        end: Epoch,
        value_if_run_failed: Option<f64>,
    ) -> Vec<f64> {
        let mut report = Vec::with_capacity(self.runs.len());
        for run in &self.runs {
            match &run.result {
                Ok(r) => {
                    for state in r.traj.every_between(step, start, end) {
                        match state.value(param) {
                            Ok(val) => report.push(val),
                            Err(e) => match value_if_run_failed {
                                Some(val) => report.push(val),
                                None => {
                                    warn!("run #{}: {}, skipping {} in report", run.index, e, param)
                                }
                            },
                        }
                    }
                }
                Err(e) => match value_if_run_failed {
                    Some(val) => report.push(val),
                    None => warn!(
                        "run #{} failed with {}, skipping {} in report",
                        run.index, e, param
                    ),
                },
            }
        }
        report
    }

    /// Returns the value of the requested state parameter for all trajectories from the start to the end of each trajectory and
    /// using the value of `value_if_run_failed` if set and skipping that run if the run failed
    pub fn every_value_of(
        &self,
        param: &StateParameter,
        step: Duration,
        value_if_run_failed: Option<f64>,
    ) -> Vec<f64> {
        let mut report = Vec::with_capacity(self.runs.len());
        for run in &self.runs {
            match &run.result {
                Ok(r) => {
                    for state in r.traj.every(step) {
                        match state.value(param) {
                            Ok(val) => report.push(val),
                            Err(e) => match value_if_run_failed {
                                Some(val) => report.push(val),
                                None => {
                                    warn!("run #{}: {}, skipping {} in report", run.index, e, param)
                                }
                            },
                        }
                    }
                }
                Err(e) => match value_if_run_failed {
                    Some(val) => report.push(val),
                    None => warn!(
                        "run #{} failed with {}, skipping {} in report",
                        run.index, e, param
                    ),
                },
            }
        }
        report
    }

    /// Returns the value of the requested state parameter for the first state and
    /// using the value of `value_if_run_failed` if set and skipping that run if the run failed
    pub fn first_values_of(
        &self,
        param: &StateParameter,
        value_if_run_failed: Option<f64>,
    ) -> Vec<f64> {
        let mut report = Vec::with_capacity(self.runs.len());
        for run in &self.runs {
            match &run.result {
                Ok(r) => match r.traj.first().value(param) {
                    Ok(val) => report.push(val),
                    Err(e) => match value_if_run_failed {
                        Some(val) => report.push(val),
                        None => {
                            warn!("run #{}: {}, skipping {} in report", run.index, e, param)
                        }
                    },
                },
                Err(e) => match value_if_run_failed {
                    Some(val) => report.push(val),
                    None => warn!(
                        "run #{} failed with {}, skipping {} in report",
                        run.index, e, param
                    ),
                },
            }
        }
        report
    }

    /// Returns the value of the requested state parameter for the first state and
    /// using the value of `value_if_run_failed` if set and skipping that run if the run failed
    pub fn last_values_of(
        &self,
        param: &StateParameter,
        value_if_run_failed: Option<f64>,
    ) -> Vec<f64> {
        let mut report = Vec::with_capacity(self.runs.len());
        for run in &self.runs {
            match &run.result {
                Ok(r) => match r.traj.last().value(param) {
                    Ok(val) => report.push(val),
                    Err(e) => match value_if_run_failed {
                        Some(val) => report.push(val),
                        None => {
                            warn!("run #{}: {}, skipping {} in report", run.index, e, param)
                        }
                    },
                },
                Err(e) => match value_if_run_failed {
                    Some(val) => report.push(val),
                    None => warn!(
                        "run #{} failed with {}, skipping {} in report",
                        run.index, e, param
                    ),
                },
            }
        }
        report
    }

    /// Returns the dispersion values of the requested state parameter
    pub fn dispersion_values_of(&self, param: &StateParameter) -> Result<Vec<f64>, NyxError> {
        let mut report = Vec::with_capacity(self.runs.len());
        'run_loop: for run in &self.runs {
            for (dparam, val) in &run.dispersed_state.actual_dispersions {
                if dparam == param {
                    report.push(*val);
                    continue 'run_loop;
                }
            }
            // Oh, this parameter was not found!
            return Err(NyxError::StateParameterUnavailable);
        }
        Ok(report)
    }
}
