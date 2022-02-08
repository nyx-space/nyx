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

/// A structure of Monte Carlo results allowing for filtering
pub struct McResults<S: InterpState>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
    <DefaultAllocator as Allocator<f64, S::VecLength>>::Buffer: Send,
{
    /// Raw data from each run
    pub data: Vec<(usize, Result<(S, Traj<S>), NyxError>)>,
}

impl<S: InterpState> McResults<S>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
    <DefaultAllocator as Allocator<f64, S::VecLength>>::Buffer: Send,
{
    /// Returns the value of the requested state parameter for all trajectories from `start` to `end` every `step` and
    /// using the value of `value_if_run_failed` if set and skipping that run if the run failed
    pub fn report_every_between(
        &self,
        param: StateParameter,
        step: Duration,
        start: Epoch,
        end: Epoch,
        value_if_run_failed: Option<f64>,
    ) -> Vec<f64> {
        let mut report = Vec::with_capacity(self.data.len());
        for (run_num, run) in &self.data {
            match run {
                Ok((_, traj)) => {
                    for state in traj.every_between(step, start, end) {
                        match state.value(&param) {
                            Ok(val) => report.push(val),
                            Err(e) => match value_if_run_failed {
                                Some(val) => report.push(val),
                                None => {
                                    warn!("run #{}: {}, skipping {} in report", run_num, e, param)
                                }
                            },
                        }
                    }
                }
                Err(e) => match value_if_run_failed {
                    Some(val) => report.push(val),
                    None => warn!(
                        "run #{} failed with {}, skipping {} in report",
                        run_num, e, param
                    ),
                },
            }
        }
        report
    }

    /// Returns the value of the requested state parameter for all trajectories from the start to the end of each trajectory and
    /// using the value of `value_if_run_failed` if set and skipping that run if the run failed
    pub fn report_every(
        &self,
        param: StateParameter,
        step: Duration,
        value_if_run_failed: Option<f64>,
    ) -> Vec<f64> {
        let mut report = Vec::with_capacity(self.data.len());
        for (run_num, run) in &self.data {
            match run {
                Ok((_, traj)) => {
                    for state in traj.every(step) {
                        match state.value(&param) {
                            Ok(val) => report.push(val),
                            Err(e) => match value_if_run_failed {
                                Some(val) => report.push(val),
                                None => {
                                    warn!("run #{}: {}, skipping {} in report", run_num, e, param)
                                }
                            },
                        }
                    }
                }
                Err(e) => match value_if_run_failed {
                    Some(val) => report.push(val),
                    None => warn!(
                        "run #{} failed with {}, skipping {} in report",
                        run_num, e, param
                    ),
                },
            }
        }
        report
    }

    /// Returns the value of the requested state parameter for the first state and
    /// using the value of `value_if_run_failed` if set and skipping that run if the run failed
    pub fn report_first(
        &self,
        param: StateParameter,
        value_if_run_failed: Option<f64>,
    ) -> Vec<f64> {
        let mut report = Vec::with_capacity(self.data.len());
        for (run_num, run) in &self.data {
            match run {
                Ok((_, traj)) => match traj.first().value(&param) {
                    Ok(val) => report.push(val),
                    Err(e) => match value_if_run_failed {
                        Some(val) => report.push(val),
                        None => {
                            warn!("run #{}: {}, skipping {} in report", run_num, e, param)
                        }
                    },
                },
                Err(e) => match value_if_run_failed {
                    Some(val) => report.push(val),
                    None => warn!(
                        "run #{} failed with {}, skipping {} in report",
                        run_num, e, param
                    ),
                },
            }
        }
        report
    }

    /// Returns the value of the requested state parameter for the first state and
    /// using the value of `value_if_run_failed` if set and skipping that run if the run failed
    pub fn report_last(&self, param: StateParameter, value_if_run_failed: Option<f64>) -> Vec<f64> {
        let mut report = Vec::with_capacity(self.data.len());
        for (run_num, run) in &self.data {
            match run {
                Ok((_, traj)) => match traj.last().value(&param) {
                    Ok(val) => report.push(val),
                    Err(e) => match value_if_run_failed {
                        Some(val) => report.push(val),
                        None => {
                            warn!("run #{}: {}, skipping {} in report", run_num, e, param)
                        }
                    },
                },
                Err(e) => match value_if_run_failed {
                    Some(val) => report.push(val),
                    None => warn!(
                        "run #{} failed with {}, skipping {} in report",
                        run_num, e, param
                    ),
                },
            }
        }
        report
    }
}
