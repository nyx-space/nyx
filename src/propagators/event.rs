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

use super::PropagationError;
use crate::dynamics::Dynamics;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::md::trajectory::{Interpolatable, Traj};
use crate::propagators::{PropAlmanacSnafu, PropAnalysisSnafu, TrajectoryEventSnafu};
use crate::time::{Duration, Epoch};
use crate::State;
use anise::analysis::event::Event;
use anise::analysis::{brent_solver, AnalysisError};
use anise::frames::Frame;
use log::info;
use rayon::iter::ParallelBridge;
use rayon::prelude::ParallelIterator;
use snafu::ResultExt;
use std::f64;
use std::sync::mpsc::channel;

use super::PropInstance;

impl<D: Dynamics> PropInstance<'_, D>
where
    DefaultAllocator: Allocator<<D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::VecLength>,
{
    /// Propagates the dynamics until the specified event has occurred one, or until `max_duration` is reached. Refer to [until_nth_event] for details.
    pub fn until_event(
        &mut self,
        max_duration: Duration,
        event: &Event,
        event_frame: Option<Frame>,
    ) -> Result<(D::StateType, Traj<D::StateType>), PropagationError>
    where
        <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
        D::StateType: Interpolatable,
    {
        self.until_nth_event(max_duration, event, event_frame, 1)
    }

    /// Propagates the dynamics until the specified event has occurred `trigger` times, or until `max_duration` is reached.
    ///
    /// This method monitors the provided `event` during propagation. Once the event condition is met
    /// `trigger` number of times (e.g., set `trigger` to 1 for the first occurrence), the propagation stops
    /// at the end of that integration step.
    ///
    /// A root-finding algorithm (Brent's method) is then used to locate the exact time of the event
    /// within the final integration step. The returned state corresponds to this precise event time,
    /// interpolated from the trajectory.
    ///
    /// # Arguments
    ///
    /// * `max_duration` - The maximum duration to propagate if the event is not triggered the requested number of times.
    /// * `event` - The event definition (scalar expression and condition) to monitor.
    /// * `trigger` - The 1-based index of the event occurrence to stop at (e.g. 1 for the first crossing, 2 for the second).
    ///
    /// # Returns
    ///
    /// A tuple containing:
    /// 1. The interpolated state exactly at the moment the $n$-th event occurred.
    /// 2. The full trajectory recorded up to the end of the propagation step where the event occurred.
    ///
    /// # Errors
    ///
    /// * `PropagationError::NthEventError`: Returned if `max_duration` is reached before the event was triggered `trigger` times.
    /// * `PropagationError::TrajectoryEvent`: Returned if the interpolation of the event state fails.
    /// * `PropagationError::Analysis`: Returned if the event evaluation fails during the search.
    pub fn until_nth_event(
        &mut self,
        max_duration: Duration,
        event: &Event,
        event_frame: Option<Frame>,
        trigger: usize,
    ) -> Result<(D::StateType, Traj<D::StateType>), PropagationError>
    where
        <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
        D::StateType: Interpolatable,
    {
        info!("Propagating until {event} or {max_duration}");

        let mut crossing_counts = 0;
        let closure_almanac = self.almanac.clone();
        let orbit = if let Some(observer_frame) = event_frame {
            self.almanac
                .transform_to(self.state.orbit(), observer_frame, None)
                .context(PropAlmanacSnafu)?
        } else {
            self.state.orbit()
        };

        let mut y_prev = event
            .eval(orbit, &self.almanac)
            .context(PropAnalysisSnafu)?;

        let enough_crossings = |next_state: D::StateType| -> Result<bool, PropagationError> {
            let orbit = if let Some(observer_frame) = event_frame {
                closure_almanac
                    .transform_to(next_state.orbit(), observer_frame, None)
                    .context(PropAlmanacSnafu)?
            } else {
                next_state.orbit()
            };

            let y_next = event
                .eval(orbit, &closure_almanac)
                .context(PropAnalysisSnafu)?;

            let delta = (y_next - y_prev).abs();

            if event.scalar.is_angle() {
                // Atan2 is a triangular signal so a bracket exists only if y_prev is negative and y_next is positive.
                // Anything else is a fluke, and we can quickly speed through the whole trajectory.
                if y_prev.signum() != y_next.signum() && delta < 180.0 {
                    crossing_counts += 1;
                }
            } else {
                // Previous step accepted, check if there was a zero crossing.
                if y_prev * y_next < 0.0 {
                    crossing_counts += 1;
                }
            }
            y_prev = y_next;

            Ok(crossing_counts >= trigger)
        };

        // Build the trajectory during the propagation.
        let end_state;
        let mut traj = Traj::new();
        let start_state = self.state;

        let rx = {
            // Channels that have a single state for the propagator
            let (tx, rx) = channel();
            // Propagate the dynamics
            end_state = self.propagate(max_duration, Some(tx), Some(enough_crossings))?;
            // Note that the end state is also sent on the channel before the return of this function.
            rx
        };

        traj.states = rx.into_iter().par_bridge().collect();
        // Push the start state -- will be reordered in the finalize call.
        // For some reason, this must happen at the end -- can't figure out why.
        traj.states.push(start_state);

        traj.finalize();

        // We pushed at least one item to this trajectory, so the following unwrap will always succeed.
        let last_traj_state = traj.states.last().cloned().unwrap();
        if end_state == last_traj_state {
            return Err(PropagationError::NthEventError {
                nth: trigger,
                found: crossing_counts,
            });
        }

        // Add the final state to the trajectory so we can search it in full.
        traj.states.push(end_state);

        // We have one bracket to search.
        let traj_at = |epoch: Epoch| -> Result<f64, AnalysisError> {
            let state = traj
                .at(epoch)
                .map_err(|e| AnalysisError::GenericAnalysisError {
                    err: format!("{e}"),
                })?;

            event.eval(state.orbit(), &self.almanac)
        };

        // Else Brent solver between these states

        let event_epoch = brent_solver(traj_at, event, last_traj_state.epoch(), end_state.epoch())
            .context(PropAnalysisSnafu)?;

        Ok((traj.at(event_epoch).context(TrajectoryEventSnafu)?, traj))
    }
}
