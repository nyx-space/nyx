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

use crate::errors::NyxError;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::md::prelude::{Interpolatable, Traj};
use crate::md::EventEvaluator;
use crate::time::Duration;
use core::fmt;

/// Enumerates the possible edges of an event in a trajectory.
///
/// `EventEdge` is used to describe the nature of a trajectory event, particularly in terms of its temporal dynamics relative to a specified condition or threshold. This enum helps in distinguishing whether the event is occurring at a rising edge, a falling edge, or if the edge is unclear due to insufficient data or ambiguous conditions.
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum EventEdge {
    /// Represents a rising edge of the event. This indicates that the event is transitioning from a lower to a higher evaluation of the event. For example, in the context of elevation, a rising edge would indicate an increase in elevation from a lower angle.
    Rising,
    /// Represents a falling edge of the event. This is the opposite of the rising edge, indicating a transition from a higher to a lower value of the event evaluator. For example, if tracking the elevation of an object, a falling edge would signify a
    Falling,
    /// If the edge cannot be clearly defined, it will be marked as unclear. This happens if the event is at a saddle point and the epoch precision is too large to find the exact slope.
    Unclear,
}

/// Represents the details of an event occurring along a trajectory.
///
/// `EventDetails` encapsulates the state at which a particular event occurs in a trajectory, along with additional information about the nature of the event. This struct is particularly useful for understanding the dynamics of the event, such as whether it represents a rising or falling edge, or if the edge is unclear.
///
/// # Generics
/// S: Interpolatable - A type that represents the state of the trajectory. This type must implement the `Interpolatable` trait, ensuring that it can be interpolated and manipulated according to the trajectory's requirements.
#[derive(Clone, Debug, PartialEq)]
pub struct EventDetails<S: Interpolatable>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    /// The state of the trajectory at the found event.
    pub state: S,
    /// Indicates whether the event is a rising edge, falling edge, or unclear. This helps in understanding the direction of change at the event point.
    pub edge: EventEdge,
    /// Numerical evaluation of the event condition, e.g. if seeking the apoapsis, this returns the near zero
    pub value: f64,
    /// Numertical evaluation of the event condition one epoch step before the found event (used to compute the rising/falling edge).
    pub prev_value: Option<f64>,
    /// Numertical evaluation of the event condition one epoch step after the found event (used to compute the rising/falling edge).
    pub next_value: Option<f64>,
    /// Precision of the epoch for this value
    pub pm_duration: Duration,
    // Store the representation of this event as a string because we can't move or clone the event reference
    pub repr: String,
}

impl<S: Interpolatable> EventDetails<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    /// Generates detailed information about an event at a specific epoch in a trajectory.
    ///
    /// This takes an `Epoch` as an input and returns a `Result<Self, NyxError>`.
    /// It is designed to determine the state of a trajectory at a given epoch, evaluate the specific event at that state, and ascertain the nature of the event (rising, falling, or unclear).
    /// The initialization intelligently determines the edge type of the event by comparing the event's value at the current, previous, and next epochs.
    /// It ensures robust event characterization in trajectories.
    ///
    /// # Returns
    /// - `Ok(EventDetails<S>)` if the state at the given epoch can be determined and the event details are successfully evaluated.
    /// - `Err(NyxError)` if there is an error in retrieving the state at the specified epoch.
    ///
    pub fn new<E: EventEvaluator<S>>(
        state: S,
        value: f64,
        event: &E,
        traj: &Traj<S>,
    ) -> Result<Self, NyxError> {
        let epoch = state.epoch();
        let prev_value = if let Ok(state) = traj.at(epoch - event.epoch_precision()) {
            Some(event.eval(&state))
        } else {
            None
        };

        let next_value = if let Ok(state) = traj.at(epoch + event.epoch_precision()) {
            Some(event.eval(&state))
        } else {
            None
        };

        let edge = if let Some(prev_value) = prev_value {
            if let Some(next_value) = next_value {
                if prev_value > value && value > next_value {
                    EventEdge::Falling
                } else if prev_value < value && value < next_value {
                    EventEdge::Rising
                } else {
                    warn!("could not determine edge of {} at {}", event, state.epoch(),);
                    EventEdge::Unclear
                }
            } else if prev_value > value {
                EventEdge::Falling
            } else {
                EventEdge::Rising
            }
        } else if let Some(next_value) = next_value {
            if next_value > value {
                EventEdge::Rising
            } else {
                EventEdge::Falling
            }
        } else {
            warn!(
                "could not determine edge of {} because trajectory could be queried around {}",
                event,
                state.epoch()
            );
            EventEdge::Unclear
        };

        Ok(EventDetails {
            edge,
            state,
            value,
            prev_value,
            next_value,
            pm_duration: event.epoch_precision(),
            repr: event.eval_string(&state).to_string(),
        })
    }
}

impl<S: Interpolatable> fmt::Display for EventDetails<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let prev_fmt = match self.prev_value {
            Some(value) => format!("{value:.6}"),
            None => "".to_string(),
        };

        let next_fmt = match self.next_value {
            Some(value) => format!("{value:.6}"),
            None => "".to_string(),
        };

        write!(
            f,
            "{} and is {:?} (roots with {} intervals: {}, {:.6}, {})",
            self.repr, self.edge, self.pm_duration, prev_fmt, self.value, next_fmt
        )
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct EventArc<S: Interpolatable>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    pub rise: EventDetails<S>,
    pub fall: EventDetails<S>,
}

impl<S: Interpolatable> fmt::Display for EventArc<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} until {} (lasts {})",
            self.rise,
            self.fall.state.epoch(),
            self.fall.state.epoch() - self.rise.state.epoch()
        )
    }
}
