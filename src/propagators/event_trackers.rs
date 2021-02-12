pub use crate::md::events::{Event, EventKind, OrbitalEvent, SCEvent};
use crate::time::{Duration, Epoch, TimeUnit};
use std::fmt;

/// A tracker for events during the propagation. Attach it directly to the propagator.
#[derive(Debug)]
pub struct EventTrackers<S: Copy> {
    pub events: Vec<Box<dyn Event<StateType = S>>>,
    pub found_bounds: Vec<Vec<(Epoch, Epoch)>>,
    prev_values: Vec<S>,
}

impl<S: Copy> EventTrackers<S> {
    /// Used to initialize no event trackers. Should not be needed publicly.
    pub fn none() -> Self {
        Self {
            events: Vec::with_capacity(0),
            prev_values: Vec::with_capacity(0),
            found_bounds: Vec::with_capacity(0),
        }
    }

    /// Track only one event
    pub fn from_event(event: Box<dyn Event<StateType = S>>) -> Self {
        Self {
            events: vec![event],
            prev_values: Vec::with_capacity(1),
            found_bounds: vec![Vec::new()],
        }
    }

    /// Track several events
    pub fn from_events(events: Vec<Box<dyn Event<StateType = S>>>) -> Self {
        let len = events.len();
        let mut found_bounds = Vec::new();
        for _ in 0..len {
            found_bounds.push(Vec::new());
        }
        Self {
            events,
            prev_values: Vec::with_capacity(len),
            found_bounds,
        }
    }

    /// Evaluate whether we have crossed the boundary of an event
    pub fn eval_and_save(&mut self, prev_time: Epoch, next_time: Epoch, state: &S) {
        for event_no in 0..self.events.len() {
            if self.prev_values.len() > event_no {
                // Evaluate the event crossing
                if self.events[event_no].eval_crossing(&self.prev_values[event_no], &state) {
                    // Append the crossing times
                    self.found_bounds[event_no].push((prev_time, next_time));
                }
                self.prev_values[event_no] = *state;
            } else {
                self.prev_values.push(*state);
            }
        }
    }

    pub fn reset(&mut self) {
        for event_no in 0..self.events.len() {
            while !self.found_bounds[event_no].is_empty() {
                self.found_bounds[event_no].remove(0);
            }
        }
    }
}

impl<S: Copy> fmt::Display for EventTrackers<S> {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for event_no in 0..self.events.len() {
            if event_no > 0 {
                writeln!(f)?;
            }
            if !self.found_bounds[event_no].is_empty() {
                let last_e = self.found_bounds[event_no][self.found_bounds[event_no].len() - 1];
                write!(
                    f,
                    "[ OK  ] Event {:?} converged on ({}, {})",
                    self.events[event_no], last_e.0, last_e.1,
                )?;
            } else {
                write!(
                    f,
                    "[ERROR] Event {:?} did NOT converge",
                    self.events[event_no]
                )?;
            }
        }
        Ok(())
    }
}

/// A condition to stop a propagator.
/// Note: min_step of propagator options will guide how precise the solution can be!
#[derive(Debug)]
pub struct StopCondition<S: Copy> {
    /// Set to a negative number to search backward
    pub max_prop_time: Duration,
    /// The event which should be the stopping condition
    pub event: Box<dyn Event<StateType = S>>,
    /// The number of times the event must be hit prior to stopping (should be at least 1)
    pub trigger: usize,
    /// Maximum number of iterations of the Brent solver.
    pub max_iter: usize,
    /// Maximum error in the event, used as time convergence criteria, defaults to one second
    pub epsilon: Duration,
    /// Maximum error in the evaluation of the event (e.g. 0.1 )
    pub epsilon_eval: f64,
}

#[allow(clippy::identity_op)]
impl<S: Copy> StopCondition<S> {
    /// Finds the closest time at which this condition is met. Stops on first occurence.
    pub fn new(event: Box<dyn Event<StateType = S>>, prop_time: Duration, epsilon: f64) -> Self {
        Self {
            max_prop_time: prop_time,
            event,
            trigger: 1,
            max_iter: 50,
            epsilon: 1 * TimeUnit::Second,
            epsilon_eval: epsilon,
        }
    }

    /// Finds the closest time at which this condition is met. Stops on `hits` occurence (must be strictly greater than 1)
    pub fn after_hits(
        event: Box<dyn Event<StateType = S>>,
        hits: usize,
        prop_time: Duration,
        epsilon: f64,
    ) -> Self {
        assert!(hits >= 1, "cannot stop on zero-th event passing");
        Self {
            max_prop_time: prop_time,
            event,
            trigger: hits,
            max_iter: 50,
            epsilon: 1 * TimeUnit::Second,
            epsilon_eval: epsilon,
        }
    }
}
