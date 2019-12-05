use crate::celestia::{Geoid, State};
use crate::utils::between_pm_180;
use std::fmt;

pub trait Event: fmt::Debug {
    /// Defines the type which will be accepted by the condition
    type StateType: Copy;

    // Evaluation of the event, must return whether the condition happened between between both states.
    fn eval_crossing(&self, prev_state: &Self::StateType, next_state: &Self::StateType) -> bool;
}

#[derive(Debug)]
pub struct EventTrackers<S: Copy> {
    pub events: Vec<Box<dyn Event<StateType = S>>>,
    pub found_bounds: Vec<Vec<(f64, f64)>>,
    prev_values: Vec<S>,
}

impl<S: Copy> EventTrackers<S> {
    pub fn none() -> Self {
        Self {
            events: Vec::with_capacity(0),
            prev_values: Vec::with_capacity(0),
            found_bounds: Vec::with_capacity(0),
        }
    }

    pub fn from_event(event: Box<dyn Event<StateType = S>>) -> Self {
        Self {
            events: vec![event],
            prev_values: Vec::with_capacity(1),
            found_bounds: Vec::with_capacity(1),
        }
    }

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

    pub fn eval_and_save(&mut self, prev_time: f64, next_time: f64, state: &S) {
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
}

impl<S: Copy> fmt::Display for EventTrackers<S> {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for event_no in 0..self.events.len() {
            if event_no > 0 {
                write!(f, "\n")?;
            }
            write!(
                f,
                "For {:?}, found times: {:?}",
                self.events[event_no], self.found_bounds[event_no]
            )?;
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug)]
pub enum StateEventKind {
    Periapse,
    Apoapse,
    TA(f64),
}

impl fmt::Display for StateEventKind {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Debug)]
pub struct OrbitalEvent {
    pub kind: StateEventKind,
}

impl Event for OrbitalEvent {
    type StateType = State<Geoid>;

    fn eval_crossing(&self, prev_state: &Self::StateType, next_state: &Self::StateType) -> bool {
        match self.kind {
            StateEventKind::TA(angle) => prev_state.ta() <= angle && next_state.ta() > angle,
            StateEventKind::Apoapse => {
                between_pm_180(prev_state.ta()) > between_pm_180(next_state.ta())
            }
            StateEventKind::Periapse => prev_state.ta() > next_state.ta(),
        }
    }
}
