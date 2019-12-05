use crate::celestia::{Geoid, State};
use crate::utils::between_pm_180;
use std::fmt::Debug;

pub trait Event: Debug {
    /// Defines the type which will be accepted by the condition
    type StateType;

    // Evaluation of the event, must return whether the event value at the provided state.
    fn eval(&self, state: &Self::StateType) -> f64;
}

#[derive(Debug)]
pub struct EventTrackers<S> {
    pub events: Vec<Box<dyn Event<StateType = S>>>,
    pub found_bounds: Vec<Vec<(f64, f64)>>,
    prev_values: Vec<f64>,
}

impl<S> EventTrackers<S> {
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
        Self {
            events,
            prev_values: Vec::with_capacity(len),
            found_bounds: Vec::with_capacity(len),
        }
    }

    pub fn eval_and_save(&mut self, prev_time: f64, next_time: f64, state: &S) {
        for event_no in 0..self.events.len() {
            // Evaluate the event
            let val = self.events[event_no].eval(state);
            if !self.prev_values.is_empty() {
                // Check if we've changed the sign of the event
                if self.prev_values[event_no].signum() * val < 0.0 {
                    // The previous value and the new value have opposite signs
                    // Let's store this as an event passage
                    self.found_bounds[event_no].push((prev_time, next_time));
                }
                self.prev_values[event_no] = val;
            } else {
                self.prev_values.push(val);
            }
        }
    }
}

#[derive(Debug)]
pub enum StateEventKind {
    TA(f64),
}

#[derive(Debug)]
pub struct StateEvent {
    pub kind: StateEventKind,
}

impl Event for StateEvent {
    type StateType = State<Geoid>;

    fn eval(&self, state: &Self::StateType) -> f64 {
        match self.kind {
            StateEventKind::TA(angle) => between_pm_180(state.ta()) - angle,
        }
    }
}
