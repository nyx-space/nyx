use crate::celestia::{Cosm, Frame, FrameInfo, State};
use crate::dynamics::spacecraft::SpacecraftState;
use crate::utils::between_pm_180;
use std::fmt;

/// A general Event
pub trait Event: fmt::Debug {
    /// Defines the type which will be accepted by the condition
    type StateType: Copy;

    // Evaluation of event crossing, must return whether the condition happened between between both states.
    fn eval_crossing(&self, prev_state: &Self::StateType, next_state: &Self::StateType) -> bool;

    // Evaluation of the event, must return a value corresponding to whether the state is before or after the event
    fn eval(&self, state: &Self::StateType) -> f64;
}

/// A tracker for events during the propagation. Attach it directly to the propagator.
#[derive(Debug)]
pub struct EventTrackers<S: Copy> {
    pub events: Vec<Box<dyn Event<StateType = S>>>,
    pub found_bounds: Vec<Vec<(f64, f64)>>,
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

/// Built-in events, will likely be expanded as development continues.
#[derive(Clone, Copy, Debug)]
pub enum EventKind {
    Sma(f64),
    Ecc(f64),
    Inc(f64),
    Raan(f64),
    Aop(f64),
    TA(f64),
    Periapse,
    Apoapse,
    Fuel(f64),
}

impl fmt::Display for EventKind {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// An orbital event, in the same frame or in another frame.
#[derive(Debug)]
pub struct OrbitalEvent<'a> {
    pub kind: EventKind,
    pub tgt: Option<FrameInfo>,
    pub cosm: Option<&'a Cosm<'a>>,
}

impl<'a> OrbitalEvent<'a> {
    pub fn new(kind: EventKind) -> Box<Self> {
        Box::new(OrbitalEvent {
            kind,
            tgt: None,
            cosm: None,
        })
    }
    pub fn in_frame(kind: EventKind, tgt: FrameInfo, cosm: &'a Cosm<'a>) -> Box<Self> {
        Box::new(OrbitalEvent {
            kind,
            tgt: Some(tgt),
            cosm: Some(cosm),
        })
    }
}

impl<'a> Event for OrbitalEvent<'a> {
    type StateType = State;

    fn eval(&self, state: &Self::StateType) -> f64 {
        let state = match self.tgt {
            Some(tgt) => self.cosm.unwrap().frame_chg(state, tgt),
            None => *state,
        };

        match self.kind {
            EventKind::Sma(sma) => state.sma() - sma,
            EventKind::Ecc(ecc) => state.ecc() - ecc,
            EventKind::Inc(inc) => state.inc() - inc,
            EventKind::Raan(raan) => state.raan() - raan,
            EventKind::Aop(aop) => state.aop() - aop,
            EventKind::TA(angle) => state.ta() - angle,
            EventKind::Periapse => state.ta(),
            EventKind::Apoapse => between_pm_180(state.ta()),
            _ => panic!("event {:?} not supported", self.kind),
        }
    }

    fn eval_crossing(&self, prev_state: &Self::StateType, next_state: &Self::StateType) -> bool {
        let prev_val = self.eval(prev_state);
        let next_val = self.eval(next_state);
        match self.kind {
            // XXX: Should this condition be applied to all angles?
            EventKind::Periapse | EventKind::Apoapse => prev_val > next_val,
            _ => self.eval(prev_state) * self.eval(next_state) <= 0.0,
        }
    }
}

#[derive(Debug)]
pub struct SCEvent<'a> {
    pub kind: EventKind,
    pub orbital: Option<Box<OrbitalEvent<'a>>>,
}

impl<'a> SCEvent<'a> {
    pub fn fuel_mass(mass: f64) -> Box<Self> {
        Box::new(Self {
            kind: EventKind::Fuel(mass),
            orbital: None,
        })
    }
    pub fn orbital(event: Box<OrbitalEvent<'a>>) -> Box<Self> {
        Box::new(Self {
            kind: event.kind,
            orbital: Some(event),
        })
    }
}

impl<'a> Event for SCEvent<'a> {
    type StateType = SpacecraftState;

    fn eval(&self, state: &Self::StateType) -> f64 {
        match self.kind {
            EventKind::Fuel(mass) => state.fuel_mass - mass,
            _ => self.orbital.as_ref().unwrap().eval(&state.orbit),
        }
    }

    fn eval_crossing(&self, prev_state: &Self::StateType, next_state: &Self::StateType) -> bool {
        match self.kind {
            EventKind::Fuel(mass) => prev_state.fuel_mass <= mass && next_state.fuel_mass > mass,
            _ => self
                .orbital
                .as_ref()
                .unwrap()
                .eval_crossing(&prev_state.orbit, &next_state.orbit),
        }
    }
}

/// A condition to stop a propagator.
/// Note: min_step of propagator options will guide how precise the solution can be!
#[derive(Debug)]
pub struct StopCondition<S: Copy> {
    /// Set to a negative number to search backward
    pub max_prop_time: f64,
    /// The event which should be the stopping condition
    pub event: Box<dyn Event<StateType = S>>,
    /// The number of times the event must be hit prior to stopping (should be at least 1)
    pub trigger: usize,
    /// Maximum number of iterations of the Brent solver.
    pub max_iter: usize,
    /// Maximum error in the event, used as convergence criteria.
    pub epsilon: f64,
}

impl<S: Copy> StopCondition<S> {
    /// Finds the closest time at which this condition is met. Stops on first occurence.
    pub fn new(event: Box<dyn Event<StateType = S>>, prop_time: f64, epsilon: f64) -> Self {
        Self {
            max_prop_time: prop_time,
            event,
            trigger: 1,
            max_iter: 50,
            epsilon,
        }
    }

    /// Finds the closest time at which this condition is met. Stops on `hits` occurence (must be strictly greater than 1)
    pub fn after_hits(
        event: Box<dyn Event<StateType = S>>,
        hits: usize,
        prop_time: f64,
        epsilon: f64,
    ) -> Self {
        assert!(hits >= 1, "cannot stop on zero-th event passing");
        Self {
            max_prop_time: prop_time,
            event,
            trigger: hits,
            max_iter: 50,
            epsilon,
        }
    }
}

/// Built-in events, will likely be expanded as development continues.
#[derive(Clone, Copy, Debug)]
pub enum ConvergenceError {
    /// Event not triggered in max prop time
    NeverTriggered,
    /// Event not hit enough times (requested, found)
    UnsufficientTriggers(usize, usize),
    /// Value corresponds to the number of iterations used
    MaxIterReached(usize),
}
