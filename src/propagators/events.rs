use crate::celestia::{Geoid, State};
use crate::dynamics::spacecraft::SpacecraftState;
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
            found_bounds: vec![Vec::new()],
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

#[derive(Debug)]
pub struct OrbitalEvent {
    pub kind: EventKind,
}

impl Event for OrbitalEvent {
    type StateType = State<Geoid>;

    fn eval_crossing(&self, prev_state: &Self::StateType, next_state: &Self::StateType) -> bool {
        match self.kind {
            EventKind::Sma(sma) => prev_state.sma() <= sma && next_state.sma() > sma,
            EventKind::Ecc(ecc) => prev_state.ecc() <= ecc && next_state.ecc() > ecc,
            EventKind::Inc(inc) => prev_state.inc() <= inc && next_state.inc() > inc,
            EventKind::Raan(raan) => prev_state.raan() <= raan && next_state.raan() > raan,
            EventKind::Aop(aop) => prev_state.aop() <= aop && next_state.aop() > aop,
            EventKind::TA(angle) => prev_state.ta() <= angle && next_state.ta() > angle,
            EventKind::Periapse => prev_state.ta() > next_state.ta(),
            EventKind::Apoapse => between_pm_180(prev_state.ta()) > between_pm_180(next_state.ta()),
            _ => panic!("event {:?} not supported"),
        }
    }
}

#[derive(Debug)]
pub struct SCEvent {
    pub kind: EventKind,
}

impl Event for SCEvent {
    type StateType = SpacecraftState;

    fn eval_crossing(&self, prev_state: &Self::StateType, next_state: &Self::StateType) -> bool {
        match self.kind {
            EventKind::Sma(sma) => prev_state.orbit.sma() <= sma && next_state.orbit.sma() > sma,
            EventKind::Ecc(ecc) => prev_state.orbit.ecc() <= ecc && next_state.orbit.ecc() > ecc,
            EventKind::Inc(inc) => prev_state.orbit.inc() <= inc && next_state.orbit.inc() > inc,
            EventKind::Raan(raan) => {
                prev_state.orbit.raan() <= raan && next_state.orbit.raan() > raan
            }
            EventKind::Aop(aop) => prev_state.orbit.aop() <= aop && next_state.orbit.aop() > aop,
            EventKind::TA(angle) => prev_state.orbit.ta() <= angle && next_state.orbit.ta() > angle,
            EventKind::Periapse => prev_state.orbit.ta() > next_state.orbit.ta(),
            EventKind::Apoapse => {
                between_pm_180(prev_state.orbit.ta()) > between_pm_180(next_state.orbit.ta())
            }
            EventKind::Fuel(mass) => prev_state.fuel_mass <= mass && next_state.fuel_mass > mass,
        }
    }
}
