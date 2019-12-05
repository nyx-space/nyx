use crate::celestia::{Cosm, Geoid, State};
use crate::dynamics::spacecraft::SpacecraftState;
use crate::utils::between_pm_180;
use std::fmt;

/// A general Event
pub trait Event: fmt::Debug {
    /// Defines the type which will be accepted by the condition
    type StateType: Copy;

    // Evaluation of the event, must return whether the condition happened between between both states.
    fn eval_crossing(&self, prev_state: &Self::StateType, next_state: &Self::StateType) -> bool;
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
}

impl<S: Copy> fmt::Display for EventTrackers<S> {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for event_no in 0..self.events.len() {
            if event_no > 0 {
                writeln!(f)?;
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
    pub tgt: Option<Geoid>,
    pub cosm: Option<&'a Cosm>,
}

impl<'a> OrbitalEvent<'a> {
    pub fn new(kind: EventKind) -> Box<Self> {
        Box::new(OrbitalEvent {
            kind,
            tgt: None,
            cosm: None,
        })
    }
    pub fn in_frame(kind: EventKind, tgt: Geoid, cosm: &'a Cosm) -> Box<Self> {
        Box::new(OrbitalEvent {
            kind,
            tgt: Some(tgt),
            cosm: Some(cosm),
        })
    }
}

impl<'a> Event for OrbitalEvent<'a> {
    type StateType = State<Geoid>;

    fn eval_crossing(&self, prev_state: &Self::StateType, next_state: &Self::StateType) -> bool {
        let (prev, next) = match self.tgt {
            Some(tgt) => (
                self.cosm.unwrap().frame_chg(prev_state, tgt),
                self.cosm.unwrap().frame_chg(next_state, tgt),
            ),
            None => (*prev_state, *next_state),
        };
        match self.kind {
            EventKind::Sma(sma) => prev.sma() <= sma && next.sma() > sma,
            EventKind::Ecc(ecc) => prev.ecc() <= ecc && next.ecc() > ecc,
            EventKind::Inc(inc) => prev.inc() <= inc && next.inc() > inc,
            EventKind::Raan(raan) => prev.raan() <= raan && next.raan() > raan,
            EventKind::Aop(aop) => prev.aop() <= aop && next.aop() > aop,
            EventKind::TA(angle) => prev.ta() <= angle && next.ta() > angle,
            EventKind::Periapse => prev.ta() > next.ta(),
            EventKind::Apoapse => between_pm_180(prev.ta()) > between_pm_180(next.ta()),
            _ => panic!("event {:?} not supported", self.kind),
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

/// A condition to stop a propagator
#[derive(Debug)]
pub struct StopCondition<S: Copy> {
    pub event: Box<dyn Event<StateType = S>>,
    pub trigger: usize,
    pub max_iter: usize,
    pub epsilon: f64,
}
