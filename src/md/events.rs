/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use crate::celestia::{Cosm, Frame, Orbit};
use crate::utils::between_pm_180;
use crate::SpacecraftState;
use std::fmt;

/// A general Event
pub trait Event: Send + Sync + fmt::Debug {
    /// Defines the type which will be accepted by the condition
    type StateType: Copy;

    // Evaluation of event crossing, must return whether the condition happened between between both states.
    fn eval_crossing(&self, prev_state: &Self::StateType, next_state: &Self::StateType) -> bool;

    // Evaluation of the event, must return a value corresponding to whether the state is before or after the event
    fn eval(&self, state: &Self::StateType) -> f64;
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
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

/// An orbital event, in the same frame or in another frame.
#[derive(Debug)]
pub struct OrbitalEvent<'a> {
    pub kind: EventKind,
    pub tgt: Option<Frame>,
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
    pub fn in_frame(kind: EventKind, tgt: Frame, cosm: &'a Cosm) -> Box<Self> {
        Box::new(OrbitalEvent {
            kind,
            tgt: Some(tgt),
            cosm: Some(cosm),
        })
    }
}

impl<'a> Event for OrbitalEvent<'a> {
    type StateType = Orbit;

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
            EventKind::Periapse => between_pm_180(state.ta()),
            EventKind::Apoapse => {
                // We use the sign change in flight path angle to determine that we have crossed the apoapse
                between_pm_180(state.ta()) - 180.0
            }
            _ => panic!("event {:?} not supported", self.kind),
        }
    }

    fn eval_crossing(&self, prev_state: &Self::StateType, next_state: &Self::StateType) -> bool {
        let prev_val = self.eval(prev_state);
        let next_val = self.eval(next_state);
        match self.kind {
            // XXX: Should this condition be applied to all angles?
            EventKind::Periapse => prev_val < 0.0 && next_val >= 0.0,
            EventKind::Apoapse => prev_val > 0.0 && next_val <= 0.0,
            _ => prev_val * next_val <= 0.0,
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
            EventKind::Fuel(mass) => state.fuel_mass_kg - mass,
            _ => self.orbital.as_ref().unwrap().eval(&state.orbit),
        }
    }

    fn eval_crossing(&self, prev_state: &Self::StateType, next_state: &Self::StateType) -> bool {
        match self.kind {
            EventKind::Fuel(mass) => {
                prev_state.fuel_mass_kg <= mass && next_state.fuel_mass_kg > mass
            }
            _ => self
                .orbital
                .as_ref()
                .unwrap()
                .eval_crossing(&prev_state.orbit, &next_state.orbit),
        }
    }
}
