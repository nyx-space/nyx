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

use anise::almanac::Almanac;
use anise::astro::Occultation;
use anise::constants::frames::{EARTH_J2000, MOON_J2000, SUN_J2000};
use anise::errors::AlmanacResult;
use snafu::ResultExt;

pub use super::{Frame, Orbit, Spacecraft};
use crate::errors::{EventAlmanacSnafu, EventError};
use crate::md::EventEvaluator;
use crate::time::{Duration, Unit};
use std::fmt;
use std::sync::Arc;

#[derive(Clone)]
pub struct EclipseLocator {
    pub light_source: Frame,
    pub shadow_bodies: Vec<Frame>,
}

impl fmt::Display for EclipseLocator {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let shadow_bodies: Vec<String> = self
            .shadow_bodies
            .iter()
            .map(|b| format!("{b:x}"))
            .collect();
        write!(
            f,
            "light-source: {:x}, shadows casted by: {}",
            self.light_source,
            shadow_bodies.join(", ")
        )
    }
}

impl EclipseLocator {
    /// Creates a new typical eclipse locator.
    /// The light source is the Sun, and the shadow bodies are the Earth and the Moon.
    pub fn cislunar(almanac: Arc<Almanac>) -> Self {
        let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
        let moon_j2k = almanac.frame_from_uid(MOON_J2000).unwrap();
        Self {
            light_source: almanac.frame_from_uid(SUN_J2000).unwrap(),
            shadow_bodies: vec![eme2k, moon_j2k],
        }
    }

    /// Compute the visibility/eclipse between an observer and an observed state
    pub fn compute(&self, observer: Orbit, almanac: Arc<Almanac>) -> AlmanacResult<Occultation> {
        let mut state = Occultation {
            epoch: observer.epoch,
            back_frame: SUN_J2000,
            front_frame: observer.frame,
            percentage: 0.0,
        };
        for eclipsing_body in &self.shadow_bodies {
            let this_state = almanac.solar_eclipsing(*eclipsing_body, observer, None)?;
            if this_state.percentage > state.percentage {
                state = this_state;
            }
        }
        Ok(state)
    }

    /// Creates an umbra event from this eclipse locator.
    /// Evaluation of the event, returns 0.0 for umbra, 1.0 for visibility (no shadow) and some value in between for penumbra
    pub fn to_umbra_event(&self) -> UmbraEvent {
        UmbraEvent {
            e_loc: self.clone(),
        }
    }

    /// Creates a penumbra event from this eclipse locator
    // Evaluation of the event, returns 0.0 for umbra, 1.0 for visibility (no shadow) and some value in between for penumbra
    pub fn to_penumbra_event(&self) -> PenumbraEvent {
        PenumbraEvent {
            e_loc: self.clone(),
        }
    }
}

/// An event to find the darkest eclipse state (more than 98% shadow)
pub struct UmbraEvent {
    e_loc: EclipseLocator,
}

impl fmt::Display for UmbraEvent {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "umbra event {}", self.e_loc)
    }
}

impl EventEvaluator<Spacecraft> for UmbraEvent {
    // Evaluation of the event
    fn eval(&self, sc: &Spacecraft, almanac: Arc<Almanac>) -> Result<f64, EventError> {
        Ok(self
            .e_loc
            .compute(sc.orbit, almanac)
            .context(EventAlmanacSnafu)?
            .factor())
        // match self
        //     .e_loc
        //     .compute(sc.orbit, almanac)
        //     .context(EventAlmanacSnafu)?
        // {
        //     EclipseState::Umbra => Ok(0.0),
        //     EclipseState::Visibilis => Ok(1.0),
        //     EclipseState::Penumbra(val) => Ok(val),
        // }
    }

    /// Stop searching when the time has converged to less than 0.1 seconds
    fn epoch_precision(&self) -> Duration {
        0.1 * Unit::Second
    }
    /// Finds the darkest part of an eclipse within 2% of penumbra (i.e. 98% in shadow)
    fn value_precision(&self) -> f64 {
        0.02
    }
    fn eval_string(&self, state: &Spacecraft, almanac: Arc<Almanac>) -> Result<String, EventError> {
        Ok(format!(
            "{}",
            self.e_loc
                .compute(state.orbit, almanac)
                .context(EventAlmanacSnafu)?
        ))
    }
}

/// An event to find the start of a penumbra
pub struct PenumbraEvent {
    e_loc: EclipseLocator,
}

impl fmt::Display for PenumbraEvent {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "penumbra event {}", self.e_loc)
    }
}

impl EventEvaluator<Spacecraft> for PenumbraEvent {
    fn eval(&self, sc: &Spacecraft, almanac: Arc<Almanac>) -> Result<f64, EventError> {
        Ok(self
            .e_loc
            .compute(sc.orbit, almanac)
            .context(EventAlmanacSnafu)?
            .factor())
    }

    /// Stop searching when the time has converged to less than 0.1 seconds
    fn epoch_precision(&self) -> Duration {
        0.1 * Unit::Second
    }
    /// Finds the slightest penumbra within 2% (i.e. 98% in visibility)
    fn value_precision(&self) -> f64 {
        0.02
    }

    fn eval_string(&self, state: &Spacecraft, almanac: Arc<Almanac>) -> Result<String, EventError> {
        Ok(format!(
            "{}",
            self.e_loc
                .compute(state.orbit, almanac)
                .context(EventAlmanacSnafu)?
        ))
    }
}
