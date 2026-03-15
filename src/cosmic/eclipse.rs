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
use anise::analysis::prelude::Event;
use anise::astro::Occultation;
use anise::constants::frames::{EARTH_J2000, MOON_J2000, SUN_J2000};
use anise::errors::AlmanacResult;

pub use super::{Frame, Orbit, Spacecraft};
use std::fmt;
use std::sync::Arc;

#[derive(Clone, Debug)]
pub struct ShadowModel {
    pub light_source: Frame,
    pub shadow_bodies: Vec<Frame>,
}

impl fmt::Display for ShadowModel {
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

impl ShadowModel {
    /// Creates a new typical eclipse locator.
    /// The light source is the Sun, and the shadow bodies are the Earth and the Moon.
    pub fn cislunar(almanac: Arc<Almanac>) -> Self {
        let eme2k = almanac.frame_info(EARTH_J2000).unwrap();
        let moon_j2k = almanac.frame_info(MOON_J2000).unwrap();
        Self {
            light_source: almanac.frame_info(SUN_J2000).unwrap(),
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
    pub fn to_umbra_events(&self) -> Vec<Event> {
        self.shadow_bodies
            .iter()
            .copied()
            .map(Event::total_eclipse)
            .collect()
    }

    /// Creates a penumbra event from this eclipse locator
    // Evaluation of the event, returns 0.0 for umbra, 1.0 for visibility (no shadow) and some value in between for penumbra
    pub fn to_penumbra_events(&self) -> Vec<Event> {
        self.shadow_bodies
            .iter()
            .copied()
            .map(Event::eclipse)
            .collect()
    }
}
