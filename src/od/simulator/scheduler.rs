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

pub use crate::dynamics::{Dynamics, NyxError};
use crate::io::{duration_from_str, duration_to_str};
pub use crate::{cosmic::Cosm, State, TimeTagged};
use hifitime::Duration;
use serde::Deserialize;
use serde_derive::Serialize;
use std::fmt::Debug;
use typed_builder::TypedBuilder;

/// Defines the handoff from a current ground station to the next one that is visible to prevent overlapping of measurements
#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
pub enum Handoff {
    /// If a new station is in visibility of the spacecraft, the "Eager" station will immediately stop tracking and switch over (default)
    Eager,
    /// If a new station is in visibility of the spacecraft, the "Greedy" station will continue to tracking until the vehicle is below its elevation mask
    Greedy,
    /// If a new station is in visibility of the spacecraft, the "Overlap" station will continue tracking, and so will the other one
    Overlap,
}

impl Default for Handoff {
    fn default() -> Self {
        Self::Eager
    }
}

/// A scheduler allows building a scheduling of spaceraft tracking for a set of ground stations.
#[derive(Copy, Clone, Debug, Default, Deserialize, PartialEq, Serialize, TypedBuilder)]
pub struct Scheduler {
    pub handoff: Handoff,
    pub cadence: Cadence,
}

/// Determines whether tracking is continuous or intermittent.
#[derive(Copy, Clone, Deserialize, PartialEq, Serialize)]
pub enum Cadence {
    Continuous,
    /// An intermittent schedule has On and Off durations.
    Intermittent {
        #[serde(
            serialize_with = "duration_to_str",
            deserialize_with = "duration_from_str"
        )]
        on: Duration,
        #[serde(
            serialize_with = "duration_to_str",
            deserialize_with = "duration_from_str"
        )]
        off: Duration,
    },
}

impl Default for Cadence {
    fn default() -> Self {
        Self::Continuous
    }
}

impl Debug for Cadence {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Continuous => write!(f, "Continuous"),
            Self::Intermittent { on, off } => f
                .debug_struct("Intermittent")
                .field("on", &format!("{on}"))
                .field("off", &format!("{off}"))
                .finish(),
        }
    }
}

#[cfg(test)]
mod scheduler_ut {
    use process::simulator::scheduler::Handoff;

    use crate::od::prelude::*;

    use super::Scheduler;

    #[test]
    fn serde_cadence() {
        use hifitime::TimeUnits;
        use serde_yaml;

        let cont: Cadence = serde_yaml::from_str("!Continuous").unwrap();
        assert_eq!(cont, Cadence::Continuous);

        let int: Cadence =
            serde_yaml::from_str("!Intermittent {on: 1 h 35 min, off: 15 h 02 min 3 s}").unwrap();
        assert_eq!(
            int,
            Cadence::Intermittent {
                on: 1.hours() + 35.0.minutes(),
                off: 15.hours() + 2.minutes() + 3.seconds()
            }
        );
        assert_eq!(
            format!("{int:?}"),
            r#"Intermittent { on: "1 h 35 min", off: "15 h 2 min 3 s" }"#
        );

        let serialized = serde_yaml::to_string(&int).unwrap();
        let deserd: Cadence = serde_yaml::from_str(&serialized).unwrap();
        assert_eq!(deserd, int);
    }

    #[test]
    fn api_and_serde_scheduler() {
        use hifitime::TimeUnits;
        use serde_yaml;

        let scheduler = Scheduler::default();
        let serialized = serde_yaml::to_string(&scheduler).unwrap();
        assert_eq!(serialized, "handoff: Eager\ncadence: Continuous\n");
        let deserd: Scheduler = serde_yaml::from_str(&serialized).unwrap();
        assert_eq!(deserd, scheduler);

        let scheduler = Scheduler::builder()
            .handoff(Handoff::Eager)
            .cadence(Cadence::Intermittent {
                on: 0.2.hours(),
                off: 17.hours() + 5.minutes(),
            })
            .build();

        let serialized = serde_yaml::to_string(&scheduler).unwrap();
        assert_eq!(
            serialized,
            "handoff: Eager\ncadence: !Intermittent\n  on: 12 min\n  off: 17 h 5 min\n"
        );
        let deserd: Scheduler = serde_yaml::from_str(&serialized).unwrap();
        assert_eq!(deserd, scheduler);
    }

    #[test]
    fn defaults() {
        let sched = Scheduler::default();

        assert_eq!(sched.cadence, Cadence::Continuous);

        assert_eq!(sched.handoff, Handoff::Eager);
    }
}
