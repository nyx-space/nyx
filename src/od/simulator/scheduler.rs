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

use crate::io::{
    duration_from_str, duration_to_str, maybe_duration_from_str, maybe_duration_to_str,
};
pub use crate::State;
use hifitime::{Duration, Unit};
use serde::Deserialize;
use serde_derive::Serialize;
use std::fmt::Debug;
use typed_builder::TypedBuilder;

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Defines the handoff from a current ground station to the next one that is visible to prevent overlapping of measurements
#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize)]
#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(module = "nyx_space.orbit_determination"))]
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
#[builder(doc)]
#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(module = "nyx_space.orbit_determination"))]
pub struct Scheduler {
    /// Handoff strategy if two trackers see the vehicle at the same time
    #[builder(default)]
    pub handoff: Handoff,
    /// On/off cadence of this scheduler
    #[builder(default)]
    pub cadence: Cadence,
    /// Minimum number of samples for a valid arc, i.e. if there are less than this many samples during a pass, the strand is discarded.
    #[builder(default = 10)]
    pub min_samples: u32,
    /// Round the time of the samples to the provided duration. For example, if the vehicle is above the horizon at 01:02:03.456 and the alignment
    /// is set to 01 seconds, then this will cause the tracking to start at 01:02:03 as it is rounded to the nearest second.
    #[builder(default = Some(Unit::Second * 1.0), setter(strip_option))]
    #[serde(
        serialize_with = "maybe_duration_to_str",
        deserialize_with = "maybe_duration_from_str"
    )]
    pub sample_alignment: Option<Duration>,
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
        assert_eq!(
            serialized,
            "handoff: Eager\ncadence: Continuous\nmin_samples: 0\nsample_alignment: null\n"
        );
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
            "handoff: Eager\ncadence: !Intermittent\n  on: 12 min\n  off: 17 h 5 min\nmin_samples: 10\nsample_alignment: 1 s\n"
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
