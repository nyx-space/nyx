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

/// Determines whether a tracking schedule is continuous or intermittent.
#[derive(Copy, Clone, Serialize, Deserialize, PartialEq)]
pub enum Schedule {
    Continuous,
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

impl Default for Schedule {
    fn default() -> Self {
        Self::Continuous
    }
}

impl Debug for Schedule {
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

#[test]
fn serde_schedule() {
    use hifitime::TimeUnits;
    use serde_yaml;

    let cont: Schedule = serde_yaml::from_str("!Continuous").unwrap();
    assert_eq!(cont, Schedule::Continuous);

    let int: Schedule =
        serde_yaml::from_str("!Intermittent {on: 1 h 35 min, off: 15 h 02 min 3 s}").unwrap();
    assert_eq!(
        int,
        Schedule::Intermittent {
            on: 1.hours() + 35.0.minutes(),
            off: 15.hours() + 2.minutes() + 3.seconds()
        }
    );
    assert_eq!(
        format!("{int:?}"),
        r#"Intermittent { on: "1 h 35 min", off: "15 h 2 min 3 s" }"#
    );

    let serialized = serde_yaml::to_string(&int).unwrap();
    let deserd: Schedule = serde_yaml::from_str(&serialized).unwrap();
    assert_eq!(deserd, int);
}
