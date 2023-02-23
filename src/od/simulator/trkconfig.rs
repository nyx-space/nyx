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
use hifitime::TimeUnits;
use serde::Deserialize;
use serde_derive::Serialize;
use std::fmt::Debug;

use super::schedule::Schedule;
use super::StartMode;

/// Stores a tracking configuration, there is one per tracking data simulator (e.g. one for ground station #1 and another for #2)
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct TrkConfig {
    #[serde(default)]
    pub start_mode: StartMode,
    #[serde(default)]
    pub schedule: Schedule,
    #[serde(
        serialize_with = "duration_to_str",
        deserialize_with = "duration_from_str"
    )]
    pub sampling: Duration,
}

impl Default for TrkConfig {
    /// The default configuration is to generate a measurement every minute (continuously) while the vehicle is visible
    fn default() -> Self {
        Self {
            start_mode: StartMode::Visible,
            schedule: Schedule::Continuous,
            sampling: 1.minutes(),
        }
    }
}

#[test]
fn serde_trkconfig() {
    use hifitime::{Epoch, TimeScale};
    use serde_yaml;

    // Test the default config
    let cfg = TrkConfig::default();
    let serialized = serde_yaml::to_string(&cfg).unwrap();
    let deserd: TrkConfig = serde_yaml::from_str(&serialized).unwrap();
    assert_eq!(deserd, cfg);

    // Specify an intermittent schedule and a specific start epoch.
    let cfg = TrkConfig {
        start_mode: StartMode::Epoch(Epoch::from_gregorian_at_midnight(
            2023,
            2,
            22,
            TimeScale::TAI,
        )),
        schedule: Schedule::Intermittent {
            on: 23.1.hours(),
            off: 0.9.hours(),
        },
        sampling: 45.2.seconds(),
    };
    let serialized = serde_yaml::to_string(&cfg).unwrap();
    println!("{serialized}");
    let deserd: TrkConfig = serde_yaml::from_str(&serialized).unwrap();
    assert_eq!(deserd, cfg);
}
