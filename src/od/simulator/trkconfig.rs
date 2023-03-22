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
use crate::io::{duration_from_str, duration_to_str, epoch_from_str, epoch_to_str};
use crate::io::{ConfigRepr, Configurable};
pub use crate::{cosmic::Cosm, State, TimeTagged};
use hifitime::TimeUnits;
use hifitime::{Duration, Epoch};
#[cfg(feature = "python")]
use pyo3::prelude::*;
use serde::Deserialize;
use serde_derive::Serialize;
use std::fmt::Debug;
use std::sync::Arc;

use super::schedule::Schedule;
use super::Availability;

/// Stores a tracking configuration, there is one per tracking data simulator (e.g. one for ground station #1 and another for #2)
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "python", pyclass)]
pub struct TrkConfig {
    /// Availability configuration to start the tracking arc
    #[serde(default)]
    pub start: Availability,
    /// Availability configuration to end the tracking arc
    #[serde(default)]
    pub end: Availability,
    #[serde(default)]
    pub schedule: Schedule,
    #[serde(
        serialize_with = "duration_to_str",
        deserialize_with = "duration_from_str"
    )]
    pub sampling: Duration,
    /// List of epoch ranges to include
    #[serde(rename = "inclusion epochs")]
    pub inclusion_epochs: Option<Vec<EpochRanges>>,
    /// List of epoch ranges to exclude
    #[serde(rename = "exclusion epochs")]
    pub exclusion_epochs: Option<Vec<EpochRanges>>,
}

impl ConfigRepr for TrkConfig {}

impl Configurable for TrkConfig {
    type IntermediateRepr = Self;

    fn from_config(
        cfg: Self::IntermediateRepr,
        _cosm: Arc<Cosm>,
    ) -> Result<Self, crate::io::ConfigError>
    where
        Self: Sized,
    {
        Ok(cfg)
    }

    fn to_config(&self) -> Result<Self::IntermediateRepr, crate::io::ConfigError> {
        Ok(self.clone())
    }
}

impl TrkConfig {
    /// Initialize a default TrkConfig providing only the sample rate
    pub fn from_sample_rate(sampling: Duration) -> Self {
        Self {
            sampling,
            ..Default::default()
        }
    }
}

impl Default for TrkConfig {
    /// The default configuration is to generate a measurement every minute (continuously) while the vehicle is visible
    fn default() -> Self {
        Self {
            start: Availability::Visible,
            end: Availability::Visible,
            schedule: Schedule::Continuous,
            sampling: 1.minutes(),
            inclusion_epochs: None,
            exclusion_epochs: None,
        }
    }
}

/// Stores an epoch range for tracking.
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "python", pyclass)]
pub struct EpochRanges {
    #[serde(serialize_with = "epoch_to_str", deserialize_with = "epoch_from_str")]
    pub start: Epoch,
    #[serde(serialize_with = "epoch_to_str", deserialize_with = "epoch_from_str")]
    pub end: Epoch,
}

impl EpochRanges {
    /// Returns whether the provided epoch is within the range
    pub fn contains(&self, epoch: Epoch) -> bool {
        (self.start..=self.end).contains(&epoch)
    }
}

#[test]
fn serde_trkconfig() {
    use hifitime::{Epoch, TimeScale};
    use serde_yaml;

    // Test the default config
    let cfg = TrkConfig::default();
    let serialized = serde_yaml::to_string(&cfg).unwrap();
    println!("{serialized}");
    let deserd: TrkConfig = serde_yaml::from_str(&serialized).unwrap();
    assert_eq!(deserd, cfg);

    // Specify an intermittent schedule and a specific start epoch.
    let cfg = TrkConfig {
        start: Availability::Epoch(Epoch::from_gregorian_at_midnight(
            2023,
            2,
            22,
            TimeScale::TAI,
        )),
        end: Availability::Visible,
        schedule: Schedule::Intermittent {
            on: 23.1.hours(),
            off: 0.9.hours(),
        },
        sampling: 45.2.seconds(),
        ..Default::default()
    };
    let serialized = serde_yaml::to_string(&cfg).unwrap();
    println!("{serialized}");
    let deserd: TrkConfig = serde_yaml::from_str(&serialized).unwrap();
    assert_eq!(deserd, cfg);
}

#[test]
fn deserialize_from_file() {
    use std::collections::HashMap;
    use std::env;
    use std::path::PathBuf;

    // Load the tracking configuration from the test data.
    let trkconfg_yaml: PathBuf = [
        &env::var("CARGO_MANIFEST_DIR").unwrap(),
        "data",
        "tests",
        "config",
        "tracking_cfg.yaml",
    ]
    .iter()
    .collect();

    let configs: HashMap<String, TrkConfig> = TrkConfig::load_named_yaml(trkconfg_yaml).unwrap();
    dbg!(configs);
}
