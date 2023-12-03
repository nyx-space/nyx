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
use crate::io::{duration_from_str, duration_to_str, epoch_from_str, epoch_to_str, ConfigError};
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

use super::scheduler::Scheduler;

/// Stores a tracking configuration, there is one per tracking data simulator (e.g. one for ground station #1 and another for #2).
/// By default, the tracking configuration is continuous and the tracking arc is from the beginning of the simulation to the end.
/// In Python, any value that is set to None at initialization will use the default values.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(module = "nyx_space.orbit_determination"))]
pub struct TrkConfig {
    /// Set to automatically build a tracking schedule based on some criteria
    #[serde(default)]
    pub scheduler: Option<Scheduler>,
    #[serde(
        serialize_with = "duration_to_str",
        deserialize_with = "duration_from_str"
    )]
    /// Sampling rate once tracking has started
    pub sampling: Duration,
    /// List of tracking strands during which the given tracker will be tracking
    pub strands: Option<Vec<EpochRanges>>,
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
            scheduler: Some(Scheduler::default()),
            ..Default::default()
        }
    }

    /// Check that the configuration is valid: a successful call means that either we have a set of tracking strands or we have a valid scheduler
    pub(crate) fn sanity_check(&self) -> Result<(), ConfigError> {
        if self.strands.is_some() && self.scheduler.is_some() {
            return Err(ConfigError::InvalidConfig(
                "Both tracking strands and a scheduler are configured, must be one or the other"
                    .to_string(),
            ));
        } else if let Some(strands) = &self.strands {
            if strands.is_empty() {
                return Err(ConfigError::InvalidConfig(
                    "Provided tracking strands is empty (set to None to use scheduler)".to_string(),
                ));
            }
            for (ii, strand) in strands.iter().enumerate() {
                if strand.duration() < self.sampling {
                    return Err(ConfigError::InvalidConfig(format!(
                        "Strand #{ii} is shorter than sampling time"
                    )));
                }
                if strand.duration().is_negative() {
                    return Err(ConfigError::InvalidConfig(format!(
                        "Strand #{ii} is anti-chronological"
                    )));
                }
            }
        } else if self.strands.is_none() && self.scheduler.is_none() {
            return Err(ConfigError::InvalidConfig(
                "Provided tracking strands is empty (set to None to use scheduler)".to_string(),
            ));
        }

        Ok(())
    }
}

impl Default for TrkConfig {
    /// The default configuration is to generate a measurement every minute (continuously) while the vehicle is visible
    fn default() -> Self {
        Self {
            scheduler: Some(Scheduler::default()),
            sampling: 1.minutes(),
            strands: None,
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

    pub fn duration(&self) -> Duration {
        self.end - self.start
    }
}

#[cfg(test)]
mod trkconfig_ut {
    use crate::io::ConfigRepr;
    use crate::od::prelude::*;

    #[test]
    fn sanity_checks() {
        let mut cfg = TrkConfig::default();
        assert!(cfg.sanity_check().is_ok(), "default config should be sane");

        cfg.scheduler = None;
        assert!(
            cfg.sanity_check().is_err(),
            "no scheduler should mark this insane"
        );

        cfg.strands = Some(Vec::new());
        assert!(
            cfg.sanity_check().is_err(),
            "no scheduler and empty strands should mark this insane"
        );

        let start = Epoch::now().unwrap();
        let end = start + 10.seconds();
        cfg.strands = Some(vec![EpochRanges { start, end }]);
        assert!(
            cfg.sanity_check().is_err(),
            "strand of too short of a duration should mark this insane"
        );

        let end = start + cfg.sampling;
        cfg.strands = Some(vec![EpochRanges { start, end }]);
        assert!(
            cfg.sanity_check().is_ok(),
            "strand allowing for a single measurement should be OK"
        );

        // An anti-chronological strand should be invalid
        cfg.strands = Some(vec![EpochRanges {
            start: end,
            end: start,
        }]);
        assert!(
            cfg.sanity_check().is_err(),
            "anti chronological strand should be insane"
        );
    }

    #[test]
    fn serde_trkconfig() {
        use serde_yaml;

        // Test the default config
        let cfg = TrkConfig::default();
        let serialized = serde_yaml::to_string(&cfg).unwrap();
        println!("{serialized}");
        let deserd: TrkConfig = serde_yaml::from_str(&serialized).unwrap();
        assert_eq!(deserd, cfg);
        assert_eq!(cfg.scheduler.unwrap(), Scheduler::default());
        assert!(cfg.strands.is_none());

        // Specify an intermittent schedule and a specific start epoch.
        let cfg = TrkConfig {
            scheduler: Some(Scheduler {
                cadence: Cadence::Intermittent {
                    on: 23.1.hours(),
                    off: 0.9.hours(),
                },
                handoff: Handoff::Eager,
            }),
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

        let configs: HashMap<String, TrkConfig> = TrkConfig::load_named(trkconfg_yaml).unwrap();
        dbg!(configs);
    }
}
