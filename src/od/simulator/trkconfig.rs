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

use super::scheduler::Scheduler;
use crate::cosmic::Cosm;
use crate::io::{duration_from_str, duration_to_str, epoch_from_str, epoch_to_str, ConfigError};
use crate::io::{ConfigRepr, Configurable};
use hifitime::TimeUnits;
use hifitime::{Duration, Epoch};
#[cfg(feature = "python")]
use pyo3::prelude::*;
use serde::Deserialize;
use serde_derive::Serialize;
use std::fmt::Debug;
use std::sync::Arc;
use typed_builder::TypedBuilder;

/// Stores a tracking configuration, there is one per tracking data simulator (e.g. one for ground station #1 and another for #2).
/// By default, the tracking configuration is continuous and the tracking arc is from the beginning of the simulation to the end.
/// In Python, any value that is set to None at initialization will use the default values.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, TypedBuilder)]
#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(module = "nyx_space.orbit_determination"))]
#[cfg_attr(feature = "python", pyo3(get_all, set_all))]
#[builder(doc)]
pub struct TrkConfig {
    /// Set to automatically build a tracking schedule based on some criteria
    #[serde(default)]
    #[builder(default, setter(strip_option))]
    pub scheduler: Option<Scheduler>,
    #[serde(
        serialize_with = "duration_to_str",
        deserialize_with = "duration_from_str"
    )]
    /// Sampling rate once tracking has started
    #[builder(default = 1.minutes())]
    pub sampling: Duration,
    /// List of tracking strands during which the given tracker will be tracking
    #[builder(default, setter(strip_option))]
    pub strands: Option<Vec<Strand>>,
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
    /// Initialize a default TrkConfig providing only the sample rate.
    /// Note: this will also set the sample alignment time to the provided duration.
    pub fn from_sample_rate(sampling: Duration) -> Self {
        Self {
            sampling,
            scheduler: Some(Scheduler::builder().sample_alignment(sampling).build()),
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
            if strands.is_empty() && self.scheduler.is_none() {
                return Err(ConfigError::InvalidConfig(
                    "Provided tracking strands is empty and no scheduler is defined".to_string(),
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
                "Neither tracking strands not a scheduler is provided".to_string(),
            ));
        }

        Ok(())
    }
}

impl Default for TrkConfig {
    /// The default configuration is to generate a measurement every minute (continuously) while the vehicle is visible
    fn default() -> Self {
        Self {
            // Allows calling the builder's defaults
            scheduler: Some(Scheduler::builder().build()),
            sampling: 1.minutes(),
            strands: None,
        }
    }
}

/// Stores a tracking strand with a start and end epoch
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(module = "nyx_space.orbit_determination"))]
#[cfg_attr(feature = "python", pyo3(get_all, set_all))]
pub struct Strand {
    #[serde(serialize_with = "epoch_to_str", deserialize_with = "epoch_from_str")]
    pub start: Epoch,
    #[serde(serialize_with = "epoch_to_str", deserialize_with = "epoch_from_str")]
    pub end: Epoch,
}

#[cfg_attr(feature = "python", pymethods)]
impl Strand {
    /// Returns whether the provided epoch is within the range
    pub fn contains(&self, epoch: Epoch) -> bool {
        (self.start..=self.end).contains(&epoch)
    }

    /// Returns the duration of this tracking strand
    pub fn duration(&self) -> Duration {
        self.end - self.start
    }

    /// Builds a new strand with the start and end epochs of this tracking strand.
    #[cfg(feature = "python")]
    #[new]
    fn py_new(start: Epoch, end: Epoch) -> Self {
        Self { start, end }
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
        cfg.strands = Some(vec![Strand { start, end }]);
        assert!(
            cfg.sanity_check().is_err(),
            "strand of too short of a duration should mark this insane"
        );

        let end = start + cfg.sampling;
        cfg.strands = Some(vec![Strand { start, end }]);
        assert!(
            cfg.sanity_check().is_ok(),
            "strand allowing for a single measurement should be OK"
        );

        // An anti-chronological strand should be invalid
        cfg.strands = Some(vec![Strand {
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
        assert_eq!(
            cfg.scheduler.unwrap(),
            Scheduler::builder().min_samples(10).build()
        );
        assert!(cfg.strands.is_none());

        // Specify an intermittent schedule and a specific start epoch.
        let cfg = TrkConfig {
            scheduler: Some(Scheduler {
                cadence: Cadence::Intermittent {
                    on: 23.1.hours(),
                    off: 0.9.hours(),
                },
                handoff: Handoff::Eager,
                min_samples: 10,
                ..Default::default()
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
        use std::collections::BTreeMap;
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

        let configs: BTreeMap<String, TrkConfig> = TrkConfig::load_named(trkconfg_yaml).unwrap();
        dbg!(configs);
    }

    #[test]
    fn api_trk_config() {
        use serde_yaml;

        let cfg = TrkConfig::builder()
            .sampling(15.seconds())
            .scheduler(Scheduler::builder().handoff(Handoff::Overlap).build())
            .build();

        let serialized = serde_yaml::to_string(&cfg).unwrap();
        println!("{serialized}");
        let deserd: TrkConfig = serde_yaml::from_str(&serialized).unwrap();
        assert_eq!(deserd, cfg);

        let cfg = TrkConfig::builder()
            .scheduler(Scheduler::builder().handoff(Handoff::Overlap).build())
            .build();

        assert_eq!(cfg.sampling, 60.seconds());
    }
}
