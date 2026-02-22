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

use super::scheduler::Scheduler;
use crate::io::ConfigRepr;
use crate::io::{duration_from_str, duration_to_str, epoch_from_str, epoch_to_str, ConfigError};
use der::{Decode, Encode, Reader};
use hifitime::TimeUnits;
use hifitime::{Duration, Epoch, TimeScale};

#[cfg(feature = "python")]
use pyo3::{exceptions::PyValueError, prelude::*, types::PyBytes, types::PyType};

use serde::Deserialize;
use serde_derive::Serialize;
use std::fmt::Debug;
use std::str::FromStr;
use typed_builder::TypedBuilder;

/// Stores a tracking configuration, there is one per tracking data simulator (e.g. one for ground station #1 and another for #2).
/// By default, the tracking configuration is continuous and the tracking arc is from the beginning of the simulation to the end.
/// In Python, any value that is set to None at initialization will use the default values.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, TypedBuilder)]
#[cfg_attr(feature = "python", pyclass)]
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

impl<'a> Decode<'a> for TrkConfig {
    fn decode<R: Reader<'a>>(decoder: &mut R) -> der::Result<Self> {
        let scheduler = if decoder.decode::<bool>()? {
            Some(decoder.decode()?)
        } else {
            None
        };
        let sampling_ns = decoder.decode::<i128>()?;
        let strands = if decoder.decode::<bool>()? {
            Some(decoder.decode()?)
        } else {
            None
        };

        Ok(Self {
            scheduler,
            sampling: Duration::from_total_nanoseconds(sampling_ns),
            strands,
        })
    }
}

impl Encode for TrkConfig {
    fn encoded_len(&self) -> der::Result<der::Length> {
        let mut len = self.scheduler.is_some().encoded_len()?;
        if let Some(sched) = &self.scheduler {
            len = (len + sched.encoded_len()?)?;
        }
        len = (len + self.sampling.total_nanoseconds().encoded_len()?)?;
        len = (len + self.strands.is_some().encoded_len()?)?;
        if let Some(strands) = &self.strands {
            len = (len + strands.encoded_len()?)?;
        }
        Ok(len)
    }

    fn encode(&self, encoder: &mut impl der::Writer) -> der::Result<()> {
        if let Some(sched) = &self.scheduler {
            true.encode(encoder)?;
            sched.encode(encoder)?;
        } else {
            false.encode(encoder)?;
        }
        self.sampling.total_nanoseconds().encode(encoder)?;
        if let Some(strands) = &self.strands {
            true.encode(encoder)?;
            strands.encode(encoder)?;
        } else {
            false.encode(encoder)?;
        }
        Ok(())
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl TrkConfig {
    #[new]
    #[pyo3(signature = (scheduler=None, sampling=1.minutes(), strands=None))]
    fn py_new(
        scheduler: Option<Scheduler>,
        sampling: Duration,
        strands: Option<Vec<Strand>>,
    ) -> Self {
        Self {
            scheduler,
            sampling,
            strands,
        }
    }

    #[getter]
    fn get_scheduler(&self) -> Option<Scheduler> {
        self.scheduler
    }

    #[setter]
    fn set_scheduler(&mut self, scheduler: Option<Scheduler>) {
        self.scheduler = scheduler;
    }

    #[getter]
    fn get_sampling(&self) -> Duration {
        self.sampling
    }

    #[setter]
    fn set_sampling(&mut self, sampling: Duration) {
        self.sampling = sampling;
    }

    #[getter]
    fn get_strands(&self) -> Option<Vec<Strand>> {
        self.strands.clone()
    }

    #[setter]
    fn set_strands(&mut self, strands: Option<Vec<Strand>>) {
        self.strands = strands;
    }

    fn __repr__(&self) -> String {
        format!("{self:?}")
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    /// Decodes an ASN.1 DER encoded byte array into a TrkConfig object.
    ///
    /// :type data: bytes
    /// :rtype: TrkConfig
    #[classmethod]
    pub fn from_asn1(_cls: &Bound<'_, PyType>, data: &[u8]) -> PyResult<Self> {
        match Self::from_der(data) {
            Ok(obj) => Ok(obj),
            Err(e) => Err(PyValueError::new_err(format!("ASN.1 decoding error: {e}"))),
        }
    }

    /// Encodes this TrkConfig object into an ASN.1 DER encoded byte array.
    ///
    /// :rtype: bytes
    pub fn to_asn1<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let mut buf = Vec::new();
        match self.encode_to_vec(&mut buf) {
            Ok(_) => Ok(PyBytes::new(py, &buf)),
            Err(e) => Err(PyValueError::new_err(format!("ASN.1 encoding error: {e}"))),
        }
    }
}

impl ConfigRepr for TrkConfig {}

impl FromStr for TrkConfig {
    type Err = ConfigError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        serde_yml::from_str(s).map_err(|source| ConfigError::ParseError { source })
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
            return Err(ConfigError::InvalidConfig {
                msg:
                    "Both tracking strands and a scheduler are configured, must be one or the other"
                        .to_string(),
            });
        } else if let Some(strands) = &self.strands {
            if strands.is_empty() && self.scheduler.is_none() {
                return Err(ConfigError::InvalidConfig {
                    msg: "Provided tracking strands is empty and no scheduler is defined"
                        .to_string(),
                });
            }
            for (ii, strand) in strands.iter().enumerate() {
                if strand.duration() < self.sampling {
                    return Err(ConfigError::InvalidConfig {
                        msg: format!(
                            "Strand #{ii} lasts {} which is shorter than sampling time of {}",
                            strand.duration(),
                            self.sampling
                        ),
                    });
                }
                if strand.duration().is_negative() {
                    return Err(ConfigError::InvalidConfig {
                        msg: format!("Strand #{ii} is anti-chronological"),
                    });
                }
            }
        } else if self.strands.is_none() && self.scheduler.is_none() {
            return Err(ConfigError::InvalidConfig {
                msg: "Neither tracking strands not a scheduler is provided".to_string(),
            });
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
pub struct Strand {
    #[serde(serialize_with = "epoch_to_str", deserialize_with = "epoch_from_str")]
    pub start: Epoch,
    #[serde(serialize_with = "epoch_to_str", deserialize_with = "epoch_from_str")]
    pub end: Epoch,
}

impl<'a> Decode<'a> for Strand {
    fn decode<R: Reader<'a>>(decoder: &mut R) -> der::Result<Self> {
        let start_ns = decoder.decode::<i128>()?;
        let start_ts_u8 = decoder.decode::<u8>()?;
        let start_ts = TimeScale::from(start_ts_u8);

        let end_ns = decoder.decode::<i128>()?;
        let end_ts_u8 = decoder.decode::<u8>()?;
        let end_ts = TimeScale::from(end_ts_u8);

        Ok(Self {
            start: Epoch::from_duration(Duration::from_total_nanoseconds(start_ns), start_ts),
            end: Epoch::from_duration(Duration::from_total_nanoseconds(end_ns), end_ts),
        })
    }
}

impl Encode for Strand {
    fn encoded_len(&self) -> der::Result<der::Length> {
        let ts_len = 1u8.encoded_len()?;
        let start_len = (self.start.duration.total_nanoseconds().encoded_len()? + ts_len)?;
        let end_len = (self.end.duration.total_nanoseconds().encoded_len()? + ts_len)?;
        start_len + end_len
    }

    fn encode(&self, encoder: &mut impl der::Writer) -> der::Result<()> {
        self.start.duration.total_nanoseconds().encode(encoder)?;
        (self.start.time_scale as u8).encode(encoder)?;

        self.end.duration.total_nanoseconds().encode(encoder)?;
        (self.end.time_scale as u8).encode(encoder)
    }
}

impl Strand {
    pub fn new(start: Epoch, end: Epoch) -> Self {
        Self { start, end }
    }

    /// Returns whether the provided epoch is within the range
    pub fn contains(&self, epoch: Epoch) -> bool {
        (self.start..=self.end).contains(&epoch)
    }

    /// Returns the duration of this tracking strand
    pub fn duration(&self) -> Duration {
        self.end - self.start
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl Strand {
    #[new]
    fn py_new(start: Epoch, end: Epoch) -> Self {
        Self::new(start, end)
    }

    #[getter]
    fn get_start(&self) -> Epoch {
        self.start
    }

    #[setter]
    fn set_start(&mut self, start: Epoch) {
        self.start = start;
    }

    #[getter]
    fn get_end(&self) -> Epoch {
        self.end
    }

    #[setter]
    fn set_end(&mut self, end: Epoch) {
        self.end = end;
    }

    fn __repr__(&self) -> String {
        format!("{self:?}")
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    /// Decodes an ASN.1 DER encoded byte array into a Strand object.
    ///
    /// :type data: bytes
    /// :rtype: Strand
    #[classmethod]
    pub fn from_asn1(_cls: &Bound<'_, PyType>, data: &[u8]) -> PyResult<Self> {
        match Self::from_der(data) {
            Ok(obj) => Ok(obj),
            Err(e) => Err(PyValueError::new_err(format!("ASN.1 decoding error: {e}"))),
        }
    }

    /// Encodes this Strand object into an ASN.1 DER encoded byte array.
    ///
    /// :rtype: bytes
    pub fn to_asn1<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let mut buf = Vec::new();
        match self.encode_to_vec(&mut buf) {
            Ok(_) => Ok(PyBytes::new(py, &buf)),
            Err(e) => Err(PyValueError::new_err(format!("ASN.1 encoding error: {e}"))),
        }
    }
}

#[cfg(test)]
mod trkconfig_ut {
    use crate::io::ConfigRepr;
    use crate::od::simulator::{Cadence, Handoff, Scheduler, Strand, TrkConfig};
    use der::{Decode, Encode};
    use hifitime::{Epoch, TimeUnits};

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
        use serde_yml;

        // Test the default config
        let cfg = TrkConfig::default();
        let serialized = serde_yml::to_string(&cfg).unwrap();
        println!("{serialized}");
        let deserd: TrkConfig = serde_yml::from_str(&serialized).unwrap();
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
        let serialized = serde_yml::to_string(&cfg).unwrap();
        println!("{serialized}");
        let deserd: TrkConfig = serde_yml::from_str(&serialized).unwrap();
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
            "03_tests",
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
        use serde_yml;

        let cfg = TrkConfig::builder()
            .sampling(15.seconds())
            .scheduler(Scheduler::builder().handoff(Handoff::Overlap).build())
            .build();

        let serialized = serde_yml::to_string(&cfg).unwrap();
        println!("{serialized}");
        let deserd: TrkConfig = serde_yml::from_str(&serialized).unwrap();
        assert_eq!(deserd, cfg);

        let cfg = TrkConfig::builder()
            .scheduler(Scheduler::builder().handoff(Handoff::Overlap).build())
            .build();

        assert_eq!(cfg.sampling, 60.seconds());
    }

    #[test]
    fn test_handoff_asn1() {
        let h = Handoff::Greedy;
        let mut buf = Vec::new();
        h.encode_to_vec(&mut buf).unwrap();
        let h2 = Handoff::from_der(&buf).unwrap();
        assert_eq!(h, h2);
    }

    #[test]
    fn test_cadence_asn1() {
        let c = Cadence::Intermittent {
            on: 1.0.hours(),
            off: 0.5.hours(),
        };
        let mut buf = Vec::new();
        c.encode_to_vec(&mut buf).unwrap();
        let c2 = Cadence::from_der(&buf).unwrap();
        assert_eq!(c, c2);

        let c = Cadence::Continuous;
        let mut buf = Vec::new();
        c.encode_to_vec(&mut buf).unwrap();
        let c2 = Cadence::from_der(&buf).unwrap();
        assert_eq!(c, c2);
    }

    #[test]
    fn test_scheduler_asn1() {
        let s = Scheduler::builder()
            .handoff(Handoff::Overlap)
            .cadence(Cadence::Intermittent {
                on: 10.0.minutes(),
                off: 5.0.minutes(),
            })
            .min_samples(5)
            .sample_alignment(1.0.seconds())
            .build();

        let mut buf = Vec::new();
        s.encode_to_vec(&mut buf).unwrap();
        let s2 = Scheduler::from_der(&buf).unwrap();
        assert_eq!(s, s2);
    }

    #[test]
    fn test_strand_asn1() {
        let epoch = Epoch::from_gregorian_utc_at_midnight(2023, 1, 1);
        let s = Strand {
            start: epoch,
            end: epoch + 1.0.hours(),
        };

        let mut buf = Vec::new();
        s.encode_to_vec(&mut buf).unwrap();
        let s2 = Strand::from_der(&buf).unwrap();

        assert_eq!(s, s2);

        // Test TAI explicitly
        let epoch_tai = Epoch::from_gregorian_utc_at_midnight(2023, 1, 1);
        let s = Strand {
            start: epoch_tai,
            end: epoch_tai + 1.0.hours(),
        };

        let mut buf = Vec::new();
        s.encode_to_vec(&mut buf).unwrap();
        let s2 = Strand::from_der(&buf).unwrap();

        assert_eq!(s, s2);
    }

    #[test]
    fn test_trkconfig_asn1() {
        // Encode one in UTC and the other in TAI
        let epoch = Epoch::from_gregorian_utc_at_midnight(2023, 1, 1);
        let strand = Strand {
            start: epoch,
            end: (epoch + 1.0.hours()).to_time_scale(hifitime::TimeScale::TAI),
        };

        let cfg = TrkConfig::builder()
            .sampling(10.0.seconds())
            .strands(vec![strand])
            .build();

        let mut buf = Vec::new();
        cfg.encode_to_vec(&mut buf).unwrap();
        let cfg2 = TrkConfig::from_der(&buf).unwrap();

        assert_eq!(cfg, cfg2);
    }
}
