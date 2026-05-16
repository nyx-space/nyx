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

use crate::io::{
    duration_from_str, duration_to_str, maybe_duration_from_str, maybe_duration_to_str,
};
pub use crate::State;
use der::{Decode, Encode, Enumerated, Reader};
use hifitime::{Duration, Unit};
use serde::Deserialize;
use serde::Serialize;
use std::fmt::Debug;
use typed_builder::TypedBuilder;

#[cfg(feature = "python")]
use pyo3::{exceptions::PyValueError, prelude::*, types::PyBytes, types::PyType};

/// Defines the handoff from a current ground station to the next one that is visible to prevent overlapping of measurements
#[derive(Copy, Clone, Debug, Deserialize, PartialEq, Serialize, Default, Enumerated)]
#[repr(u8)]
#[cfg_attr(feature = "python", pyclass(eq, eq_int))]
pub enum Handoff {
    /// If a new station is in visibility of the spacecraft, the "Eager" station will immediately stop tracking and switch over (default)
    #[default]
    Eager = 0,
    /// If a new station is in visibility of the spacecraft, the "Greedy" station will continue to tracking until the vehicle is below its elevation mask
    Greedy = 1,
    /// If a new station is in visibility of the spacecraft, the "Overlap" station will continue tracking, and so will the other one
    Overlap = 2,
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl Handoff {
    fn __repr__(&self) -> String {
        format!("{self:?}")
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }
}

/// A scheduler allows building a scheduling of spaceraft tracking for a set of ground stations.
#[derive(Copy, Clone, Debug, Default, Deserialize, PartialEq, Serialize, TypedBuilder)]
#[cfg_attr(feature = "python", pyclass)]
#[builder(doc)]
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
#[derive(Copy, Clone, Deserialize, PartialEq, Serialize, Default)]
pub enum Cadence {
    #[default]
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

#[cfg(feature = "python")]
#[pyclass(name = "Cadence")]
#[derive(Clone, Debug)]
pub struct PyCadence {
    pub inner: Cadence,
}

impl<'a> Decode<'a> for Cadence {
    fn decode<R: Reader<'a>>(decoder: &mut R) -> der::Result<Self> {
        let tag = decoder.decode::<u8>()?;
        match tag {
            0 => Ok(Self::Continuous),
            1 => {
                let on_ns = decoder.decode::<i128>()?;
                let off_ns = decoder.decode::<i128>()?;
                Ok(Self::Intermittent {
                    on: Duration::from_total_nanoseconds(on_ns),
                    off: Duration::from_total_nanoseconds(off_ns),
                })
            }
            _ => Err(der::ErrorKind::Value {
                tag: der::Tag::Integer,
            }
            .into()),
        }
    }
}

impl Encode for Cadence {
    fn encoded_len(&self) -> der::Result<der::Length> {
        match self {
            Self::Continuous => 0u8.encoded_len(),
            Self::Intermittent { on, off } => {
                1u8.encoded_len()?
                    + on.total_nanoseconds().encoded_len()?
                    + off.total_nanoseconds().encoded_len()?
            }
        }
    }

    fn encode(&self, encoder: &mut impl der::Writer) -> der::Result<()> {
        match self {
            Self::Continuous => 0u8.encode(encoder),
            Self::Intermittent { on, off } => {
                1u8.encode(encoder)?;
                on.total_nanoseconds().encode(encoder)?;
                off.total_nanoseconds().encode(encoder)
            }
        }
    }
}

#[cfg(feature = "python")]
#[pymethods]
impl PyCadence {
    #[staticmethod]
    fn continuous() -> Self {
        Self {
            inner: Cadence::Continuous,
        }
    }

    #[staticmethod]
    fn intermittent(on: Duration, off: Duration) -> Self {
        Self {
            inner: Cadence::Intermittent { on, off },
        }
    }

    fn __repr__(&self) -> String {
        format!("{:?}", self.inner)
    }

    fn __str__(&self) -> String {
        format!("{:?}", self.inner)
    }

    /// Decodes an ASN.1 DER encoded byte array into a Cadence object.
    ///
    /// :type data: bytes
    /// :rtype: Cadence
    #[classmethod]
    pub fn from_asn1(_cls: &Bound<'_, PyType>, data: &[u8]) -> PyResult<Self> {
        match Cadence::from_der(data) {
            Ok(obj) => Ok(Self { inner: obj }),
            Err(e) => Err(PyValueError::new_err(format!("ASN.1 decoding error: {e}"))),
        }
    }

    /// Encodes this Cadence object into an ASN.1 DER encoded byte array.
    ///
    /// :rtype: bytes
    pub fn to_asn1<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let mut buf = Vec::new();
        match self.inner.encode_to_vec(&mut buf) {
            Ok(_) => Ok(PyBytes::new(py, &buf)),
            Err(e) => Err(PyValueError::new_err(format!("ASN.1 encoding error: {e}"))),
        }
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

impl<'a> Decode<'a> for Scheduler {
    fn decode<R: Reader<'a>>(decoder: &mut R) -> der::Result<Self> {
        let handoff = decoder.decode()?;
        let cadence = decoder.decode()?;
        let min_samples = decoder.decode()?;
        let sample_alignment_ns = if decoder.decode::<bool>()? {
            Some(decoder.decode::<i128>()?)
        } else {
            None
        };

        Ok(Self {
            handoff,
            cadence,
            min_samples,
            sample_alignment: sample_alignment_ns.map(Duration::from_total_nanoseconds),
        })
    }
}

impl Encode for Scheduler {
    fn encoded_len(&self) -> der::Result<der::Length> {
        let mut len = (self.handoff.encoded_len()?
            + self.cadence.encoded_len()?
            + self.min_samples.encoded_len()?
            + self.sample_alignment.is_some().encoded_len()?)?;

        if let Some(sa) = self.sample_alignment {
            len = (len + sa.total_nanoseconds().encoded_len()?)?;
        }
        Ok(len)
    }

    fn encode(&self, encoder: &mut impl der::Writer) -> der::Result<()> {
        self.handoff.encode(encoder)?;
        self.cadence.encode(encoder)?;
        self.min_samples.encode(encoder)?;
        if let Some(sa) = self.sample_alignment {
            true.encode(encoder)?;
            sa.total_nanoseconds().encode(encoder)?;
        } else {
            false.encode(encoder)?;
        }
        Ok(())
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl Scheduler {
    #[new]
    #[pyo3(signature = (handoff=Handoff::Eager, cadence=None, min_samples=10, sample_alignment=None))]
    fn py_new(
        handoff: Handoff,
        cadence: Option<PyCadence>,
        min_samples: u32,
        sample_alignment: Option<Duration>,
    ) -> Self {
        Self {
            handoff,
            cadence: cadence.map(|c| c.inner).unwrap_or_default(),
            min_samples,
            sample_alignment,
        }
    }

    #[getter]
    fn get_handoff(&self) -> Handoff {
        self.handoff
    }

    #[setter]
    fn set_handoff(&mut self, handoff: Handoff) {
        self.handoff = handoff;
    }

    #[getter]
    fn get_cadence(&self) -> PyCadence {
        PyCadence {
            inner: self.cadence,
        }
    }

    #[setter]
    fn set_cadence(&mut self, cadence: PyCadence) {
        self.cadence = cadence.inner;
    }

    #[getter]
    fn get_min_samples(&self) -> u32 {
        self.min_samples
    }

    #[setter]
    fn set_min_samples(&mut self, min_samples: u32) {
        self.min_samples = min_samples;
    }

    #[getter]
    fn get_sample_alignment(&self) -> Option<Duration> {
        self.sample_alignment
    }

    #[setter]
    fn set_sample_alignment(&mut self, sample_alignment: Option<Duration>) {
        self.sample_alignment = sample_alignment;
    }

    fn __repr__(&self) -> String {
        format!("{self:?}")
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    /// Decodes an ASN.1 DER encoded byte array into a Scheduler object.
    ///
    /// :type data: bytes
    /// :rtype: Scheduler
    #[classmethod]
    pub fn from_asn1(_cls: &Bound<'_, PyType>, data: &[u8]) -> PyResult<Self> {
        match Self::from_der(data) {
            Ok(obj) => Ok(obj),
            Err(e) => Err(PyValueError::new_err(format!("ASN.1 decoding error: {e}"))),
        }
    }

    /// Encodes this Scheduler object into an ASN.1 DER encoded byte array.
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
mod scheduler_ut {
    use process::simulator::scheduler::Handoff;

    use crate::od::prelude::*;

    use super::Scheduler;

    #[test]
    fn serde_cadence() {
        use hifitime::TimeUnits;
        use serde_yml;

        let cont: Cadence = serde_yml::from_str("!Continuous").unwrap();
        assert_eq!(cont, Cadence::Continuous);

        let int: Cadence =
            serde_yml::from_str("!Intermittent {on: 1 h 35 min, off: 15 h 02 min 3 s}").unwrap();
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

        let serialized = serde_yml::to_string(&int).unwrap();
        let deserd: Cadence = serde_yml::from_str(&serialized).unwrap();
        assert_eq!(deserd, int);
    }

    #[test]
    fn api_and_serde_scheduler() {
        use hifitime::TimeUnits;
        use serde_yml;

        let scheduler = Scheduler::default();
        let serialized = serde_yml::to_string(&scheduler).unwrap();
        assert_eq!(
            serialized,
            "handoff: Eager\ncadence: Continuous\nmin_samples: 0\nsample_alignment: null\n"
        );
        let deserd: Scheduler = serde_yml::from_str(&serialized).unwrap();
        assert_eq!(deserd, scheduler);

        let scheduler = Scheduler::builder()
            .handoff(Handoff::Eager)
            .cadence(Cadence::Intermittent {
                on: 0.2.hours(),
                off: 17.hours() + 5.minutes(),
            })
            .build();

        let serialized = serde_yml::to_string(&scheduler).unwrap();
        assert_eq!(
            serialized,
            "handoff: Eager\ncadence: !Intermittent\n  'on': '12 min'\n  'off': '17 h 5 min'\nmin_samples: 10\nsample_alignment: '1 s'\n"
        );
        let deserd: Scheduler = serde_yml::from_str(&serialized).unwrap();
        assert_eq!(deserd, scheduler);
    }

    #[test]
    fn defaults() {
        let sched = Scheduler::default();

        assert_eq!(sched.cadence, Cadence::Continuous);

        assert_eq!(sched.handoff, Handoff::Eager);
    }
}
