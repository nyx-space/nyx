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

use crate::errors::NyxError;
use crate::md::StateParameter;
use crate::time::Epoch;
use crate::Orbit;
pub(crate) mod watermark;
use hifitime::prelude::{Format, Formatter};
use hifitime::Duration;
use serde::de::DeserializeOwned;
use serde::ser::SerializeSeq;
use serde::{Deserialize, Deserializer};
use serde::{Serialize, Serializer};
use serde_yaml::Error as YamlError;
use std::collections::{BTreeMap, HashMap};
use std::convert::From;
use std::fmt::Debug;
use std::fs::File;
use std::io::BufReader;
use std::io::Error as IoError;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::Arc;
use typed_builder::TypedBuilder;

use self::orbit::OrbitSerde;
use crate::cosmic::{Cosm, Frame};

/// Handles writing to an XYZV file
pub mod cosmo;
pub mod dynamics;
pub mod estimate;
/// Handles reading from frames defined in input files
pub mod frame_serde;
/// Handles loading of gravity models using files of NASA PDS and GMAT COF. Several gunzipped files are provided with nyx.
pub mod gravity;
pub mod matrices;
pub mod orbit;
pub mod tracking_data;
pub mod trajectory_data;

use std::io;
use thiserror::Error;

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Configuration for exporting a trajectory to parquet.
#[derive(Clone, Default, Serialize, Deserialize, TypedBuilder)]
#[cfg_attr(feature = "python", pyclass)]
#[builder(doc)]
pub struct ExportCfg {
    /// Fields to export, if unset, defaults to all possible fields.
    #[builder(default, setter(strip_option))]
    pub fields: Option<Vec<StateParameter>>,
    /// Start epoch to export, defaults to the start of the trajectory
    #[builder(default, setter(strip_option))]
    pub start_epoch: Option<Epoch>,
    /// End epoch to export, defaults to the end of the trajectory
    #[builder(default, setter(strip_option))]
    pub end_epoch: Option<Epoch>,
    /// An optional step, defaults to every state in the trajectory (which likely isn't equidistant)
    #[builder(default, setter(strip_option))]
    pub step: Option<Duration>,
    /// Additional metadata to store in the Parquet metadata
    #[builder(default, setter(strip_option))]
    pub metadata: Option<HashMap<String, String>>,
    /// Set to true to append the timestamp to the filename
    #[builder(default)]
    pub timestamp: bool,
}

impl ExportCfg {
    /// Initialize a new configuration with the given metadata entries.
    pub fn from_metadata(metadata: Vec<(String, String)>) -> Self {
        let mut me = ExportCfg {
            metadata: Some(HashMap::new()),
            ..Default::default()
        };
        for (k, v) in metadata {
            me.metadata.as_mut().unwrap().insert(k, v);
        }
        me
    }

    /// Initialize a new default configuration but timestamp the filename.
    pub fn timestamped() -> Self {
        Self {
            timestamp: true,
            ..Default::default()
        }
    }

    pub fn append_field(&mut self, field: StateParameter) {
        if let Some(fields) = self.fields.as_mut() {
            fields.push(field);
        } else {
            self.fields = Some(vec![field]);
        }
    }

    /// Modifies the provided path to include the timestamp if required.
    pub(crate) fn actual_path<P: AsRef<Path>>(&self, path: P) -> PathBuf {
        let mut path_buf = path.as_ref().to_path_buf();
        if self.timestamp {
            if let Some(file_name) = path_buf.file_name() {
                if let Some(file_name_str) = file_name.to_str() {
                    if let Some(extension) = path_buf.extension() {
                        let stamp = Formatter::new(
                            Epoch::now().unwrap(),
                            Format::from_str("%Y-%m-%dT%H-%M-%S").unwrap(),
                        );
                        let ext = extension.to_str().unwrap();
                        let file_name = file_name_str.replace(&format!(".{ext}"), "");
                        let new_file_name = format!("{file_name}-{stamp}.{}", ext);
                        path_buf.set_file_name(new_file_name);
                    }
                }
            }
        };
        path_buf
    }
}

#[cfg(feature = "python")]
#[pymethods]
impl ExportCfg {
    #[new]
    #[pyo3(
        text_signature = "(timestamp=None, fields=None, start_epoch=None, step=None, end_epoch=None, metadata=None)"
    )]
    fn py_new(
        timestamp: Option<bool>,
        fields: Option<Vec<StateParameter>>,
        start_epoch: Option<Epoch>,
        end_epoch: Option<Epoch>,
        metadata: Option<HashMap<String, String>>,
    ) -> Self {
        Self {
            timestamp: timestamp.unwrap_or(false),
            fields,
            start_epoch,
            end_epoch,
            metadata,
            ..Default::default()
        }
    }
}

#[derive(Error, Debug)]
pub enum ConfigError {
    #[error("Failed to read configuration file: {0}")]
    ReadError(#[from] io::Error),

    #[error("Failed to parse YAML configuration file: {0}")]
    ParseError(#[source] serde_yaml::Error),

    #[error("Invalid configuration: {0}")]
    InvalidConfig(String),
}

impl PartialEq for ConfigError {
    /// No two configuration errors match
    fn eq(&self, _other: &Self) -> bool {
        false
    }
}

pub trait ConfigRepr: Debug + Sized + Serialize + DeserializeOwned {
    /// Builds the configuration representation from the path to a yaml
    fn load<P>(path: P) -> Result<Self, ConfigError>
    where
        P: AsRef<Path>,
    {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        serde_yaml::from_reader(reader).map_err(ConfigError::ParseError)
    }

    /// Builds a sequence of "Selves" from the provided path to a yaml
    fn load_many<P>(path: P) -> Result<Vec<Self>, ConfigError>
    where
        P: AsRef<Path>,
    {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        serde_yaml::from_reader(reader).map_err(ConfigError::ParseError)
    }

    /// Builds a map of names to "selves" from the provided path to a yaml
    fn load_named<P>(path: P) -> Result<BTreeMap<String, Self>, ConfigError>
    where
        P: AsRef<Path>,
    {
        let file = File::open(path)?;
        let reader = BufReader::new(file);

        serde_yaml::from_reader(reader).map_err(ConfigError::ParseError)
    }

    /// Builds a sequence of "Selves" from the provided string of a yaml
    fn loads_many(data: &str) -> Result<Vec<Self>, ConfigError> {
        debug!("Loading YAML:\n{data}");
        serde_yaml::from_str(data).map_err(ConfigError::ParseError)
    }
}

/// Trait to specify that a structure can be configured from a file, either in TOML, YAML, JSON, INI, etc.
pub trait Configurable
where
    Self: Sized,
{
    /// The intermediate representation needed to create `Self` or to serialize Self.
    type IntermediateRepr: ConfigRepr;

    fn from_yaml<P: AsRef<Path>>(path: P, cosm: Arc<Cosm>) -> Result<Self, ConfigError> {
        Self::from_config(Self::IntermediateRepr::load(path)?, cosm)
    }

    /// Creates a new instance of `self` from the configuration.
    fn from_config(cfg: Self::IntermediateRepr, cosm: Arc<Cosm>) -> Result<Self, ConfigError>
    where
        Self: Sized;

    /// Converts self into the intermediate representation which is serializable.
    fn to_config(&self) -> Result<Self::IntermediateRepr, ConfigError>;
}

pub(crate) fn epoch_to_str<S>(epoch: &Epoch, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.serialize_str(&format!("{epoch}"))
}

/// A deserializer from Epoch string
pub(crate) fn epoch_from_str<'de, D>(deserializer: D) -> Result<Epoch, D::Error>
where
    D: Deserializer<'de>,
{
    // implementation of the custom deserialization function
    let s = String::deserialize(deserializer)?;
    Epoch::from_str(&s).map_err(serde::de::Error::custom)
}

pub(crate) fn duration_to_str<S>(duration: &Duration, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.serialize_str(&format!("{duration}"))
}

/// A deserializer from Duration string
pub(crate) fn duration_from_str<'de, D>(deserializer: D) -> Result<Duration, D::Error>
where
    D: Deserializer<'de>,
{
    // implementation of the custom deserialization function
    let s = String::deserialize(deserializer)?;
    Duration::from_str(&s).map_err(serde::de::Error::custom)
}

pub(crate) fn frame_to_str<S>(frame: &Frame, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    serializer.serialize_str(&format!("{frame}"))
}

pub(crate) fn frame_from_str<'de, D>(deserializer: D) -> Result<Frame, D::Error>
where
    D: Deserializer<'de>,
{
    let s = String::deserialize(deserializer)?;
    // TODO: Figure out how to use DeserializeSeed here, but I'm not sure it would work. -- https://github.com/nyx-space/nyx/issues/86
    let cosm = Cosm::de438();
    cosm.try_frame(&s).map_err(serde::de::Error::custom)
}

pub(crate) fn frames_to_str<S>(frames: &Vec<Frame>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    let mut seq = serializer.serialize_seq(Some(frames.len()))?;
    for frame in frames {
        seq.serialize_element(&format!("{frame}"))?;
    }
    seq.end()
}

pub(crate) fn frames_from_str<'de, D>(deserializer: D) -> Result<Vec<Frame>, D::Error>
where
    D: Deserializer<'de>,
{
    let frame_names: Vec<String> = Vec::deserialize(deserializer)?;
    let cosm = Cosm::de438();
    let mut frames = Vec::new();
    for name in frame_names {
        frames.push(cosm.try_frame(&name).map_err(serde::de::Error::custom)?)
    }
    Ok(frames)
}

/// A deserializer from Epoch string
pub(crate) fn orbit_from_str<'de, D>(deserializer: D) -> Result<Orbit, D::Error>
where
    D: Deserializer<'de>,
{
    let orbit_serde: OrbitSerde = Deserialize::deserialize(deserializer)?;
    Ok(orbit_serde.into())
}

pub(crate) fn maybe_duration_to_str<S>(
    duration: &Option<Duration>,
    serializer: S,
) -> Result<S::Ok, S::Error>
where
    S: Serializer,
{
    if let Some(duration) = duration {
        duration_to_str(duration, serializer)
    } else {
        serializer.serialize_none()
    }
}

pub(crate) fn maybe_duration_from_str<'de, D>(deserializer: D) -> Result<Option<Duration>, D::Error>
where
    D: Deserializer<'de>,
{
    if let Ok(s) = String::deserialize(deserializer) {
        if let Ok(duration) = Duration::from_str(&s) {
            Ok(Some(duration))
        } else {
            Ok(None)
        }
    } else {
        Ok(None)
    }
}

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug)]
pub enum ParsingError {
    MD(String),
    OD(String),
    UseOdInstead,
    UseMdInstead,
    FileNotFound(String),
    FileNotUTF8(String),
    SequenceNotFound(String),
    LoadingError(String),
    PropagatorNotFound(String),
    Duration(String),
    Quantity(String),
    Distance(String),
    Velocity(String),
    IllDefined(String),
    ExecutionError(NyxError),
    IoError(IoError),
    Yaml(YamlError),
}

impl From<NyxError> for ParsingError {
    fn from(error: NyxError) -> Self {
        Self::ExecutionError(error)
    }
}
