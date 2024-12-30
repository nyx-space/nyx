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

use crate::errors::NyxError;
use crate::md::StateParameter;
use crate::time::Epoch;

use arrow::error::ArrowError;
use parquet::errors::ParquetError;
use snafu::prelude::*;
pub(crate) mod watermark;
use hifitime::prelude::{Format, Formatter};
use hifitime::Duration;
use serde::de::DeserializeOwned;
use serde::{Deserialize, Deserializer};
use serde::{Serialize, Serializer};
use serde_yml::Error as YamlError;
use std::collections::{BTreeMap, HashMap};
use std::convert::From;
use std::fmt::Debug;
use std::fs::File;
use std::io::BufReader;
use std::io::Error as IoError;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use typed_builder::TypedBuilder;

/// Handles loading of gravity models using files of NASA PDS and GMAT COF. Several gunzipped files are provided with nyx.
pub mod gravity;

use std::io;

/// Configuration for exporting a trajectory to parquet.
#[derive(Clone, Default, Serialize, Deserialize, TypedBuilder)]
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

#[derive(Debug, Snafu)]
#[snafu(visibility(pub(crate)))]
pub enum ConfigError {
    #[snafu(display("failed to read configuration file: {source}"))]
    ReadError { source: io::Error },

    #[snafu(display("failed to parse YAML configuration file: {source}"))]
    ParseError { source: serde_yml::Error },

    #[snafu(display("of invalid configuration: {msg}"))]
    InvalidConfig { msg: String },
}

impl PartialEq for ConfigError {
    /// No two configuration errors match
    fn eq(&self, _other: &Self) -> bool {
        false
    }
}

#[derive(Debug, Snafu)]
#[snafu(visibility(pub(crate)))]
pub enum InputOutputError {
    #[snafu(display("{action} encountered i/o error: {source}"))]
    StdIOError {
        source: io::Error,
        action: &'static str,
    },
    #[snafu(display("missing required data {which}"))]
    MissingData { which: String },
    #[snafu(display("unknown data `{which}`"))]
    UnsupportedData { which: String },
    #[snafu(display("{action} encountered a Parquet error: {source}"))]
    ParquetError {
        source: ParquetError,
        action: &'static str,
    },
    #[snafu(display("inconsistency detected: {msg}"))]
    Inconsistency { msg: String },
    #[snafu(display("{action} encountered an Arrow error: {source}"))]
    ArrowError {
        source: ArrowError,
        action: &'static str,
    },
    #[snafu(display("error parsing `{data}` as Dhall config: {err}"))]
    ParseDhall { data: String, err: String },
    #[snafu(display("error serializing {what} to Dhall: {err}"))]
    SerializeDhall { what: String, err: String },
    #[snafu(display("empty dataset error when (de)serializing for {action}"))]
    EmptyDataset { action: &'static str },
}

impl PartialEq for InputOutputError {
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
        let file = File::open(path).context(ReadSnafu)?;
        let reader = BufReader::new(file);

        serde_yml::from_reader(reader).context(ParseSnafu)
    }

    /// Builds a sequence of "Selves" from the provided path to a yaml
    fn load_many<P>(path: P) -> Result<Vec<Self>, ConfigError>
    where
        P: AsRef<Path>,
    {
        let file = File::open(path).context(ReadSnafu)?;
        let reader = BufReader::new(file);

        serde_yml::from_reader(reader).context(ParseSnafu)
    }

    /// Builds a map of names to "selves" from the provided path to a yaml
    fn load_named<P>(path: P) -> Result<BTreeMap<String, Self>, ConfigError>
    where
        P: AsRef<Path>,
    {
        let file = File::open(path).context(ReadSnafu)?;
        let reader = BufReader::new(file);

        serde_yml::from_reader(reader).context(ParseSnafu)
    }

    /// Builds a sequence of "Selves" from the provided string of a yaml
    fn loads_many(data: &str) -> Result<Vec<Self>, ConfigError> {
        debug!("Loading YAML:\n{data}");
        serde_yml::from_str(data).context(ParseSnafu)
    }

    /// Builds a sequence of "Selves" from the provided string of a yaml
    fn loads_named(data: &str) -> Result<BTreeMap<String, Self>, ConfigError> {
        debug!("Loading YAML:\n{data}");
        serde_yml::from_str(data).context(ParseSnafu)
    }
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
