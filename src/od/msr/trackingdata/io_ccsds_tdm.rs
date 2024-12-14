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

use crate::io::ExportCfg;
use crate::io::{InputOutputError, StdIOSnafu};
use crate::od::msr::{Measurement, MeasurementType};
use hifitime::prelude::Epoch;
use hifitime::TimeScale;
use snafu::ResultExt;
use std::collections::{BTreeMap, HashMap};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use super::TrackingDataArc;

impl TrackingDataArc {
    /// Loads a tracking arc from its serialization in CCSDS TDM.
    pub fn from_tdm<P: AsRef<Path>>(
        path: P,
        aliases: Option<HashMap<String, String>>,
    ) -> Result<Self, InputOutputError> {
        let file = File::open(&path).context(StdIOSnafu {
            action: "opening CCSDS TDM file for tracking arc",
        })?;

        let mut measurements = BTreeMap::new();

        let reader = BufReader::new(file);

        let mut in_data_section = false;
        let mut current_tracker = String::new();
        let mut time_system = TimeScale::UTC;

        for line in reader.lines() {
            let line = line.context(StdIOSnafu {
                action: "reading CCSDS TDM file",
            })?;
            let line = line.trim();

            if line == "DATA_START" {
                in_data_section = true;
                continue;
            } else if line == "DATA_STOP" {
                in_data_section = false;
            }

            if !in_data_section {
                if line.starts_with("PARTICIPANT_1") {
                    current_tracker = line.split('=').nth(1).unwrap_or("").trim().to_string();
                    // If aliases are provided, try to map them.
                    if let Some(aliases) = &aliases {
                        if let Some(alias) = aliases.get(&current_tracker) {
                            current_tracker = alias.clone();
                        }
                    }
                }
                if line.starts_with("TIME_SYSTEM") {
                    let ts = line.split('=').nth(1).unwrap_or("UTC").trim();
                    match ts {
                        "UTC" => time_system = TimeScale::UTC,
                        "TAI" => time_system = TimeScale::TAI,
                        "GPS" => time_system = TimeScale::GPST,
                        _ => {
                            return Err(InputOutputError::UnsupportedData {
                                which: format!("time scale `{ts}` not supported"),
                            })
                        }
                    }
                }
                continue;
            }

            if let Some((mtype, epoch, value)) = parse_measurement_line(line, time_system)? {
                measurements
                    .entry(epoch)
                    .or_insert_with(|| Measurement {
                        tracker: current_tracker.clone(),
                        epoch,
                        data: HashMap::new(),
                    })
                    .data
                    .insert(mtype, value);
            }
        }

        Ok(Self {
            measurements,
            source: Some(path.as_ref().to_path_buf().display().to_string()),
        })
    }

    /// Store this tracking arc to a CCSDS TDM file, with optional metadata and a timestamp appended to the filename.
    pub fn to_tdm<P: AsRef<Path>>(
        &self,
        path: P,
        cfg: ExportCfg,
    ) -> Result<PathBuf, Box<dyn Error>> {
        todo!()
    }
}

fn parse_measurement_line(
    line: &str,
    time_system: TimeScale,
) -> Result<Option<(MeasurementType, Epoch, f64)>, InputOutputError> {
    let parts: Vec<&str> = line.split('=').collect();
    if parts.len() != 2 {
        return Ok(None);
    }

    let (mtype_str, data) = (parts[0].trim(), parts[1].trim());
    let mtype = match mtype_str {
        "RANGE" => MeasurementType::Range,
        "DOPPLER_INSTANTANEOUS" | "DOPPLER_INTEGRATED" => MeasurementType::Doppler,
        "ANGLE_1" => MeasurementType::Azimuth,
        "ANGLE_2" => MeasurementType::Elevation,
        _ => {
            return Err(InputOutputError::UnsupportedData {
                which: mtype_str.to_string(),
            })
        }
    };

    let data_parts: Vec<&str> = data.split_whitespace().collect();
    if data_parts.len() != 2 {
        return Ok(None);
    }

    let epoch =
        Epoch::from_gregorian_str(&format!("{} {time_system}", data_parts[0])).map_err(|e| {
            InputOutputError::Inconsistency {
                msg: format!("{e} when parsing epoch"),
            }
        })?;

    let value = data_parts[1]
        .parse::<f64>()
        .map_err(|e| InputOutputError::UnsupportedData {
            which: format!("`{}` is not a float: {e}", data_parts[1]),
        })?;

    Ok(Some((mtype, epoch, value)))
}

#[cfg(test)]
mod ut_tdm {
    extern crate pretty_env_logger as pel;
    use hifitime::TimeUnits;
    use std::{collections::HashMap, path::PathBuf};

    use crate::od::msr::TrackingDataArc;

    #[test]
    fn test_tdm_no_alias() {
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "output_data",
            "demo_tracking_arc.tdm",
        ]
        .iter()
        .collect();

        let tdm = TrackingDataArc::from_tdm(path, None).unwrap();
        println!("{tdm}");
    }

    #[test]
    fn test_tdm_with_alias() {
        pel::init();
        let mut aliases = HashMap::new();
        aliases.insert("DSS-65 Madrid".to_string(), "DSN Madrid".to_string());

        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "output_data",
            "demo_tracking_arc.tdm",
        ]
        .iter()
        .collect();

        let tdm = TrackingDataArc::from_tdm(path, Some(aliases)).unwrap();
        println!("{tdm}");

        let tdm_failed_downsample = tdm.clone().downsample(10.seconds());
        assert_eq!(
            tdm_failed_downsample.len(),
            tdm.len(),
            "downsampling should have failed because it's upsampling"
        );

        let tdm_downsample = tdm.clone().downsample(10.minutes());
        println!("{tdm_downsample}");
        assert_eq!(
            tdm_downsample.len(),
            tdm.len() / 10 + 1,
            "downsampling has wrong sample count"
        );
    }
}
