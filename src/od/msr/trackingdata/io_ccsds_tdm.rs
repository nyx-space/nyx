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

use crate::io::watermark::prj_name_ver;
use crate::io::ExportCfg;
use crate::io::{InputOutputError, StdIOSnafu};
use crate::od::msr::{Measurement, MeasurementType};
use hifitime::efmt::{Format, Formatter};
use hifitime::prelude::Epoch;
use hifitime::TimeScale;
use indexmap::IndexMap;
use snafu::ResultExt;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::str::FromStr;

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
                        data: IndexMap::new(),
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
    pub fn to_tdm_file<P: AsRef<Path>>(
        mut self,
        spacecraft_name: String,
        aliases: Option<HashMap<String, String>>,
        path: P,
        cfg: ExportCfg,
    ) -> Result<PathBuf, InputOutputError> {
        if self.is_empty() {
            return Err(InputOutputError::MissingData {
                which: " - empty tracking data cannot be exported to TDM".to_string(),
            });
        }

        // Filter epochs if needed.
        if cfg.start_epoch.is_some() && cfg.end_epoch.is_some() {
            self = self.filter_by_epoch(cfg.start_epoch.unwrap()..cfg.end_epoch.unwrap());
        } else if cfg.start_epoch.is_some() {
            self = self.filter_by_epoch(cfg.start_epoch.unwrap()..);
        } else if cfg.end_epoch.is_some() {
            self = self.filter_by_epoch(..cfg.end_epoch.unwrap());
        }

        let tick = Epoch::now().unwrap();
        info!("Exporting tracking data to CCSDS TDM file...");

        // Grab the path here before we move stuff.
        let path_buf = cfg.actual_path(path);

        let metadata = cfg.metadata.unwrap_or_default();

        let file = File::create(&path_buf).context(StdIOSnafu {
            action: "creating CCSDS TDM file for tracking arc",
        })?;
        let mut writer = BufWriter::new(file);

        let err_hdlr = |source| InputOutputError::StdIOError {
            source,
            action: "writing data to TDM file",
        };

        // Epoch formmatter.
        let iso8601_no_ts = Format::from_str("%Y-%m-%dT%H:%M:%S.%f").unwrap();

        // Write mandatory metadata
        writeln!(writer, "CCSDS_TDM_VERS = 2.0").map_err(err_hdlr)?;
        writeln!(
            writer,
            "\nCOMMENT Build by {} -- https://nyxspace.com",
            prj_name_ver()
        )
        .map_err(err_hdlr)?;
        writeln!(
            writer,
            "COMMENT Nyx Space provided under the AGPL v3 open source license -- https://nyxspace.com/pricing\n"
        )
        .map_err(err_hdlr)?;
        writeln!(
            writer,
            "CREATION_DATE = {}",
            Formatter::new(Epoch::now().unwrap(), iso8601_no_ts)
        )
        .map_err(err_hdlr)?;
        writeln!(
            writer,
            "ORIGINATOR = {}\n",
            metadata
                .get("originator")
                .unwrap_or(&"Nyx Space".to_string())
        )
        .map_err(err_hdlr)?;

        // Create a new meta section for each tracker
        // Get unique trackers and process each one separately
        let trackers = self.unique_aliases();

        for tracker in trackers {
            let tracker_data = self.clone().filter_by_tracker(tracker.clone());

            writeln!(writer, "META_START").map_err(err_hdlr)?;
            writeln!(writer, "\tTIME_SYSTEM = UTC").map_err(err_hdlr)?;
            writeln!(
                writer,
                "\tSTART_TIME = {}",
                Formatter::new(tracker_data.start_epoch().unwrap(), iso8601_no_ts)
            )
            .map_err(err_hdlr)?;
            writeln!(
                writer,
                "\tSTOP_TIME = {}",
                Formatter::new(tracker_data.end_epoch().unwrap(), iso8601_no_ts)
            )
            .map_err(err_hdlr)?;

            writeln!(
                writer,
                "\tPARTICIPANT_1 = {}",
                if let Some(aliases) = &aliases {
                    if let Some(alias) = aliases.get(&tracker) {
                        alias
                    } else {
                        &tracker
                    }
                } else {
                    &tracker
                }
            )
            .map_err(err_hdlr)?;

            writeln!(writer, "\tPARTICIPANT_2 = {spacecraft_name}").map_err(err_hdlr)?;

            // Add additional metadata, could include timetag ref for example.
            for (k, v) in &metadata {
                if k != "originator" {
                    writeln!(writer, "\t{k} = {v}").map_err(err_hdlr)?;
                }
            }

            if self.unique_types().contains(&MeasurementType::Range) {
                writeln!(writer, "\tRANGE_UNITS = km").map_err(err_hdlr)?;
            }

            if self.unique_types().contains(&MeasurementType::Azimuth)
                || self.unique_types().contains(&MeasurementType::Elevation)
            {
                writeln!(writer, "\tANGLE_TYPE = AZEL").map_err(err_hdlr)?;
            }

            writeln!(writer, "META_STOP\n").map_err(err_hdlr)?;

            // Write the data section
            writeln!(writer, "DATA_START").map_err(err_hdlr)?;

            // Process measurements for this tracker
            for (epoch, measurement) in tracker_data.measurements {
                for (mtype, value) in &measurement.data {
                    let type_str = match mtype {
                        MeasurementType::Range => "RANGE",
                        MeasurementType::Doppler => "DOPPLER_INTEGRATED",
                        MeasurementType::Azimuth => "ANGLE_1",
                        MeasurementType::Elevation => "ANGLE_2",
                    };

                    writeln!(
                        writer,
                        "\t{:<20} = {:<23}\t{:.12}",
                        type_str,
                        Formatter::new(epoch, iso8601_no_ts),
                        value
                    )
                    .map_err(err_hdlr)?;
                }
            }

            writeln!(writer, "DATA_STOP\n").map_err(err_hdlr)?;
        }

        #[allow(clippy::writeln_empty_string)]
        writeln!(writer, "").map_err(err_hdlr)?;

        // Return the path this was written to
        let tock_time = Epoch::now().unwrap() - tick;
        info!("CCSDS TDM written to {} in {tock_time}", path_buf.display());
        Ok(path_buf)
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
