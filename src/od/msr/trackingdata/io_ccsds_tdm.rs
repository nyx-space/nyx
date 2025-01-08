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
use anise::constants::SPEED_OF_LIGHT_KM_S;
use hifitime::efmt::{Format, Formatter};
use hifitime::prelude::Epoch;
use hifitime::TimeScale;
use indexmap::{IndexMap, IndexSet};
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

        let source = path.as_ref().to_path_buf().display().to_string();
        info!("parsing CCSDS TDM {source}");

        let mut measurements = BTreeMap::new();
        let mut metadata = HashMap::new();

        let reader = BufReader::new(file);

        let mut in_data_section = false;
        let mut current_tracker = String::new();
        let mut time_system = TimeScale::UTC;
        let mut has_freq_data = false;
        let mut msr_divider = 1.0;

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
                } else if line.starts_with("TIME_SYSTEM") {
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
                } else if line.starts_with("PATH") {
                    match line.split(",").count() {
                        2 => msr_divider = 1.0,
                        3 => msr_divider = 2.0,
                        cnt => {
                            return Err(InputOutputError::UnsupportedData {
                                which: format!(
                                    "found {cnt} paths in TDM, only 1 or 2 are supported"
                                ),
                            })
                        }
                    }
                }

                let mut splt = line.split('=');
                if let Some(keyword) = splt.nth(0) {
                    // Get the zeroth item again since we've consumed the first zeroth one.
                    if let Some(value) = splt.nth(0) {
                        metadata.insert(keyword.trim().to_string(), value.trim().to_string());
                    }
                }

                continue;
            }

            if let Some((mtype, epoch, value)) = parse_measurement_line(line, time_system)? {
                if [
                    MeasurementType::ReceiveFrequency,
                    MeasurementType::TransmitFrequency,
                ]
                .contains(&mtype)
                {
                    has_freq_data = true;
                    // Don't modify the values.
                    msr_divider = 1.0;
                }
                measurements
                    .entry(epoch)
                    .or_insert_with(|| Measurement {
                        tracker: current_tracker.clone(),
                        epoch,
                        data: IndexMap::new(),
                    })
                    .data
                    .insert(mtype, value / msr_divider);
            }
        }

        let mut turnaround_ratio = None;
        let drop_freq_data;
        if has_freq_data {
            // If there is any frequency measurement, compute the turn-around ratio.
            if let Some(ta_num_str) = metadata.get("TURNAROUND_NUMERATOR") {
                if let Some(ta_denom_str) = metadata.get("TURNAROUND_DENOMINATOR") {
                    if let Ok(ta_num) = ta_num_str.parse::<i32>() {
                        if let Ok(ta_denom) = ta_denom_str.parse::<i32>() {
                            // turn-around ratio is set.
                            turnaround_ratio = Some(f64::from(ta_num) / f64::from(ta_denom));
                            info!("turn-around ratio is {ta_num}/{ta_denom}");
                            drop_freq_data = false;
                        } else {
                            error!("turn-around denominator `{ta_denom_str}` is not a valid double precision float");
                            drop_freq_data = true;
                        }
                    } else {
                        error!("turn-around numerator `{ta_num_str}` is not a valid double precision float");
                        drop_freq_data = true;
                    }
                } else {
                    error!("required turn-around denominator missing from metadata -- dropping ALL RECEIVE/TRANSMIT data");
                    drop_freq_data = true;
                }
            } else {
                error!("required turn-around numerator missing from metadata -- dropping ALL RECEIVE/TRANSMIT data");
                drop_freq_data = true;
            }
        } else {
            drop_freq_data = true;
        }

        // Now, let's convert the receive and transmit frequencies to Doppler measurements in velocity units.
        // We expect the transmit and receive frequencies to have the exact same timestamp.
        let mut freq_types = IndexSet::new();
        freq_types.insert(MeasurementType::ReceiveFrequency);
        freq_types.insert(MeasurementType::TransmitFrequency);
        let mut latest_transmit_freq = None;
        for (epoch, measurement) in measurements.iter_mut() {
            if drop_freq_data {
                for freq in &freq_types {
                    measurement.data.swap_remove(freq);
                }
                continue;
            }

            let avail = measurement.availability(&freq_types);
            let use_prev_transmit_freq;
            let num_freq_msr = avail
                .iter()
                .copied()
                .map(|v| if v { 1 } else { 0 })
                .sum::<u8>();
            if num_freq_msr == 0 {
                // No frequency measurements
                continue;
            } else if num_freq_msr == 1 {
                // avail[0] means that Receive Freq is available
                // avail[1] means that Transmit Freq is available
                // We can only compute Doppler data from one data point if that data point
                // if the receive frequency and the transmit frequency was previously set.
                if latest_transmit_freq.is_some() && avail[0] {
                    use_prev_transmit_freq = true;
                    warn!(
                        "no transmit frequency at {epoch}, using previous value of {} Hz",
                        latest_transmit_freq.unwrap()
                    );
                } else {
                    warn!("only one of receive or transmit frequencies found at {epoch}, ignoring");
                    for freq in &freq_types {
                        measurement.data.swap_remove(freq);
                    }
                    continue;
                }
            } else {
                use_prev_transmit_freq = false;
            }

            if !use_prev_transmit_freq {
                // Update the latest transmit frequency since it's set.
                latest_transmit_freq = Some(
                    *measurement
                        .data
                        .get(&MeasurementType::TransmitFrequency)
                        .unwrap(),
                );
            }

            let transmit_freq_hz = latest_transmit_freq.unwrap();
            let receive_freq_hz = *measurement
                .data
                .get(&MeasurementType::ReceiveFrequency)
                .unwrap();

            // Compute the Doppler shift, equation from section 3.5.2.8.2 of CCSDS TDM v2 specs
            let doppler_shift_hz = transmit_freq_hz * turnaround_ratio.unwrap() - receive_freq_hz;
            // Compute the expected Doppler measurement as range-rate.
            let rho_dot_km_s = (doppler_shift_hz * SPEED_OF_LIGHT_KM_S)
                / (2.0 * transmit_freq_hz * turnaround_ratio.unwrap());

            // Finally, replace the frequency data with a Doppler measurement.
            for freq in &freq_types {
                measurement.data.swap_remove(freq);
            }
            measurement
                .data
                .insert(MeasurementType::Doppler, rho_dot_km_s);
        }

        let moduli = if let Some(range_modulus) = metadata.get("RANGE_MODULUS") {
            if let Ok(value) = range_modulus.parse::<f64>() {
                let mut modulos = IndexMap::new();
                modulos.insert(MeasurementType::Range, value);
                // Only range modulus exists in TDM files.
                Some(modulos)
            } else {
                warn!("could not parse RANGE_MODULO of `{range_modulus}` as a double");
                None
            }
        } else {
            None
        };

        let trk = Self {
            measurements,
            source: Some(source),
            moduli,
        };

        if trk.unique_types().is_empty() {
            Err(InputOutputError::EmptyDataset {
                action: "CCSDS TDM file",
            })
        } else {
            Ok(trk)
        }
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

        // Create a new meta section for each tracker and for each measurement type that is one or two way.
        // Get unique trackers and process each one separately
        let trackers = self.unique_aliases();

        for tracker in trackers {
            let tracker_data = self.clone().filter_by_tracker(tracker.clone());

            let types = tracker_data.unique_types();

            let two_way_types = types
                .iter()
                .filter(|msr_type| msr_type.may_be_two_way())
                .copied()
                .collect::<Vec<_>>();

            let one_way_types = types
                .iter()
                .filter(|msr_type| !msr_type.may_be_two_way())
                .copied()
                .collect::<Vec<_>>();

            // Add the two-way data first.
            for (tno, types) in [two_way_types, one_way_types].iter().enumerate() {
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

                let multiplier = if tno == 0 {
                    writeln!(writer, "\tPATH = 1,2,1").map_err(err_hdlr)?;
                    2.0
                } else {
                    writeln!(writer, "\tPATH = 1,2").map_err(err_hdlr)?;
                    1.0
                };

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

                writeln!(writer, "\tMODE = SEQUENTIAL").map_err(err_hdlr)?;

                // Add additional metadata, could include timetag ref for example.
                for (k, v) in &metadata {
                    if k != "originator" {
                        writeln!(writer, "\t{k} = {v}").map_err(err_hdlr)?;
                    }
                }

                if types.contains(&MeasurementType::Range) {
                    writeln!(writer, "\tRANGE_UNITS = km").map_err(err_hdlr)?;

                    if let Some(moduli) = &self.moduli {
                        if let Some(range_modulus) = moduli.get(&MeasurementType::Range) {
                            writeln!(writer, "\tRANGE_MODULUS = {range_modulus:E}")
                                .map_err(err_hdlr)?;
                        }
                    }
                }

                if types.contains(&MeasurementType::Azimuth)
                    || types.contains(&MeasurementType::Elevation)
                {
                    writeln!(writer, "\tANGLE_TYPE = AZEL").map_err(err_hdlr)?;
                }

                writeln!(writer, "META_STOP\n").map_err(err_hdlr)?;

                // Write the data section
                writeln!(writer, "DATA_START").map_err(err_hdlr)?;

                // Process measurements for this tracker
                for (epoch, measurement) in &tracker_data.measurements {
                    for (mtype, value) in &measurement.data {
                        if !types.contains(mtype) {
                            continue;
                        }
                        let type_str = match mtype {
                            MeasurementType::Range => "RANGE",
                            MeasurementType::Doppler => "DOPPLER_INTEGRATED",
                            MeasurementType::Azimuth => "ANGLE_1",
                            MeasurementType::Elevation => "ANGLE_2",
                            MeasurementType::ReceiveFrequency => "RECEIVE_FREQ",
                            MeasurementType::TransmitFrequency => "TRANSMIT_FREQ",
                        };

                        writeln!(
                            writer,
                            "\t{:<20} = {:<23}\t{:.12}",
                            type_str,
                            Formatter::new(*epoch, iso8601_no_ts),
                            value * multiplier
                        )
                        .map_err(err_hdlr)?;
                    }
                }

                writeln!(writer, "DATA_STOP\n").map_err(err_hdlr)?;
            }
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
        "RECEIVE_FREQ" | "RECEIVE_FREQ_1" | "RECEIVE_FREQ_2" | "RECEIVE_FREQ_3"
        | "RECEIVE_FREQ_4" | "RECEIVE_FREQ_5" => MeasurementType::ReceiveFrequency,
        "TRANSMIT_FREQ" | "TRANSMIT_FREQ_1" | "TRANSMIT_FREQ_2" | "TRANSMIT_FREQ_3"
        | "TRANSMIT_FREQ_4" | "TRANSMIT_FREQ_5" => MeasurementType::TransmitFrequency,
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
