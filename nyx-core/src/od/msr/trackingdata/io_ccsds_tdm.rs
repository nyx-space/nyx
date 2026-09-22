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
use crate::io::watermark::prj_name_ver;
use crate::io::{InputOutputError, StdIOSnafu};
use crate::od::ground_station::DopplerConfig;
use crate::od::msr::{IntegrationRef, Measurement, MeasurementType};
use anise::constants::SPEED_OF_LIGHT_KM_S;
use hifitime::efmt::{Format, Formatter};
use hifitime::{Duration, Epoch, TimeScale, Unit};
use indexmap::{IndexMap, IndexSet};
use log::{error, info, warn};
use snafu::ResultExt;
use std::collections::HashMap;
use std::fs::File;
use std::io::Write;
use std::io::{BufRead, BufReader, BufWriter};
use std::path::{Path, PathBuf};
use std::str::FromStr;

use super::TrackingDataArc;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum TdmParserState {
    Header,
    Metadata,
    Data,
}

#[allow(clippy::too_many_arguments)]
fn finish_segment(
    measurements: &mut Vec<Measurement>,
    segment_measurements: &mut Vec<Measurement>,
    segment_metadata: &HashMap<String, String>,
    msr_divider: f64,
    has_freq_data: bool,
    all_applied_corrections: &mut IndexSet<MeasurementType>,
    moduli: &mut Option<IndexMap<MeasurementType, f64>>,
) -> Result<(), InputOutputError> {
    if segment_measurements.is_empty() {
        return Ok(());
    }

    let mut turnaround_ratio = None;
    let drop_freq_data;
    if has_freq_data {
        // If there is any frequency measurement, compute the turn-around ratio.
        if let Some(ta_num_str) = segment_metadata.get("TURNAROUND_NUMERATOR") {
            if let Some(ta_denom_str) = segment_metadata.get("TURNAROUND_DENOMINATOR") {
                if let Ok(ta_num) = ta_num_str.parse::<i32>() {
                    if let Ok(ta_denom) = ta_denom_str.parse::<i32>() {
                        // turn-around ratio is set.
                        turnaround_ratio = Some(f64::from(ta_num) / f64::from(ta_denom));
                        info!("turn-around ratio is {ta_num}/{ta_denom}");
                        drop_freq_data = false;
                    } else {
                        error!("turn-around denominator `{ta_denom_str}` is not a valid integer");
                        drop_freq_data = true;
                    }
                } else {
                    error!("turn-around numerator `{ta_num_str}` is not a valid integer");
                    drop_freq_data = true;
                }
            } else {
                error!(
                    "required turn-around denominator missing from metadata -- dropping ALL RECEIVE/TRANSMIT data"
                );
                drop_freq_data = true;
            }
        } else {
            error!(
                "required turn-around numerator missing from metadata -- dropping ALL RECEIVE/TRANSMIT data"
            );
            drop_freq_data = true;
        }
    } else {
        drop_freq_data = true;
    }

    let corrections_applied = if let Some(corr_flag) = segment_metadata.get("CORRECTIONS_APPLIED") {
        match corr_flag.trim().to_lowercase().as_str() {
            "no" => false,
            "yes" => true,
            _ => {
                warn!("invalid CORRECTIONS_APPLIED `{corr_flag}`");
                false
            }
        }
    } else {
        false
    };

    // Now, let's convert the receive and transmit frequencies to Doppler measurements in velocity units.
    // We expect the transmit and receive frequencies to have the exact same timestamp.
    let mut freq_types = IndexSet::new();
    freq_types.insert(MeasurementType::ReceiveFrequency);
    freq_types.insert(MeasurementType::TransmitFrequency);
    freq_types.insert(MeasurementType::TransmitFrequencyRate);

    let mut latest_transmit_freq = None;
    let mut latest_transmit_epoch = None;
    let mut latest_transmit_rate = 0.0;

    for measurement in segment_measurements.iter_mut() {
        let epoch = measurement.epoch;
        // Apply corrections if any
        if !corrections_applied {
            for msr_type in [
                MeasurementType::Range,
                MeasurementType::Doppler,
                MeasurementType::Azimuth,
                MeasurementType::Elevation,
                MeasurementType::ReceiveFrequency,
                MeasurementType::TransmitFrequency,
                MeasurementType::TransmitFrequencyRate,
            ] {
                let kws = match msr_type {
                    MeasurementType::Doppler => vec![
                        "CORRECTION_DOPPLER".to_string(),
                        "CORRECTION_DOPPLER_INTEGRATED".to_string(),
                        "CORRECTION_DOPPLER_INSTANTANEOUS".to_string(),
                    ],
                    _ => vec![format!("CORRECTION_{}", msr_type.ccsds_tdm_name())],
                };

                for kw in kws {
                    if let Some(correction_str) = segment_metadata.get(&kw) {
                        if let Ok(correction) = correction_str.parse::<f64>() {
                            let scaled_correction = match msr_type {
                                MeasurementType::Range => {
                                    if let Some(range_units) = segment_metadata.get("RANGE_UNITS") {
                                        match convert_range_units(
                                            correction,
                                            range_units,
                                            msr_divider,
                                        ) {
                                            Ok(sc) => sc,
                                            Err(e) => {
                                                warn!("failed to convert CORRECTION_RANGE: {e}");
                                                continue;
                                            }
                                        }
                                    } else {
                                        warn!(
                                            "RANGE_UNITS missing when converting CORRECTION_RANGE"
                                        );
                                        correction / msr_divider
                                    }
                                }
                                MeasurementType::Doppler => correction / msr_divider,
                                _ => correction,
                            };

                            measurement.correct(msr_type, scaled_correction);
                            all_applied_corrections.insert(msr_type);
                        } else {
                            warn!("invalid correction value for {kw}");
                        }
                    }
                }
            }
        }

        if drop_freq_data {
            for freq in &freq_types {
                measurement.data.swap_remove(freq);
            }
            continue;
        }

        // Update the transmit frequency and rate if they are set.
        if let Some(rate) = measurement
            .data
            .get(&MeasurementType::TransmitFrequencyRate)
        {
            if let (Some(last_f), Some(last_e)) = (latest_transmit_freq, latest_transmit_epoch) {
                let dt: Duration = epoch - last_e;
                latest_transmit_freq = Some(last_f + latest_transmit_rate * dt.to_seconds());
            }
            latest_transmit_epoch = Some(epoch);
            latest_transmit_rate = *rate;
        }

        if let Some(freq) = measurement.data.get(&MeasurementType::TransmitFrequency) {
            latest_transmit_freq = Some(*freq);
            latest_transmit_epoch = Some(epoch);
        }

        if !measurement
            .data
            .contains_key(&MeasurementType::ReceiveFrequency)
        {
            // If there's no receive frequency, we just continue (having updated the transmit freq rate)
            // but we must remove the transmit freq rate from the measurement.
            for freq in &freq_types {
                measurement.data.swap_remove(freq);
            }
            continue;
        }

        // There is a receive frequency
        if latest_transmit_freq.is_none() {
            warn!(
                "receive frequency found at {epoch} but no transmit frequency was ever set, ignoring"
            );
            for freq in &freq_types {
                measurement.data.swap_remove(freq);
            }
            continue;
        }

        let dt: Duration = epoch - latest_transmit_epoch.unwrap();
        let transmit_freq_hz =
            latest_transmit_freq.unwrap() + latest_transmit_rate * dt.to_seconds();

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

    if let Some(range_modulus) = segment_metadata.get("RANGE_MODULUS") {
        if let Ok(value) = range_modulus.parse::<f64>() {
            if value > 0.0 {
                let map = moduli.get_or_insert_with(IndexMap::new);
                map.insert(MeasurementType::Range, value);
            }
        } else {
            warn!("could not parse RANGE_MODULUS of `{range_modulus}` as a double");
        }
    }

    // Remove measurements that have no data left after our processing.
    segment_measurements.retain(|m| !m.data.is_empty());

    measurements.append(segment_measurements);

    Ok(())
}

impl TrackingDataArc {
    /// Loads a tracking arc from its serialization in CCSDS TDM.
    ///
    /// # Support level
    ///
    /// - Only the KVN format is supported.
    /// - Support is limited to orbit determination in "xGEO", i.e. cislunar and deep space missions.
    /// - Supports multiple metadata and data sections per file.
    ///
    /// ## Data types
    ///
    /// Fully supported:
    ///     - RANGE
    ///     - DOPPLER_INSTANTANEOUS, DOPPLER_INTEGRATED
    ///     - ANGLE_1 / ANGLE_2, as azimuth/elevation only
    ///
    /// Partially supported:
    ///     - TRANSMIT_FREQ / RECEIVE_FREQ : these will be converted to Doppler measurements using the TURNAROUND_NUMERATOR and TURNAROUND_DENOMINATOR in the TDM. The freq rate is _not_ supported.
    ///
    /// ## Metadata support
    ///
    /// ### Mode
    ///
    /// Only the MODE = SEQUENTIAL is supported.
    ///
    /// ### Time systems / time scales
    ///
    /// All timescales supported by hifitime are supported here. This includes: UTC, TAI, GPS, TT, TDB, TAI, GST, QZSST, TL, TCL.
    ///
    /// ### Path
    ///
    /// Only one way or two way data is supported, i.e. path must be either `PATH n,m,n` or `PATH n,m`.
    ///
    /// Note that the actual indexes of the path are ignored.
    ///
    /// ### Participants
    ///
    /// `PARTICIPANT_1` must be the ground station / tracker.
    /// The second participant is ignored: the user must ensure that the Orbit Determination Process is properly configured and the proper arc is given.
    ///
    /// ### Turnaround ratio
    ///
    /// The turnaround ratio is only accounted for when the data contains RECEIVE_FREQ and TRANSMIT_FREQ data.
    ///
    /// ### Range and modulus
    ///
    /// Only kilometers are supported in range units. Range modulus is accounted for to compute range ambiguity.
    ///
    pub fn from_tdm<P: AsRef<Path>>(
        path: P,
        aliases: Option<HashMap<String, String>>,
    ) -> Result<Self, InputOutputError> {
        let file = File::open(&path).context(StdIOSnafu {
            action: "opening CCSDS TDM file for tracking arc",
        })?;

        let source = path.as_ref().to_path_buf().display().to_string();
        info!("parsing CCSDS TDM {source}");

        let reader = BufReader::new(file);

        let mut parser_state = TdmParserState::Header;
        let mut has_metadata_for_segment = false;

        let mut measurements = Vec::new();
        let mut segment_measurements = Vec::new();
        let mut segment_metadata = HashMap::new();
        let mut all_applied_corrections = IndexSet::new();
        let mut moduli = None;

        let mut current_tracker = String::new();
        let mut time_system = TimeScale::UTC;
        let mut has_freq_data = false;
        let mut msr_divider = 1.0;
        let mut integration_ref = None;
        let mut integration_time = None;

        let parse_one_val =
            |lno: usize, line: &str, err: &str| -> Result<String, InputOutputError> {
                match line.split_once('=') {
                    Some((_, val_str)) => Ok(val_str.trim().to_string()),
                    None => Err(InputOutputError::TDMError {
                        msg: format!("line {lno}: {err}"),
                    }),
                }
            };

        for (lno, line) in reader.lines().enumerate() {
            let line = line.context(StdIOSnafu {
                action: "reading CCSDS TDM file",
            })?;
            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            if line.starts_with("META_START") {
                match parser_state {
                    TdmParserState::Metadata => {
                        return Err(InputOutputError::TDMError {
                            msg: format!("line {lno}: nested META_START is not allowed"),
                        });
                    }
                    TdmParserState::Data => {
                        return Err(InputOutputError::TDMError {
                            msg: format!(
                                "line {lno}: META_START cannot appear inside an open DATA block"
                            ),
                        });
                    }
                    TdmParserState::Header => {}
                }
                parser_state = TdmParserState::Metadata;
                has_metadata_for_segment = false;
                current_tracker.clear();
                time_system = TimeScale::UTC;
                has_freq_data = false;
                msr_divider = 1.0;
                integration_ref = None;
                integration_time = None;
                segment_metadata.clear();
                segment_measurements.clear();
                continue;
            }

            if line.starts_with("META_STOP") {
                if parser_state != TdmParserState::Metadata {
                    return Err(InputOutputError::TDMError {
                        msg: format!("line {lno}: META_STOP without META_START"),
                    });
                }
                parser_state = TdmParserState::Header;
                has_metadata_for_segment = true;
                continue;
            }

            if line.starts_with("DATA_START") {
                match parser_state {
                    TdmParserState::Metadata => {
                        return Err(InputOutputError::TDMError {
                            msg: format!(
                                "line {lno}: DATA_START cannot appear inside metadata block"
                            ),
                        });
                    }
                    TdmParserState::Data => {
                        return Err(InputOutputError::TDMError {
                            msg: format!("line {lno}: nested DATA_START is not allowed"),
                        });
                    }
                    TdmParserState::Header => {
                        if !has_metadata_for_segment {
                            return Err(InputOutputError::TDMError {
                                msg: format!("line {lno}: DATA_START without prior metadata block"),
                            });
                        }
                    }
                }
                parser_state = TdmParserState::Data;
                continue;
            }

            if line.starts_with("DATA_STOP") {
                if parser_state != TdmParserState::Data {
                    return Err(InputOutputError::TDMError {
                        msg: format!("line {lno}: DATA_STOP without DATA_START"),
                    });
                }
                finish_segment(
                    &mut measurements,
                    &mut segment_measurements,
                    &segment_metadata,
                    msr_divider,
                    has_freq_data,
                    &mut all_applied_corrections,
                    &mut moduli,
                )?;
                parser_state = TdmParserState::Header;
                has_metadata_for_segment = false;
                continue;
            }

            if line.starts_with("COMMENT") {
                continue;
            }

            // Validate metadata keys appearing outside META_START / META_STOP
            let metadata_key = line.split_once('=').map(|(key, _)| key.trim());
            if matches!(
                metadata_key,
                Some(
                    "TIME_SYSTEM"
                        | "START_TIME"
                        | "STOP_TIME"
                        | "PARTICIPANT_1"
                        | "PARTICIPANT_2"
                        | "PARTICIPANT_3"
                        | "PARTICIPANT_4"
                        | "PARTICIPANT_5"
                        | "MODE"
                        | "PATH"
                        | "PATH_1"
                        | "PATH_2"
                        | "TRANSMIT_BAND"
                        | "RECEIVE_BAND"
                        | "INTEGRATION_INTERVAL"
                        | "INTEGRATION_TIME"
                        | "INTEGRATION_REF"
                        | "FREQ_OFFSET"
                        | "RANGE_MODE"
                        | "RANGE_MODULUS"
                        | "RANGE_UNITS"
                        | "ANGLE_TYPE"
                        | "DATA_QUALITY"
                        | "CORRECTIONS_APPLIED"
                        | "CORRECTION_RANGE"
                        | "CORRECTION_DOPPLER"
                        | "CORRECTION_DOPPLER_INTEGRATED"
                        | "CORRECTION_DOPPLER_INSTANTANEOUS"
                        | "CORRECTION_ANGLE_1"
                        | "CORRECTION_ANGLE_2"
                        | "CORRECTION_AZIMUTH"
                        | "CORRECTION_ELEVATION"
                        | "CORRECTION_RECEIVE_FREQ"
                        | "CORRECTION_TRANSMIT_FREQ"
                        | "CORRECTION_TRANSMIT_FREQ_RATE"
                        | "TURNAROUND_NUMERATOR"
                        | "TURNAROUND_DENOMINATOR"
                        | "INTERPOLATION"
                        | "INTERPOLATION_DEGREE"
                )
            ) && parser_state != TdmParserState::Metadata
            {
                return Err(InputOutputError::TDMError {
                    msg: format!(
                        "metadata field `{}` appears outside META_START/META_STOP (line {lno})",
                        metadata_key.unwrap()
                    ),
                });
            }

            if parser_state == TdmParserState::Header {
                if line.starts_with("CCSDS_TDM_VERS") {
                    let version_str = parse_one_val(lno, line, "no value for CCSDS_TDM_VERS")?;
                    match version_str.parse::<f32>() {
                        Ok(version_val) => match version_val as i16 {
                            1..=3 => {}
                            _ => {
                                return Err(InputOutputError::UnsupportedData {
                                    which: format!(
                                        "CCSDS TDM version {version_val} not supported (line {lno})"
                                    ),
                                });
                            }
                        },
                        Err(_) => {
                            return Err(InputOutputError::TDMError {
                                msg: format!(
                                    "could not parse TDM version `{version_str}` (line {lno})"
                                ),
                            });
                        }
                    }
                }
                continue;
            }

            if parser_state == TdmParserState::Metadata {
                if line.starts_with("PARTICIPANT_1") {
                    current_tracker = parse_one_val(lno, line, "no value for PARTICIPANT_1")?;
                    if let Some(aliases) = &aliases
                        && let Some(alias) = aliases.get(&current_tracker)
                    {
                        current_tracker = alias.clone();
                    }
                } else if line.starts_with("TIME_SYSTEM") {
                    let ts = parse_one_val(lno, line, "no value for TIME_SYSTEM")?;
                    if let Ok(ts_scale) = TimeScale::from_str(&ts) {
                        time_system = ts_scale;
                    } else {
                        return Err(InputOutputError::UnsupportedData {
                            which: format!("time scale `{ts}` not supported"),
                        });
                    }
                } else if line.starts_with("PATH") {
                    let path_val = parse_one_val(lno, line, "no value for PATH")?;
                    match path_val.split(',').count() {
                        2 => msr_divider = 1.0,
                        3 => msr_divider = 2.0,
                        cnt => {
                            return Err(InputOutputError::UnsupportedData {
                                which: format!(
                                    "found {cnt} paths in TDM, only 1 or 2 are supported"
                                ),
                            });
                        }
                    }
                } else if line.starts_with("INTEGRATION_REF") {
                    let value = parse_one_val(lno, line, "no value for INTEGRATION_REF")?;
                    integration_ref = Some(IntegrationRef::from_str(&value)?);
                } else if line.starts_with("INTEGRATION_TIME")
                    || line.starts_with("INTEGRATION_INTERVAL")
                {
                    let value = parse_one_val(lno, line, "no value for INTEGRATION_TIME/INTERVAL")?;
                    let dur = if let Ok(val) = value.parse::<f64>() {
                        Unit::Second * val
                    } else if let Ok(dur) = Duration::from_str(&value) {
                        dur
                    } else {
                        return Err(InputOutputError::UnsupportedData {
                            which: format!("invalid integration time `{value}`"),
                        });
                    };
                    integration_time = Some(dur);
                }

                if let Some((keyword, value)) = line.split_once('=') {
                    segment_metadata.insert(keyword.trim().to_string(), value.trim().to_string());
                }
                continue;
            }

            if parser_state == TdmParserState::Data {
                if let Some((mtype, epoch, value)) = parse_measurement_line(line, time_system)? {
                    let effective_divider = if mtype.may_be_two_way() {
                        msr_divider
                    } else {
                        if [
                            MeasurementType::ReceiveFrequency,
                            MeasurementType::TransmitFrequency,
                            MeasurementType::TransmitFrequencyRate,
                        ]
                        .contains(&mtype)
                        {
                            has_freq_data = true;
                        }
                        1.0
                    };

                    let mut scaled_value = value;
                    if mtype == MeasurementType::Range {
                        if let Some(range_units) = segment_metadata.get("RANGE_UNITS") {
                            scaled_value = convert_range_units(
                                value,
                                range_units.as_str(),
                                effective_divider,
                            )?;
                        } else {
                            return Err(InputOutputError::MissingData {
                                which:
                                    "RANGE_UNITS not specified in metadata for RANGE measurement"
                                        .to_string(),
                            });
                        }
                    } else {
                        scaled_value /= effective_divider;
                    }

                    let is_concurrent =
                        segment_measurements
                            .last()
                            .is_some_and(|last: &Measurement| {
                                last.epoch == epoch && last.tracker == current_tracker
                            });

                    let doppler_config = match (integration_time, integration_ref) {
                        (Some(time), Some(reference)) => Some(DopplerConfig {
                            integration_time: time,
                            integration_ref: reference,
                        }),
                        (Some(time), None) => Some(DopplerConfig {
                            integration_time: time,
                            integration_ref: IntegrationRef::default(),
                        }),
                        (None, Some(reference)) => Some(DopplerConfig {
                            integration_time: DopplerConfig::default().integration_time,
                            integration_ref: reference,
                        }),
                        (None, None) => None,
                    };

                    if is_concurrent {
                        let last = segment_measurements.last_mut().unwrap();
                        last.data.insert(mtype, scaled_value);
                        if last.doppler_config.is_none() {
                            last.doppler_config = doppler_config;
                        }
                    } else {
                        let mut data = IndexMap::new();
                        data.insert(mtype, scaled_value);

                        segment_measurements.push(Measurement {
                            tracker: current_tracker.clone(),
                            epoch,
                            data,
                            rejected: false,
                            doppler_config,
                        });
                    }
                }
            }
        }

        if parser_state == TdmParserState::Metadata {
            return Err(InputOutputError::TDMError {
                msg: "unterminated META_START section at end of file".to_string(),
            });
        }
        if parser_state == TdmParserState::Data {
            return Err(InputOutputError::TDMError {
                msg: "unterminated DATA_START section at end of file".to_string(),
            });
        }

        if !all_applied_corrections.is_empty() {
            info!("applied corrections for {all_applied_corrections:?}");
        }

        let mut trk = Self {
            measurements,
            source: Some(source),
            moduli,
            force_reject: false,
        };

        // Ensure data is sorted (TDM spec requires that, but you never know).
        trk.sort();

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
        if let Some(start_epoch) = cfg.start_epoch {
            if let Some(end_epoch) = cfg.end_epoch {
                self = self.filter_by_epoch(start_epoch..end_epoch);
            } else {
                self = self.filter_by_epoch(start_epoch..);
            }
        } else if let Some(end_epoch) = cfg.end_epoch {
            self = self.filter_by_epoch(..end_epoch);
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
                if types.is_empty() {
                    continue;
                }
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
                    let k_upper = k.to_uppercase();
                    if k != "originator"
                        && (!types.contains(&MeasurementType::Doppler)
                            || (k_upper != "INTEGRATION_INTERVAL" && k_upper != "INTEGRATION_REF"))
                    {
                        writeln!(writer, "\t{k} = {v}").map_err(err_hdlr)?;
                    }
                }

                if types.contains(&MeasurementType::Doppler)
                    && let Some(doppler_cfg) = tracker_data
                        .measurements
                        .iter()
                        .find_map(|m| m.doppler_config)
                {
                    writeln!(
                        writer,
                        "\tINTEGRATION_INTERVAL = {:.6}",
                        doppler_cfg.integration_time.to_seconds()
                    )
                    .map_err(err_hdlr)?;
                    let ref_str = match doppler_cfg.integration_ref {
                        IntegrationRef::Start => "START",
                        IntegrationRef::Middle => "MIDDLE",
                        IntegrationRef::End => "END",
                    };
                    writeln!(writer, "\tINTEGRATION_REF = {ref_str}").map_err(err_hdlr)?;
                }

                if types.contains(&MeasurementType::Range) {
                    writeln!(writer, "\tRANGE_UNITS = km").map_err(err_hdlr)?;

                    if let Some(moduli) = &self.moduli
                        && let Some(range_modulus) = moduli.get(&MeasurementType::Range)
                    {
                        writeln!(writer, "\tRANGE_MODULUS = {range_modulus:E}")
                            .map_err(err_hdlr)?;
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
                for m in &tracker_data.measurements {
                    for (mtype, value) in &m.data {
                        if !types.contains(mtype) {
                            continue;
                        }

                        writeln!(
                            writer,
                            "\t{:<20} = {:<23}\t{:.12}",
                            mtype.ccsds_tdm_name(),
                            Formatter::new(m.epoch, iso8601_no_ts),
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

fn convert_range_units(
    value: f64,
    range_units: &str,
    divider: f64,
) -> Result<f64, InputOutputError> {
    match range_units {
        "km" => Ok(value / divider),
        "RU" => Err(InputOutputError::UnsupportedData {
            which: "RANGE_UNITS `RU` requires mission-specific conversion and is not currently supported".to_string(),
        }),
        "s" => Ok((value * SPEED_OF_LIGHT_KM_S) / divider),
        "m" => {
            warn!(
                "RANGE_UNITS in TDM file is `m`, which is not CCSDS compliant. Proceeding with conversion to km."
            );
            Ok((value / 1000.0) / divider)
        }
        "ms" => {
            warn!(
                "RANGE_UNITS in TDM file is `ms`, which is not CCSDS compliant. Proceeding with conversion to km."
            );
            Ok((value * 1e-3 * SPEED_OF_LIGHT_KM_S) / divider)
        }
        "us" => {
            warn!(
                "RANGE_UNITS in TDM file is `us`, which is not CCSDS compliant. Proceeding with conversion to km."
            );
            Ok((value * 1e-6 * SPEED_OF_LIGHT_KM_S) / divider)
        }
        "ns" | "NANOSEC" => {
            warn!(
                "RANGE_UNITS in TDM file is `ns`, which is not CCSDS compliant. Proceeding with conversion to km."
            );
            Ok((value * 1e-9 * SPEED_OF_LIGHT_KM_S) / divider)
        }
        _ => Err(InputOutputError::UnsupportedData {
            which: format!("unsupported RANGE_UNITS `{range_units}`"),
        }),
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
        "TRANSMIT_FREQ_RATE"
        | "TRANSMIT_FREQ_RATE_1"
        | "TRANSMIT_FREQ_RATE_2"
        | "TRANSMIT_FREQ_RATE_3"
        | "TRANSMIT_FREQ_RATE_4"
        | "TRANSMIT_FREQ_RATE_5" => MeasurementType::TransmitFrequencyRate,
        _ => {
            return Err(InputOutputError::UnsupportedData {
                which: mtype_str.to_string(),
            });
        }
    };

    let data_parts: Vec<&str> = data.split_whitespace().collect();
    if data_parts.len() != 2 {
        return Ok(None);
    }

    let epoch =
        Epoch::from_gregorian_str(&format!("{} {time_system}", data_parts[0])).map_err(|e| {
            InputOutputError::TDMError {
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
