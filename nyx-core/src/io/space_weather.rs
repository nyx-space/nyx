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

use crate::io::InputOutputError;
use flate2::read::GzDecoder;
use hifitime::prelude::*;
use serde::Deserialize;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::str::FromStr;

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use std::path::PathBuf;

/// Comprehensive representation of a single daily record in the CelesTrak Space Weather CSV.
#[derive(Debug, Clone, Deserialize)]
#[cfg_attr(feature = "python", pyclass(from_py_object, get_all))]
pub struct RawSpaceWeatherRow {
    #[serde(rename = "DATE")]
    pub date: String,
    #[serde(rename = "BSRN")]
    pub bsrn: u32,
    #[serde(rename = "ND")]
    pub nd: u32,

    // Kp Planetary Indices (Encoded as Kp * 10 in CelesTrak CSV)
    #[serde(rename = "KP1")]
    pub kp1: f64,
    #[serde(rename = "KP2")]
    pub kp2: f64,
    #[serde(rename = "KP3")]
    pub kp3: f64,
    #[serde(rename = "KP4")]
    pub kp4: f64,
    #[serde(rename = "KP5")]
    pub kp5: f64,
    #[serde(rename = "KP6")]
    pub kp6: f64,
    #[serde(rename = "KP7")]
    pub kp7: f64,
    #[serde(rename = "KP8")]
    pub kp8: f64,
    #[serde(rename = "KP_SUM")]
    pub kp_sum: f64,

    // Ap Linear Indices
    #[serde(rename = "AP1")]
    pub ap1: f64,
    #[serde(rename = "AP2")]
    pub ap2: f64,
    #[serde(rename = "AP3")]
    pub ap3: f64,
    #[serde(rename = "AP4")]
    pub ap4: f64,
    #[serde(rename = "AP5")]
    pub ap5: f64,
    #[serde(rename = "AP6")]
    pub ap6: f64,
    #[serde(rename = "AP7")]
    pub ap7: f64,
    #[serde(rename = "AP8")]
    pub ap8: f64,
    #[serde(rename = "AP_AVG")]
    pub ap_avg: f64,

    // Geophysical & Solar Indicators
    #[serde(rename = "CP")]
    pub cp: f64,
    #[serde(rename = "C9")]
    pub c9: u8,
    #[serde(rename = "ISN")]
    pub isn: u32,

    // Solar Radio Flux (10.7 cm) - Observed & Adjusted
    #[serde(rename = "F10.7_OBS")]
    pub f107_obs: f64,
    #[serde(rename = "F10.7_ADJ")]
    pub f107_adj: f64,
    #[serde(rename = "F10.7_DATA_TYPE")]
    pub f107_data_type: String,
    #[serde(rename = "F10.7_OBS_CENTER81")]
    pub f107_obs_center81: f64,
    #[serde(rename = "F10.7_OBS_LAST81")]
    pub f107_obs_last81: f64,
    #[serde(rename = "F10.7_ADJ_CENTER81")]
    pub f107_adj_center81: f64,
    #[serde(rename = "F10.7_ADJ_LAST81")]
    pub f107_adj_last81: f64,
}

impl RawSpaceWeatherRow {
    /// Returns the eight 3-hour Kp values rescaled to standard floating-point bounds [0.0, 9.0].
    #[inline]
    pub fn kp_bins(&self) -> [f64; 8] {
        [
            self.kp1 / 10.0,
            self.kp2 / 10.0,
            self.kp3 / 10.0,
            self.kp4 / 10.0,
            self.kp5 / 10.0,
            self.kp6 / 10.0,
            self.kp7 / 10.0,
            self.kp8 / 10.0,
        ]
    }

    /// Returns the eight 3-hour linear Ap values.
    #[inline]
    pub fn ap_bins(&self) -> [f64; 8] {
        [
            self.ap1, self.ap2, self.ap3, self.ap4, self.ap5, self.ap6, self.ap7, self.ap8,
        ]
    }
}

#[cfg_attr(feature = "python", pyclass(from_py_object, get_all))]
#[derive(Clone, Debug)]
pub struct SpaceWeatherData {
    records: BTreeMap<Epoch, RawSpaceWeatherRow>,
}

impl SpaceWeatherData {
    /// Ingests a complete CelesTrak Space Weather CSV file into an Epoch-indexed map.
    pub fn from_csv_file<P: AsRef<Path>>(path: P) -> Result<Self, InputOutputError> {
        let path_ref = path.as_ref();
        let file = File::open(path_ref).map_err(|e| InputOutputError::StdIOError {
            source: e,
            action: "reading space weather file",
        })?;

        let mut buf_reader = BufReader::new(file);

        // Peek at the first 2 bytes to detect Gzip magic header (0x1F, 0x8B) without consuming the buffer
        let is_gzipped = match buf_reader.fill_buf() {
            Ok(header) => header.len() >= 2 && header[0] == 0x1f && header[1] == 0x8b,
            Err(source) => {
                return Err(InputOutputError::StdIOError {
                    source,
                    action: "reading header of CSV file",
                })
            }
        };

        let stream: Box<dyn Read> = if is_gzipped {
            Box::new(GzDecoder::new(buf_reader))
        } else {
            Box::new(buf_reader)
        };

        let mut rdr = csv::ReaderBuilder::new()
            .trim(csv::Trim::All)
            .from_reader(stream);

        let mut records = BTreeMap::new();

        for result in rdr.deserialize() {
            let record: RawSpaceWeatherRow =
                result.map_err(|source| InputOutputError::CsvData {
                    source,
                    action: "reading space weather",
                })?;
            if let Ok(epoch) = Epoch::from_str(&format!("{}T00:00:00 UTC", record.date)) {
                records.insert(epoch, record);
            }
        }

        Ok(Self { records })
    }

    /// Returns a reference to the raw, unparsed daily row for an exact midnight UTC epoch.
    pub fn raw_daily_record(&self, midnight_epoch: Epoch) -> Option<&RawSpaceWeatherRow> {
        self.records.get(&midnight_epoch)
    }

    /// Evaluates the space weather state at `epoch` and constructs the `Msise00DailyWeather` payload.
    pub fn msise_weather(&self, epoch: Epoch) -> Option<Msise00DailyWeather> {
        let target_midnight = epoch.with_hms(0, 0, 0);
        let current_day = self.records.get(&target_midnight)?;

        let seconds_into_day = (epoch - target_midnight).to_seconds();
        // Bins are 3 hours large
        let bin_idx = ((seconds_into_day / (Unit::Hour * 3).to_seconds()).floor() as usize).min(7);

        let ap_history = self.build_ap_history(target_midnight, bin_idx)?;

        Some(Msise00DailyWeather {
            f107_daily_sfu: current_day.f107_obs,
            f107_avg_sfu: current_day.f107_obs_center81,
            ap_daily: current_day.ap_avg,
            ap_3hour_history: ap_history,
        })
    }

    /// Assembles the 7-element Ap array spanning current bin back 57 hours across 4 calendar days.
    fn build_ap_history(&self, midnight: Epoch, bin_idx: usize) -> Option<[f64; 7]> {
        let one_day = Unit::Day * 1.0;

        let day_0 = self.records.get(&midnight)?;
        let day_m1 = self.records.get(&(midnight - one_day))?;
        let day_m2 = self.records.get(&(midnight - one_day * 2.0))?;
        let day_m3 = self.records.get(&(midnight - one_day * 3.0))?;

        let mut continuous_ap = [0.0; 32];
        continuous_ap[0..8].copy_from_slice(&day_m3.ap_bins());
        continuous_ap[8..16].copy_from_slice(&day_m2.ap_bins());
        continuous_ap[16..24].copy_from_slice(&day_m1.ap_bins());
        continuous_ap[24..32].copy_from_slice(&day_0.ap_bins());

        let idx = 24 + bin_idx;

        let avg_slice = |start: usize, end: usize| -> f64 {
            let slice = &continuous_ap[start..=end];
            slice.iter().sum::<f64>() / slice.len() as f64
        };

        Some([
            day_0.ap_avg,                  // ap_3hour_history[0]: Daily Ap
            continuous_ap[idx],            // ap_3hour_history[1]: Ap at target epoch
            continuous_ap[idx - 1],        // ap_3hour_history[2]: Ap at T - 3h
            continuous_ap[idx - 2],        // ap_3hour_history[3]: Ap at T - 6h
            continuous_ap[idx - 3],        // ap_3hour_history[4]: Ap at T - 9h
            avg_slice(idx - 11, idx - 4),  // ap_3hour_history[5]: Average Ap from T-12h to T-33h
            avg_slice(idx - 19, idx - 12), // ap_3hour_history[6]: Average Ap from T-36h to T-57h
        ])
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl SpaceWeatherData {
    #[new]
    fn py_new(path: PathBuf) -> Result<Self, InputOutputError> {
        Self::from_csv_file(path)
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}

// Define the model-specific weather data extracted from the space weather

/// Target weather payload required by the NRLMSISE-00 density model.
#[derive(Debug, Clone, Copy, Default)]
pub struct Msise00DailyWeather {
    /// Daily F10.7 solar radio flux: $\text{ SFU} = 10^{-22} \text{ W}\cdot\text{m}^{-2}\cdot\text{Hz}^{-1}$.
    pub f107_daily_sfu: f64,
    /// 81-day centered average F10.7 solar radio flux [SFU].
    pub f107_avg_sfu: f64,
    /// Daily mean planetary Ap index.
    pub ap_daily: f64,
    /// 7-element Ap historical array covering the 57-hour lookback window.
    pub ap_3hour_history: [f64; 7],
}
