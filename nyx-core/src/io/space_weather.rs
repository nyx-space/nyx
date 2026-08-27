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
use serde::{Deserialize, Serialize};
use serde_dhall::{SimpleType, StaticType};
use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;
use std::str::FromStr;

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use std::path::PathBuf;

/// Strategy for resolving missing predictive space weather parameters (F10.7, Ap, Kp).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "python", pyclass(from_py_object, get_all))]
pub enum StaticSpaceWeather {
    /// Solar minimum conditions (F10.7 = 65.0 SFU, Ap = 4.0, Kp = 1.0).
    SolarMinimum(),
    /// 11-year solar cycle average (F10.7 = 130.0 SFU, Ap = 15.0, Kp = 3.0). Standard default.
    SolarAverage(),
    /// Sustained solar maximum conditions (F10.7 = 200.0 SFU, Ap = 30.0, Kp = 4.3).
    SolarMaximum(),
    /// Custom operator-defined fallback values for missing data.
    Custom { f107: f64, ap: f64, kp: f64 },
}

impl Default for StaticSpaceWeather {
    fn default() -> Self {
        Self::SolarAverage()
    }
}

impl StaticSpaceWeather {
    /// Resolves the effective F10.7 solar flux [SFU].
    pub fn resolve_f107(&self, value: Option<f64>) -> f64 {
        value.unwrap_or(match self {
            Self::SolarMinimum() => 65.0,
            Self::SolarAverage() => 130.0,
            Self::SolarMaximum() => 200.0,
            Self::Custom { f107, .. } => *f107,
        })
    }

    /// Resolves the effective Ap index.
    pub fn resolve_ap(&self, value: Option<f64>) -> f64 {
        value.unwrap_or(match self {
            Self::SolarMinimum() => 4.0,
            Self::SolarAverage() => 15.0,
            Self::SolarMaximum() => 30.0,
            Self::Custom { ap, .. } => *ap,
        })
    }

    /// Resolves the effective Kp index (unscaled float, e.g., 3.0).
    pub fn resolve_kp(&self, value: Option<f64>) -> f64 {
        value.unwrap_or(match self {
            Self::SolarMinimum() => 1.0,
            Self::SolarAverage() => 3.0,
            Self::SolarMaximum() => 4.3,
            Self::Custom { kp, .. } => *kp,
        })
    }
}

/// Comprehensive representation of a single daily record in the CelesTrak Space Weather CSV.
#[derive(Debug, Clone, Deserialize, Serialize, StaticType, PartialEq)]
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
    pub kp1: Option<f64>,
    #[serde(rename = "KP2")]
    pub kp2: Option<f64>,
    #[serde(rename = "KP3")]
    pub kp3: Option<f64>,
    #[serde(rename = "KP4")]
    pub kp4: Option<f64>,
    #[serde(rename = "KP5")]
    pub kp5: Option<f64>,
    #[serde(rename = "KP6")]
    pub kp6: Option<f64>,
    #[serde(rename = "KP7")]
    pub kp7: Option<f64>,
    #[serde(rename = "KP8")]
    pub kp8: Option<f64>,
    #[serde(rename = "KP_SUM")]
    pub kp_sum: Option<f64>,

    // Ap Linear Indices
    #[serde(rename = "AP1")]
    pub ap1: Option<f64>,
    #[serde(rename = "AP2")]
    pub ap2: Option<f64>,
    #[serde(rename = "AP3")]
    pub ap3: Option<f64>,
    #[serde(rename = "AP4")]
    pub ap4: Option<f64>,
    #[serde(rename = "AP5")]
    pub ap5: Option<f64>,
    #[serde(rename = "AP6")]
    pub ap6: Option<f64>,
    #[serde(rename = "AP7")]
    pub ap7: Option<f64>,
    #[serde(rename = "AP8")]
    pub ap8: Option<f64>,
    #[serde(rename = "AP_AVG")]
    pub ap_avg: Option<f64>,

    // Geophysical & Solar Indicators
    #[serde(rename = "CP")]
    pub cp: Option<f64>,
    #[serde(rename = "C9")]
    pub c9: Option<u16>,
    #[serde(rename = "ISN")]
    pub isn: Option<u32>,

    // Solar Radio Flux (10.7 cm) - Observed & Adjusted
    #[serde(rename = "F10.7_OBS")]
    pub f107_obs: f64,
    #[serde(rename = "F10.7_ADJ")]
    pub f107_adj: f64,
    #[serde(rename = "F10.7_DATA_TYPE")]
    pub f107_data_type: String,
    #[serde(rename = "F10.7_OBS_CENTER81")]
    pub f107_obs_center81: Option<f64>,
    #[serde(rename = "F10.7_OBS_LAST81")]
    pub f107_obs_last81: Option<f64>,
    #[serde(rename = "F10.7_ADJ_CENTER81")]
    pub f107_adj_center81: Option<f64>,
    #[serde(rename = "F10.7_ADJ_LAST81")]
    pub f107_adj_last81: Option<f64>,
}

impl RawSpaceWeatherRow {
    /// Returns the eight 3-hour Kp values rescaled to standard floating-point bounds [0.0, 9.0].
    ///
    /// Missing bins default first to the row's mean daily Kp (derived from `KP_SUM`),
    /// and secondarily to the provided `SpaceWeatherFallback` policy.
    #[inline]
    pub fn kp_bins(&self, fallback: StaticSpaceWeather) -> [f64; 8] {
        // KP_SUM in CelesTrak CSV is the sum of the eight 3-hour $K_p \times 10$ values.
        // Dividing KP_SUM by $80.0$ yields the mean 3-hour $K_p$ index in standard $0.0\text{--}9.0$ scale
        let daily_mean_kp = self
            .kp_sum
            .map(|sum| sum / 80.0)
            .unwrap_or_else(|| fallback.resolve_kp(None));

        let resolve = |bin: Option<f64>| bin.map(|v| v / 10.0).unwrap_or(daily_mean_kp);

        [
            resolve(self.kp1),
            resolve(self.kp2),
            resolve(self.kp3),
            resolve(self.kp4),
            resolve(self.kp5),
            resolve(self.kp6),
            resolve(self.kp7),
            resolve(self.kp8),
        ]
    }

    /// Returns the eight 3-hour linear Ap values.
    ///
    /// Missing bins default first to the row's `AP_AVG`, and secondarily to the
    /// provided `SpaceWeatherFallback` policy.
    #[inline]
    pub fn ap_bins(&self, fallback: StaticSpaceWeather) -> [f64; 8] {
        let daily_mean_ap = self.ap_avg.unwrap_or_else(|| fallback.resolve_ap(None));

        let resolve = |bin: Option<f64>| bin.unwrap_or(daily_mean_ap);

        [
            resolve(self.ap1),
            resolve(self.ap2),
            resolve(self.ap3),
            resolve(self.ap4),
            resolve(self.ap5),
            resolve(self.ap6),
            resolve(self.ap7),
            resolve(self.ap8),
        ]
    }
}

/// Stores SpaceWeather data as provided by [CelesTrak](https://celestrak.org/SpaceData/).
/// Data may be provided either as original CSV or in a compressed (non-archived) gunzip (gz) format.
///
/// :type path: str | None
/// :type fallback: StaticSpaceWeather | None
#[cfg_attr(feature = "python", pyclass(from_py_object))]
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct SpaceWeatherData {
    #[serde(with = "as_vec")]
    pub records: BTreeMap<Epoch, RawSpaceWeatherRow>,
    pub fallback: StaticSpaceWeather,
}

impl SpaceWeatherData {
    /// Initialize new space weather by provided fixed values only.
    pub fn from_static_weather(weather: StaticSpaceWeather) -> Self {
        Self {
            records: BTreeMap::new(),
            fallback: weather,
        }
    }

    /// Ingests a complete CelesTrak Space Weather CSV file into an Epoch-indexed map.
    pub fn from_csv_file<P: AsRef<Path>>(
        path: P,
        fallback: StaticSpaceWeather,
    ) -> Result<Self, InputOutputError> {
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
                });
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

        Ok(Self { records, fallback })
    }

    /// Returns a reference to the raw, unparsed daily row for an exact midnight UTC epoch.
    pub fn raw_daily_record(&self, midnight_epoch: Epoch) -> Option<&RawSpaceWeatherRow> {
        self.records.get(&midnight_epoch)
    }
}

#[cfg_attr(feature = "python", pymethods)]
impl SpaceWeatherData {
    /// Evaluates the space weather state at `epoch` and constructs the `Msise00DailyWeather` payload.
    ///
    /// Missing daily records or unforecasted fields are resolved using `SpaceWeatherFallback`.
    ///
    /// :type epoch: Epoch
    /// :rtype: Msise00DailyWeather
    pub fn msise_weather(&self, epoch: Epoch) -> Msise00DailyWeather {
        let target_midnight = epoch.with_hms_strict(0, 0, 0) - Unit::Day * 1;
        let current_day = self.records.get(&target_midnight);

        let seconds_into_day = (epoch - target_midnight).to_seconds();
        // Bins are 3 hours large (0..7)
        let bin_idx = ((seconds_into_day / (Unit::Hour * 3).to_seconds()).floor() as usize).min(7);

        let ap_history = self.build_ap_history(target_midnight, bin_idx);

        // 1. Daily F10.7: Prefer observed, fall back to adjusted, then global fallback
        let f107_daily = self.fallback.resolve_f107(current_day.map(|r| r.f107_obs));

        // 2. 81-day Centered Mean F10.7: Prefer observed 81d, then adjusted 81d,
        // fall back to resolved daily F10.7 before applying static global fallback
        let f107_avg = current_day
            .and_then(|r| r.f107_obs_center81.or(r.f107_adj_center81))
            .unwrap_or(f107_daily);

        // 3. Daily Ap: Prefer recorded ap_avg, fall back to fallback policy
        let ap_daily = self.fallback.resolve_ap(current_day.and_then(|r| r.ap_avg));

        Msise00DailyWeather {
            f107_daily_sfu: f107_daily,
            f107_avg_sfu: f107_avg,
            ap_daily,
            ap_3hour_history: ap_history,
        }
    }

    /// Assembles the 7-element Ap array spanning current bin back 57 hours across 4 calendar days.
    ///
    /// Missing daily records or unforecasted bins are populated using the configured `SpaceWeatherFallback`.
    ///
    /// :type midnight: Epoch
    /// :type bin_idx: int
    /// :rtype: list[float]
    fn build_ap_history(&self, midnight: Epoch, bin_idx: usize) -> [f64; 7] {
        let one_day = Unit::Day * 1.0;

        // Helper to retrieve or synthesize a 8-bin 3-hour Ap slice for a given day offset.
        let get_ap_bins = |offset_days: f64| -> [f64; 8] {
            let target_epoch = midnight - one_day * offset_days;
            match self.records.get(&target_epoch) {
                Some(row) => row.ap_bins(self.fallback),
                None => [self.fallback.resolve_ap(None); 8],
            }
        };

        // Extract day 0 metadata and bins
        let day_0_row = self.records.get(&midnight);
        let daily_ap = self.fallback.resolve_ap(day_0_row.and_then(|r| r.ap_avg));
        let day_0_bins = match day_0_row {
            Some(row) => row.ap_bins(self.fallback),
            None => [self.fallback.resolve_ap(None); 8],
        };

        let mut continuous_ap = [0.0; 32];
        continuous_ap[0..8].copy_from_slice(&get_ap_bins(3.0));
        continuous_ap[8..16].copy_from_slice(&get_ap_bins(2.0));
        continuous_ap[16..24].copy_from_slice(&get_ap_bins(1.0));
        continuous_ap[24..32].copy_from_slice(&day_0_bins);

        let idx = 24 + bin_idx;

        let avg_slice = |start: usize, end: usize| -> f64 {
            let slice = &continuous_ap[start..=end];
            slice.iter().sum::<f64>() / slice.len() as f64
        };

        [
            daily_ap,                      // ap_3hour_history[0]: Daily Ap
            continuous_ap[idx],            // ap_3hour_history[1]: Ap at target epoch
            continuous_ap[idx - 1],        // ap_3hour_history[2]: Ap at T - 3h
            continuous_ap[idx - 2],        // ap_3hour_history[3]: Ap at T - 6h
            continuous_ap[idx - 3],        // ap_3hour_history[4]: Ap at T - 9h
            avg_slice(idx - 11, idx - 4),  // ap_3hour_history[5]: Average Ap from T-12h to T-33h
            avg_slice(idx - 19, idx - 12), // ap_3hour_history[6]: Average Ap from T-36h to T-57h
        ]
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl SpaceWeatherData {
    #[new]
    fn py_new(
        path: Option<PathBuf>,
        fallback: Option<StaticSpaceWeather>,
    ) -> Result<Self, InputOutputError> {
        if let Some(path) = path {
            Self::from_csv_file(path, fallback.unwrap_or_default())
        } else if let Some(weather) = fallback {
            Ok(Self::from_static_weather(weather))
        } else {
            Err(InputOutputError::MissingData {
                which:
                    "must provide at least either a path to a weather file or a fallback, or both"
                        .to_string(),
            })
        }
    }

    fn __str__(&self) -> String {
        format!("{self}")
    }

    fn __repr__(&self) -> String {
        format!("{self} @ {self:p}")
    }
}

impl fmt::Display for SpaceWeatherData {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.records.is_empty() {
            write!(f, "empty SpaceWeatherData")
        } else {
            write!(
                f,
                "SpaceWeatherData from {} to {} ({:?})",
                self.records.first_key_value().unwrap().0,
                self.records.last_key_value().unwrap().0,
                self.fallback
            )
        }
    }
}

// Implement StaticType manually by mapping BTreeMap to a List of key-value pairs
impl StaticType for SpaceWeatherData {
    fn static_type() -> SimpleType {
        let mut rcrd = HashMap::new();
        rcrd.insert("epoch".to_string(), String::static_type());
        rcrd.insert("raw_weather".to_string(), RawSpaceWeatherRow::static_type());

        SimpleType::List(Box::new(SimpleType::Record(rcrd)))
    }
}

/// Serde helper module to serialize BTreeMap as a vector of pairs
/// so it matches Dhall's List of records representation.
mod as_vec {
    use super::*;
    use serde::{Deserializer, Serializer};

    #[derive(Serialize, Deserialize)]
    struct WeatherEntry {
        epoch: Epoch,
        raw_weather: RawSpaceWeatherRow,
    }

    pub fn serialize<S>(
        map: &BTreeMap<Epoch, RawSpaceWeatherRow>,
        serializer: S,
    ) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let vec: Vec<WeatherEntry> = map
            .iter()
            .map(|(epoch, raw_weather)| WeatherEntry {
                epoch: *epoch,
                raw_weather: raw_weather.clone(),
            })
            .collect();
        vec.serialize(serializer)
    }

    pub fn deserialize<'de, D>(
        deserializer: D,
    ) -> Result<BTreeMap<Epoch, RawSpaceWeatherRow>, D::Error>
    where
        D: Deserializer<'de>,
    {
        use serde::Deserialize;
        let vec: Vec<WeatherEntry> = Vec::deserialize(deserializer)?;
        let mut rcrd = BTreeMap::new();
        for entry in vec {
            rcrd.insert(entry.epoch, entry.raw_weather);
        }
        Ok(rcrd)
    }
}

// Define the model-specific weather data extracted from the space weather

/// Target weather payload required by the NRLMSISE-00 density model.
#[derive(Debug, Clone, Copy, Default)]
#[cfg_attr(feature = "python", pyclass(from_py_object))]
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

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl Msise00DailyWeather {
    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}
