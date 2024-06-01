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

use crate::{
    io::{MissingDataSnafu, ParquetSnafu},
    linalg::{allocator::Allocator, DefaultAllocator, OVector},
    od::{msr::TrackingArc, Measurement},
};

#[cfg(feature = "python")]
use crate::NyxError;

use arrow::{
    array::{Float64Array, StringArray},
    record_batch::RecordBatchReader,
};
use hifitime::Epoch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use snafu::prelude::*;
use std::fs::File;
use std::{collections::HashMap, error::Error, fmt::Display, path::Path};

#[cfg(feature = "python")]
use pyo3::prelude::*;

use super::{InputOutputError, StdIOSnafu};

/// A dynamic tracking arc allows loading a set of measurements from a parquet file and converting them
/// to the concrete measurement type when desired.
#[cfg_attr(feature = "python", pyclass)]
pub struct DynamicTrackingArc {
    pub device_cfg: String,
    pub path: String,
    metadata: HashMap<String, String>,
}

impl DynamicTrackingArc {
    pub fn from_parquet<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn Error>> {
        let file = File::open(&path)?;

        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();

        let mut metadata = HashMap::new();
        let mut device_cfg = String::new();
        // Build the custom metadata
        if let Some(file_metadata) = builder.metadata().file_metadata().key_value_metadata() {
            for key_value in file_metadata {
                if !key_value.key.starts_with("ARROW:") {
                    metadata.insert(
                        key_value.key.clone(),
                        key_value.value.clone().unwrap_or("[unset]".to_string()),
                    );
                    if key_value.key == "devices" {
                        device_cfg = key_value.value.clone().unwrap_or("[unset]".to_string());
                    }
                }
            }
        }

        let me = Self {
            path: path.as_ref().to_string_lossy().to_string(),
            metadata,
            device_cfg,
        };

        for item in me.repr() {
            info!("{item}");
        }

        Ok(me)
    }

    /// Reads through the loaded parquet file and attempts to convert to the provided tracking arc.
    pub fn to_tracking_arc<Msr>(&self) -> Result<TrackingArc<Msr>, InputOutputError>
    where
        Msr: Measurement,
        DefaultAllocator: Allocator<f64, Msr::MeasurementSize>,
    {
        // Read the file since we closed it earlier
        let file = File::open(&self.path).context(StdIOSnafu {
            action: "opening file for tracking arc",
        })?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();

        let reader = builder.build().context(ParquetSnafu {
            action: "reading tracking arc",
        })?;

        // Check the schema
        let mut has_epoch = false;
        let mut has_tracking_dev = false;
        let mut range_avail = false;
        let mut rate_avail = false;
        for field in &reader.schema().fields {
            match field.name().as_str() {
                "Epoch:TAI (s)" => has_epoch = true,
                "Tracking device" => has_tracking_dev = true,
                "Range (km)" => range_avail = true,
                "Doppler (km/s)" => rate_avail = true,
                _ => {}
            }
        }

        ensure!(
            has_epoch,
            MissingDataSnafu {
                which: "Epoch: TAI (s)"
            }
        );

        ensure!(
            has_tracking_dev,
            MissingDataSnafu {
                which: "Tracking device"
            }
        );

        ensure!(
            range_avail || rate_avail,
            MissingDataSnafu {
                which: "`Range (km)` or `Doppler (km/s)`"
            }
        );

        let expected_type = std::any::type_name::<Msr>().split("::").last().unwrap();

        // Only check that the file contains the data we need
        match expected_type {
            "RangeDoppler" => {
                if !range_avail || !rate_avail {
                    return Err(InputOutputError::MissingData {
                        which: "`Range (km)` and `Doppler (km/s)`".to_string(),
                    });
                }
            }
            "RangeMsr" => {
                if !range_avail {
                    return Err(InputOutputError::MissingData {
                        which: "`Range (km)`".to_string(),
                    });
                }
            }
            "RangeRate" => {
                return Err(InputOutputError::MissingData {
                    which: "`Doppler (km/s)`".to_string(),
                });
            }
            _ => {
                return Err(InputOutputError::UnsupportedData {
                    which: expected_type.to_string(),
                });
            }
        }

        // At this stage, we know that the measurement is valid and the conversion is supported.
        let mut arc = TrackingArc {
            device_cfg: self.device_cfg.clone(),
            measurements: Vec::new(),
        };

        // Now convert each batch on the fly
        for maybe_batch in reader {
            let batch = maybe_batch.unwrap();

            let tracking_device = batch
                .column_by_name("Tracking device")
                .unwrap()
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();

            let epochs = batch
                .column_by_name("Epoch:TAI (s)")
                .unwrap()
                .as_any()
                .downcast_ref::<Float64Array>()
                .unwrap();

            // Now read the data depending on what we're deserializing as
            match expected_type {
                "RangeDoppler" => {
                    let range_data = batch
                        .column_by_name("Range (km)")
                        .unwrap()
                        .as_any()
                        .downcast_ref::<Float64Array>()
                        .unwrap();

                    let rate_data = batch
                        .column_by_name("Doppler (km/s)")
                        .unwrap()
                        .as_any()
                        .downcast_ref::<Float64Array>()
                        .unwrap();

                    // Set the measurements in the tracking arc
                    for i in 0..batch.num_rows() {
                        arc.measurements.push((
                            tracking_device.value(i).to_string(),
                            Msr::from_observation(
                                Epoch::from_tai_seconds(epochs.value(i)),
                                OVector::<f64, Msr::MeasurementSize>::from_iterator([
                                    range_data.value(i),
                                    rate_data.value(i),
                                ]),
                            ),
                        ));
                    }
                }
                "RangeMsr" => {
                    let range_data = batch
                        .column_by_name("Range (km)")
                        .unwrap()
                        .as_any()
                        .downcast_ref::<Float64Array>()
                        .unwrap();

                    // Set the measurements in the tracking arc
                    for i in 0..batch.num_rows() {
                        arc.measurements.push((
                            tracking_device.value(i).to_string(),
                            Msr::from_observation(
                                Epoch::from_tdb_seconds(epochs.value(i)),
                                OVector::<f64, Msr::MeasurementSize>::from_iterator([
                                    range_data.value(i)
                                ]),
                            ),
                        ));
                    }
                }
                "RangeRate" => {
                    let rate_data = batch
                        .column_by_name("Doppler (km/s)")
                        .unwrap()
                        .as_any()
                        .downcast_ref::<Float64Array>()
                        .unwrap();

                    // Set the measurements in the tracking arc
                    for i in 0..batch.num_rows() {
                        arc.measurements.push((
                            tracking_device.value(i).to_string(),
                            Msr::from_observation(
                                Epoch::from_tdb_seconds(epochs.value(i)),
                                OVector::<f64, Msr::MeasurementSize>::from_iterator([
                                    rate_data.value(i)
                                ]),
                            ),
                        ));
                    }
                }
                _ => unreachable!("should have errored earlier"),
            }
        }

        Ok(arc)
    }

    fn repr(&self) -> Vec<String> {
        let mut r = Vec::new();
        r.push(format!("File: {}", self.path));
        for (k, v) in &self.metadata {
            if k != "devices" {
                r.push(format!("{k}: {v}"));
            }
        }
        r
    }
}

impl Display for DynamicTrackingArc {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for item in self.repr() {
            writeln!(f, "{item}")?;
        }
        Ok(())
    }
}

#[cfg(feature = "python")]
#[pymethods]
impl DynamicTrackingArc {
    /// Initializes a new dynamic tracking arc from the provided parquet file
    #[new]
    fn new(path: String) -> Result<Self, NyxError> {
        Self::from_parquet(path).map_err(|e| NyxError::CustomError { msg: e.to_string() })
    }

    fn __repr__(&self) -> String {
        format!("{self}")
    }
}
