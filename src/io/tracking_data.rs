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
    linalg::{allocator::Allocator, DefaultAllocator, DimName, OVector},
    od::{msr::arc::TrackingArc, Measurement},
    NyxError, State,
};
use arrow::{
    array::{Float64Array, StringArray},
    record_batch::RecordBatchReader,
};
use hifitime::Epoch;
use parquet::arrow::arrow_reader::{ParquetRecordBatchReader, ParquetRecordBatchReaderBuilder};
use serde::{Deserialize, Serialize};
use std::{collections::HashMap, error::Error, fmt::Display, path::Path};
use std::{fs::File, sync::Arc};

use crate::od::SimMeasurement;

/// A dynamic measurement is only used when reading tracking data file
#[derive(Debug, Serialize, Deserialize)]
pub struct DynamicMeasurement {
    measurement_type: String,
    epoch: Epoch,
    data: Vec<f64>,
}

impl DynamicMeasurement {
    /// Initializes a new dynamic measurement from a concrete measurement.
    /// TODO: Is this ever needed? The measurements will be serialized directly, so I'm not sure I'll ever need this function.
    pub fn new<M: SimMeasurement>(measurement: M) -> DynamicMeasurement
    where
        Self: Sized,
        DefaultAllocator: Allocator<f64, M::MeasurementSize>
            + Allocator<f64, <M::State as State>::Size>
            + Allocator<f64, <M::State as State>::Size, <M::State as State>::Size>
            + Allocator<f64, <M::State as State>::VecLength>
            + Allocator<f64, M::MeasurementSize, <M::State as State>::Size>,
    {
        let mut data = Vec::new();
        for obs in measurement.observation().iter() {
            data.push(*obs);
        }
        // for data in measurement.
        DynamicMeasurement {
            measurement_type: std::any::type_name::<M>().to_string(),
            data,
            epoch: measurement.epoch(),
        }
    }

    /// Attempts to convert this dynamic measurement into the concrete turbofish measurement
    pub fn try_into<M: SimMeasurement>(self) -> Result<M, Box<dyn std::error::Error>>
    where
        Self: Sized,
        DefaultAllocator: Allocator<f64, M::MeasurementSize>
            + Allocator<f64, <M::State as State>::Size>
            + Allocator<f64, <M::State as State>::Size, <M::State as State>::Size>
            + Allocator<f64, <M::State as State>::VecLength>
            + Allocator<f64, M::MeasurementSize, <M::State as State>::Size>,
    {
        let expected_type = std::any::type_name::<M>().to_string();
        if self.measurement_type != expected_type {
            return Err(format!(
                "Expected measurement of type {}, but got {}",
                expected_type, self.measurement_type
            )
            .into());
        }
        // Try to build the vector of the observation
        if self.data.len() != M::MeasurementSize::USIZE {
            panic!();
        }
        let obs = OVector::<f64, M::MeasurementSize>::from_iterator(self.data);
        let measurement = M::from_observation(self.epoch, obs);
        Ok(measurement)
    }
}

pub struct DynamicTrackingArc {
    pub device_cfg: String,
    pub path: String,
    metadata: HashMap<String, String>,
    reader: ParquetRecordBatchReader,
}

impl DynamicTrackingArc {
    pub fn from_parquet<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn Error>> {
        let file = File::open(&path)?;

        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();

        let mut metadata = HashMap::new();
        // Build the custom metadata
        if let Some(file_metadata) = builder.metadata().file_metadata().key_value_metadata() {
            for key_value in file_metadata {
                if !key_value.key.starts_with("ARROW:") {
                    metadata.insert(
                        key_value.key.clone(),
                        key_value.value.clone().unwrap_or("[unset]".to_string()),
                    );
                }
            }
        }

        let me = Self {
            path: path.as_ref().to_string_lossy().to_string(),
            reader: builder.build()?,
            metadata,
            device_cfg: "".to_string(),
        };

        info!("{me}");

        Ok(me)

        // let record_batch = reader.next().unwrap().unwrap();

        // println!("Read {} records.", record_batch.num_rows());
    }

    /// Reads through the loaded parquet file and attempts to convert to the provided tracking arc.
    pub fn to_tracking_arc<Msr>(mut self) -> Result<TrackingArc<Msr>, NyxError>
    where
        Msr: Measurement,
        DefaultAllocator: Allocator<f64, Msr::MeasurementSize>,
    {
        // Check the schema
        let mut has_epoch = false;
        let mut has_tracking_dev = false;
        let mut range_avail = false;
        let mut rate_avail = false;
        for field in &self.reader.schema().fields {
            match field.name().as_str() {
                "Epoch TDB (s)" => has_epoch = true,
                "Tracking device" => has_tracking_dev = true,
                "Range (km)" => range_avail = true,
                "Range rate (km/s)" => rate_avail = true,
                _ => {}
            }
        }

        if !has_epoch {
            return Err(NyxError::FileUnreadable(
                "Missing `Epoch TDB (s)` field".to_string(),
            ));
        } else if !has_tracking_dev {
            return Err(NyxError::FileUnreadable(
                "Missing `Tracking device` field".to_string(),
            ));
        } else if !range_avail && !rate_avail {
            return Err(NyxError::FileUnreadable(
                "Missing both range and range rate data".to_string(),
            ));
        }

        let expected_type = std::any::type_name::<Msr>().to_string();
        match expected_type.as_str() {
            "StdMeasurement" => {
                if !range_avail || !rate_avail {
                    return Err(NyxError::FileUnreadable(
                        "Cannot convert to simultaneous range and range rate, missing one of them"
                            .to_string(),
                    ));
                }
            }
            "RangeMsr" => {
                if !range_avail {
                    return Err(NyxError::FileUnreadable(
                        "Cannot convert to range measurement: data missing".to_string(),
                    ));
                }
            }
            "RangeRate" => {
                if !rate_avail {
                    return Err(NyxError::FileUnreadable(
                        "Cannot convert to range rate measurement: data missing".to_string(),
                    ));
                }
            }
            _ => {
                return Err(NyxError::FileUnreadable(format!(
                    "Yet unsupported measurement {expected_type}"
                )));
            }
        }

        // At this stage, we know that the measurement is valid and the conversion is supported.

        let mut arc = TrackingArc {
            device_cfg: self.metadata["devices"].clone(),
            measurements: Vec::new(),
        };
        // Now convert each data on the fly

        let batch = self.reader.next().unwrap().unwrap();

        let tracking_device = batch
            .column_by_name("Tracking device")
            .unwrap()
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();

        let epochs = batch
            .column_by_name("Epoch TDB (s)")
            .unwrap()
            .as_any()
            .downcast_ref::<Float64Array>()
            .unwrap();

        let expected_type = std::any::type_name::<Msr>().to_string();
        match expected_type.as_str() {
            "StdMeasurement" => {
                let range_data = batch
                    .column_by_name("Range (km)")
                    .unwrap()
                    .as_any()
                    .downcast_ref::<Float64Array>()
                    .unwrap();

                let rate_data = batch
                    .column_by_name("Range rate (km/s)")
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
                                range_data.value(i),
                                rate_data.value(i),
                            ]),
                        ),
                    ));
                }
            }
            "RangeMsr" => {
                if !range_avail {
                    return Err(NyxError::FileUnreadable(
                        "Cannot convert to range measurement: data missing".to_string(),
                    ));
                }
            }
            "RangeRate" => {
                if !rate_avail {
                    return Err(NyxError::FileUnreadable(
                        "Cannot convert to range rate measurement: data missing".to_string(),
                    ));
                }
            }
            _ => {
                return Err(NyxError::FileUnreadable(format!(
                    "Yet unsupported measurement {expected_type}"
                )));
            }
        }

        Ok(arc)
    }
}

impl Display for DynamicTrackingArc {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "File: {}", self.path)?;
        for (k, v) in &self.metadata {
            if k != "devices" {
                write!(f, "{k}: {v}\n")?;
            }
        }
        Ok(())
    }
}
