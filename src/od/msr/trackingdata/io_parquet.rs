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
use crate::io::watermark::pq_writer;
use crate::io::{ArrowSnafu, InputOutputError, MissingDataSnafu, ParquetSnafu, StdIOSnafu};
use crate::io::{EmptyDatasetSnafu, ExportCfg};
use crate::od::msr::{Measurement, MeasurementType};
use arrow::array::{Array, Float64Builder, StringBuilder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use arrow::{
    array::{Float64Array, PrimitiveArray, StringArray},
    datatypes,
    record_batch::RecordBatchReader,
};
use hifitime::prelude::Epoch;
use hifitime::TimeScale;
use indexmap::IndexMap;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;
use snafu::{ensure, ResultExt};
use std::collections::{BTreeMap, HashMap};
use std::error::Error;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use super::TrackingDataArc;

impl TrackingDataArc {
    /// Loads a tracking arc from its serialization in parquet.
    pub fn from_parquet<P: AsRef<Path>>(path: P) -> Result<Self, InputOutputError> {
        let file = File::open(&path).context(StdIOSnafu {
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
        let mut doppler_avail = false;
        let mut az_avail = false;
        let mut el_avail = false;
        for field in &reader.schema().fields {
            match field.name().as_str() {
                "Epoch (UTC)" => has_epoch = true,
                "Tracking device" => has_tracking_dev = true,
                "Range (km)" => range_avail = true,
                "Doppler (km/s)" => doppler_avail = true,
                "Azimuth (deg)" => az_avail = true,
                "Elevation (deg)" => el_avail = true,
                _ => {}
            }
        }

        ensure!(
            has_epoch,
            MissingDataSnafu {
                which: "Epoch (UTC)"
            }
        );

        ensure!(
            has_tracking_dev,
            MissingDataSnafu {
                which: "Tracking device"
            }
        );

        ensure!(
            range_avail || doppler_avail || az_avail || el_avail,
            MissingDataSnafu {
                which: "`Range (km)` or `Doppler (km/s)` or `Azimuth (deg)` or `Elevation (deg)`"
            }
        );

        let mut measurements = BTreeMap::new();

        // We can safely unwrap the columns since we've checked for their existance just before.
        for maybe_batch in reader {
            let batch = maybe_batch.context(ArrowSnafu {
                action: "reading batch of tracking data",
            })?;

            let tracking_device = batch
                .column_by_name("Tracking device")
                .unwrap()
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();

            let epochs = batch
                .column_by_name("Epoch (UTC)")
                .unwrap()
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();

            let range_data: Option<&PrimitiveArray<datatypes::Float64Type>> = if range_avail {
                Some(
                    batch
                        .column_by_name("Range (km)")
                        .unwrap()
                        .as_any()
                        .downcast_ref::<Float64Array>()
                        .unwrap(),
                )
            } else {
                None
            };

            let doppler_data: Option<&PrimitiveArray<datatypes::Float64Type>> = if doppler_avail {
                Some(
                    batch
                        .column_by_name("Doppler (km/s)")
                        .unwrap()
                        .as_any()
                        .downcast_ref::<Float64Array>()
                        .unwrap(),
                )
            } else {
                None
            };

            let azimuth_data: Option<&PrimitiveArray<datatypes::Float64Type>> = if az_avail {
                Some(
                    batch
                        .column_by_name("Azimuth (deg)")
                        .unwrap()
                        .as_any()
                        .downcast_ref::<Float64Array>()
                        .unwrap(),
                )
            } else {
                None
            };

            let elevation_data: Option<&PrimitiveArray<datatypes::Float64Type>> = if el_avail {
                Some(
                    batch
                        .column_by_name("Elevation (deg)")
                        .unwrap()
                        .as_any()
                        .downcast_ref::<Float64Array>()
                        .unwrap(),
                )
            } else {
                None
            };

            // Set the measurements in the tracking arc
            for i in 0..batch.num_rows() {
                let epoch = Epoch::from_gregorian_str(epochs.value(i)).map_err(|e| {
                    InputOutputError::Inconsistency {
                        msg: format!("{e} when parsing epoch"),
                    }
                })?;

                let mut measurement = Measurement {
                    epoch,
                    tracker: tracking_device.value(i).to_string(),
                    data: IndexMap::new(),
                };

                if range_avail {
                    measurement
                        .data
                        .insert(MeasurementType::Range, range_data.unwrap().value(i));
                }

                if doppler_avail {
                    measurement
                        .data
                        .insert(MeasurementType::Doppler, doppler_data.unwrap().value(i));
                }

                if az_avail {
                    measurement
                        .data
                        .insert(MeasurementType::Azimuth, azimuth_data.unwrap().value(i));
                }

                if el_avail {
                    measurement
                        .data
                        .insert(MeasurementType::Elevation, elevation_data.unwrap().value(i));
                }

                measurements.insert(epoch, measurement);
            }
        }

        Ok(Self {
            measurements,
            source: Some(path.as_ref().to_path_buf().display().to_string()),
        })
    }
    /// Store this tracking arc to a parquet file.
    pub fn to_parquet_simple<P: AsRef<Path>>(&self, path: P) -> Result<PathBuf, Box<dyn Error>> {
        self.to_parquet(path, ExportCfg::default())
    }

    /// Store this tracking arc to a parquet file, with optional metadata and a timestamp appended to the filename.
    pub fn to_parquet<P: AsRef<Path>>(
        &self,
        path: P,
        cfg: ExportCfg,
    ) -> Result<PathBuf, Box<dyn Error>> {
        ensure!(
            !self.is_empty(),
            EmptyDatasetSnafu {
                action: "exporting tracking data arc"
            }
        );

        let path_buf = cfg.actual_path(path);

        if cfg.step.is_some() {
            warn!("The `step` parameter in the export is not supported for tracking arcs.");
        }

        if cfg.fields.is_some() {
            warn!("The `fields` parameter in the export is not supported for tracking arcs.");
        }

        // Build the schema
        let mut hdrs = vec![
            Field::new("Epoch (UTC)", DataType::Utf8, false),
            Field::new("Tracking device", DataType::Utf8, false),
        ];

        let msr_types = self.unique_types();
        let mut msr_fields = msr_types
            .iter()
            .map(|msr_type| msr_type.to_field())
            .collect::<Vec<Field>>();

        hdrs.append(&mut msr_fields);

        // Build the schema
        let schema = Arc::new(Schema::new(hdrs));
        let mut record: Vec<Arc<dyn Array>> = Vec::new();

        // Build the measurement iterator

        let measurements =
            if cfg.start_epoch.is_some() || cfg.end_epoch.is_some() || cfg.step.is_some() {
                let start = cfg
                    .start_epoch
                    .unwrap_or_else(|| self.start_epoch().unwrap());
                let end = cfg.end_epoch.unwrap_or_else(|| self.end_epoch().unwrap());

                info!("Exporting measurements from {start} to {end}.");

                self.measurements
                    .range(start..end)
                    .map(|(k, v)| (*k, v.clone()))
                    .collect()
            } else {
                self.measurements.clone()
            };

        // Build all of the records

        // Epochs
        let mut utc_epoch = StringBuilder::new();
        for epoch in measurements.keys() {
            utc_epoch.append_value(epoch.to_time_scale(TimeScale::UTC).to_isoformat());
        }
        record.push(Arc::new(utc_epoch.finish()));

        // Device names
        let mut device_names = StringBuilder::new();
        for m in measurements.values() {
            device_names.append_value(m.tracker.clone());
        }
        record.push(Arc::new(device_names.finish()));

        // Measurement data, column by column
        for msr_type in msr_types {
            let mut data_builder = Float64Builder::new();

            for m in measurements.values() {
                match m.data.get(&msr_type) {
                    Some(value) => data_builder.append_value(*value),
                    None => data_builder.append_null(),
                };
            }
            record.push(Arc::new(data_builder.finish()));
        }

        // Serialize all of the devices and add that to the parquet file too.
        let mut metadata = HashMap::new();
        metadata.insert("Purpose".to_string(), "Tracking Arc Data".to_string());
        if let Some(add_meta) = cfg.metadata {
            for (k, v) in add_meta {
                metadata.insert(k, v);
            }
        }

        let props = pq_writer(Some(metadata));

        let file = File::create(&path_buf)?;

        let mut writer = ArrowWriter::try_new(file, schema.clone(), props).unwrap();

        let batch = RecordBatch::try_new(schema, record)?;
        writer.write(&batch)?;
        writer.close()?;

        info!("Serialized {self} to {}", path_buf.display());

        // Return the path this was written to
        Ok(path_buf)
    }
}
