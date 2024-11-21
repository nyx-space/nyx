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
use arrow::{
    array::{Float64Array, PrimitiveArray, StringArray},
    datatypes,
    record_batch::RecordBatchReader,
};
use core::fmt;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::File;
use std::path::Path;

use hifitime::Epoch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use snafu::{ensure, ResultExt};

use crate::io::{ArrowSnafu, InputOutputError, MissingDataSnafu, ParquetSnafu, StdIOSnafu};

use super::{measurement::Measurement, MeasurementType};

/// Tracking data storing all of measurements as a B-Tree.
pub struct TrackingDataArc {
    /// All measurements in this data arc
    pub measurements: BTreeMap<Epoch, Measurement>,
    /// Source file if loaded from a file
    pub source: Option<String>,
}

impl TrackingDataArc {
    /// Loads a tracking arc from its serialization in parquet.
    pub fn from_parquet<P: AsRef<Path> + fmt::Display>(path: P) -> Result<Self, InputOutputError> {
        // Read the file since we closed it earlier
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
        let mut rate_avail = false;
        for field in &reader.schema().fields {
            match field.name().as_str() {
                "Epoch (UTC)" => has_epoch = true,
                "Tracking device" => has_tracking_dev = true,
                "Range (km)" => range_avail = true,
                "Doppler (km/s)" => rate_avail = true,
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
            range_avail || rate_avail,
            MissingDataSnafu {
                which: "`Range (km)` or `Doppler (km/s)`"
            }
        );

        let mut measurements = BTreeMap::new();

        // Now convert each batch on the fly
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

            let doppler_data: Option<&PrimitiveArray<datatypes::Float64Type>> = if rate_avail {
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
                    data: HashMap::new(),
                };

                if range_avail {
                    measurement
                        .data
                        .insert(MeasurementType::Range, range_data.unwrap().value(i));
                }

                if rate_avail {
                    measurement
                        .data
                        .insert(MeasurementType::Doppler, doppler_data.unwrap().value(i));
                }

                measurements.insert(epoch, measurement);
            }
        }

        Ok(Self {
            measurements,
            source: Some(path.to_string()),
        })
    }

    /// Compute the list of unique aliases in this tracking arc
    pub fn unique_aliases(&self) -> HashSet<String> {
        let mut aliases = HashSet::new();
        for msr in self.measurements.values() {
            aliases.insert(msr.tracker.clone());
        }
        aliases
    }

    /// Returns the start epoch of this tracking arc
    pub fn start_epoch(&self) -> Option<Epoch> {
        self.measurements.first_key_value().map(|(k, _)| *k)
    }

    /// Returns the end epoch of this tracking arc
    pub fn end_epoch(&self) -> Option<Epoch> {
        self.measurements.last_key_value().map(|(k, _)| *k)
    }

    /// Returns the number of measurements in this data arc
    pub fn len(&self) -> usize {
        self.measurements.len()
    }

    /// Returns whether this arc has no measurements.
    pub fn is_empty(&self) -> bool {
        self.measurements.is_empty()
    }
}

impl fmt::Display for TrackingDataArc {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            write!(f, "Empty tracking arc")
        } else {
            let start = self.start_epoch().unwrap();
            let end = self.end_epoch().unwrap();
            let src = match &self.source {
                Some(src) => format!(" (source: {src})"),
                None => String::new(),
            };
            write!(
                f,
                "Tracking arc with {} measurements over {} (from {start} to {end}) with trackers {:?}{src}",
                self.len(),
                end - start,
                self.unique_aliases()
            )
        }
    }
}

#[cfg(test)]
mod ut_tracker {
    use super::TrackingDataArc;

    #[test]
    fn test_lro_data() {
        let trk = TrackingDataArc::from_parquet("04_lro_simulated_tracking.parquet").unwrap();
        println!("{trk}");
        // TODO: Add test nulls in specific sets, and missing columns.
    }
}
