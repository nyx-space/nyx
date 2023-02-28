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
    linalg::{allocator::Allocator, DefaultAllocator},
    md::{trajectory::InterpState, ui::Traj, StateParameter},
    NyxError,
};
use arrow::{array::Float64Array, record_batch::RecordBatchReader};
use hifitime::Epoch;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use std::fs::File;
use std::{collections::HashMap, error::Error, fmt::Display, path::Path};

#[cfg(feature = "python")]
use pyo3::prelude::*;

use super::odp::Cosm;

/// A dynamic trajectory allows loading a trajectory Parquet file and converting it
/// to the concrete trajectory state type when desired.
/// TODO: Trajectory export might benefit from storing the dynamics etc. as YAML to reload that?
#[cfg_attr(feature = "python", pyclass)]
#[derive(Clone)]
pub struct DynamicTrajectory {
    pub path: String,
    metadata: HashMap<String, String>,
}

impl DynamicTrajectory {
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
            metadata,
        };

        for item in me.repr() {
            info!("{item}");
        }

        Ok(me)
    }

    /// Reads through the loaded parquet file and attempts to convert to the provided concrete state.
    ///
    /// # Design limitations
    /// For Python compatibility, the file is actually re-read here, although it was read and closed during initialization.
    /// This is required because the parquet file reader is not clonable.
    pub fn to_traj<S>(&self) -> Result<Traj<S>, Box<dyn Error>>
    where
        S: InterpState,
        DefaultAllocator: Allocator<f64, S::VecLength>
            + Allocator<f64, S::Size>
            + Allocator<f64, S::Size, S::Size>,
    {
        // Check the schema
        let mut has_epoch = false; // Required
        let mut frame = None;

        let mut found_fields = [
            (StateParameter::X, false),
            (StateParameter::Y, false),
            (StateParameter::Z, false),
            (StateParameter::VX, false),
            (StateParameter::VY, false),
            (StateParameter::VZ, false),
            (StateParameter::FuelMass, false),
        ];

        let file = File::open(&self.path)?;

        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();

        let reader = builder.build()?;

        for field in &reader.schema().fields {
            if field.name().as_str() == "Epoch:TAI (s)" {
                has_epoch = true;
            } else {
                for potential_field in &mut found_fields {
                    if field.name() == potential_field.0.to_field(None).name() {
                        potential_field.1 = true;
                        if potential_field.0 != StateParameter::FuelMass {
                            if let Some(frame_info) = field.metadata().get("Frame") {
                                if frame.is_none() {
                                    frame = Some(frame_info.to_string())
                                } else if frame.as_ref().unwrap().as_str() != frame_info {
                                    return Err(Box::new(NyxError::FileUnreadable(format!(
                                        "Frame previous set to `{}` but set to `{}` in field `{}`",
                                        frame.unwrap(),
                                        frame_info,
                                        field.name()
                                    ))));
                                }
                            }
                        }
                        break;
                    }
                }
            }
        }

        if !has_epoch {
            return Err(Box::new(NyxError::FileUnreadable(
                "Missing `Epoch:TAI (s)` field".to_string(),
            )));
        } else if frame.is_none() {
            return Err(Box::new(NyxError::FileUnreadable(
                "Missing `Frame` in column metadata".to_string(),
            )));
        }

        for (field, exists) in found_fields.iter().take(found_fields.len() - 1) {
            if !exists {
                return Err(Box::new(NyxError::FileUnreadable(format!(
                    "Missing `{}` field",
                    field.to_field(None).name()
                ))));
            }
        }

        let sc_compat = found_fields.last().unwrap().1;

        let expected_type = std::any::type_name::<S>().split("::").last().unwrap();

        if expected_type == "Spacecraft" && !sc_compat {
            // Oops, missing fuel data (using the full call in case the field name changes in the future)
            return Err(Box::new(NyxError::FileUnreadable(format!(
                "Missing `{}` field",
                found_fields.last().unwrap().0.to_field(None).name()
            ))));
        }

        // At this stage, we know that the measurement is valid and the conversion is supported.
        let mut traj = Traj::default();

        // Now convert each batch on the fly
        for maybe_batch in reader {
            let batch = maybe_batch.unwrap();

            let epochs = batch
                .column_by_name("Epoch:TAI (s)")
                .unwrap()
                .as_any()
                .downcast_ref::<Float64Array>()
                .unwrap();

            let mut shared_data = vec![];

            for (field, _) in found_fields.iter().take(found_fields.len() - 1) {
                shared_data.push(
                    batch
                        .column_by_name(field.to_field(None).name())
                        .unwrap()
                        .as_any()
                        .downcast_ref::<Float64Array>()
                        .unwrap(),
                );
            }

            if sc_compat {
                // Read the fuel
                shared_data.push(
                    batch
                        .column_by_name("fuel (kg)")
                        .unwrap()
                        .as_any()
                        .downcast_ref::<Float64Array>()
                        .unwrap(),
                );
            }

            // Grab the frame
            // TODO: This is ugly and wrong because we might be using another frame info!
            // So arguably, the cosm to load should be in the metadata or passed to this function!
            let cosm = Cosm::de438();
            let frame = cosm.try_frame(frame.as_ref().unwrap().as_str())?;

            // Build the states
            for i in 0..batch.num_rows() {
                let mut state = S::zeros();
                state.set_epoch(Epoch::from_tai_seconds(epochs.value(i)));
                state.set_frame(frame);
                state.unset_stm(); // We don't have any STM data, so let's unset this.

                for (j, (param, exists)) in found_fields.iter().enumerate() {
                    if *exists {
                        state.set_value(param, shared_data[j].value(i))?;
                    }
                }

                traj.states.push(state);
            }
        }

        Ok(traj)
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

impl Display for DynamicTrajectory {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for item in self.repr() {
            writeln!(f, "{item}")?;
        }
        Ok(())
    }
}

#[cfg(feature = "python")]
#[pymethods]
impl DynamicTrajectory {
    /// Initializes a new dynamic trajectory from the provided parquet file
    #[new]
    fn new(path: String) -> Result<Self, NyxError> {
        Self::from_parquet(path).map_err(|e| NyxError::CustomError(e.to_string()))
    }

    fn __repr__(&self) -> String {
        format!("{self}")
    }
}
