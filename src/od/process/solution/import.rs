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

// Add this code within the impl block for ODSolution,
// potentially in a new file like `src/od/process/solution/import.rs`
// and ensure necessary imports are present.

use crate::io::{ArrowSnafu, InputOutputError, MissingDataSnafu, ParquetSnafu, StdIOSnafu};
use crate::linalg::allocator::Allocator;
use crate::linalg::{Const, DefaultAllocator, DimName, OMatrix, OVector, SMatrix};
use crate::od::estimate::*;
use crate::od::msr::MeasurementType;
use crate::od::*;
use crate::Spacecraft;
use anise::frames::Frame;
use anise::prelude::Orbit;
use anise::structure::spacecraft::{DragData, Mass, SRPData};
use arrow::array::RecordBatchReader;
use arrow::array::{Array, BooleanArray, Float64Array, StringArray};
use hifitime::Epoch;
use indexmap::IndexSet;
use msr::sensitivity::TrackerSensitivity;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use snafu::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::path::Path;
use std::str::FromStr;
use std::sync::Arc;

use super::ODSolution; // Needed for StateType bounds

// --- Function Definition ---

impl<MsrSize, Trk> ODSolution<Spacecraft, KfEstimate<Spacecraft>, MsrSize, Trk>
where
    MsrSize: DimName + std::fmt::Debug + Clone, // Added Debug+Clone for error messages/vec construction
    Trk: TrackerSensitivity<Spacecraft, Spacecraft> + Clone, // Added Clone for devices
    // Bounds needed for KfEstimate and Spacecraft
    DefaultAllocator: Allocator<MsrSize>
        + Allocator<MsrSize, <Spacecraft as State>::Size>
        + Allocator<Const<1>, MsrSize>
        + Allocator<<Spacecraft as State>::Size>
        + Allocator<<Spacecraft as State>::Size, <Spacecraft as State>::Size>
        + Allocator<MsrSize, MsrSize>
        + Allocator<<Spacecraft as State>::Size, MsrSize>
        + Allocator<<Spacecraft as State>::VecLength>,
    <DefaultAllocator as Allocator<<Spacecraft as State>::VecLength>>::Buffer<f64>: Send,
    <DefaultAllocator as Allocator<<Spacecraft as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<Spacecraft as State>::Size, <Spacecraft as State>::Size>>::Buffer<f64>: Copy,

{
    /// Loads an OD solution from a Parquet file created by `ODSolution::to_parquet`.
    ///
    /// The `devices` map must be provided by the caller as it contains potentially complex
    /// state (like Almanac references) that isn't serialized in the Parquet file.
    ///
    /// Note: This function currently assumes the StateType is `Spacecraft` and the
    /// estimate type is `KfEstimate<Spacecraft>`.
    pub fn from_parquet<P: AsRef<Path>>(
        path: P,
        devices: BTreeMap<String, Trk>,
    ) -> Result<Self, InputOutputError> {


     let file = File::open(&path).context(StdIOSnafu {
          action: "opening OD solution file",
      })?;

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

      // Check the schema
      let mut has_epoch = false; // Required
      let mut frame = None;
      let mut srp_area_m2 = None;
      let mut drag_area_m2 = None;

      let schema = builder.schema().clone();

      let reader = builder.build().context(ParquetSnafu {
          action: "building Parquet reader for OD results",
      })?;

      for field in &reader.schema().fields {
          if field.name().as_str() == "Epoch (UTC)" {
              has_epoch = true;
          } else {
            if let Some(frame_info) = field.metadata().get("Frame") {
                // Frame is expected to be serialized as Dhall.
                match serde_dhall::from_str(frame_info).parse::<Frame>() {
                    Err(e) => {
                        return Err(InputOutputError::ParseDhall {
                            data: frame_info.to_string(),
                            err: format!("{e}"),
                        })
                    }
                    Ok(deser_frame) => frame = Some(deser_frame),
                };
            }
            if let Some(info) = field.metadata().get("SRP Area (m2)") {
                srp_area_m2 = Some(info.parse::<f64>().unwrap_or(0.0));
            }
            if let Some(info) = field.metadata().get("Drag Area (m2)"){
                drag_area_m2 = Some(info.parse::<f64>().unwrap_or(0.0));
            }
        }
      }

      ensure!(
          has_epoch,
          MissingDataSnafu {
              which: "Epoch (UTC)"
          }
      );

      ensure!(
          frame.is_some(),
          MissingDataSnafu {
              which: "Frame in metadata"
          }
      );


        let mut estimates: Vec<KfEstimate<Spacecraft>> = Vec::new();
        let mut residuals: Vec<Option<Residual<MsrSize>>> = Vec::new();
        let mut gains: Vec<Option<OMatrix<f64, <Spacecraft as State>::Size, MsrSize>>> = Vec::new();
        let mut filter_smoother_ratios: Vec<Option<OVector<f64, <Spacecraft as State>::Size>>> =
            Vec::new();
        let mut measurement_types_found = IndexSet::new();

        let state_size = <Spacecraft as State>::Size::USIZE;

        // State item names used in column naming
        let state_items = ["X", "Y", "Z", "Vx", "Vy", "Vz", "Cr", "Cd", "Mass"];
        let mut cov_units = vec![];

        for i in 0..state_items.len() {
            for j in i..state_items.len() {
                let cov_unit = if i < 3 {
                    if j < 3 {
                        "km^2"
                    } else if (3..6).contains(&j) {
                        "km^2/s"
                    } else if j == 8 {
                        "km*kg"
                    } else {
                        "km"
                    }
                } else if (3..6).contains(&i) {
                    if (3..6).contains(&j) {
                        "km^2/s^2"
                    } else if j == 8 {
                        "km/s*kg"
                    } else {
                        "km/s"
                    }
                } else if i == 8 || j == 8 {
                    "kg^2"
                } else {
                    "unitless"
                };

                cov_units.push(cov_unit);
            }
        }

        // --- Pre-parse Measurement Types from Schema ---
        // Infer measurement types from residual column names
        for field in schema.fields() {
             if let Some(msr_type_str) = field.name().strip_prefix("Prefit residual: ") {
                 if let Some(bracket_pos) = msr_type_str.find(" (") {
                     let type_name = &msr_type_str[..bracket_pos];
                     if let Ok(msr_type) = MeasurementType::from_str(type_name) {
                          measurement_types_found.insert(msr_type);
                     } else {
                         warn!("Could not parse measurement type from column: {}", field.name());
                     }
                 }
             }
        }
        if measurement_types_found.is_empty() {
             warn!("Could not automatically detect any measurement types from Parquet column names. Residuals may not be loaded correctly.");
        } else {
             info!("Detected measurement types: {measurement_types_found:?}");
        }



        // while let Some(record_batch) = reader.next() {
        for record_batch in reader {
            let batch = record_batch.context(ArrowSnafu {
                action: "reading record batch from OD results",
            })?;

            let num_rows = batch.num_rows();

            // --- Helper to get column data ---
            let get_col = |name: &str| -> Result<Arc<dyn Array>, InputOutputError> {
                batch
                    .column_by_name(name)
                    .ok_or_else(|| InputOutputError::MissingData {
                        which: format!("column '{name}' in OD results"),
                    })
                    .cloned() // Clone the Arc<dyn Array>
            };

            // --- Extract Columns (handle potential errors) ---

            let epoch_col = get_col("Epoch (UTC)")?
                .as_any()
                .downcast_ref::<StringArray>()
                .ok_or_else(|| InputOutputError::ArrowError {
                     action: "downcasting Epoch column",
                     source: arrow::error::ArrowError::CastError("Could not cast Epoch to StringArray".to_string()),
                 })?.clone(); // Clone the concrete array

            // State component columns
            let x_col = get_col("x (km)")?.as_any().downcast_ref::<Float64Array>().ok_or_else(|| InputOutputError::ArrowError{action: "downcasting X", source: arrow::error::ArrowError::CastError("".to_string())})?.clone();
            let y_col = get_col("y (km)")?.as_any().downcast_ref::<Float64Array>().ok_or_else(|| InputOutputError::ArrowError{action: "downcasting Y", source: arrow::error::ArrowError::CastError("".to_string())})?.clone();
            let z_col = get_col("z (km)")?.as_any().downcast_ref::<Float64Array>().ok_or_else(|| InputOutputError::ArrowError{action: "downcasting Z", source: arrow::error::ArrowError::CastError("".to_string())})?.clone();
            let vx_col = get_col("vx (km/s)")?.as_any().downcast_ref::<Float64Array>().ok_or_else(|| InputOutputError::ArrowError{action: "downcasting VX", source: arrow::error::ArrowError::CastError("".to_string())})?.clone();
            let vy_col = get_col("vy (km/s)")?.as_any().downcast_ref::<Float64Array>().ok_or_else(|| InputOutputError::ArrowError{action: "downcasting VY", source: arrow::error::ArrowError::CastError("".to_string())})?.clone();
            let vz_col = get_col("vz (km/s)")?.as_any().downcast_ref::<Float64Array>().ok_or_else(|| InputOutputError::ArrowError{action: "downcasting VZ", source: arrow::error::ArrowError::CastError("".to_string())})?.clone();

            let cr_col = get_col("cr")?.as_any().downcast_ref::<Float64Array>().ok_or_else(|| InputOutputError::ArrowError{action: "downcasting Cr (unitless)", source: arrow::error::ArrowError::CastError("".to_string())})?.clone();
            let cd_col = get_col("cd")?.as_any().downcast_ref::<Float64Array>().ok_or_else(|| InputOutputError::ArrowError{action: "downcasting Cd (unitless)", source: arrow::error::ArrowError::CastError("".to_string())})?.clone();

            let dry_mass_col = get_col("dry_mass (kg)")?.as_any().downcast_ref::<Float64Array>().ok_or_else(|| InputOutputError::ArrowError{action: "downcasting dry_mass (kg)", source: arrow::error::ArrowError::CastError("".to_string())})?.clone();
            let prop_mass_col = get_col("prop_mass (kg)")?.as_any().downcast_ref::<Float64Array>().ok_or_else(|| InputOutputError::ArrowError{action: "downcasting prop_mass (kg)", source: arrow::error::ArrowError::CastError("".to_string())})?.clone();

            // Covariance columns (store them for iteration)
            let mut cov_cols = Vec::new();
            for i in 0..state_size {
                for j in i..state_size {
                     // Column names need frame and units, which were part of the export but hard to reconstruct perfectly here.
                     // We'll guess the base name format. Robust parsing would require metadata storage.
                     let base_name = format!("Covariance {}*{}", state_items[i], state_items[j]);
                     // Find the actual column name (it has frame/units appended)
                     let col_name = schema.fields().iter()
                         .find(|f| f.name().starts_with(&base_name))
                         .map(|f| f.name().as_str())
                         .ok_or_else(|| InputOutputError::ParquetError {
                              action: "seeking covariance column",
                              source: parquet::errors::ParquetError::General("Column not found".to_string()),
                          })?;
                     cov_cols.push(get_col(col_name)?.as_any().downcast_ref::<Float64Array>().ok_or_else(|| InputOutputError::ArrowError{action: "downcasting covariance column", source: arrow::error::ArrowError::CastError("".to_string())})?.clone());
                }
            }


            // Residual related columns
            let rejected_col = get_col("Residual Rejected").ok().and_then(|arr| arr.as_any().downcast_ref::<BooleanArray>().cloned());
            let tracker_col = get_col("Tracker").ok().and_then(|arr| arr.as_any().downcast_ref::<StringArray>().cloned());
            let ratio_col = get_col("Residual ratio").ok().and_then(|arr| arr.as_any().downcast_ref::<Float64Array>().cloned());

            let mut residual_data_cols: HashMap<MeasurementType, BTreeMap<String, Float64Array>> = HashMap::new();
            for msr_type in &measurement_types_found {
                 let mut type_cols = BTreeMap::new();
                 let prefixes = ["Prefit residual", "Postfit residual", "Measurement noise", "Real observation", "Computed observation"];
                 for prefix in prefixes {
                      // Again, guessing column name format
                      let base_name = format!("{}: {:?} ({})", prefix, msr_type, msr_type.unit());
                      if let Ok(col) = get_col(&base_name) {
                            if let Some(arr) = col.as_any().downcast_ref::<Float64Array>() {
                                 type_cols.insert(prefix.to_string(), arr.clone());
                            }
                      }
                 }
                 residual_data_cols.insert(*msr_type, type_cols);
            }

            // Gain columns (store for iteration)
            let mut gain_cols: Vec<Option<Float64Array>> = Vec::new();
            let mut gain_available = true;
            for i in 0..state_size {
                for f in &measurement_types_found {
                    // Guessing format - needs robust parsing or metadata
                    let base_name = format!(
                        "Gain {}*{f:?} ({}*{})",
                        state_items[i],
                        cov_units[i],
                        f.unit()
                    );
                    let col_name = schema.fields().iter()
                         .find(|f| f.name().starts_with(&base_name))
                         .map(|f| f.name().as_str());

                    if let Some(name) = col_name {
                          if let Ok(col) = get_col(name) {
                               gain_cols.push(col.as_any().downcast_ref::<Float64Array>().cloned());
                          } else {
                               gain_cols.push(None); // Column missing
                               gain_available = false;
                          }
                    } else {
                          // If *any* gain column is missing, assume no gains were stored (e.g., smoother run)
                          gain_available = false;
                          break; // No need to check further gain columns
                    }
                }
                if !gain_available { break; }
            }
            if !gain_available { gain_cols.clear(); } // Ensure empty if incomplete


            // FSR columns (store for iteration)
            let mut fsr_cols: Vec<Option<Float64Array>> = Vec::new();
            let mut fsr_available = true;
            // for i in 0..state_size {
            for state_item in state_items.iter().take(state_size) {
                 // Guessing format
                 let base_name = format!("Filter-smoother ratio {state_item}");
                 let col_name = schema.fields().iter()
                     .find(|f| f.name().starts_with(&base_name))
                     .map(|f| f.name().as_str());

                 if let Some(name) = col_name {
                      if let Ok(col) = get_col(name) {
                            fsr_cols.push(col.as_any().downcast_ref::<Float64Array>().cloned());
                      } else {
                            fsr_cols.push(None);
                            fsr_available = false;
                      }
                 } else {
                      fsr_available = false;
                      break;
                 }
            }
             if !fsr_available { fsr_cols.clear(); }


            // --- Iterate through rows in the batch ---
            for i in 0..num_rows {

                let epoch = Epoch::from_gregorian_str(epoch_col.value(i)).map_err(|e| {
                    InputOutputError::Inconsistency {
                        msg: format!("{e} when parsing epoch"),
                    }
                })?;

                // Reconstruct spacecraft
                let nominal_state = Spacecraft::builder().orbit(
                    Orbit::cartesian(
                        x_col.value(i),
                        y_col.value(i),
                        z_col.value(i),
                        vx_col.value(i),
                        vy_col.value(i),
                        vz_col.value(i),
                        epoch,
                        frame.expect("somehow frame isn't set")
                    )).mass(
                        Mass::from_dry_and_prop_masses(
                            dry_mass_col.value(i),
                            prop_mass_col.value(i))
                    ).srp(SRPData {
                        area_m2: srp_area_m2.expect("somehow srp area isn't set"),
                        coeff_reflectivity: cr_col.value(i)
                    }).drag(DragData{
                        area_m2: drag_area_m2.expect("somehow dragarea isn't set"),
                        coeff_drag: cd_col.value(i)
                    }).build();

                // Reconstruct Covariance
                let mut covar = SMatrix::<f64, 9, 9>::zeros();
                let mut cov_col_idx = 0;
                for row in 0..state_size {
                    for col in row..state_size {
                        let val = cov_cols[cov_col_idx].value(i);
                        covar[(row, col)] = val;
                        if row != col {
                            covar[(col, row)] = val; // Symmetric
                        }
                        cov_col_idx += 1;
                    }
                }

                // Reconstruct KfEstimate
                let estimate = KfEstimate {
                    nominal_state,
                    state_deviation: OVector::<f64, Const<9>>::zeros(), // Deviation not stored
                    covar,
                    covar_bar: covar, // Not stored, use covar
                    stm: OMatrix::<f64, Const<9>, Const<9>>::identity(), // Not stored
                    predicted: false, // Not stored
                };
                estimates.push(estimate);

                // Reconstruct Residual (if applicable)
                let is_rejected_opt = rejected_col.as_ref().and_then(|col| if col.is_valid(i) { Some(col.value(i)) } else { None });
                let tracker_opt = tracker_col.as_ref().and_then(|col| if col.is_valid(i) { Some(col.value(i).to_string()) } else { None });
                let ratio_opt = ratio_col.as_ref().and_then(|col| if col.is_valid(i) { Some(col.value(i)) } else { None });

                let current_residual: Option<Residual<MsrSize>> = if let (Some(is_rejected), Some(tracker), Some(ratio)) = (is_rejected_opt, tracker_opt.clone(), ratio_opt) {
                     // It's a measurement update
                     let mut prefit_vec = OVector::<f64, MsrSize>::zeros();
                     let mut postfit_vec = OVector::<f64, MsrSize>::zeros();
                     let mut noise_vec = OVector::<f64, MsrSize>::zeros();
                     let mut real_obs_vec = OVector::<f64, MsrSize>::zeros();
                     let mut comp_obs_vec = OVector::<f64, MsrSize>::zeros();
                     let mut current_msr_types = IndexSet::with_capacity(MsrSize::USIZE);

                     let mut msr_idx = 0;
                     for (msr_type, type_cols) in &residual_data_cols {
                           if msr_idx >= MsrSize::USIZE { break; } // Should not happen if MsrSize matches data

                           // Check if data exists for this type *at this row*
                           let prefit_val = type_cols.get("Prefit residual").and_then(|col| if col.is_valid(i) { Some(col.value(i)) } else { None });
                           let postfit_val = type_cols.get("Postfit residual").and_then(|col| if col.is_valid(i) { Some(col.value(i)) } else { None });
                           let noise_val = type_cols.get("Measurement noise").and_then(|col| if col.is_valid(i) { Some(col.value(i)) } else { None });
                           let real_val = type_cols.get("Real observation").and_then(|col| if col.is_valid(i) { Some(col.value(i)) } else { None });
                           let comp_val = type_cols.get("Computed observation").and_then(|col| if col.is_valid(i) { Some(col.value(i)) } else { None });

                           // Only include if *at least one* value is present for this type in this row
                           if prefit_val.is_some() || postfit_val.is_some() || noise_val.is_some() || real_val.is_some() || comp_val.is_some() {
                                prefit_vec[msr_idx] = prefit_val.unwrap_or(f64::NAN); // Or handle differently
                                postfit_vec[msr_idx] = postfit_val.unwrap_or(f64::NAN);
                                noise_vec[msr_idx] = noise_val.unwrap_or(f64::NAN);
                                real_obs_vec[msr_idx] = real_val.unwrap_or(f64::NAN);
                                comp_obs_vec[msr_idx] = comp_val.unwrap_or(f64::NAN);
                                current_msr_types.insert(*msr_type);
                                msr_idx += 1;
                           }
                     }

                     // Resize vectors if fewer than MsrSize types were found for this row
                     // This part is tricky and depends on how multi-type residuals were stored.
                     // Assuming vectors should always have MsrSize, filled potentially with NaN.

                     let resid = Residual {
                          epoch,
                          prefit: prefit_vec,
                          postfit: postfit_vec,
                          tracker_msr_noise: noise_vec,
                          ratio,
                          real_obs: real_obs_vec,
                          computed_obs: comp_obs_vec,
                          msr_types: current_msr_types, // Store types found for this row
                          rejected: is_rejected,
                          tracker: Some(tracker),
                     };
                     Some(resid)

                } else {
                    // Not all parts of a residual were present, assume it was a time update
                    None
                };
                residuals.push(current_residual);


                // Reconstruct Gain (if available)
                let current_gain: Option<OMatrix<f64, <Spacecraft as State>::Size, MsrSize>> = if gain_available && !gain_cols.is_empty() {
                     let mut gain_mat = OMatrix::<f64, <Spacecraft as State>::Size, MsrSize>::zeros();
                     let mut all_valid = true;
                     let mut col_idx = 0;
                     'gain_outer: for row in 0..state_size {
                          for col in 0..MsrSize::USIZE {
                               if let Some(gain_col) = &gain_cols[col_idx] {
                                    if gain_col.is_valid(i) {
                                         gain_mat[(row, col)] = gain_col.value(i);
                                    } else {
                                         all_valid = false;
                                         break 'gain_outer; // Null found, entire matrix is None
                                    }
                               } else {
                                    // Should not happen if gain_available is true, but safeguard
                                    all_valid = false;
                                    break 'gain_outer;
                               }
                               col_idx += 1;
                          }
                     }
                     if all_valid { Some(gain_mat) } else { None }
                } else { None };
                gains.push(current_gain);

                // Reconstruct FSR (if available)
                let current_fsr: Option<OVector<f64, <Spacecraft as State>::Size>> = if fsr_available && !fsr_cols.is_empty() {
                      let mut fsr_vec = OVector::<f64, <Spacecraft as State>::Size>::zeros();
                      let mut all_valid = true;
                      for k in 0..state_size {
                           if let Some(fsr_col) = &fsr_cols[k] {
                                if fsr_col.is_valid(i) {
                                     fsr_vec[k] = fsr_col.value(i);
                                } else {
                                     all_valid = false;
                                     break;
                                }
                           } else {
                                all_valid = false;
                                break;
                           }
                      }
                      if all_valid { Some(fsr_vec) } else { None }
                 } else { None };
                 filter_smoother_ratios.push(current_fsr);

            } // End row loop
        } // End batch loop


        // --- Final Construction ---
        Ok(ODSolution {
            estimates,
            residuals,
            gains,
            filter_smoother_ratios,
            devices, // Provided by user
            measurement_types: measurement_types_found, // Determined from columns
        })
    }
}
