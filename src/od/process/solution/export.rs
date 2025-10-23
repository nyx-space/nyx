/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

use crate::io::watermark::pq_writer;
use crate::io::{ArrowSnafu, ExportCfg, ParquetSnafu, StdIOSnafu};
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName};
use crate::md::trajectory::Interpolatable;
use crate::md::StateParameter;
use crate::od::estimate::*;
use crate::State;
use crate::{od::*, Spacecraft};
use arrow::array::{Array, BooleanBuilder, Float64Builder, StringBuilder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use hifitime::TimeScale;
use log::{info, warn};
use msr::sensitivity::TrackerSensitivity;
use nalgebra::Const;
use parquet::arrow::ArrowWriter;
use snafu::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};

use super::ODSolution;

impl<MsrSize: DimName, Trk: TrackerSensitivity<Spacecraft, Spacecraft>>
    ODSolution<Spacecraft, KfEstimate<Spacecraft>, MsrSize, Trk>
where
    DefaultAllocator: Allocator<MsrSize>
        + Allocator<MsrSize, <Spacecraft as State>::Size>
        + Allocator<Const<1>, MsrSize>
        + Allocator<<Spacecraft as State>::Size>
        + Allocator<<Spacecraft as State>::Size, <Spacecraft as State>::Size>
        + Allocator<MsrSize, MsrSize>
        + Allocator<MsrSize, <Spacecraft as State>::Size>
        + Allocator<<Spacecraft as State>::Size, MsrSize>
        + Allocator<<Spacecraft as State>::Size>
        + Allocator<<Spacecraft as State>::VecLength>
        + Allocator<<Spacecraft as State>::Size, <Spacecraft as State>::Size>,
{
    /// Store the estimates and residuals in a parquet file
    pub fn to_parquet<P: AsRef<Path>>(&self, path: P, cfg: ExportCfg) -> Result<PathBuf, ODError> {
        ensure!(
            !self.estimates.is_empty(),
            TooFewMeasurementsSnafu {
                need: 1_usize,
                action: "exporting OD results"
            }
        );

        if self.estimates.len() != self.residuals.len() {
            return Err(ODError::ODConfigError {
                source: ConfigError::InvalidConfig {
                    msg: format!(
                        "Estimates ({}) and residuals ({}) are not aligned.",
                        self.estimates.len(),
                        self.residuals.len()
                    ),
                },
            });
        }

        if self.estimates.len() != self.gains.len() {
            return Err(ODError::ODConfigError {
                source: ConfigError::InvalidConfig {
                    msg: format!(
                        "Estimates ({}) and filter gains ({}) are not aligned.",
                        self.estimates.len(),
                        self.gains.len()
                    ),
                },
            });
        }

        if self.estimates.len() != self.filter_smoother_ratios.len() {
            return Err(ODError::ODConfigError {
                source: ConfigError::InvalidConfig {
                    msg: format!(
                        "Estimates ({}) and filter-smoother ratios ({}) are not aligned.",
                        self.estimates.len(),
                        self.filter_smoother_ratios.len()
                    ),
                },
            });
        }

        let tick = Epoch::now().unwrap();
        info!("Exporting orbit determination result to parquet file...");

        if cfg.step.is_some() {
            warn!("The `step` parameter in the export is not supported for orbit determination exports.");
        }

        // Grab the path here before we move stuff.
        let path_buf = cfg.actual_path(path);

        // Build the schema
        let mut hdrs = vec![Field::new("Epoch (UTC)", DataType::Utf8, false)];

        let frame = self.estimates[0].state().frame();

        let more_meta = Some(vec![
            (
                "Frame".to_string(),
                serde_dhall::serialize(&frame)
                    .static_type_annotation()
                    .to_string()
                    .map_err(|e| ODError::ODIOError {
                        source: InputOutputError::SerializeDhall {
                            what: format!("frame `{frame}`"),
                            err: e.to_string(),
                        },
                    })?,
            ),
            (
                "SRP Area (m2)".to_string(),
                self.estimates
                    .first()
                    .as_ref()
                    .unwrap()
                    .nominal_state()
                    .srp
                    .area_m2
                    .to_string(),
            ),
            (
                "Drag Area (m2)".to_string(),
                self.estimates
                    .first()
                    .as_ref()
                    .unwrap()
                    .nominal_state()
                    .drag
                    .area_m2
                    .to_string(),
            ),
        ]);

        let mut fields = match cfg.fields {
            Some(fields) => fields,
            None => Spacecraft::export_params(),
        };

        // Check that we can retrieve this information
        fields.retain(|param| match self.estimates[0].state().value(*param) {
            Ok(_) => param != &StateParameter::GuidanceMode,
            Err(_) => false,
        });

        for field in &fields {
            hdrs.push(field.to_field(more_meta.clone()));
        }

        let mut sigma_fields = fields.clone();
        // Check that we can retrieve this information
        sigma_fields.retain(|param| {
            !matches!(
                param,
                &StateParameter::X
                    | &StateParameter::Y
                    | &StateParameter::Z
                    | &StateParameter::VX
                    | &StateParameter::VY
                    | &StateParameter::VZ
            ) && self.estimates[0].sigma_for(*param).is_ok()
        });

        for field in &sigma_fields {
            hdrs.push(field.to_cov_field(more_meta.clone()));
        }

        let state_items = ["X", "Y", "Z", "Vx", "Vy", "Vz", "Cr", "Cd", "Mass"];
        let state_units = [
            "km", "km", "km", "km/s", "km/s", "km/s", "unitless", "unitless", "kg",
        ];
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

        let est_size = <Spacecraft as State>::Size::dim();

        let mut idx = 0;
        for i in 0..state_items.len() {
            for j in i..state_items.len() {
                hdrs.push(Field::new(
                    format!(
                        "Covariance {}*{} ({frame:x}) ({})",
                        state_items[i], state_items[j], cov_units[idx]
                    ),
                    DataType::Float64,
                    false,
                ));
                idx += 1;
            }
        }

        // Add the uncertainty in the integration frame
        for (i, coord) in state_items.iter().enumerate() {
            hdrs.push(Field::new(
                format!("Sigma {coord} ({frame:x}) ({})", state_units[i]),
                DataType::Float64,
                false,
            ));
        }

        // Add the position and velocity uncertainty in the RIC frame
        for (i, coord) in state_items.iter().enumerate().take(6) {
            hdrs.push(Field::new(
                format!("Sigma {coord} (RIC) ({})", state_units[i]),
                DataType::Float64,
                false,
            ));
        }

        // Add the fields of the residuals
        let mut msr_fields = Vec::new();
        for f in &self.measurement_types {
            msr_fields.push(
                f.to_field()
                    .with_nullable(true)
                    .with_name(format!("Prefit residual: {f:?} ({})", f.unit())),
            );
        }
        for f in &self.measurement_types {
            msr_fields.push(
                f.to_field()
                    .with_nullable(true)
                    .with_name(format!("Postfit residual: {f:?} ({})", f.unit())),
            );
        }
        for f in &self.measurement_types {
            msr_fields.push(
                f.to_field()
                    .with_nullable(true)
                    .with_name(format!("Measurement noise: {f:?} ({})", f.unit())),
            );
        }
        for f in &self.measurement_types {
            msr_fields.push(
                f.to_field()
                    .with_nullable(true)
                    .with_name(format!("Real observation: {f:?} ({})", f.unit())),
            );
        }
        for f in &self.measurement_types {
            msr_fields.push(
                f.to_field()
                    .with_nullable(true)
                    .with_name(format!("Computed observation: {f:?} ({})", f.unit())),
            );
        }

        msr_fields.push(Field::new("Residual ratio", DataType::Float64, true));
        msr_fields.push(Field::new("Residual Rejected", DataType::Boolean, true));
        msr_fields.push(Field::new("Tracker", DataType::Utf8, true));

        hdrs.append(&mut msr_fields);

        // Add the filter gain columns
        if self.measurement_types.len() == MsrSize::DIM {
            for i in 0..state_items.len() {
                for f in &self.measurement_types {
                    hdrs.push(Field::new(
                        format!(
                            "Gain {}*{f:?} ({}*{})",
                            state_items[i],
                            cov_units[i],
                            f.unit()
                        ),
                        DataType::Float64,
                        true,
                    ));
                }
            }
        } else {
            for state_item in &state_items {
                for j in 0..MsrSize::DIM {
                    hdrs.push(Field::new(
                        format!("Gain {state_item}*[{j}]"),
                        DataType::Float64,
                        true,
                    ));
                }
            }
        }

        // Add the filter-smoother ratio columns
        for i in 0..state_items.len() {
            hdrs.push(Field::new(
                format!(
                    "Filter-smoother ratio {} ({})",
                    state_items[i], cov_units[i],
                ),
                DataType::Float64,
                true,
            ));
        }

        // Build the schema
        let schema = Arc::new(Schema::new(hdrs));
        let mut record: Vec<Arc<dyn Array>> = Vec::new();

        // Build the states iterator -- this does require copying the current states but I can't either get a reference or a copy of all the states.
        let (estimates, residuals) =
            if cfg.start_epoch.is_some() || cfg.end_epoch.is_some() || cfg.step.is_some() {
                // Must interpolate the data!
                let start = cfg
                    .start_epoch
                    .unwrap_or_else(|| self.estimates.first().unwrap().state().epoch());
                let end = cfg
                    .end_epoch
                    .unwrap_or_else(|| self.estimates.last().unwrap().state().epoch());

                let mut residuals: Vec<Option<Residual<MsrSize>>> =
                    Vec::with_capacity(self.residuals.len());
                let mut estimates = Vec::with_capacity(self.estimates.len());

                for (estimate, residual) in self.estimates.iter().zip(self.residuals.iter()) {
                    if estimate.epoch() >= start && estimate.epoch() <= end {
                        estimates.push(*estimate);
                        residuals.push(residual.clone());
                    }
                }

                (estimates, residuals)
            } else {
                (self.estimates.to_vec(), self.residuals.to_vec())
            };

        // Build all of the records

        // Epochs
        let mut utc_epoch = StringBuilder::new();
        for s in &estimates {
            utc_epoch.append_value(s.epoch().to_time_scale(TimeScale::UTC).to_isoformat());
        }
        record.push(Arc::new(utc_epoch.finish()));

        // Add all of the fields
        for field in fields {
            let mut data = Float64Builder::new();
            for s in &estimates {
                data.append_value(
                    s.state()
                        .value(field)
                        .context(ODStateSnafu { action: "export" })?,
                );
            }
            record.push(Arc::new(data.finish()));
        }

        // Add all of the 1-sigma uncertainties
        for field in sigma_fields {
            let mut data = Float64Builder::new();
            for s in &estimates {
                data.append_value(s.sigma_for(field).unwrap());
            }
            record.push(Arc::new(data.finish()));
        }

        // Add the 1-sigma covariance in the integration frame
        for i in 0..est_size {
            for j in i..est_size {
                let mut data = Float64Builder::new();
                for s in &estimates {
                    data.append_value(s.covar()[(i, j)]);
                }
                record.push(Arc::new(data.finish()));
            }
        }

        // Add the sigma/uncertainty in the integration frame
        for i in 0..est_size {
            let mut data = Float64Builder::new();
            for s in &estimates {
                data.append_value(s.covar()[(i, i)].sqrt());
            }
            record.push(Arc::new(data.finish()));
        }

        // Add the sigma/uncertainty covariance in the RIC frame
        let mut ric_covariances = Vec::new();

        for s in &estimates {
            let dcm_ric2inertial = s
                .state()
                .orbit()
                .dcm_from_ric_to_inertial()
                .unwrap()
                .state_dcm();

            // Build the matrix view of the orbit part of the covariance.
            let cov = s.covar();
            let orbit_cov = cov.fixed_view::<6, 6>(0, 0);

            // Rotate back into the RIC frame
            let ric_covar = dcm_ric2inertial * orbit_cov * dcm_ric2inertial.transpose();
            ric_covariances.push(ric_covar);
        }

        // Now store the RIC covariance data.
        for i in 0..6 {
            let mut data = Float64Builder::new();
            for cov in ric_covariances.iter().take(estimates.len()) {
                data.append_value(cov[(i, i)].sqrt());
            }
            record.push(Arc::new(data.finish()));
        }

        // Finally, add the residuals.
        // Prefits
        for msr_type in &self.measurement_types {
            let mut data = Float64Builder::new();
            for resid_opt in &residuals {
                if let Some(resid) = resid_opt {
                    match resid.prefit(*msr_type) {
                        Some(prefit) => data.append_value(prefit),
                        None => data.append_null(),
                    };
                } else {
                    data.append_null();
                }
            }
            record.push(Arc::new(data.finish()));
        }
        // Postfit
        for msr_type in &self.measurement_types {
            let mut data = Float64Builder::new();
            for resid_opt in &residuals {
                if let Some(resid) = resid_opt {
                    match resid.postfit(*msr_type) {
                        Some(postfit) => data.append_value(postfit),
                        None => data.append_null(),
                    };
                } else {
                    data.append_null();
                }
            }
            record.push(Arc::new(data.finish()));
        }
        // Measurement noise
        for msr_type in &self.measurement_types {
            let mut data = Float64Builder::new();
            for resid_opt in &residuals {
                if let Some(resid) = resid_opt {
                    match resid.trk_noise(*msr_type) {
                        Some(noise) => data.append_value(noise),
                        None => data.append_null(),
                    };
                } else {
                    data.append_null();
                }
            }
            record.push(Arc::new(data.finish()));
        }
        // Real observation
        for msr_type in &self.measurement_types {
            let mut data = Float64Builder::new();
            for resid_opt in &residuals {
                if let Some(resid) = resid_opt {
                    match resid.real_obs(*msr_type) {
                        Some(postfit) => data.append_value(postfit),
                        None => data.append_null(),
                    };
                } else {
                    data.append_null();
                }
            }
            record.push(Arc::new(data.finish()));
        }
        // Computed observation
        for msr_type in &self.measurement_types {
            let mut data = Float64Builder::new();
            for resid_opt in &residuals {
                if let Some(resid) = resid_opt {
                    match resid.computed_obs(*msr_type) {
                        Some(postfit) => data.append_value(postfit),
                        None => data.append_null(),
                    };
                } else {
                    data.append_null();
                }
            }
            record.push(Arc::new(data.finish()));
        }
        // Residual ratio (unique entry regardless of the size)
        let mut data = Float64Builder::new();
        for resid_opt in &residuals {
            if let Some(resid) = resid_opt {
                data.append_value(resid.ratio);
            } else {
                data.append_null();
            }
        }
        record.push(Arc::new(data.finish()));

        // Residual acceptance (unique entry regardless of the size)
        let mut data = BooleanBuilder::new();
        for resid_opt in &residuals {
            if let Some(resid) = resid_opt {
                data.append_value(resid.rejected);
            } else {
                data.append_null();
            }
        }
        record.push(Arc::new(data.finish()));

        // Residual tracker (unique entry regardless of the size)
        let mut data = StringBuilder::new();
        for resid_opt in &residuals {
            if let Some(resid) = resid_opt {
                data.append_value(
                    resid
                        .tracker
                        .clone()
                        .unwrap_or("Undefined tracker".to_string()),
                );
            } else {
                data.append_null();
            }
        }
        record.push(Arc::new(data.finish()));

        // Add the filter gains
        for i in 0..est_size {
            for j in 0..MsrSize::DIM {
                let mut data = Float64Builder::new();
                for opt_k in &self.gains {
                    if let Some(k) = opt_k {
                        data.append_value(k[(i, j)]);
                    } else {
                        data.append_null();
                    }
                }
                record.push(Arc::new(data.finish()));
            }
        }

        // Add the filter-smoother consistency ratios
        for i in 0..est_size {
            let mut data = Float64Builder::new();
            for opt_fsr in &self.filter_smoother_ratios {
                if let Some(fsr) = opt_fsr {
                    data.append_value(fsr[i]);
                } else {
                    data.append_null();
                }
            }
            record.push(Arc::new(data.finish()));
        }

        info!("Serialized {} estimates and residuals", estimates.len());

        // Serialize all of the devices and add that to the parquet file too.
        let mut metadata = HashMap::new();
        metadata.insert(
            "Purpose".to_string(),
            "Orbit determination results".to_string(),
        );
        if let Some(add_meta) = cfg.metadata {
            for (k, v) in add_meta {
                metadata.insert(k, v);
            }
        }

        let props = pq_writer(Some(metadata));

        let file = File::create(&path_buf)
            .context(StdIOSnafu {
                action: "creating OD results file",
            })
            .context(ODIOSnafu)?;

        let mut writer = ArrowWriter::try_new(file, schema.clone(), props)
            .context(ParquetSnafu {
                action: "exporting OD results",
            })
            .context(ODIOSnafu)?;

        let batch = RecordBatch::try_new(schema, record)
            .context(ArrowSnafu {
                action: "writing OD results (building batch record)",
            })
            .context(ODIOSnafu)?;

        writer
            .write(&batch)
            .context(ParquetSnafu {
                action: "writing OD results",
            })
            .context(ODIOSnafu)?;

        writer
            .close()
            .context(ParquetSnafu {
                action: "closing OD results file",
            })
            .context(ODIOSnafu)?;

        // Return the path this was written to
        let tock_time = Epoch::now().unwrap() - tick;
        info!(
            "Orbit determination results written to {} in {tock_time}",
            path_buf.display()
        );
        Ok(path_buf)
    }
}
