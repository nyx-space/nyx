/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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
use crate::io::ExportCfg;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName};
use crate::md::prelude::Frame;
use crate::md::trajectory::Interpolatable;
pub use crate::od::estimate::*;
pub use crate::od::ground_station::*;
pub use crate::od::snc::*;
pub use crate::od::*;
use crate::propagators::error_ctrl::ErrorCtrl;
pub use crate::time::{Duration, Unit};
use crate::State;
use arrow::array::{Array, Float64Builder, StringBuilder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::ops::Add;
use std::path::{Path, PathBuf};

use super::ODProcess;

impl<
        'a,
        D: Dynamics,
        E: ErrorCtrl,
        Msr: Measurement,
        A: DimName,
        S: EstimateFrom<D::StateType, Msr> + Interpolatable,
        K: Filter<S, A, Msr::MeasurementSize>,
    > ODProcess<'a, D, E, Msr, A, S, K>
where
    D::StateType: Interpolatable + Add<OVector<f64, <S as State>::Size>, Output = D::StateType>,
    <DefaultAllocator as Allocator<f64, <D::StateType as State>::VecLength>>::Buffer: Send,
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, S::Size>
        + Allocator<f64, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, <D::StateType as State>::Size>
        + Allocator<f64, Msr::MeasurementSize, <S as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, Msr::MeasurementSize>
        + Allocator<f64, <S as State>::Size, Msr::MeasurementSize>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<usize, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>
        + Allocator<f64, A>
        + Allocator<f64, A, A>
        + Allocator<f64, <D::StateType as State>::Size, A>
        + Allocator<f64, A, <D::StateType as State>::Size>
        + Allocator<f64, <S as State>::Size>
        + Allocator<f64, <S as State>::VecLength>
        + Allocator<f64, <S as State>::Size, <S as State>::Size>
        + Allocator<f64, <S as State>::Size, A>
        + Allocator<f64, A, <S as State>::Size>,
{
    /// Store the estimates and residuals in a parquet file
    pub fn to_parquet<P: AsRef<Path>>(
        &self,
        path: P,
        cfg: ExportCfg,
    ) -> Result<PathBuf, Box<dyn Error>> {
        if self.estimates.is_empty() {
            return Err(Box::new(NyxError::CustomError(
                "No data: run the ODProcess before exporting it.".to_string(),
            )));
        } else if self.estimates.len() != self.residuals.len() {
            return Err(Box::new(NyxError::CustomError(
                "Estimates and residuals are not aligned.".to_string(),
            )));
        }

        let tick = Epoch::now().unwrap();
        info!("Exporting orbit determination result to parquet file...");

        if cfg.step.is_some() {
            warn!("The `step` parameter in the export is not supported for orbit determination exports.");
        }

        // Grab the path here before we move stuff.
        let path_buf = cfg.actual_path(path);

        // Build the schema
        let mut hdrs = vec![
            Field::new("Epoch:Gregorian UTC", DataType::Utf8, false),
            Field::new("Epoch:Gregorian TAI", DataType::Utf8, false),
            Field::new("Epoch:TAI (s)", DataType::Float64, false),
        ];

        let frame_name = self.estimates[0].state().frame();

        let more_meta = Some(vec![("Frame".to_string(), format!("{}", frame_name))]);

        let mut fields = match cfg.fields {
            Some(fields) => fields,
            None => S::export_params(),
        };

        // Check that we can retrieve this information
        fields.retain(|param| match self.estimates[0].state().value(*param) {
            Ok(_) => true,
            Err(_) => {
                warn!("Removed unavailable field `{param}` from orbit determination export",);
                false
            }
        });

        for field in &fields {
            hdrs.push(field.to_field(more_meta.clone()));
        }

        let cov_hdrs = match <S as State>::Size::dim() {
            6 => {
                // Add orbit 1-sigma covariance info, plotting to perform computations as desired
                vec![
                    "Covariance XX",
                    "Covariance YZ",
                    "Covariance YY",
                    "Covariance ZX",
                    "Covariance ZY",
                    "Covariance ZZ",
                    "Covariance VxX",
                    "Covariance VxY",
                    "Covariance VxZ",
                    "Covariance VxVx",
                    "Covariance VyX",
                    "Covariance VyY",
                    "Covariance VyZ",
                    "Covariance VyVx",
                    "Covariance VyVy",
                    "Covariance VzX",
                    "Covariance VzY",
                    "Covariance VzZ",
                    "Covariance VzVx",
                    "Covariance VzVy",
                    "Covariance VzVz",
                ]
            }
            _ => todo!(
                "exporting a state of size {} is not yet supported",
                <S as State>::Size::dim()
            ),
        };

        // Add the covariance in the integration frame
        for hdr in &cov_hdrs {
            hdrs.push(Field::new(
                format!("{hdr} ({frame_name})"),
                DataType::Float64,
                false,
            ));
        }

        // Add the covariance in the RIC frame
        for hdr in &cov_hdrs {
            hdrs.push(Field::new(format!("{hdr} (RIC)"), DataType::Float64, false));
        }

        // Add the fields of the residuals
        let mut msr_fields = Vec::new();
        for f in Msr::fields() {
            msr_fields.push(
                f.clone()
                    .with_nullable(true)
                    .with_name(format!("Prefit residual: {}", f.name())),
            );
        }
        for f in Msr::fields() {
            msr_fields.push(
                f.clone()
                    .with_nullable(true)
                    .with_name(format!("Postfit residual: {}", f.name())),
            );
        }

        msr_fields.push(Field::new("Residual ratio", DataType::Float64, true));

        hdrs.append(&mut msr_fields);

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

                let mut residuals: Vec<Option<Residual<Msr::MeasurementSize>>> =
                    Vec::with_capacity(self.residuals.len());
                let mut estimates = Vec::with_capacity(self.estimates.len());

                for (estimate, residual) in self.estimates.iter().zip(self.residuals.iter()) {
                    if estimate.epoch() >= start && estimate.epoch() <= end {
                        estimates.push(estimate.clone());
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
        let mut tai_epoch = StringBuilder::new();
        let mut tai_s = Float64Builder::new();
        for s in &estimates {
            utc_epoch.append_value(format!("{}", s.epoch()));
            tai_epoch.append_value(format!("{:x}", s.epoch()));
            tai_s.append_value(s.epoch().to_tai_seconds());
        }
        record.push(Arc::new(utc_epoch.finish()));
        record.push(Arc::new(tai_epoch.finish()));
        record.push(Arc::new(tai_s.finish()));

        // Add all of the fields
        for field in fields {
            let mut data = Float64Builder::new();
            for s in &estimates {
                data.append_value(s.state().value(field).unwrap());
            }
            record.push(Arc::new(data.finish()));
        }
        // Add the 1-sigma covariance in the integration frame
        for i in 0..<S as State>::Size::dim() {
            for j in 0..<S as State>::Size::dim() {
                let mut data = Float64Builder::new();
                for s in &estimates {
                    data.append_value(s.covar()[(i, j)]);
                }
                record.push(Arc::new(data.finish()));
            }
        }
        // Add the 1-sigma covariance in the RIC frame
        let mut ric_covariances = Vec::new();

        for s in &estimates {
            let dcm6x6 = s
                .state()
                .orbit()
                .dcm6x6_from_traj_frame(Frame::RIC)
                .unwrap();
            // Create a full DCM and only rotate the orbit part of it.
            let mut dcm = OMatrix::<f64, S::Size, S::Size>::identity();
            for i in 0..6 {
                for j in 0..6 {
                    dcm[(i, j)] = dcm6x6[(i, j)];
                }
            }
            let ric_covar = &dcm * s.covar() * &dcm.transpose();

            ric_covariances.push(ric_covar);
        }

        // Now store the RIC covariance data.
        for i in 0..<S as State>::Size::dim() {
            for j in 0..<S as State>::Size::dim() {
                let mut data = Float64Builder::new();
                for k in 0..estimates.len() {
                    data.append_value(ric_covariances[k][(i, j)]);
                }
                record.push(Arc::new(data.finish()));
            }
        }

        // Finally, add the residuals.
        // Prefits
        for i in 0..Msr::MeasurementSize::dim() {
            let mut data = Float64Builder::new();
            for resid_opt in &residuals {
                if let Some(resid) = resid_opt {
                    data.append_value(resid.prefit[i]);
                } else {
                    data.append_null();
                }
            }
        }
        // Postfit
        for i in 0..Msr::MeasurementSize::dim() {
            let mut data = Float64Builder::new();
            for resid_opt in &residuals {
                if let Some(resid) = resid_opt {
                    data.append_value(resid.postfit[i]);
                } else {
                    data.append_null();
                }
            }
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

        let file = File::create(&path_buf)?;
        let mut writer = ArrowWriter::try_new(file, schema.clone(), props).unwrap();

        let batch = RecordBatch::try_new(schema, record)?;
        writer.write(&batch)?;
        writer.close()?;

        // Return the path this was written to
        let tock_time = Epoch::now().unwrap() - tick;
        info!(
            "Orbit determination results written to {} in {tock_time}",
            path_buf.display()
        );
        Ok(path_buf)
    }
}
