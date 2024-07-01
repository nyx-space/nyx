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

use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::errors::{MonteCarloError, NoSuccessfulRunsSnafu, StateError};
use crate::io::watermark::pq_writer;
use crate::io::{ExportCfg, InputOutputError};
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::md::prelude::GuidanceMode;
use crate::md::trajectory::{Interpolatable, Traj};
use crate::md::{EventEvaluator, StateParameter};
use crate::propagators::PropagationError;
use crate::time::{Duration, Epoch, TimeUnits};
use anise::almanac::Almanac;
use anise::constants::frames::EARTH_J2000;
use arrow::array::{Array, Float64Builder, Int32Builder, StringBuilder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
pub use rstats::Stats;
use snafu::ensure;

use super::DispersedState;

/// A structure storing the result of a single Monte Carlo run
pub struct Run<S: Interpolatable, R>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
    <DefaultAllocator as Allocator<f64, S::VecLength>>::Buffer: Send,
{
    /// The index of this run
    pub index: usize,
    /// The original dispersed state
    pub dispersed_state: DispersedState<S>,
    /// The result from this run
    pub result: Result<R, PropagationError>,
}

/// A structure of Monte Carlo results
pub struct Results<S: Interpolatable, R>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
    <DefaultAllocator as Allocator<f64, S::VecLength>>::Buffer: Send,
{
    /// Raw data from each run, sorted by run index for O(1) access to each run
    pub runs: Vec<Run<S, R>>,
    /// Name of this scenario
    pub scenario: String,
}

/// A structure that stores the result of a propagation segment of a Monte Carlo.
pub struct PropResult<S: Interpolatable>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
    <DefaultAllocator as Allocator<f64, S::VecLength>>::Buffer: Send,
{
    pub state: S,
    pub traj: Traj<S>,
}

impl<S: Interpolatable> Results<S, PropResult<S>>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<usize, S::Size, S::Size>
        + Allocator<f64, S::VecLength>,
    <DefaultAllocator as Allocator<f64, S::VecLength>>::Buffer: Send,
{
    /// Returns the value of the requested state parameter for all trajectories from `start` to `end` every `step` and
    /// using the value of `value_if_run_failed` if set and skipping that run if the run failed
    pub fn every_value_of_between(
        &self,
        param: StateParameter,
        step: Duration,
        start: Epoch,
        end: Epoch,
        value_if_run_failed: Option<f64>,
    ) -> Vec<f64> {
        let mut report = Vec::with_capacity(self.runs.len());
        for run in &self.runs {
            match &run.result {
                Ok(r) => {
                    for state in r.traj.every_between(step, start, end) {
                        match state.value(param) {
                            Ok(val) => report.push(val),
                            Err(e) => match value_if_run_failed {
                                Some(val) => report.push(val),
                                None => {
                                    warn!("run #{}: {}, skipping {} in report", run.index, e, param)
                                }
                            },
                        }
                    }
                }
                Err(e) => match value_if_run_failed {
                    Some(val) => report.push(val),
                    None => warn!(
                        "run #{} failed with {}, skipping {} in report",
                        run.index, e, param
                    ),
                },
            }
        }
        report
    }

    /// Returns the value of the requested state parameter for all trajectories from the start to the end of each trajectory and
    /// using the value of `value_if_run_failed` if set and skipping that run if the run failed
    pub fn every_value_of(
        &self,
        param: StateParameter,
        step: Duration,
        value_if_run_failed: Option<f64>,
    ) -> Vec<f64> {
        let mut report = Vec::with_capacity(self.runs.len());
        for run in &self.runs {
            match &run.result {
                Ok(r) => {
                    for state in r.traj.every(step) {
                        match state.value(param) {
                            Ok(val) => report.push(val),
                            Err(e) => match value_if_run_failed {
                                Some(val) => report.push(val),
                                None => {
                                    warn!("run #{}: {}, skipping {} in report", run.index, e, param)
                                }
                            },
                        }
                    }
                }
                Err(e) => match value_if_run_failed {
                    Some(val) => report.push(val),
                    None => warn!(
                        "run #{} failed with {}, skipping {} in report",
                        run.index, e, param
                    ),
                },
            }
        }
        report
    }

    /// Returns the value of the requested state parameter for the first state and
    /// using the value of `value_if_run_failed` if set and skipping that run if the run failed
    pub fn first_values_of(
        &self,
        param: StateParameter,
        value_if_run_failed: Option<f64>,
    ) -> Vec<f64> {
        let mut report = Vec::with_capacity(self.runs.len());
        for run in &self.runs {
            match &run.result {
                Ok(r) => match r.traj.first().value(param) {
                    Ok(val) => report.push(val),
                    Err(e) => match value_if_run_failed {
                        Some(val) => report.push(val),
                        None => {
                            warn!("run #{}: {}, skipping {} in report", run.index, e, param)
                        }
                    },
                },
                Err(e) => match value_if_run_failed {
                    Some(val) => report.push(val),
                    None => warn!(
                        "run #{} failed with {}, skipping {} in report",
                        run.index, e, param
                    ),
                },
            }
        }
        report
    }

    /// Returns the value of the requested state parameter for the first state and
    /// using the value of `value_if_run_failed` if set and skipping that run if the run failed
    pub fn last_values_of(
        &self,
        param: StateParameter,
        value_if_run_failed: Option<f64>,
    ) -> Vec<f64> {
        let mut report = Vec::with_capacity(self.runs.len());
        for run in &self.runs {
            match &run.result {
                Ok(r) => match r.traj.last().value(param) {
                    Ok(val) => report.push(val),
                    Err(e) => match value_if_run_failed {
                        Some(val) => report.push(val),
                        None => {
                            warn!("run #{}: {}, skipping {} in report", run.index, e, param)
                        }
                    },
                },
                Err(e) => match value_if_run_failed {
                    Some(val) => report.push(val),
                    None => warn!(
                        "run #{} failed with {}, skipping {} in report",
                        run.index, e, param
                    ),
                },
            }
        }
        report
    }

    /// Returns the dispersion values of the requested state parameter
    pub fn dispersion_values_of(&self, param: StateParameter) -> Result<Vec<f64>, MonteCarloError> {
        let mut report = Vec::with_capacity(self.runs.len());
        'run_loop: for run in &self.runs {
            for (dparam, val) in &run.dispersed_state.actual_dispersions {
                if dparam == &param {
                    report.push(*val);
                    continue 'run_loop;
                }
            }
            // Oh, this parameter was not found!
            return Err(MonteCarloError::StateError {
                source: StateError::Unavailable { param },
            });
        }
        Ok(report)
    }

    pub fn to_parquet<P: AsRef<Path>>(
        &self,
        path: P,
        events: Option<Vec<&dyn EventEvaluator<S>>>,
        cfg: ExportCfg,
        almanac: Arc<Almanac>,
    ) -> Result<PathBuf, Box<dyn Error>> {
        let tick = Epoch::now().unwrap();
        info!("Exporting Monte Carlo results to parquet file...");

        // Grab the path here before we move stuff.
        let path_buf = cfg.actual_path(path);

        // Build the schema
        let mut hdrs = vec![
            Field::new("Epoch:Gregorian UTC", DataType::Utf8, false),
            Field::new("Epoch:Gregorian TAI", DataType::Utf8, false),
            Field::new("Epoch:TAI (s)", DataType::Float64, false),
            Field::new("Monte Carlo Run Index", DataType::Int32, false),
        ];

        // Use the first successful run to build up some data shared for all
        let mut frame = EARTH_J2000;
        let mut fields = match cfg.fields {
            Some(fields) => fields,
            None => S::export_params(),
        };

        let mut start = None;
        let mut end = None;

        // Literally all of the states of all the successful runs.
        let mut all_states: Vec<S> = vec![];
        let mut run_indexes: Vec<i32> = vec![];

        for run in &self.runs {
            if let Ok(success) = &run.result {
                if start.is_none() {
                    // No need to check other states.
                    frame = success.state.frame();

                    // Check that we can retrieve this information
                    fields.retain(|param| match success.state.value(*param) {
                        Ok(_) => true,
                        Err(_) => false,
                    });

                    start = Some(success.traj.first().epoch());
                    end = Some(success.state.epoch());
                }

                // Build the states iterator.
                let states =
                    if cfg.start_epoch.is_some() || cfg.end_epoch.is_some() || cfg.step.is_some() {
                        // Must interpolate the data!
                        let start = cfg.start_epoch.unwrap_or_else(|| start.unwrap());
                        let end = cfg.end_epoch.unwrap_or_else(|| end.unwrap());
                        let step = cfg.step.unwrap_or_else(|| 1.minutes());
                        success
                            .traj
                            .every_between(step, start, end)
                            .collect::<Vec<S>>()
                    } else {
                        success.traj.states.to_vec()
                    };
                all_states.extend(states.iter());
                run_indexes.push(run.index as i32);
            }
        }

        ensure!(
            start.is_some(),
            NoSuccessfulRunsSnafu {
                action: "export",
                num_runs: self.runs.len()
            }
        );

        let more_meta = Some(vec![(
            "Frame".to_string(),
            serde_dhall::serialize(&frame).to_string().map_err(|e| {
                Box::new(InputOutputError::SerializeDhall {
                    what: format!("frame `{frame}`"),
                    err: e.to_string(),
                })
            })?,
        )]);

        for field in &fields {
            hdrs.push(field.to_field(more_meta.clone()));
        }

        if let Some(events) = events.as_ref() {
            for event in events {
                let field = Field::new(format!("{event}"), DataType::Float64, false);
                hdrs.push(field);
            }
        }

        // Build the schema
        let schema = Arc::new(Schema::new(hdrs));
        let mut record: Vec<Arc<dyn Array>> = Vec::new();

        // Build all of the records

        // Epochs
        let mut utc_epoch = StringBuilder::new();
        let mut tai_epoch = StringBuilder::new();
        let mut tai_s = Float64Builder::new();
        let mut idx_col = Int32Builder::new();
        for (sno, s) in all_states.iter().enumerate() {
            utc_epoch.append_value(format!("{}", s.epoch()));
            tai_epoch.append_value(format!("{:x}", s.epoch()));
            tai_s.append_value(s.epoch().to_tai_seconds());
            // Copy this a bunch of times because all columns must have the same length
            // TODO: I need to keep track of when a new run actually start here!
            idx_col.append_value(run_indexes[sno]);
        }
        record.push(Arc::new(utc_epoch.finish()));
        record.push(Arc::new(tai_epoch.finish()));
        record.push(Arc::new(tai_s.finish()));
        record.push(Arc::new(idx_col.finish()));

        // Add all of the fields
        for field in fields {
            if field == StateParameter::GuidanceMode {
                let mut guid_mode = StringBuilder::new();
                for s in &all_states {
                    guid_mode
                        .append_value(format!("{:?}", GuidanceMode::from(s.value(field).unwrap())));
                }
                record.push(Arc::new(guid_mode.finish()));
            } else {
                let mut data = Float64Builder::new();
                for s in &all_states {
                    data.append_value(s.value(field).unwrap());
                }
                record.push(Arc::new(data.finish()));
            }
        }

        info!(
            "Serialized {} states from {} to {}",
            all_states.len(),
            start.unwrap(),
            end.unwrap()
        );

        // Add all of the evaluated events
        if let Some(events) = events {
            info!("Evaluating {} event(s)", events.len());
            for event in events {
                let mut data = Float64Builder::new();
                for s in &all_states {
                    data.append_value(event.eval(s, almanac.clone()).map_err(Box::new)?);
                }
                record.push(Arc::new(data.finish()));
            }
        }

        // Serialize all of the devices and add that to the parquet file too.
        let mut metadata = HashMap::new();
        metadata.insert(
            "Purpose".to_string(),
            "Monte Carlo Trajectory data".to_string(),
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
            "Trajectory written to {} in {tock_time}",
            path_buf.display()
        );
        Ok(path_buf)
    }
}
