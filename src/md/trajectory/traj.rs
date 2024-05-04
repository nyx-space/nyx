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

use super::traj_it::TrajIterator;
use super::{ExportCfg, INTERPOLATION_SAMPLES};
use super::{Interpolatable, TrajError};
use crate::errors::NyxError;
use crate::io::watermark::pq_writer;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::md::prelude::{GuidanceMode, StateParameter};
use crate::md::EventEvaluator;
use crate::time::{Duration, Epoch, TimeSeries, TimeUnits};
use arrow::array::{Array, Float64Builder, StringBuilder};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::iter::Iterator;
use std::ops;
use std::path::{Path, PathBuf};
use std::sync::Arc;

/// Store a trajectory of any State.
#[derive(Clone, PartialEq)]
pub struct Traj<S: Interpolatable>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    /// Optionally name this trajectory
    pub name: Option<String>,
    /// We use a vector because we know that the states are produced in a chronological manner (the direction does not matter).
    pub states: Vec<S>,
}

impl<S: Interpolatable> Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    pub fn new() -> Self {
        Self {
            name: None,
            states: Vec::new(),
        }
    }
    /// Orders the states, can be used to store the states out of order
    pub fn finalize(&mut self) {
        // Remove duplicate epochs
        self.states.dedup_by(|a, b| a.epoch().eq(&b.epoch()));
        // And sort
        self.states.sort_by_key(|a| a.epoch());
    }

    /// Evaluate the trajectory at this specific epoch.
    pub fn at(&self, epoch: Epoch) -> Result<S, TrajError> {
        if self.states.is_empty() || self.first().epoch() > epoch || self.last().epoch() < epoch {
            return Err(TrajError::NoInterpolationData { epoch });
        }
        match self
            .states
            .binary_search_by(|state| state.epoch().cmp(&epoch))
        {
            Ok(idx) => {
                // Oh wow, we actually had this exact state!
                Ok(self.states[idx])
            }
            Err(idx) => {
                if idx == 0 || idx >= self.states.len() {
                    // The binary search returns where we should insert the data, so if it's at either end of the list, then we're out of bounds.
                    // This condition should have been handled by the check at the start of this function.
                    return Err(TrajError::NoInterpolationData { epoch });
                }
                // This is the closest index, so let's grab the items around it.
                // NOTE: This is essentially the same code as in ANISE for the Hermite SPK type 13

                // We didn't find it, so let's build an interpolation here.
                let num_left = INTERPOLATION_SAMPLES / 2;

                // Ensure that we aren't fetching out of the window
                let mut first_idx = idx.saturating_sub(num_left);
                let last_idx = self.states.len().min(first_idx + INTERPOLATION_SAMPLES);

                // Check that we have enough samples
                if last_idx == self.states.len() {
                    first_idx = last_idx.saturating_sub(2 * num_left);
                }

                let mut states = Vec::with_capacity(last_idx - first_idx);
                for idx in first_idx..last_idx {
                    states.push(self.states[idx]);
                }

                Ok(self.states[idx].interpolate(epoch, &states))
            }
        }
    }

    /// Returns the first state in this ephemeris
    pub fn first(&self) -> &S {
        // This is done after we've ordered the states we received, so we can just return the first state.
        self.states.first().unwrap()
    }

    /// Returns the last state in this ephemeris
    pub fn last(&self) -> &S {
        self.states.last().unwrap()
    }

    /// Creates an iterator through the trajectory by the provided step size
    pub fn every(&self, step: Duration) -> TrajIterator<S> {
        self.every_between(step, self.first().epoch(), self.last().epoch())
    }

    /// Creates an iterator through the trajectory by the provided step size between the provided bounds
    pub fn every_between(&self, step: Duration, start: Epoch, end: Epoch) -> TrajIterator<S> {
        TrajIterator {
            time_series: TimeSeries::inclusive(
                start.max(self.first().epoch()),
                end.min(self.last().epoch()),
                step,
            ),
            traj: self,
        }
    }

    /// Store this trajectory arc to a parquet file with the default configuration (depends on the state type, search for `export_params` in the documentation for details).
    pub fn to_parquet_simple<P: AsRef<Path>>(&self, path: P) -> Result<PathBuf, Box<dyn Error>> {
        self.to_parquet(path, None, ExportCfg::default())
    }

    /// Store this trajectory arc to a parquet file with the provided configuration
    pub fn to_parquet_with_cfg<P: AsRef<Path>>(
        &self,
        path: P,
        cfg: ExportCfg,
    ) -> Result<PathBuf, Box<dyn Error>> {
        self.to_parquet(path, None, cfg)
    }

    /// Store this trajectory arc to a parquet file with the provided configuration and event evaluators
    pub fn to_parquet<P: AsRef<Path>>(
        &self,
        path: P,
        events: Option<Vec<&dyn EventEvaluator<S>>>,
        cfg: ExportCfg,
    ) -> Result<PathBuf, Box<dyn Error>> {
        let tick = Epoch::now().unwrap();
        info!("Exporting trajectory to parquet file...");

        // Grab the path here before we move stuff.
        let path_buf = cfg.actual_path(path);

        // Build the schema
        let mut hdrs = vec![
            Field::new("Epoch:Gregorian UTC", DataType::Utf8, false),
            Field::new("Epoch:Gregorian TAI", DataType::Utf8, false),
            Field::new("Epoch:TAI (s)", DataType::Float64, false),
        ];

        let more_meta = Some(vec![(
            "Frame".to_string(),
            format!("{}", self.states[0].frame()),
        )]);

        let mut fields = match cfg.fields {
            Some(fields) => fields,
            None => S::export_params(),
        };

        // Check that we can retrieve this information
        fields.retain(|param| match self.first().value(*param) {
            Ok(_) => true,
            Err(_) => {
                warn!("Removed unavailable field `{param}` from trajectory export",);
                false
            }
        });

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

        // Build the states iterator -- this does require copying the current states but I can't either get a reference or a copy of all the states.
        let states = if cfg.start_epoch.is_some() || cfg.end_epoch.is_some() || cfg.step.is_some() {
            // Must interpolate the data!
            let start = cfg.start_epoch.unwrap_or_else(|| self.first().epoch());
            let end = cfg.end_epoch.unwrap_or_else(|| self.last().epoch());
            let step = cfg.step.unwrap_or_else(|| 1.minutes());
            self.every_between(step, start, end).collect::<Vec<S>>()
        } else {
            self.states.to_vec()
        };

        // Build all of the records

        // Epochs
        let mut utc_epoch = StringBuilder::new();
        let mut tai_epoch = StringBuilder::new();
        let mut tai_s = Float64Builder::new();
        for s in &states {
            utc_epoch.append_value(format!("{}", s.epoch()));
            tai_epoch.append_value(format!("{:x}", s.epoch()));
            tai_s.append_value(s.epoch().to_tai_seconds());
        }
        record.push(Arc::new(utc_epoch.finish()));
        record.push(Arc::new(tai_epoch.finish()));
        record.push(Arc::new(tai_s.finish()));

        // Add all of the fields
        for field in fields {
            if field == StateParameter::GuidanceMode {
                let mut guid_mode = StringBuilder::new();
                for s in &states {
                    guid_mode
                        .append_value(format!("{:?}", GuidanceMode::from(s.value(field).unwrap())));
                }
                record.push(Arc::new(guid_mode.finish()));
            } else {
                let mut data = Float64Builder::new();
                for s in &states {
                    data.append_value(s.value(field).unwrap());
                }
                record.push(Arc::new(data.finish()));
            }
        }

        info!(
            "Serialized {} states from {} to {}",
            states.len(),
            states.first().unwrap().epoch(),
            states.last().unwrap().epoch()
        );

        // Add all of the evaluated events
        if let Some(_events) = events {
            unimplemented!("Removed in ANISE updated");
        }

        // Serialize all of the devices and add that to the parquet file too.
        let mut metadata = HashMap::new();
        metadata.insert("Purpose".to_string(), "Trajectory data".to_string());
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

    /// Allows resampling this trajectory at a fixed interval instead of using the propagator step size.
    /// This may lead to aliasing due to the Nyquist–Shannon sampling theorem.
    pub fn resample(&self, step: Duration) -> Result<Self, NyxError> {
        if self.states.is_empty() {
            return Err(NyxError::Trajectory {
                source: TrajError::CreationError {
                    msg: "No trajectory to convert".to_string(),
                },
            });
        }

        let mut traj = Self::new();
        for state in self.every(step) {
            traj.states.push(state);
        }

        traj.finalize();

        Ok(traj)
    }

    /// Rebuilds this trajectory with the provided epochs.
    /// This may lead to aliasing due to the Nyquist–Shannon sampling theorem.
    pub fn rebuild(&self, epochs: &[Epoch]) -> Result<Self, NyxError> {
        if self.states.is_empty() {
            return Err(NyxError::Trajectory {
                source: TrajError::CreationError {
                    msg: "No trajectory to convert".to_string(),
                },
            });
        }

        let mut traj = Self::new();
        for epoch in epochs {
            traj.states.push(self.at(*epoch)?);
        }

        traj.finalize();

        Ok(traj)
    }

    /// Export the difference in RIC from of this trajectory compare to the "other" trajectory in parquet format.
    ///
    /// # Notes
    /// + The RIC frame accounts for the transport theorem by performing a finite differencing of the RIC frame.
    pub fn ric_diff_to_parquet<P: AsRef<Path>>(
        &self,
        other: &Self,
        path: P,
        cfg: ExportCfg,
    ) -> Result<PathBuf, Box<dyn Error>> {
        let tick = Epoch::now().unwrap();
        info!("Exporting trajectory to parquet file...");

        // Grab the path here before we move stuff.
        let path_buf = cfg.actual_path(path);

        // Build the schema
        let mut hdrs = vec![
            Field::new("Epoch:Gregorian UTC", DataType::Utf8, false),
            Field::new("Epoch:Gregorian TAI", DataType::Utf8, false),
            Field::new("Epoch:TAI (s)", DataType::Float64, false),
        ];

        // Add the RIC headers
        for coord in ["x", "y", "z"] {
            let mut meta = HashMap::new();
            meta.insert("unit".to_string(), "km".to_string());

            let field = Field::new(format!("delta_{coord}_ric (km)"), DataType::Float64, false)
                .with_metadata(meta);

            hdrs.push(field);
        }

        for coord in ["x", "y", "z"] {
            let mut meta = HashMap::new();
            meta.insert("unit".to_string(), "km/s".to_string());

            let field = Field::new(
                format!("delta_v{coord}_ric (km/s)"),
                DataType::Float64,
                false,
            )
            .with_metadata(meta);

            hdrs.push(field);
        }

        let more_meta = Some(vec![(
            "Frame".to_string(),
            format!("{}", self.states[0].frame()),
        )]);

        let mut cfg = cfg;

        let mut fields = match cfg.fields {
            Some(fields) => fields,
            None => S::export_params(),
        };

        // Remove disallowed field and check that we can retrieve this information
        fields.retain(|param| {
            param != &StateParameter::GuidanceMode && self.first().value(*param).is_ok()
        });

        for field in &fields {
            hdrs.push(field.to_field(more_meta.clone()));
        }

        // Build the schema
        let schema = Arc::new(Schema::new(hdrs));
        let mut record: Vec<Arc<dyn Array>> = Vec::new();

        // Ensure the times match.
        cfg.start_epoch = if self.first().epoch() > other.first().epoch() {
            Some(self.first().epoch())
        } else {
            Some(other.first().epoch())
        };

        cfg.end_epoch = if self.last().epoch() > other.last().epoch() {
            Some(other.last().epoch())
        } else {
            Some(self.last().epoch())
        };

        // Build the states iterator
        let step = cfg.step.unwrap_or_else(|| 1.minutes());
        let self_states = self
            .every_between(step, cfg.start_epoch.unwrap(), cfg.end_epoch.unwrap())
            .collect::<Vec<S>>();

        let other_states = other
            .every_between(step, cfg.start_epoch.unwrap(), cfg.end_epoch.unwrap())
            .collect::<Vec<S>>();

        // Build an array of all the RIC differences
        let mut ric_diff = Vec::with_capacity(other_states.len());
        for (ii, other_state) in other_states.iter().enumerate() {
            let self_orbit = *self_states[ii].orbit();
            let other_orbit = *other_state.orbit();

            let this_ric_diff = self_orbit.ric_difference(&other_orbit).map_err(Box::new)?;

            ric_diff.push(this_ric_diff);
        }

        // Build all of the records

        // Epochs (both match for self and others)
        let mut utc_epoch = StringBuilder::new();
        let mut tai_epoch = StringBuilder::new();
        let mut tai_s = Float64Builder::new();
        for s in &self_states {
            utc_epoch.append_value(format!("{}", s.epoch()));
            tai_epoch.append_value(format!("{:x}", s.epoch()));
            tai_s.append_value(s.epoch().to_tai_seconds());
        }
        record.push(Arc::new(utc_epoch.finish()));
        record.push(Arc::new(tai_epoch.finish()));
        record.push(Arc::new(tai_s.finish()));

        // Add the RIC data
        for coord_no in 0..6 {
            let mut data = Float64Builder::new();
            for this_ric_dff in &ric_diff {
                data.append_value(this_ric_dff.to_cartesian_pos_vel()[coord_no]);
            }
            record.push(Arc::new(data.finish()));
        }

        // Add all of the fields
        for field in fields {
            let mut data = Float64Builder::new();
            for (ii, self_state) in self_states.iter().enumerate() {
                let self_val = self_state.value(field).unwrap();
                let other_val = other_states[ii].value(field).unwrap();
                data.append_value(self_val - other_val);
            }
            record.push(Arc::new(data.finish()));
        }

        info!("Serialized {} states differences", self_states.len());

        // Serialize all of the devices and add that to the parquet file too.
        let mut metadata = HashMap::new();
        metadata.insert(
            "Purpose".to_string(),
            "Trajectory difference data".to_string(),
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

impl<S: Interpolatable> ops::Add for Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    type Output = Result<Traj<S>, NyxError>;

    /// Add one trajectory to another. If they do not overlap to within 10ms, a warning will be printed.
    fn add(self, other: Traj<S>) -> Self::Output {
        &self + &other
    }
}

impl<S: Interpolatable> ops::Add<&Traj<S>> for &Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    type Output = Result<Traj<S>, NyxError>;

    /// Add one trajectory to another, returns an error if the frames don't match
    fn add(self, other: &Traj<S>) -> Self::Output {
        if self.first().frame() != other.first().frame() {
            Err(NyxError::Trajectory {
                source: TrajError::CreationError {
                    msg: format!(
                        "Frame mismatch in add operation: {} != {}",
                        self.first().frame(),
                        other.first().frame()
                    ),
                },
            })
        } else {
            if self.last().epoch() < other.first().epoch() {
                let gap = other.first().epoch() - self.last().epoch();
                warn!(
                    "Resulting merged trajectory will have a time-gap of {} starting at {}",
                    gap,
                    self.last().epoch()
                );
            }

            let mut me = self.clone();
            // Now start adding the other segments while correcting the index
            for state in &other
                .states
                .iter()
                .filter(|s| s.epoch() > self.last().epoch())
                .collect::<Vec<&S>>()
            {
                me.states.push(**state);
            }
            me.finalize();

            Ok(me)
        }
    }
}

impl<S: Interpolatable> ops::AddAssign<&Traj<S>> for Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    /// Attempt to add two trajectories together and assign it to `self`
    ///
    /// # Warnings
    /// 1. This will panic if the frames mismatch!
    /// 2. This is inefficient because both `self` and `rhs` are cloned.
    fn add_assign(&mut self, rhs: &Self) {
        *self = (self.clone() + rhs.clone()).unwrap();
    }
}

impl<S: Interpolatable> fmt::Display for Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.states.is_empty() {
            write!(f, "Empty Trajectory!")
        } else {
            let dur = self.last().epoch() - self.first().epoch();
            write!(
                f,
                "Trajectory {}in {} from {} to {} ({}, or {:.3} s) [{} states]",
                match &self.name {
                    Some(name) => format!("of {name}"),
                    None => String::new(),
                },
                self.first().frame(),
                self.first().epoch(),
                self.last().epoch(),
                dur,
                dur.to_seconds(),
                self.states.len()
            )
        }
    }
}

impl<S: Interpolatable> fmt::Debug for Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self}",)
    }
}

impl<S: Interpolatable> Default for Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    fn default() -> Self {
        Self::new()
    }
}
