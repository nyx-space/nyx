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
use super::{InterpState, TrajError};
use crate::cosmic::{Cosm, Frame, Orbit, Spacecraft};
use crate::errors::NyxError;
use crate::io::watermark::pq_writer;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::md::events::EventEvaluator;
use crate::md::ui::GuidanceMode;
use crate::md::StateParameter;
use crate::time::{Duration, Epoch, TimeSeries, TimeUnits, Unit};
use arrow::array::{ArrayRef, Float64Array, StringArray};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use hifitime::prelude::{Format, Formatter};
use parquet::arrow::ArrowWriter;
use rayon::prelude::*;
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::iter::Iterator;
use std::ops;
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::mpsc::channel;
use std::sync::Arc;
use std::time::Instant;

/// Store a trajectory of any State.
#[derive(Clone, PartialEq)]
pub struct Traj<S: InterpState>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    /// Optionally name this trajectory
    pub name: Option<String>,
    /// We use a vector because we know that the states are produced in a chronological manner (the direction does not matter).
    pub states: Vec<S>,
}

impl<S: InterpState> Traj<S>
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
    pub fn at(&self, epoch: Epoch) -> Result<S, NyxError> {
        if self.states.is_empty() || self.first().epoch() > epoch || self.last().epoch() < epoch {
            return Err(NyxError::Trajectory(TrajError::NoInterpolationData(epoch)));
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
                    return Err(NyxError::Trajectory(TrajError::NoInterpolationData(epoch)));
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

                self.states[idx].interpolate(epoch, &states)
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
            time_series: TimeSeries::inclusive(start, end, step),
            traj: self,
        }
    }

    /// Find the exact state where the request event happens. The event function is expected to be monotone in the provided interval because we find the event using a Brent solver.
    #[allow(clippy::identity_op)]
    pub fn find_bracketed<E>(&self, start: Epoch, end: Epoch, event: &E) -> Result<S, NyxError>
    where
        E: EventEvaluator<S>,
    {
        let max_iter = 50;

        // Helper lambdas, for f64s only
        let has_converged =
            |x1: f64, x2: f64| (x1 - x2).abs() <= event.epoch_precision().to_seconds();
        let arrange = |a: f64, ya: f64, b: f64, yb: f64| {
            if ya.abs() > yb.abs() {
                (a, ya, b, yb)
            } else {
                (b, yb, a, ya)
            }
        };

        let xa_e = start;
        let xb_e = end;

        // Search in seconds (convert to epoch just in time)
        let mut xa = 0.0;
        let mut xb = (xb_e - xa_e).to_seconds();
        // Evaluate the event at both bounds
        let mut ya = event.eval(&self.at(xa_e)?);
        let mut yb = event.eval(&self.at(xb_e)?);

        // Check if we're already at the root
        if ya.abs() <= event.value_precision().abs() {
            debug!(
                "{event} -- found with |{ya}| < {} @ {xa_e}",
                event.value_precision().abs()
            );
            return self.at(xa_e);
        } else if yb.abs() <= event.value_precision().abs() {
            debug!(
                "{event} -- found with |{yb}| < {} @ {xb_e}",
                event.value_precision().abs()
            );
            return self.at(xb_e);
        }

        debug!("{event}: eval@{xa_e} = {ya}\t eval@{xb_e} = {yb}");

        // The Brent solver, from the roots crate (sadly could not directly integrate it here)
        // Source: https://docs.rs/roots/0.0.5/src/roots/numerical/brent.rs.html#57-131

        let (mut xc, mut yc, mut xd) = (xa, ya, xa);
        let mut flag = true;

        for _ in 0..max_iter {
            if ya.abs() < event.value_precision().abs() {
                return self.at(xa_e + xa * Unit::Second);
            }
            if yb.abs() < event.value_precision().abs() {
                return self.at(xa_e + xb * Unit::Second);
            }
            if has_converged(xa, xb) {
                // The event isn't in the bracket
                return Err(NyxError::from(TrajError::EventNotFound {
                    start,
                    end,
                    event: format!("{event}"),
                }));
            }
            let mut s = if (ya - yc).abs() > f64::EPSILON && (yb - yc).abs() > f64::EPSILON {
                xa * yb * yc / ((ya - yb) * (ya - yc))
                    + xb * ya * yc / ((yb - ya) * (yb - yc))
                    + xc * ya * yb / ((yc - ya) * (yc - yb))
            } else {
                xb - yb * (xb - xa) / (yb - ya)
            };
            let cond1 = (s - xb) * (s - (3.0 * xa + xb) / 4.0) > 0.0;
            let cond2 = flag && (s - xb).abs() >= (xb - xc).abs() / 2.0;
            let cond3 = !flag && (s - xb).abs() >= (xc - xd).abs() / 2.0;
            let cond4 = flag && has_converged(xb, xc);
            let cond5 = !flag && has_converged(xc, xd);
            if cond1 || cond2 || cond3 || cond4 || cond5 {
                s = (xa + xb) / 2.0;
                flag = true;
            } else {
                flag = false;
            }
            let next_try = self.at(xa_e + s * Unit::Second)?;
            let ys = event.eval(&next_try);
            xd = xc;
            xc = xb;
            yc = yb;
            if ya * ys < 0.0 {
                // Root bracketed between a and s
                let next_try = self.at(xa_e + xa * Unit::Second)?;
                let ya_p = event.eval(&next_try);
                let (_a, _ya, _b, _yb) = arrange(xa, ya_p, s, ys);
                {
                    xa = _a;
                    ya = _ya;
                    xb = _b;
                    yb = _yb;
                }
            } else {
                // Root bracketed between s and b
                let next_try = self.at(xa_e + xb * Unit::Second)?;
                let yb_p = event.eval(&next_try);
                let (_a, _ya, _b, _yb) = arrange(s, ys, xb, yb_p);
                {
                    xa = _a;
                    ya = _ya;
                    xb = _b;
                    yb = _yb;
                }
            }
        }
        Err(NyxError::MaxIterReached(format!(
            "Brent solver failed after {max_iter} iterations",
        )))
    }

    /// Find all of the states where the event happens (usually, and with caveats).
    ///
    /// # Limitations
    /// This method uses a Brent solver. If the function that defines the event is not unimodal, the event finder may not converge correctly.
    ///
    /// # Heuristic detail
    /// The initial search step is 1% of the duration of the trajectory duration.
    /// For example, if the trajectory is 100 days long, then we split the trajectory into 100 chunks of 1 day and see whether
    /// the event is in there. If the event happens twice or more times within 1% of the trajectory duration, only the _one_ of
    /// such events will be found.
    ///
    /// If this heuristic fails to find any such events, then `find_minmax` is called on the event with a time precision of `Unit::Second`.
    /// Then we search only within the min and max bounds of the provided event.
    #[allow(clippy::identity_op)]
    pub fn find_all<E>(&self, event: &E) -> Result<Vec<S>, NyxError>
    where
        E: EventEvaluator<S>,
    {
        let start_epoch = self.first().epoch();
        let end_epoch = self.last().epoch();
        if start_epoch == end_epoch {
            return Err(NyxError::from(TrajError::EventNotFound {
                start: start_epoch,
                end: end_epoch,
                event: format!("{event}"),
            }));
        }
        let heuristic = (end_epoch - start_epoch) / 100;
        info!("Searching for {event} with initial heuristic of {heuristic}",);

        let (sender, receiver) = channel();

        let epochs: Vec<Epoch> = TimeSeries::inclusive(start_epoch, end_epoch, heuristic).collect();
        epochs.into_par_iter().for_each_with(sender, |s, epoch| {
            if let Ok(event_state) = self.find_bracketed(epoch, epoch + heuristic, event) {
                s.send(event_state).unwrap()
            };
        });

        let mut states: Vec<_> = receiver.iter().collect();

        if states.is_empty() {
            warn!("Heuristic failed to find any {event} event, using slower approach");
            // Crap, we didn't find the event.
            // Let's find the min and max of this event throughout the trajectory, and search around there.
            match self.find_minmax(event, Unit::Second) {
                Ok((min_event, max_event)) => {
                    let lower_min_epoch =
                        if min_event.epoch() - 1 * Unit::Millisecond < self.first().epoch() {
                            self.first().epoch()
                        } else {
                            min_event.epoch() - 1 * Unit::Millisecond
                        };

                    let lower_max_epoch =
                        if min_event.epoch() + 1 * Unit::Millisecond > self.last().epoch() {
                            self.last().epoch()
                        } else {
                            min_event.epoch() + 1 * Unit::Millisecond
                        };

                    let upper_min_epoch =
                        if max_event.epoch() - 1 * Unit::Millisecond < self.first().epoch() {
                            self.first().epoch()
                        } else {
                            max_event.epoch() - 1 * Unit::Millisecond
                        };

                    let upper_max_epoch =
                        if max_event.epoch() + 1 * Unit::Millisecond > self.last().epoch() {
                            self.last().epoch()
                        } else {
                            max_event.epoch() + 1 * Unit::Millisecond
                        };

                    // Search around the min event
                    if let Ok(event_state) =
                        self.find_bracketed(lower_min_epoch, lower_max_epoch, event)
                    {
                        states.push(event_state);
                    };

                    // Search around the max event
                    if let Ok(event_state) =
                        self.find_bracketed(upper_min_epoch, upper_max_epoch, event)
                    {
                        states.push(event_state);
                    };

                    // If there still isn't any match, report that the event was not found
                    if states.is_empty() {
                        return Err(NyxError::from(TrajError::EventNotFound {
                            start: start_epoch,
                            end: end_epoch,
                            event: format!("{event}"),
                        }));
                    }
                }
                Err(_) => {
                    return Err(NyxError::from(TrajError::EventNotFound {
                        start: start_epoch,
                        end: end_epoch,
                        event: format!("{event}"),
                    }));
                }
            };
        }
        // Remove duplicates and reorder
        states.sort_by(|s1, s2| s1.epoch().partial_cmp(&s2.epoch()).unwrap());
        states.dedup();
        for (cnt, event_state) in states.iter().enumerate() {
            info!(
                "{event} #{}: {} for {event_state}",
                cnt + 1,
                event.eval_string(event_state)
            );
        }
        Ok(states)
    }

    /// Find the minimum and maximum of the provided event through the trajectory
    #[allow(clippy::identity_op)]
    pub fn find_minmax<E>(&self, event: &E, precision: Unit) -> Result<(S, S), NyxError>
    where
        E: EventEvaluator<S>,
    {
        let step: Duration = 1 * precision;
        let mut min_val = std::f64::INFINITY;
        let mut max_val = std::f64::NEG_INFINITY;
        let mut min_state = S::zeros();
        let mut max_state = S::zeros();

        let (sender, receiver) = channel();

        let epochs: Vec<Epoch> =
            TimeSeries::inclusive(self.first().epoch(), self.last().epoch(), step).collect();

        epochs.into_par_iter().for_each_with(sender, |s, epoch| {
            let state = self.at(epoch).unwrap();
            let this_eval = event.eval(&state);
            s.send((this_eval, state)).unwrap();
        });

        let evald_states: Vec<_> = receiver.iter().collect();
        for (this_eval, state) in evald_states {
            if this_eval < min_val {
                min_val = this_eval;
                min_state = state;
            }
            if this_eval > max_val {
                max_val = this_eval;
                max_state = state;
            }
        }

        Ok((min_state, max_state))
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
        let mut record = Vec::new();

        // Build the states iterator

        if cfg.start_epoch.is_some() || cfg.end_epoch.is_some() || cfg.step.is_some() {
            let start = if let Some(start) = cfg.start_epoch {
                start
            } else {
                self.first().epoch()
            };

            let end = if let Some(end) = cfg.end_epoch {
                end
            } else {
                self.last().epoch()
            };

            let step = if let Some(step) = cfg.step {
                step
            } else {
                1.minutes()
            };

            // Build all of the records
            let mut data = Vec::new();
            for s in self.every_between(step, start, end) {
                data.push(format!("{}", s.epoch()));
            }
            record.push(Arc::new(StringArray::from(data)) as ArrayRef);

            // TDB epoch
            let mut data = Vec::new();
            for s in self.every_between(step, start, end) {
                data.push(format!("{:x}", s.epoch()));
            }
            record.push(Arc::new(StringArray::from(data)) as ArrayRef);

            // TAI Epoch seconds
            let mut data = Vec::new();
            for s in self.every_between(step, start, end) {
                data.push(s.epoch().to_tai_seconds());
            }
            record.push(Arc::new(Float64Array::from(data)) as ArrayRef);

            // Add all of the fields
            // This is super ugly, but I can't seem to convert the TrajIterator into an `Iter<S>`
            for field in fields {
                if field == StateParameter::GuidanceMode {
                    // This is the only string field
                    record.push(Arc::new(StringArray::from({
                        let mut data = Vec::new();
                        for s in self.every_between(step, start, end) {
                            let mode = GuidanceMode::from(s.value(field).unwrap());

                            data.push(format!("{mode:?}"));
                        }
                        data
                    })) as ArrayRef);
                } else {
                    record.push(Arc::new(Float64Array::from({
                        let mut data = Vec::new();
                        for s in self.every_between(step, start, end) {
                            data.push(s.value(field).unwrap());
                        }
                        data
                    })) as ArrayRef);
                }
            }
            // Add all of the evaluated events
            if let Some(events) = events {
                for event in events {
                    record.push(Arc::new(Float64Array::from({
                        let mut data = Vec::new();
                        for s in self.every_between(step, start, end) {
                            data.push(event.eval(&s));
                        }
                        data
                    })) as ArrayRef);
                }
            }
        } else {
            // Build all of the records
            record.push(Arc::new(StringArray::from(
                self.states
                    .iter()
                    .map(|s| format!("{}", s.epoch()))
                    .collect::<Vec<String>>(),
            )) as ArrayRef);

            // TDB epoch
            record.push(Arc::new(StringArray::from(
                self.states
                    .iter()
                    .map(|s| format!("{:x}", s.epoch()))
                    .collect::<Vec<String>>(),
            )) as ArrayRef);

            // TDB Epoch seconds
            record.push(Arc::new(Float64Array::from(
                self.states
                    .iter()
                    .map(|s| s.epoch().to_tai_seconds())
                    .collect::<Vec<f64>>(),
            )) as ArrayRef);

            // Add all of the fields
            for field in fields {
                if field == StateParameter::GuidanceMode {
                    record.push(Arc::new(StringArray::from(
                        self.states
                            .iter()
                            .map(|s| format!("{:?}", GuidanceMode::from(s.value(field).unwrap())))
                            .collect::<Vec<String>>(),
                    )) as ArrayRef);
                } else {
                    record.push(Arc::new(Float64Array::from(
                        self.states
                            .iter()
                            .map(|s| s.value(field).unwrap())
                            .collect::<Vec<f64>>(),
                    )) as ArrayRef);
                }
            }

            // Add all of the evaluated events
            if let Some(events) = events {
                for event in events {
                    record.push(Arc::new(Float64Array::from({
                        self.states
                            .iter()
                            .map(|s| event.eval(s))
                            .collect::<Vec<f64>>()
                    })) as ArrayRef);
                }
            }
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

        let mut path_buf = path.as_ref().to_path_buf();

        if cfg.timestamp {
            if let Some(file_name) = path_buf.file_name() {
                if let Some(file_name_str) = file_name.to_str() {
                    if let Some(extension) = path_buf.extension() {
                        let stamp = Formatter::new(
                            Epoch::now().unwrap(),
                            Format::from_str("%Y-%m-%dT%H-%M-%S").unwrap(),
                        );
                        let new_file_name =
                            format!("{file_name_str}-{stamp}.{}", extension.to_str().unwrap());
                        path_buf.set_file_name(new_file_name);
                    }
                }
            }
        } else {
        };

        let file = File::create(&path_buf)?;
        let mut writer = ArrowWriter::try_new(file, schema.clone(), props).unwrap();

        let batch = RecordBatch::try_new(schema, record)?;
        writer.write(&batch)?;
        writer.close()?;

        // Return the path this was written to
        Ok(path_buf)
    }
}

impl<S: InterpState> ops::Add for Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    type Output = Traj<S>;

    /// Add one trajectory to another. If they do not overlap to within 10ms, a warning will be printed.
    fn add(self, other: Traj<S>) -> Self::Output {
        self + &other
    }
}

impl<S: InterpState> ops::Add<&Traj<S>> for Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    type Output = Traj<S>;

    /// Add one trajectory to another. If they do not overlap to within 10ms, a warning will be printed.
    fn add(self, other: &Traj<S>) -> Self::Output {
        let (first, second) = if self.first().epoch() < other.first().epoch() {
            (&self, other)
        } else {
            (other, &self)
        };

        if first.last().epoch() < second.first().epoch() {
            let gap = second.first().epoch() - first.last().epoch();
            warn!(
                "Resulting merged trajectory will have a time-gap of {} starting at {}",
                gap,
                first.last().epoch()
            );
        }

        let mut me = self.clone();
        // Now start adding the other segments while correcting the index
        for state in &second.states {
            me.states.push(*state);
        }
        me.finalize();
        me
    }
}

impl<S: InterpState> ops::AddAssign for Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}

impl<S: InterpState> ops::AddAssign<&Traj<S>> for Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    fn add_assign(&mut self, rhs: &Self) {
        *self = self.clone() + rhs;
    }
}

impl Traj<Orbit> {
    /// Allows converting the source trajectory into the (almost) equivalent trajectory in another frame.
    /// This simply converts each state into the other frame and may lead to aliasing due to the Nyquist–Shannon sampling theorem.
    #[allow(clippy::map_clone)]
    pub fn to_frame(&self, new_frame: Frame, cosm: Arc<Cosm>) -> Result<Self, NyxError> {
        if self.states.is_empty() {
            return Err(NyxError::Trajectory(TrajError::CreationError(
                "No trajectory to convert".to_string(),
            )));
        }
        let start_instant = Instant::now();
        let mut traj = Self::new();
        for state in &self.states {
            traj.states.push(cosm.frame_chg(state, new_frame));
        }
        traj.finalize();

        info!(
            "Converted trajectory from {} to {} in {} ms: {traj}",
            self.first().frame,
            new_frame,
            (Instant::now() - start_instant).as_millis()
        );
        Ok(traj)
    }

    /// Exports this trajectory to the provided filename in parquet format with only the epoch, the geodetic latitude, longitude, and height at one state per minute.
    /// Must provide a body fixed frame to correctly compute the latitude and longitude.
    #[allow(clippy::identity_op)]
    pub fn to_groundtrack_parquet<P: AsRef<Path>>(
        &self,
        path: P,
        body_fixed_frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Result<(), Box<dyn Error>> {
        let traj = self.to_frame(body_fixed_frame, cosm)?;

        traj.to_parquet_with_cfg(
            path,
            ExportCfg {
                fields: Some(vec![
                    StateParameter::GeodeticLatitude,
                    StateParameter::GeodeticLongitude,
                    StateParameter::GeodeticHeight,
                ]),
                step: Some(1 * Unit::Minute),
                ..Default::default()
            },
        )?;

        Ok(())
    }
}

impl Traj<Spacecraft> {
    /// Allows converting the source trajectory into the (almost) equivalent trajectory in another frame
    #[allow(clippy::map_clone)]
    pub fn to_frame(&self, new_frame: Frame, cosm: Arc<Cosm>) -> Result<Self, NyxError> {
        if self.states.is_empty() {
            return Err(NyxError::Trajectory(TrajError::CreationError(
                "No trajectory to convert".to_string(),
            )));
        }
        let start_instant = Instant::now();
        let mut traj = Self::new();
        for state in &self.states {
            traj.states
                .push(state.with_orbit(cosm.frame_chg(&state.orbit, new_frame)));
        }
        traj.finalize();

        info!(
            "Converted trajectory from {} to {} in {} ms: {traj}",
            self.first().orbit.frame,
            new_frame,
            (Instant::now() - start_instant).as_millis()
        );
        Ok(traj)
    }

    /// A shortcut to `to_parquet_with_csv`
    pub fn to_parquet_with_step<P: AsRef<Path>>(
        &self,
        path: P,
        step: Duration,
    ) -> Result<(), Box<dyn Error>> {
        self.to_parquet_with_cfg(
            path,
            ExportCfg {
                step: Some(step),
                ..Default::default()
            },
        )?;

        Ok(())
    }

    /// Exports this trajectory to the provided filename in parquet format with only the epoch, the geodetic latitude, longitude, and height at one state per minute.
    /// Must provide a body fixed frame to correctly compute the latitude and longitude.
    #[allow(clippy::identity_op)]
    pub fn to_groundtrack_parquet<P: AsRef<Path>>(
        &self,
        path: P,
        body_fixed_frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Result<(), Box<dyn Error>> {
        let traj = self.to_frame(body_fixed_frame, cosm)?;

        traj.to_parquet_with_cfg(
            path,
            ExportCfg {
                fields: Some(vec![
                    StateParameter::GeodeticLatitude,
                    StateParameter::GeodeticLongitude,
                    StateParameter::GeodeticHeight,
                ]),
                step: Some(1 * Unit::Minute),
                ..Default::default()
            },
        )?;

        Ok(())
    }
}

impl<S: InterpState> fmt::Display for Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let dur = self.last().epoch() - self.first().epoch();
        write!(
            f,
            "Trajectory from {} to {} ({}, or {:.3} s) [{} states]",
            self.first().epoch(),
            self.last().epoch(),
            dur,
            dur.to_seconds(),
            self.states.len()
        )
    }
}

impl<S: InterpState> fmt::Debug for Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{self}",)
    }
}

impl<S: InterpState> Default for Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    fn default() -> Self {
        Self::new()
    }
}
