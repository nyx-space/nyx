/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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

extern crate bacon_sci;
extern crate crossbeam;
extern crate rayon;

use self::crossbeam::thread;
use self::rayon::prelude::*;
use super::spline::{Spline, SPLINE_DEGREE};
use super::traj_it::TrajIterator;
use super::InterpState;
use crate::cosmic::{Cosm, Frame, Orbit, Spacecraft};
use crate::errors::NyxError;
use crate::io::formatter::StateFormatter;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName};
use crate::md::{events::EventEvaluator, MdHdlr, OrbitStateOutput};
use crate::polyfit::hermite;
use crate::time::{Duration, Epoch, TimeSeries, TimeUnit};
use crate::utils::normalize;
use crate::State;
use std::collections::BTreeMap;
use std::fmt;
use std::iter::Iterator;
use std::sync::mpsc::{channel, Receiver};
use std::sync::Arc;
use std::time::Duration as StdDur;
use std::time::Instant;

/// Store a trajectory of any State.
#[derive(Clone)]
pub struct Traj<S: InterpState>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    /// Splines are organized as a binary tree map whose index is the
    /// number of seconds since the start of this ephemeris rounded down
    pub segments: BTreeMap<i32, Spline<S>>,
    /// Timeout is used to stop listening to new state
    pub timeout_ms: u64,
    pub(crate) start_state: S,
    pub(crate) max_offset: i32,
    /// Store the items per segment to convert to another frame.
    pub(crate) items_per_segments: usize,
}

impl<S: InterpState> Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    /// Creates a new trajectory with the provided starting state (used as a template) and a receiving channel.
    /// The trajectories are always generated on a separate thread.
    pub fn new(state: S, rx: Receiver<S>) -> Result<Self, NyxError> {
        // Bug? With a spacecraft, we need more interpolation windows than just an orbit.
        // I've spent 12h trying to understand why, but I can't, so screw it for it.
        Self::new_bucket_as(state, if S::VecLength::dim() == 90 { 6 } else { 32 }, rx)
    }

    /// Creates a new trajectory but specifies the number of items per segment
    pub fn new_bucket_as(
        state: S,
        items_per_segments: usize,
        rx: Receiver<S>,
    ) -> Result<Self, NyxError> {
        thread::scope(|s| {
            // Initialize the interpolator
            let mut me = Self {
                segments: BTreeMap::new(),
                start_state: state,
                timeout_ms: 100,
                max_offset: 0,
                items_per_segments,
            };

            let mut children = vec![];
            let mut window_states: Vec<S> = Vec::with_capacity(items_per_segments);
            // Push the initial state
            window_states.push(state);

            // Note that we're using the typical map+reduce pattern
            // Start receiving states on a blocking call (map)
            while let Ok(state) = rx.recv_timeout(StdDur::from_millis(me.timeout_ms)) {
                if window_states.len() == items_per_segments {
                    let this_wdn = window_states.clone();
                    children.push(
                        s.spawn(move |_| -> Result<Spline<S>, NyxError> { interpolate(this_wdn) }),
                    );
                    // Copy the last state as the first state of the next window
                    let last_wdn_state = window_states[items_per_segments - 1];
                    window_states.clear();
                    window_states.push(last_wdn_state);
                }
                window_states.push(state);
            }
            // And interpolate the remaining states too, even if the buffer is not full!
            children.push(
                s.spawn(move |_| -> Result<Spline<S>, NyxError> { interpolate(window_states) }),
            );

            // Reduce
            for child in children {
                // collect each child thread's return-value
                let segment = child.join().unwrap()?;
                me.append_spline(segment);
            }

            Ok(me)
        })
        .unwrap()
    }

    pub(crate) fn append_spline(&mut self, segment: Spline<S>) {
        // Compute the number of seconds since start of trajectory
        let offset_s = ((segment.start_epoch - self.start_state.epoch())
            .in_seconds()
            .floor()) as i32;
        self.segments.insert(offset_s, segment);
        if offset_s > self.max_offset {
            self.max_offset = offset_s;
        }
    }

    /// Evaluate the trajectory at this specific epoch.
    pub fn at(&self, epoch: Epoch) -> Result<S, NyxError> {
        let offset_s = ((epoch - self.start_state.epoch()).in_seconds().floor()) as i32;

        // Retrieve that segment
        match self.segments.range(..=offset_s).rev().next() {
            None => {
                // Let's see if this corresponds to the max offset value
                let last_item = self.segments[&self.max_offset].end_state;
                if last_item.epoch() == epoch {
                    Ok(last_item)
                } else {
                    Err(NyxError::NoInterpolationData(format!("{}", epoch)))
                }
            }
            Some((_, segment)) => segment.evaluate(self.start_state, epoch),
        }
    }

    /// Returns the first state in this ephemeris
    pub fn first(&self) -> S {
        self.start_state
    }

    /// Returns the last state in this ephemeris
    pub fn last(&self) -> S {
        self.segments[&self.max_offset].end_state
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

    /// Find the exact state where the request event happens. The event function is expected to be monotone in the provided interval.
    #[allow(clippy::identity_op)]
    pub fn find_bracketed<E>(&self, start: Epoch, end: Epoch, event: &E) -> Result<S, NyxError>
    where
        E: EventEvaluator<S>,
    {
        use std::f64::EPSILON;
        let max_iter = 50;

        // Helper lambdas, for f64s only
        let has_converged =
            |x1: f64, x2: f64| (x1 - x2).abs() <= event.epoch_precision().in_seconds();
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
        let mut xb = (xb_e - xa_e).in_seconds();
        // Evaluate the event at both bounds
        let mut ya = event.eval(&self.at(xa_e)?);
        let mut yb = event.eval(&self.at(xb_e)?);

        // Check if we're already at the root
        if ya.abs() <= event.value_precision().abs() {
            return self.at(xa_e);
        } else if yb.abs() <= event.value_precision().abs() {
            return self.at(xb_e);
        }
        // The Brent solver, from the roots crate (sadly could not directly integrate it here)
        // Source: https://docs.rs/roots/0.0.5/src/roots/numerical/brent.rs.html#57-131

        let (mut xc, mut yc, mut xd) = (xa, ya, xa);
        let mut flag = true;

        for _ in 0..max_iter {
            if ya.abs() < event.value_precision().abs() {
                return self.at(xa_e + xa * TimeUnit::Second);
            }
            if yb.abs() < event.value_precision().abs() {
                return self.at(xa_e + xb * TimeUnit::Second);
            }
            if has_converged(xa, xb) {
                // The event isn't in the bracket
                return Err(NyxError::EventNotInEpochBraket(
                    start.to_string(),
                    end.to_string(),
                ));
            }
            let mut s = if (ya - yc).abs() > EPSILON && (yb - yc).abs() > EPSILON {
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
            let next_try = self.at(xa_e + s * TimeUnit::Second)?;
            let ys = event.eval(&next_try);
            xd = xc;
            xc = xb;
            yc = yb;
            if ya * ys < 0.0 {
                // Root bracketed between a and s
                let next_try = self.at(xa_e + xa * TimeUnit::Second)?;
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
                let next_try = self.at(xa_e + xb * TimeUnit::Second)?;
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
            "Brent solver failed after {} iterations",
            max_iter
        )))
    }

    /// Find (usually) all of the states where the event happens.
    /// WARNING: The initial search step is 1% of the duration of the trajectory duration!
    /// For example, if the trajectory is 100 days long, then we split the trajectory into 100 chunks of 1 day and see whether
    /// the event is in there. If the event happens twice or more times within 1% of the trajectory duration, only the _one_ of
    /// such events will be found.
    #[allow(clippy::identity_op)]
    pub fn find_all<E>(&self, event: &E) -> Result<Vec<S>, NyxError>
    where
        E: EventEvaluator<S>,
    {
        let start_epoch = self.first().epoch();
        let end_epoch = self.last().epoch();
        if start_epoch == end_epoch {
            return Err(NyxError::NoInterpolationData(format!(
                "Trajector of {} segments has identical start and end dates -- cannot search it",
                self.segments.len()
            )));
        }
        let heuristic = (end_epoch - start_epoch) / 100;
        info!(
            "Searching for {} with initial heuristic of {}",
            event, heuristic
        );

        let (sender, receiver) = channel();

        let epochs: Vec<Epoch> = TimeSeries::inclusive(start_epoch, end_epoch, heuristic).collect();
        epochs.into_par_iter().for_each_with(sender, |s, epoch| {
            if let Ok(event_state) = self.find_bracketed(epoch, epoch + heuristic, event) {
                s.send(event_state).unwrap()
            };
        });

        let mut states: Vec<_> = receiver.iter().collect();

        if states.is_empty() {
            warn!(
                "Heuristic failed to find any {} event, using slower approach",
                event
            );
            // Crap, we didn't find the event.
            // Let's find the min and max of this event throughout the trajectory, and search around there.
            match self.find_minmax(event, TimeUnit::Second) {
                Ok((min_event, max_event)) => {
                    let lower_min_epoch =
                        if min_event.epoch() - 1 * TimeUnit::Millisecond < self.first().epoch() {
                            self.first().epoch()
                        } else {
                            min_event.epoch() - 1 * TimeUnit::Millisecond
                        };

                    let lower_max_epoch =
                        if min_event.epoch() + 1 * TimeUnit::Millisecond > self.last().epoch() {
                            self.last().epoch()
                        } else {
                            min_event.epoch() + 1 * TimeUnit::Millisecond
                        };

                    let upper_min_epoch =
                        if max_event.epoch() - 1 * TimeUnit::Millisecond < self.first().epoch() {
                            self.first().epoch()
                        } else {
                            max_event.epoch() - 1 * TimeUnit::Millisecond
                        };

                    let upper_max_epoch =
                        if max_event.epoch() + 1 * TimeUnit::Millisecond > self.last().epoch() {
                            self.last().epoch()
                        } else {
                            max_event.epoch() + 1 * TimeUnit::Millisecond
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
                        return Err(NyxError::EventNotInEpochBraket(
                            start_epoch.to_string(),
                            end_epoch.to_string(),
                        ));
                    }
                }
                Err(_) => {
                    return Err(NyxError::EventNotInEpochBraket(
                        start_epoch.to_string(),
                        end_epoch.to_string(),
                    ))
                }
            };
        }
        // Remove duplicates and reorder
        states.sort_by(|s1, s2| s1.epoch().partial_cmp(&s2.epoch()).unwrap());
        states.dedup();
        for (cnt, event_state) in states.iter().enumerate() {
            info!("{} #{}: {}", event, cnt + 1, event_state);
        }
        Ok(states)
    }

    /// Find the minimum and maximum of the provided event through the trajectory
    #[allow(clippy::identity_op)]
    pub fn find_minmax<E>(&self, event: &E, precision: TimeUnit) -> Result<(S, S), NyxError>
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
}

impl Traj<Orbit> {
    /// Allows converting the source trajectory into the (almost) equivalent trajectory in another frame
    /// This is super slow.
    pub fn to_frame(&self, new_frame: Frame, cosm: Arc<Cosm>) -> Result<Self, NyxError> {
        let start_instant = Instant::now();
        let start_state = cosm.frame_chg(&self.first(), new_frame);

        let rx = {
            // Channels that have the states in a bucket of the correct length
            let (tx, rx) = channel();

            // Nyquist–Shannon sampling theorem
            let sample_rate =
                1.0 / (((self.items_per_segments * 2 + 1) * self.segments.len()) as f64);
            let step = sample_rate * (self.last().epoch() - self.first().epoch());

            /* *** */
            /* Map: bucket the states and send on a channel */
            /* *** */

            let items_per_segments = <Orbit as InterpState>::INTERPOLATION_SAMPLES;
            let mut window_states = Vec::with_capacity(items_per_segments);
            // Push the initial state
            window_states.push(start_state);

            // Note that we're using the typical map+reduce pattern
            for original_state in self.every(step) {
                let state = cosm.frame_chg(&original_state, new_frame);
                if window_states.len() == items_per_segments {
                    let this_wdn = window_states.clone();
                    tx.send(this_wdn)
                        .map_err(|_| NyxError::TrajectoryCreationError)?;
                    // Copy the last state as the first state of the next window
                    let last_wdn_state = window_states[items_per_segments - 1];
                    window_states.clear();
                    window_states.push(last_wdn_state);
                }
                window_states.push(state);
            }
            // if window_states.len() < 3 {
            //     // Add more states or the interpolation will fail
            //     window_states.clear();
            //     for original_state in
            //         self.every_between(step, self.last().epoch() - 4 * step, self.last().epoch())
            //     {
            //         let state = cosm.frame_chg(&original_state, new_frame);
            //         window_states.push(state)
            //     }
            //     dbg!(window_states.len());
            // }
            // And interpolate the remaining states too, even if the buffer is not full!
            tx.send(window_states)
                .map_err(|_| NyxError::TrajectoryCreationError)?;

            // Return the rx channel for these buckets
            rx
        };

        /* *** */
        /* Reduce: Build an interpolation of each of the segments */
        /* *** */
        let splines: Vec<_> = rx
            .into_iter()
            .par_bridge()
            .map(|window_states| {
                let rslt = interpolate(window_states.clone());
                if rslt.is_err() {
                    for s in window_states {
                        println!("{}", s);
                    }
                }
                rslt
            })
            .collect();

        // Finally, build the whole trajectory
        let mut traj = Traj {
            segments: BTreeMap::new(),
            start_state,
            timeout_ms: 100,
            max_offset: 0,
            items_per_segments: 8,
        };

        for maybe_spline in splines {
            let spline = maybe_spline?;
            traj.append_spline(spline);
        }

        info!(
            "Converted trajectory from {} to {} in {} ms",
            self.first().frame,
            new_frame,
            (Instant::now() - start_instant).as_millis()
        );

        Ok(traj)
    }

    /// Exports this trajectory to the provided filename in CSV format with the default headers and the provided step
    pub fn to_csv_with_step(
        &self,
        filename: &str,
        step: Duration,
        cosm: Arc<Cosm>,
    ) -> Result<(), NyxError> {
        let fmtr = StateFormatter::default(filename.to_string(), cosm);
        let mut out = OrbitStateOutput::new(fmtr)?;
        for state in self.every(step) {
            out.handle(&state);
        }
        Ok(())
    }

    /// Exports this trajectory to the provided filename in CSV format with the default headers and the provided step
    pub fn to_csv_between_with_step(
        &self,
        filename: &str,
        start: Option<Epoch>,
        end: Option<Epoch>,
        step: Duration,
        cosm: Arc<Cosm>,
    ) -> Result<(), NyxError> {
        let fmtr = StateFormatter::default(filename.to_string(), cosm);
        let mut out = OrbitStateOutput::new(fmtr)?;
        let start = match start {
            Some(s) => s,
            None => self.first().epoch(),
        };
        let end = match end {
            Some(e) => e,
            None => self.last().epoch(),
        };
        for state in self.every_between(step, start, end) {
            out.handle(&state);
        }
        Ok(())
    }

    /// Exports this trajectory to the provided filename in CSV format with the default headers, one state per minute
    #[allow(clippy::identity_op)]
    pub fn to_csv(&self, filename: &str, cosm: Arc<Cosm>) -> Result<(), NyxError> {
        self.to_csv_with_step(filename, 1 * TimeUnit::Minute, cosm)
    }

    /// Exports this trajectory to the provided filename in CSV format with the default headers, one state per minute
    #[allow(clippy::identity_op)]
    pub fn to_csv_between(
        &self,
        filename: &str,
        start: Option<Epoch>,
        end: Option<Epoch>,
        cosm: Arc<Cosm>,
    ) -> Result<(), NyxError> {
        self.to_csv_between_with_step(filename, start, end, 1 * TimeUnit::Minute, cosm)
    }

    /// Exports this trajectory to the provided filename in CSV format with only the epoch, the geodetic latitude, longitude, and height at one state per minute.
    /// Must provide a body fixed frame to correctly compute the latitude and longitude.
    #[allow(clippy::identity_op)]
    pub fn to_groundtrack_csv(
        &self,
        filename: &str,
        body_fixed_frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Result<(), NyxError> {
        let fmtr = StateFormatter::from_headers(
            vec![
                "epoch",
                "geodetic_latitude",
                "geodetic_longitude",
                "geodetic_height",
            ],
            filename.to_string(),
            cosm.clone(),
        )?;
        let mut out = OrbitStateOutput::new(fmtr)?;
        for state in self
            .to_frame(body_fixed_frame, cosm)?
            .every(1 * TimeUnit::Minute)
        {
            out.handle(&state);
        }
        Ok(())
    }

    /// Exports this trajectory to the provided filename in CSV format with the provided headers and the provided step
    pub fn to_csv_custom(
        &self,
        filename: &str,
        headers: Vec<&str>,
        step: Duration,
        cosm: Arc<Cosm>,
    ) -> Result<(), NyxError> {
        let fmtr = StateFormatter::from_headers(headers, filename.to_string(), cosm)?;
        let mut out = OrbitStateOutput::new(fmtr)?;
        for state in self.every(step) {
            out.handle(&state);
        }
        Ok(())
    }
}

impl Traj<Spacecraft> {
    /// Allows converting the source trajectory into the (almost) equivalent trajectory in another frame
    /// This is super slow.
    pub fn to_frame(&self, new_frame: Frame, cosm: Arc<Cosm>) -> Result<Self, NyxError> {
        let start_instant = Instant::now();

        let mut start_state = self.first();
        let start_orbit = cosm.frame_chg(&start_state.orbit, new_frame);
        // Update the orbit component and the frame
        start_state.orbit.x = start_orbit.x;
        start_state.orbit.y = start_orbit.y;
        start_state.orbit.z = start_orbit.z;
        start_state.orbit.vx = start_orbit.vx;
        start_state.orbit.vy = start_orbit.vy;
        start_state.orbit.vz = start_orbit.vz;
        start_state.orbit.frame = start_orbit.frame;

        let rx = {
            // Channels that have the states in a bucket of the correct length
            let (tx_bucket, rx_bucket) = channel();

            // Channels that have a single state for the propagator
            let (tx, rx) = channel();
            // Nyquist–Shannon sampling theorem
            let sample_rate =
                1.0 / (((self.items_per_segments * 2 + 1) * self.segments.len()) as f64);
            let step = sample_rate * (self.last().epoch() - self.first().epoch());

            for state in self.every(step) {
                let converted_orbit = cosm.frame_chg(&state.orbit, new_frame);
                let mut converted_state = state;
                // Update the orbit component and the frame
                converted_state.orbit.x = converted_orbit.x;
                converted_state.orbit.y = converted_orbit.y;
                converted_state.orbit.z = converted_orbit.z;
                converted_state.orbit.vx = converted_orbit.vx;
                converted_state.orbit.vy = converted_orbit.vy;
                converted_state.orbit.vz = converted_orbit.vz;
                converted_state.orbit.frame = converted_orbit.frame;

                tx.send(converted_state).unwrap();
            }

            /* *** */
            /* Map: bucket the states and send on a channel */
            /* *** */

            let items_per_segments = <Orbit as InterpState>::INTERPOLATION_SAMPLES;
            let mut window_states = Vec::with_capacity(items_per_segments);
            // Push the initial state
            window_states.push(start_state);

            // Note that we're using the typical map+reduce pattern
            // Start receiving states on a blocking call (map)
            while let Ok(state) = rx.recv_timeout(StdDur::from_millis(100)) {
                if window_states.len() == items_per_segments {
                    let this_wdn = window_states.clone();
                    tx_bucket
                        .send(this_wdn)
                        .map_err(|_| NyxError::TrajectoryCreationError)?;
                    // Copy the last state as the first state of the next window
                    let last_wdn_state = window_states[items_per_segments - 1];
                    window_states.clear();
                    window_states.push(last_wdn_state);
                }
                window_states.push(state);
            }
            // And interpolate the remaining states too, even if the buffer is not full!
            tx_bucket
                .send(window_states)
                .map_err(|_| NyxError::TrajectoryCreationError)?;

            // Return the rx channel for these buckets
            rx_bucket
        };

        /* *** */
        /* Reduce: Build an interpolation of each of the segments */
        /* *** */
        let splines: Vec<_> = rx
            .into_iter()
            .par_bridge()
            .map(|window_states| interpolate(window_states))
            .collect();

        // Finally, build the whole trajectory
        let mut traj = Traj {
            segments: BTreeMap::new(),
            start_state,
            timeout_ms: 100,
            max_offset: 0,
            items_per_segments: 8,
        };

        for maybe_spline in splines {
            let spline = maybe_spline?;
            traj.append_spline(spline);
        }

        info!(
            "Converted trajectory from {} to {} in {} ms",
            self.first().orbit.frame,
            new_frame,
            (Instant::now() - start_instant).as_millis()
        );

        Ok(traj)
    }

    /// Exports this trajectory to the provided filename in CSV format with the default headers and the provided step
    pub fn to_csv_with_step(
        &self,
        filename: &str,
        step: Duration,
        cosm: Arc<Cosm>,
    ) -> Result<(), NyxError> {
        let fmtr = StateFormatter::default(filename.to_string(), cosm);
        let mut out = OrbitStateOutput::new(fmtr)?;
        for state in self.every(step) {
            out.handle(&state);
        }
        Ok(())
    }

    /// Exports this trajectory to the provided filename in CSV format with the default headers and the provided step
    pub fn to_csv_between_with_step(
        &self,
        filename: &str,
        start: Option<Epoch>,
        end: Option<Epoch>,
        step: Duration,
        cosm: Arc<Cosm>,
    ) -> Result<(), NyxError> {
        let fmtr = StateFormatter::default(filename.to_string(), cosm);
        let mut out = OrbitStateOutput::new(fmtr)?;
        let start = match start {
            Some(s) => s,
            None => self.first().epoch(),
        };
        let end = match end {
            Some(e) => e,
            None => self.last().epoch(),
        };
        for state in self.every_between(step, start, end) {
            out.handle(&state);
        }
        Ok(())
    }

    /// Exports this trajectory to the provided filename in CSV format with the default headers, one state per minute
    #[allow(clippy::identity_op)]
    pub fn to_csv(&self, filename: &str, cosm: Arc<Cosm>) -> Result<(), NyxError> {
        self.to_csv_with_step(filename, 1 * TimeUnit::Minute, cosm)
    }

    /// Exports this trajectory to the provided filename in CSV format with the default headers, one state per minute
    #[allow(clippy::identity_op)]
    pub fn to_csv_between(
        &self,
        filename: &str,
        start: Option<Epoch>,
        end: Option<Epoch>,
        cosm: Arc<Cosm>,
    ) -> Result<(), NyxError> {
        self.to_csv_between_with_step(filename, start, end, 1 * TimeUnit::Minute, cosm)
    }

    /// Exports this trajectory to the provided filename in CSV format with only the epoch, the geodetic latitude, longitude, and height at one state per minute.
    /// Must provide a body fixed frame to correctly compute the latitude and longitude.
    #[allow(clippy::identity_op)]
    pub fn to_groundtrack_csv(
        &self,
        filename: &str,
        body_fixed_frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Result<(), NyxError> {
        let fmtr = StateFormatter::from_headers(
            vec![
                "epoch",
                "geodetic_latitude",
                "geodetic_longitude",
                "geodetic_height",
            ],
            filename.to_string(),
            cosm.clone(),
        )?;
        let mut out = OrbitStateOutput::new(fmtr)?;
        for state in self
            .to_frame(body_fixed_frame, cosm)?
            .every(1 * TimeUnit::Minute)
        {
            out.handle(&state);
        }
        Ok(())
    }

    /// Exports this trajectory to the provided filename in CSV format with the provided headers and the provided step
    pub fn to_csv_custom(
        &self,
        filename: &str,
        headers: Vec<&str>,
        step: Duration,
        cosm: Arc<Cosm>,
    ) -> Result<(), NyxError> {
        let fmtr = StateFormatter::from_headers(headers, filename.to_string(), cosm)?;
        let mut out = OrbitStateOutput::new(fmtr)?;
        for state in self.every(step) {
            out.handle(&state);
        }
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
            "Trajectory from {} to {} ({}, or {:.3} s) [{} splines]",
            self.first().epoch(),
            self.last().epoch(),
            dur,
            dur.in_seconds(),
            self.segments.len()
        )
    }
}

pub(crate) fn interpolate<S: InterpState>(this_wdn: Vec<S>) -> Result<Spline<S>, NyxError>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    if this_wdn.len() == 0 {
        return Err(NyxError::NoInterpolationData(format!(
            "Cannot interpolate with only {} items",
            this_wdn.len()
        )));
    }
    // Generate interpolation and flush.
    let start_win_epoch = this_wdn[0].epoch();
    let end_win_epoch = this_wdn[this_wdn.len() - 1].epoch();
    let window_duration = end_win_epoch - start_win_epoch;
    let mut ts = Vec::with_capacity(this_wdn.len());
    let mut values = Vec::with_capacity(S::params().len());
    let mut values_dt = Vec::with_capacity(S::params().len());
    let mut polynomials = Vec::with_capacity(S::params().len());

    // Initialize the vector of values and coefficients.
    for _ in 0..S::params().len() {
        values.push(Vec::with_capacity(this_wdn.len()));
        values_dt.push(Vec::with_capacity(this_wdn.len()));
    }
    for state in &this_wdn {
        let t_prime = if this_wdn.len() == 1 {
            1.0
        } else {
            normalize(
                (state.epoch() - start_win_epoch).in_seconds(),
                0.0,
                window_duration.in_seconds(),
            )
        };
        // Deduplicate
        if let Some(latest_t) = ts.last() {
            if *latest_t == t_prime {
                continue;
            }
        }

        ts.push(t_prime);
        for (pos, param) in S::params().iter().enumerate() {
            let (value, value_dt) = state.value_and_deriv(param)?;
            values[pos].push(value);
            values_dt[pos].push(value_dt);
        }
    }

    // Generate the polynomials
    for pos in 0..values.len() {
        let poly = hermite::<SPLINE_DEGREE>(&ts, &values[pos], &values_dt[pos])?;
        polynomials.push(poly);
    }

    Ok(Spline {
        start_epoch: start_win_epoch,
        duration: window_duration,
        polynomials,
        end_state: this_wdn[this_wdn.len() - 1],
    })
}
