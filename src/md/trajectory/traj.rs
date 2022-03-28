/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

extern crate rayon;

use self::rayon::prelude::*;
use super::spline::{Spline, INTERPOLATION_SAMPLES, SPLINE_DEGREE};
use super::traj_it::TrajIterator;
use super::{InterpState, TrajError};
use crate::cosmic::{Cosm, Frame, Orbit, Spacecraft};
use crate::errors::NyxError;
use crate::io::formatter::StateFormatter;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::md::{events::EventEvaluator, MdHdlr, OrbitStateOutput};
use crate::polyfit::{hermite, Polynomial};
use crate::time::{Duration, Epoch, TimeSeries, Unit};
use crate::utils::normalize;
use crate::State;
use std::collections::BTreeMap;
use std::fmt;
use std::iter::Iterator;
use std::ops;
use std::sync::mpsc::channel;
use std::sync::Arc;
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
    pub(crate) start_state: S,
    pub(crate) backward: bool,
}

impl<S: InterpState> Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    pub(crate) fn append_spline(&mut self, segment: Spline<S>) -> Result<(), NyxError> {
        // Compute the number of seconds since start of trajectory
        let offset_s = ((100.0 * (segment.start_epoch - self.start_state.epoch()).in_seconds())
            .floor()) as i32;
        if offset_s.is_negative() {
            if !self.segments.is_empty() && !self.backward {
                return Err(NyxError::from(TrajError::CreationError(format!(
                    "Offset = {} but traj is not backward",
                    offset_s
                ))));
            }
            self.backward = true;
        }
        assert!(
            self.segments.get(&offset_s).is_none(),
            "{offset_s} already exists!"
        );
        self.segments.insert(offset_s, segment);
        Ok(())
    }

    /// Evaluate the trajectory at this specific epoch.
    pub fn at(&self, epoch: Epoch) -> Result<S, NyxError> {
        // Durations are darn precise and converting a -2.6e-23 into an i32 will be -1
        let offset_s = ((100.0 * (epoch - self.start_state.epoch()).in_seconds()).floor()) as i32;

        // Retrieve that segment
        match self.segments.range(..=offset_s).rev().next() {
            None => {
                // Let's see if this corresponds to the max offset value
                let last_item = self.last();
                if last_item.epoch() == epoch {
                    Ok(last_item)
                } else {
                    Err(NyxError::from(TrajError::NoInterpolationData(epoch)))
                }
            }
            Some((_, segment)) => segment.evaluate(self.start_state, epoch),
        }
    }

    /// Returns the first state in this ephemeris
    pub fn first(&self) -> S {
        if self.backward {
            self.segments[self.segments.keys().last().unwrap()].end_state
        } else {
            self.start_state
        }
    }

    /// Returns the last state in this ephemeris
    pub fn last(&self) -> S {
        if self.backward {
            // Note that this trajectory's "start_state" is actually the first state received, so it's chronologically the last state of the first spline.
            let spline = &self.segments[self.segments.keys().next().unwrap()];
            spline
                .evaluate(spline.end_state, spline.start_epoch)
                .unwrap()
        } else {
            self.segments[self.segments.keys().last().unwrap()].end_state
        }
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
                    event: format!("{}", event),
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
            return Err(NyxError::from(TrajError::EventNotFound {
                start: start_epoch,
                end: end_epoch,
                event: format!("{}", event),
            }));
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
                            event: format!("{}", event),
                        }));
                    }
                }
                Err(_) => {
                    return Err(NyxError::from(TrajError::EventNotFound {
                        start: start_epoch,
                        end: end_epoch,
                        event: format!("{}", event),
                    }));
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

        let mut me = Self {
            segments: first.segments.clone(),
            start_state: first.start_state,
            backward: false,
        };
        // Now start adding the other segments while correcting the index
        for spline in second.segments.values() {
            me.append_spline(spline.clone()).unwrap();
        }
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
    /// Allows converting the source trajectory into the (almost) equivalent trajectory in another frame
    #[allow(clippy::map_clone)]
    pub fn to_frame(&self, new_frame: Frame, cosm: Arc<Cosm>) -> Result<Self, NyxError> {
        let start_instant = Instant::now();
        let start_state = cosm.frame_chg(&self.first(), new_frame);

        let rx = {
            // Channels that have the states in a bucket of the correct length
            let (tx, rx) = channel();

            let items_per_segments = INTERPOLATION_SAMPLES;

            // Nyquist–Shannon sampling theorem
            let sample_rate = 1.0 / (((items_per_segments * 2 + 1) * self.segments.len()) as f64);
            let step = sample_rate * (self.last().epoch() - self.first().epoch());

            /* *** */
            /* Map: bucket the states and send on a channel */
            /* *** */

            let mut window_states = Vec::with_capacity(2 * items_per_segments);

            // Note that we're using the typical map+reduce pattern
            for original_state in self.every(step) {
                let state = cosm.frame_chg(&original_state, new_frame);
                window_states.push(state);
                if window_states.len() == 2 * items_per_segments {
                    // Publish the first items of this vector
                    let this_wdn = window_states[..items_per_segments]
                        .iter()
                        .map(|&x| x)
                        .collect::<Vec<Orbit>>();

                    tx.send(this_wdn).map_err(|_| {
                        NyxError::from(TrajError::CreationError(
                            "could not send onto channel".to_string(),
                        ))
                    })?;

                    // Now, let's remove the first states
                    for _ in 0..items_per_segments - 1 {
                        window_states.remove(0);
                    }
                }
            }

            if window_states.last().unwrap().epoch() != self.last().epoch() {
                // Our final step placed us out of the trajectory epochs, so let's add it
                let state = cosm.frame_chg(&self.last(), new_frame);
                window_states.push(state);
            }

            // And interpolate the remaining states too, even if the buffer is not full!
            let mut start_idx = 0;
            loop {
                if window_states.is_empty() {
                    break;
                }

                tx.send(
                    window_states[start_idx..start_idx + items_per_segments]
                        .iter()
                        .map(|&x| x)
                        .collect::<Vec<Orbit>>(),
                )
                .map_err(|_| {
                    NyxError::from(TrajError::CreationError(
                        "could not send final window onto channel".to_string(),
                    ))
                })?;
                if start_idx > 0 {
                    break;
                }
                start_idx = window_states.len() - items_per_segments;
                if start_idx == 0 {
                    // This means that the window states are exactly the correct size, break here
                    break;
                }
            }

            rx
        };

        /* *** */
        /* Reduce: Build an interpolation of each of the segments */
        /* *** */
        let splines: Vec<_> = rx.into_iter().par_bridge().map(interpolate).collect();

        // Finally, build the whole trajectory
        let mut traj = Traj {
            segments: BTreeMap::new(),
            start_state,
            backward: false,
        };

        for maybe_spline in splines {
            let spline = maybe_spline?;
            // This shouldn't fail because we're progressing through the trajectory in a chronological fashion
            traj.append_spline(spline).unwrap();
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
        self.to_csv_with_step(filename, 1 * Unit::Minute, cosm)
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
        self.to_csv_between_with_step(filename, start, end, 1 * Unit::Minute, cosm)
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
            .every(1 * Unit::Minute)
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
    #[allow(clippy::map_clone)]
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
            let (tx, rx) = channel();

            let items_per_segments = INTERPOLATION_SAMPLES;

            // Nyquist–Shannon sampling theorem
            let sample_rate = 1.0 / (((items_per_segments * 2 + 1) * self.segments.len()) as f64);
            let step = sample_rate * (self.last().epoch() - self.first().epoch());

            /* *** */
            /* Map: bucket the states and send on a channel */
            /* *** */

            let mut window_states = Vec::with_capacity(2 * items_per_segments);

            // Note that we're using the typical map+reduce pattern
            for original_state in self.every(step) {
                let converted_orbit = cosm.frame_chg(&original_state.orbit, new_frame);
                let mut converted_state = original_state;
                // Update the orbit component and the frame
                converted_state.orbit.x = converted_orbit.x;
                converted_state.orbit.y = converted_orbit.y;
                converted_state.orbit.z = converted_orbit.z;
                converted_state.orbit.vx = converted_orbit.vx;
                converted_state.orbit.vy = converted_orbit.vy;
                converted_state.orbit.vz = converted_orbit.vz;
                converted_state.orbit.frame = converted_orbit.frame;

                window_states.push(converted_state);
                if window_states.len() == 2 * items_per_segments {
                    // Publish the first items of this vector
                    let this_wdn = window_states[..items_per_segments]
                        .iter()
                        .map(|&x| x)
                        .collect::<Vec<Spacecraft>>();

                    tx.send(this_wdn).map_err(|_| {
                        NyxError::from(TrajError::CreationError(
                            "could not send onto channel".to_string(),
                        ))
                    })?;

                    // Now, let's remove the first states
                    for _ in 0..items_per_segments - 1 {
                        window_states.remove(0);
                    }
                }
            }
            if window_states.last().unwrap().epoch() != self.last().epoch() {
                let original_state = self.last();
                // Our final step placed us out of the trajectory epochs, so let's add it
                let converted_orbit = cosm.frame_chg(&original_state.orbit, new_frame);
                window_states.push(original_state.with_orbit(converted_orbit));
            }

            // And interpolate the remaining states too, even if the buffer is not full!
            let mut start_idx = 0;
            loop {
                if window_states.is_empty() {
                    break;
                }

                tx.send(
                    window_states[start_idx..start_idx + items_per_segments]
                        .iter()
                        .map(|&x| x)
                        .collect::<Vec<Spacecraft>>(),
                )
                .map_err(|_| {
                    NyxError::from(TrajError::CreationError(
                        "could not send final window onto channel".to_string(),
                    ))
                })?;
                if start_idx > 0 {
                    break;
                }
                start_idx = window_states.len() - items_per_segments;
                if start_idx == 0 {
                    // This means that the window states are exactly the correct size, break here
                    break;
                }
            }

            rx
        };

        /* *** */
        /* Reduce: Build an interpolation of each of the segments */
        /* *** */
        let splines: Vec<_> = rx.into_iter().par_bridge().map(interpolate).collect();

        // Finally, build the whole trajectory
        let mut traj = Traj {
            segments: BTreeMap::new(),
            start_state,
            backward: false,
        };

        for maybe_spline in splines {
            let spline = maybe_spline?;
            // This shouldn't fail because we're progressing through the trajectory in a chronological fashion
            traj.append_spline(spline).unwrap();
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
        self.to_csv_with_step(filename, 1 * Unit::Minute, cosm)
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
        self.to_csv_between_with_step(filename, start, end, 1 * Unit::Minute, cosm)
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
            .every(1 * Unit::Minute)
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
    if this_wdn.is_empty() {
        return Err(NyxError::from(TrajError::CreationError(
            "cannot interpolate with zero items".to_string(),
        )));
    }
    // Generate interpolation and flush.
    let mut start_win_epoch = this_wdn.first().unwrap().epoch();
    let mut end_win_epoch = this_wdn.last().unwrap().epoch();
    let mut end_state = this_wdn.last().unwrap();
    if end_win_epoch < start_win_epoch {
        // Backward propagation, swap times
        std::mem::swap(&mut start_win_epoch, &mut end_win_epoch);
        // Swap end states
        end_state = this_wdn.first().unwrap();
    }
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
            let delta_t: f64 = *latest_t - t_prime;
            if delta_t.abs() < f64::EPSILON {
                continue;
            }
        }

        ts.push(t_prime);

        for (pos, param) in S::params().iter().enumerate() {
            let (value, deriv) = state.value_and_deriv(param)?;
            values[pos].push(value);
            values_dt[pos].push(deriv);
        }
    }

    // Generate the polynomials
    for pos in 0..values.len() {
        let poly = if values[pos].len() != INTERPOLATION_SAMPLES {
            // Do not fit the data, just build the polynomials for these values
            let p =
                hermite::hermite::<{ 2 * SPLINE_DEGREE + 1 }>(&ts, &values[pos], &values_dt[pos])?;
            // println!("{:x}", p);
            let mut coefficients: [f64; SPLINE_DEGREE] = [0.0; SPLINE_DEGREE];
            // Rebuild the coeffs ignoring the highest power
            coefficients[..SPLINE_DEGREE].clone_from_slice(&p.coefficients[..SPLINE_DEGREE]);
            Polynomial { coefficients }
        } else {
            hermite::hermfit::<INTERPOLATION_SAMPLES, SPLINE_DEGREE>(&ts, &values[pos])?
        };

        polynomials.push(poly);
    }

    Ok(Spline {
        start_epoch: start_win_epoch,
        duration: window_duration,
        polynomials,
        end_state: *end_state,
    })
}
