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

use super::bacon_sci::interp::lagrange;
use super::bacon_sci::polynomial::Polynomial;
use super::crossbeam::thread;
use super::rayon::prelude::*;
use crate::celestia::{Cosm, Frame, Orbit, SpacecraftState};
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, DimName, VectorN};
use crate::errors::NyxError;
use crate::md::events::EventEvaluator;
use crate::time::{Duration, Epoch, TimeSeries, TimeUnit};
use crate::utils::normalize;
use crate::{State, TimeTagged};
use std::collections::BTreeMap;
use std::iter::Iterator;
use std::sync::mpsc::{channel, Receiver};
use std::sync::Arc;
use std::time::Duration as StdDur;

const INTERP_TOLERANCE: f64 = 1e-10;
const ITEMS_PER_SEGMENTS: usize = 12;

/// Stores a segment of an interpolation
#[derive(Clone)]
pub struct Segment<S: State>
where
    DefaultAllocator: Allocator<f64, S::PropVecSize> + Allocator<f64, S::Size>,
{
    start_epoch: Epoch,
    duration: Duration,
    coefficients: Vec<Vec<f64>>,
    end_state: S,
}

/// Store a trajectory of any State.
#[derive(Clone)]
pub struct Traj<S: State>
where
    DefaultAllocator: Allocator<f64, S::PropVecSize> + Allocator<f64, S::Size>,
{
    /// Segments are organized as a binary tree map whose index is the
    /// number of seconds since the start of this ephemeris rounded down
    pub segments: BTreeMap<i32, Segment<S>>,
    /// Timeout is used to stop listening to new state (default is 50ms, should work well in release and debug mode).
    pub timeout_ms: u64,
    start_state: S,
    max_offset: i32,
}

pub struct TrajIterator<'a, S: State>
where
    DefaultAllocator: Allocator<f64, S::PropVecSize> + Allocator<f64, S::Size>,
{
    pub time_series: TimeSeries,
    /// A shared pointer to the original trajectory.
    pub traj: &'a Traj<S>,
}

impl<S: State> Segment<S>
where
    DefaultAllocator: Allocator<f64, S::PropVecSize> + Allocator<f64, S::Size>,
{
    /// Evaluate a specific segment at the provided Epoch, requires an initial state as a "template"
    pub fn evaluate(&self, from: S, epoch: Epoch) -> Result<S, NyxError> {
        // Compute the normalized time
        let dur_into_window = epoch - self.start_epoch;
        if dur_into_window > self.duration {
            return Err(NyxError::OutOfInterpolationWindow(format!(
                "Requested trajectory at time {} but that is past the final interpolation window by {}",
                epoch, dur_into_window
            )));
        } else if dur_into_window.in_seconds() < -1.0 {
            // We should not be in this window, but in the next one
            // We allow for a delta of one second because of the rounding of the indexing.
            return Err(NyxError::InvalidInterpolationData(format!(
                "Bug: should be in next window: {}",
                dur_into_window
            )));
        }

        let t_prime = normalize(
            dur_into_window.in_seconds(),
            0.0,
            self.duration.in_seconds(),
        );

        let mut state = from;

        // Rebuild the polynominals
        let mut state_vec = VectorN::<f64, S::PropVecSize>::zeros();
        for (cno, coeffs) in self.coefficients.iter().enumerate() {
            state_vec[cno] = Polynomial::from_slice(coeffs).evaluate(t_prime)
        }
        state.set(epoch, &state_vec)?;

        Ok(state)
    }
}

impl<S: State> Traj<S>
where
    DefaultAllocator: Allocator<f64, S::PropVecSize> + Allocator<f64, S::Size>,
{
    /// Creates a new trajectory with the provided starting state (used as a template) and a receiving channel.
    /// The trajectories are always generated on a separate thread.
    pub fn new(state: S, rx: Receiver<S>) -> Result<Self, NyxError> {
        thread::scope(|s| {
            // Initialize the interpolator
            let mut me = Self {
                segments: BTreeMap::new(),
                start_state: state,
                timeout_ms: 100,
                max_offset: 0,
            };

            let mut children = vec![];
            let mut window_states: Vec<S> = Vec::with_capacity(ITEMS_PER_SEGMENTS);
            // Push the initial state
            window_states.push(state);

            // Note that we're using the typical map+reduce pattern
            // Start receiving states on a blocking call (map)
            while let Ok(state) = rx.recv_timeout(StdDur::from_millis(me.timeout_ms)) {
                if window_states.len() == ITEMS_PER_SEGMENTS {
                    let this_wdn = window_states.clone();
                    children.push(
                        s.spawn(move |_| -> Result<Segment<S>, NyxError> { interpolate(this_wdn) }),
                    );
                    // Copy the last state as the first state of the next window
                    let last_wdn_state = window_states[ITEMS_PER_SEGMENTS - 1];
                    window_states.clear();
                    window_states.push(last_wdn_state);
                }
                window_states.push(state);
            }
            // And interpolate the remaining states too, even if the buffer is not full!
            children.push(
                s.spawn(move |_| -> Result<Segment<S>, NyxError> { interpolate(window_states) }),
            );

            // Reduce
            for child in children {
                // collect each child thread's return-value
                let segment = child.join().unwrap()?;
                me.append_segment(segment);
            }

            Ok(me)
        })
        .unwrap()
    }

    fn append_segment(&mut self, segment: Segment<S>) {
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
    pub fn evaluate(&self, epoch: Epoch) -> Result<S, NyxError> {
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
        TrajIterator {
            time_series: TimeSeries::inclusive(self.first().epoch(), self.last().epoch(), step),
            traj: &self,
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
        let mut ya = event.eval(&self.evaluate(xa_e)?);
        let mut yb = event.eval(&self.evaluate(xb_e)?);

        // Check if we're already at the root
        if ya.abs() <= event.value_precision().abs() {
            return Ok(self.evaluate(xa_e)?);
        } else if yb.abs() <= event.value_precision().abs() {
            return Ok(self.evaluate(xb_e)?);
        }
        // The Brent solver, from the roots crate (sadly could not directly integrate it here)
        // Source: https://docs.rs/roots/0.0.5/src/roots/numerical/brent.rs.html#57-131

        let (mut xc, mut yc, mut xd) = (xa, ya, xa);
        let mut flag = true;

        for _ in 0..max_iter {
            if ya.abs() < event.value_precision().abs() {
                return Ok(self.evaluate(xa_e + xa * TimeUnit::Second)?);
            }
            if yb.abs() < event.value_precision().abs() {
                return Ok(self.evaluate(xa_e + xb * TimeUnit::Second)?);
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
            let next_try = self.evaluate(xa_e + s * TimeUnit::Second)?;
            let ys = event.eval(&next_try);
            xd = xc;
            xc = xb;
            yc = yb;
            if ya * ys < 0.0 {
                // Root bracketed between a and s
                let next_try = self.evaluate(xa_e + xa * TimeUnit::Second)?;
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
                let next_try = self.evaluate(xa_e + xb * TimeUnit::Second)?;
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
        Err(NyxError::MaxIterReached(max_iter))
    }

    /// Find (usually) all of the states where the event happens.
    /// WARNING: The initial search step is 1% of the duration of the trajectory duration!
    /// For example, if the trajectory is 100 days long, then we split the trajectory into 100 chunks of 1 day and see whether
    /// the event is in there. If the event happens twice or more times within 1% of the trajectory duration, only the _one_ of
    /// such events will be found.
    pub fn find_all<E>(&self, event: &E) -> Result<Vec<S>, NyxError>
    where
        E: EventEvaluator<S>,
    {
        let start_epoch = self.first().epoch();
        let end_epoch = self.last().epoch();
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
            return Err(NyxError::EventNotInEpochBraket(
                start_epoch.to_string(),
                end_epoch.to_string(),
            ));
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
    pub fn find_minmax<E>(&self, event: E, precision: TimeUnit) -> Result<(S, S), NyxError>
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
        println!("search over {}", epochs.len());
        epochs.into_par_iter().for_each_with(sender, |s, epoch| {
            let state = self.evaluate(epoch).unwrap();
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
        thread::scope(|s| {
            let (tx, rx) = channel();

            let start_state = cosm.frame_chg(&self.first(), new_frame);
            // The trajectory must always be generated on its own thread.
            let traj_thread = s.spawn(move |_| Self::new(start_state, rx));
            // And start sampling, converting, and flushing into the new trajectory.

            // Nyquist–Shannon sampling theorem
            let sample_rate = 1.0 / (((32 * 2 + 1) * self.segments.len()) as f64);
            let step = sample_rate * (self.last().epoch() - self.first().epoch());

            for state in self.every(step) {
                let converted_state = cosm.frame_chg(&state, new_frame);
                tx.send(converted_state).unwrap();
            }

            let traj = traj_thread.join().unwrap_or_else(|_| {
                Err(NyxError::NoInterpolationData(
                    "Could not generate trajectory".to_string(),
                ))
            })?;

            Ok(traj)
        })
        .unwrap()
    }
}

impl Traj<SpacecraftState> {
    /// Allows converting the source trajectory into the (almost) equivalent trajectory in another frame
    /// This is super slow.
    pub fn to_frame(&self, new_frame: Frame, cosm: Arc<Cosm>) -> Result<Self, NyxError> {
        thread::scope(|s| {
            let (tx, rx) = channel();

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

            // The trajectory must always be generated on its own thread.
            let traj_thread = s.spawn(move |_| Self::new(start_state, rx));
            // And start sampling, converting, and flushing into the new trajectory.

            // Nyquist–Shannon sampling theorem
            let sample_rate = 1.0 / (((16 * 2 + 1) * self.segments.len()) as f64);
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

            let traj = traj_thread.join().unwrap_or_else(|_| {
                Err(NyxError::NoInterpolationData(
                    "Could not generate trajectory".to_string(),
                ))
            })?;

            Ok(traj)
        })
        .unwrap()
    }
}

impl<S: State> Iterator for TrajIterator<'_, S>
where
    DefaultAllocator: Allocator<f64, S::PropVecSize> + Allocator<f64, S::Size>,
{
    type Item = S;

    fn next(&mut self) -> Option<Self::Item> {
        match self.time_series.next() {
            Some(next_epoch) => match self.traj.evaluate(next_epoch) {
                Ok(item) => Some(item),
                _ => None,
            },
            None => None,
        }
    }
}

fn interpolate<S: State>(this_wdn: Vec<S>) -> Result<Segment<S>, NyxError>
where
    DefaultAllocator: Allocator<f64, S::PropVecSize> + Allocator<f64, S::Size>,
{
    // Generate interpolation and flush.
    let start_win_epoch = this_wdn[0].epoch();
    let end_win_epoch = this_wdn[this_wdn.len() - 1].epoch();
    let window_duration = end_win_epoch - start_win_epoch;
    let mut ts = Vec::new();
    let mut values = Vec::with_capacity(S::PropVecSize::dim());
    let mut coefficients = Vec::new();
    // Initialize the vector of values and coefficients.
    for _ in 0..S::PropVecSize::dim() {
        values.push(Vec::with_capacity(this_wdn.len()));
        coefficients.push(Vec::with_capacity(this_wdn.len()));
    }
    for state in &this_wdn {
        let t_prime = normalize(
            (state.epoch() - start_win_epoch).in_seconds(),
            0.0,
            window_duration.in_seconds(),
        );
        ts.push(t_prime);
        for (pos, val) in state.as_vector().as_ref().unwrap().iter().enumerate() {
            values[pos].push(*val);
        }
    }

    // Generate the polynomials
    for (pos, these_values) in values.iter().enumerate() {
        match lagrange(&ts, these_values, INTERP_TOLERANCE) {
            Ok(polyn) => coefficients[pos] = polyn.get_coefficients(),
            Err(e) => {
                return Err(NyxError::InvalidInterpolationData(format!(
                    "Interpolation failed: {}",
                    e
                )))
            }
        };
    }

    Ok(Segment {
        start_epoch: start_win_epoch,
        duration: window_duration,
        coefficients,
        end_state: this_wdn[this_wdn.len() - 1],
    })
}
