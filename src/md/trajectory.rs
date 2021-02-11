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
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, DimName, VectorN};
use crate::errors::NyxError;
use crate::propagators::events::Event;
use crate::time::{Duration, Epoch, TimeUnit};
use crate::State;
use std::collections::BTreeMap;
use std::sync::mpsc::Receiver;
use std::time::Duration as StdDur;

const INTERP_TOLERANCE: f64 = 1e-10;

/// Stores a segment of an interpolation
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
pub struct Traj<S: State>
where
    DefaultAllocator: Allocator<f64, S::PropVecSize> + Allocator<f64, S::Size>,
{
    /// Segments are organized as a binary tree map whose index is the
    /// number of seconds since the start of this ephemeris rounded down
    pub segments: BTreeMap<u32, Segment<S>>,
    /// Timeout is used to stop listening to new state (default is 50ms, should work well in release and debug mode).
    pub timeout_ms: u64,
    start_state: S,
    max_offset: u32,
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
        } else if dur_into_window.in_seconds() < 0.0 {
            // We should not be in this window, but in the next one
            println!("Oh no: {}", dur_into_window.in_seconds());
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

impl<S: State + 'static> Traj<S>
where
    DefaultAllocator: Allocator<f64, S::PropVecSize> + Allocator<f64, S::Size>,
{
    /// Creates a new trajectory with the provided starting state (used as a template) and a receiving channel.
    /// The trajectories are always generated on a separate thread.
    pub fn new(state: S, rx: Receiver<S>) -> Result<Self, NyxError> {
        // Initialize the interpolator
        let mut me = Self {
            segments: BTreeMap::new(),
            start_state: state,
            timeout_ms: 100,
            max_offset: 0,
        };

        // Bug? With a spacecraft, we need more interpolation windows than just an orbit.
        // I've spent 12h trying to understand why, but I can't, so screw it for it.
        let items_per_segments = if S::PropVecSize::dim() == 43 { 16 } else { 32 };

        let mut children = vec![];
        let mut window_states: Vec<S> = Vec::with_capacity(items_per_segments);
        // Push the initial state
        window_states.push(state);

        // Note that we're using the typical map+reduce pattern
        // Start receiving states on a blocking call (map)
        while let Ok(state) = rx.recv_timeout(StdDur::from_millis(me.timeout_ms)) {
            if window_states.len() == items_per_segments {
                let this_wdn = window_states.clone();
                children.push(std::thread::spawn(
                    move || -> Result<Segment<S>, NyxError> { interpolate(this_wdn) },
                ));

                window_states.clear();
            }
            window_states.push(state);
        }
        // And interpolate the remaining states too, even if the buffer is not full!
        children.push(std::thread::spawn(
            move || -> Result<Segment<S>, NyxError> { interpolate(window_states) },
        ));

        // Reduce
        for child in children {
            // collect each child thread's return-value
            let segment = child.join().unwrap()?;
            me.append_segment(segment);
        }

        Ok(me)
    }

    fn append_segment(&mut self, segment: Segment<S>) {
        // Compute the number of seconds since start of trajectory
        let offset_s = ((segment.start_epoch - self.start_state.epoch())
            .in_seconds()
            .floor()) as u32;
        self.segments.insert(offset_s, segment);
        if offset_s > self.max_offset {
            self.max_offset = offset_s;
        }
    }

    /// Evaluate the trajectory at this specific epoch.
    pub fn evaluate(&self, epoch: Epoch) -> Result<S, NyxError> {
        let offset_s = ((epoch - self.start_state.epoch()).in_seconds().floor()) as u32;
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

    /// Find the exact state where the request event happens. The event function is expected to be monotone in the provided interval.
    pub fn find(
        &self,
        start: Epoch,
        end: Epoch,
        event: Box<dyn Event<StateType = S>>,
        epsilon: Duration,
        epsilon_eval: f64,
    ) -> Result<S, NyxError> {
        use std::f64::EPSILON;
        let xa_e = start;
        let xb_e = end;
        // let (xa_e, xb_e) = self.event_trackers.found_bounds[0][condition.trigger - 1];
        let mut xa = 0.0;
        let mut xb = (xb_e - xa_e).in_seconds();
        // Evaluate the event at both bounds
        let mut ya = event.eval(&self.evaluate(xa_e)?);
        let mut yb = event.eval(&self.evaluate(xb_e)?);
        dbg!(ya, yb);
        // The Brent solver, from the roots crate (sadly could not directly integrate it here)
        // Source: https://docs.rs/roots/0.0.5/src/roots/numerical/brent.rs.html#57-131

        // Helper lambdas, for f64s only
        let has_converged = |x1: f64, x2: f64| (x1 - x2).abs() <= epsilon.in_seconds();
        let arrange = |a: f64, ya: f64, b: f64, yb: f64| {
            if ya.abs() > yb.abs() {
                dbg!(a, ya, b, yb)
            } else {
                dbg!(b, yb, a, ya)
            }
        };

        let (mut c, mut yc, mut d) = (xa, ya, xa);
        let mut flag = true;
        let mut iter = 0;
        let closest_t;
        loop {
            if ya.abs() < epsilon_eval.abs() {
                closest_t = xa;
                break;
            }
            if yb.abs() < epsilon_eval.abs() {
                closest_t = xb;
                break;
            }
            if has_converged(xa, xb) {
                closest_t = c;
                break;
            }
            let mut s = if (ya - yc).abs() > EPSILON && (yb - yc).abs() > EPSILON {
                xa * yb * yc / ((ya - yb) * (ya - yc))
                    + xb * ya * yc / ((yb - ya) * (yb - yc))
                    + c * ya * yb / ((yc - ya) * (yc - yb))
            } else {
                xb - yb * (xb - xa) / (yb - ya)
            };
            let cond1 = (s - xb) * (s - (3.0 * xa + xb) / 4.0) > 0.0;
            let cond2 = flag && (s - xb).abs() >= (xb - c).abs() / 2.0;
            let cond3 = !flag && (s - xb).abs() >= (c - d).abs() / 2.0;
            let cond4 = flag && has_converged(xb, c);
            let cond5 = !flag && has_converged(c, d);
            if cond1 || cond2 || cond3 || cond4 || cond5 {
                s = (xa + xb) / 2.0;
                flag = true;
            } else {
                flag = false;
            }
            // Propagate until time s
            // self.for_duration(s * TimeUnit::Second)?;
            // let ys = self.event_trackers.events[0].eval(&self.state);
            println!("240 Fetching state at time {}", xa_e + s * TimeUnit::Second);
            let ys = event.eval(&self.evaluate(xa_e + s * TimeUnit::Second)?);
            d = c;
            c = xb;
            yc = yb;
            if ya * ys < 0.0 {
                // Root bracketed between a and s
                // Propagate until time xa
                // self.for_duration(xa * TimeUnit::Second)?;
                // let ya_p = self.event_trackers.events[0].eval(&self.state);
                println!(
                    "250 Fetching state at time {}",
                    xa_e + xa * TimeUnit::Second
                );
                let ya_p = event.eval(&self.evaluate(xa_e + xa * TimeUnit::Second)?);
                let (_a, _ya, _b, _yb) = arrange(xa, ya_p, s, ys);
                {
                    xa = _a;
                    ya = _ya;
                    xb = _b;
                    yb = _yb;
                }
            } else {
                // Root bracketed between s and b
                // Propagate until time xb
                // self.for_duration(xb * TimeUnit::Second)?;
                // let yb_p = self.event_trackers.events[0].eval(&self.state);
                println!(
                    "264 Fetching state at time {}",
                    xa_e + xb * TimeUnit::Second
                );
                let yb_p = event.eval(&self.evaluate(xa_e + xb * TimeUnit::Second)?);
                let (_a, _ya, _b, _yb) = arrange(s, ys, xb, yb_p);
                {
                    xa = _a;
                    ya = _ya;
                    xb = _b;
                    yb = _yb;
                }
            }
            iter += 1;
            if iter >= 50 {
                // condition.max_iter
                return Err(NyxError::MaxIterReached(iter));
            }
        }

        // Now that we have the time at which the condition is matched, let's propagate until then
        let converged_time = xa_e - (closest_t * TimeUnit::Second);
        println!("{}", closest_t);

        self.evaluate(converged_time)
    }
}

// Normalize between -1.0 and 1.0
fn normalize(x: f64, min_x: f64, max_x: f64) -> f64 {
    2.0 * (x - min_x) / (max_x - min_x) - 1.0
}

// Denormalize between -1.0 and 1.0
fn _denormalize(xp: f64, min_x: f64, max_x: f64) -> f64 {
    (max_x - min_x) * (xp + 1.0) / 2.0 + min_x
}

fn interpolate<S: State + 'static>(this_wdn: Vec<S>) -> Result<Segment<S>, NyxError>
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
