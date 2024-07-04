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

use super::details::{EventArc, EventDetails, EventEdge};
use crate::errors::{EventError, EventTrajSnafu};
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::md::prelude::{Interpolatable, Traj};
use crate::md::EventEvaluator;
use crate::time::{Duration, Epoch, TimeSeries, Unit};
use anise::almanac::Almanac;
use rayon::prelude::*;
use snafu::ResultExt;
use std::iter::Iterator;
use std::sync::mpsc::channel;
use std::sync::Arc;

impl<S: Interpolatable> Traj<S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    /// Find the exact state where the request event happens. The event function is expected to be monotone in the provided interval because we find the event using a Brent solver.
    #[allow(clippy::identity_op)]
    pub fn find_bracketed<E>(
        &self,
        start: Epoch,
        end: Epoch,
        event: &E,
        almanac: Arc<Almanac>,
    ) -> Result<EventDetails<S>, EventError>
    where
        E: EventEvaluator<S>,
    {
        let max_iter = 50;

        // Helper lambdas, for f64s only
        let has_converged =
            |xa: f64, xb: f64| (xa - xb).abs() <= event.epoch_precision().to_seconds();
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
        let ya_state = self.at(xa_e).context(EventTrajSnafu {})?;
        let yb_state = self.at(xb_e).context(EventTrajSnafu {})?;
        let mut ya = event.eval(&ya_state, almanac.clone())?;
        let mut yb = event.eval(&yb_state, almanac.clone())?;

        // Check if we're already at the root
        if ya.abs() <= event.value_precision().abs() {
            debug!(
                "{event} -- found with |{ya}| < {} @ {xa_e}",
                event.value_precision().abs()
            );
            return EventDetails::new(ya_state, ya, event, self, almanac.clone());
        } else if yb.abs() <= event.value_precision().abs() {
            debug!(
                "{event} -- found with |{yb}| < {} @ {xb_e}",
                event.value_precision().abs()
            );
            return EventDetails::new(yb_state, yb, event, self, almanac.clone());
        }

        // The Brent solver, from the roots crate (sadly could not directly integrate it here)
        // Source: https://docs.rs/roots/0.0.5/src/roots/numerical/brent.rs.html#57-131

        let (mut xc, mut yc, mut xd) = (xa, ya, xa);
        let mut flag = true;

        for _ in 0..max_iter {
            if ya.abs() < event.value_precision().abs() {
                // Can't fail, we got it earlier
                let state = self.at(xa_e + xa * Unit::Second).unwrap();
                debug!(
                    "{event} -- found with |{ya}| < {} @ {}",
                    event.value_precision().abs(),
                    state.epoch(),
                );
                return EventDetails::new(state, ya, event, self, almanac.clone());
            }
            if yb.abs() < event.value_precision().abs() {
                // Can't fail, we got it earlier
                let state = self.at(xa_e + xb * Unit::Second).unwrap();
                debug!(
                    "{event} -- found with |{yb}| < {} @ {}",
                    event.value_precision().abs(),
                    state.epoch()
                );
                return EventDetails::new(state, yb, event, self, almanac.clone());
            }
            if has_converged(xa, xb) {
                // The event isn't in the bracket
                return Err(EventError::NotFound {
                    start,
                    end,
                    event: format!("{event}"),
                });
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
            let next_try = self
                .at(xa_e + s * Unit::Second)
                .context(EventTrajSnafu {})?;
            let ys = event.eval(&next_try, almanac.clone())?;
            xd = xc;
            xc = xb;
            yc = yb;
            if ya * ys < 0.0 {
                // Root bracketed between a and s
                let next_try = self
                    .at(xa_e + xa * Unit::Second)
                    .context(EventTrajSnafu {})?;
                let ya_p = event.eval(&next_try, almanac.clone())?;
                let (_a, _ya, _b, _yb) = arrange(xa, ya_p, s, ys);
                {
                    xa = _a;
                    ya = _ya;
                    xb = _b;
                    yb = _yb;
                }
            } else {
                // Root bracketed between s and b
                let next_try = self
                    .at(xa_e + xb * Unit::Second)
                    .context(EventTrajSnafu {})?;
                let yb_p = event.eval(&next_try, almanac.clone())?;
                let (_a, _ya, _b, _yb) = arrange(s, ys, xb, yb_p);
                {
                    xa = _a;
                    ya = _ya;
                    xb = _b;
                    yb = _yb;
                }
            }
        }
        error!("Brent solver failed after {max_iter} iterations");
        Err(EventError::NotFound {
            start,
            end,
            event: format!("{event}"),
        })
    }

    /// Find all of the states where the event happens
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
    pub fn find<E>(
        &self,
        event: &E,
        almanac: Arc<Almanac>,
    ) -> Result<Vec<EventDetails<S>>, EventError>
    where
        E: EventEvaluator<S>,
    {
        let start_epoch = self.first().epoch();
        let end_epoch = self.last().epoch();
        if start_epoch == end_epoch {
            return Err(EventError::NotFound {
                start: start_epoch,
                end: end_epoch,
                event: format!("{event}"),
            });
        }
        let heuristic = (end_epoch - start_epoch) / 100;
        info!("Searching for {event} with initial heuristic of {heuristic}");

        let (sender, receiver) = channel();

        let epochs: Vec<Epoch> = TimeSeries::inclusive(start_epoch, end_epoch, heuristic).collect();
        epochs.into_par_iter().for_each_with(sender, |s, epoch| {
            if let Ok(event_state) =
                self.find_bracketed(epoch, epoch + heuristic, event, almanac.clone())
            {
                s.send(event_state).unwrap()
            };
        });

        let mut states: Vec<_> = receiver.iter().collect();

        if states.is_empty() {
            warn!("Heuristic failed to find any {event} event, using slower approach");
            // Crap, we didn't find the event.
            // Let's find the min and max of this event throughout the trajectory, and search around there.
            match self.find_minmax(event, Unit::Second, almanac.clone()) {
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
                    if let Ok(event_state) = self.find_bracketed(
                        lower_min_epoch,
                        lower_max_epoch,
                        event,
                        almanac.clone(),
                    ) {
                        states.push(event_state);
                    };

                    // Search around the max event
                    if let Ok(event_state) = self.find_bracketed(
                        upper_min_epoch,
                        upper_max_epoch,
                        event,
                        almanac.clone(),
                    ) {
                        states.push(event_state);
                    };

                    // If there still isn't any match, report that the event was not found
                    if states.is_empty() {
                        return Err(EventError::NotFound {
                            start: start_epoch,
                            end: end_epoch,
                            event: format!("{event}"),
                        });
                    }
                }
                Err(_) => {
                    return Err(EventError::NotFound {
                        start: start_epoch,
                        end: end_epoch,
                        event: format!("{event}"),
                    });
                }
            };
        }
        // Remove duplicates and reorder
        states.sort_by(|s1, s2| s1.state.epoch().partial_cmp(&s2.state.epoch()).unwrap());
        states.dedup();

        match states.len() {
            0 => info!("Event {event} not found"),
            1 => info!("Event {event} found once on {}", states[0].state.epoch()),
            _ => {
                info!(
                    "Event {event} found {} times from {} until {}",
                    states.len(),
                    states.first().unwrap().state.epoch(),
                    states.last().unwrap().state.epoch()
                )
            }
        };

        Ok(states)
    }

    /// Find the minimum and maximum of the provided event through the trajectory
    #[allow(clippy::identity_op)]
    pub fn find_minmax<E>(
        &self,
        event: &E,
        precision: Unit,
        almanac: Arc<Almanac>,
    ) -> Result<(S, S), EventError>
    where
        E: EventEvaluator<S>,
    {
        let step: Duration = 1 * precision;
        let mut min_val = f64::INFINITY;
        let mut max_val = f64::NEG_INFINITY;
        let mut min_state = S::zeros();
        let mut max_state = S::zeros();

        let (sender, receiver) = channel();

        let epochs: Vec<Epoch> =
            TimeSeries::inclusive(self.first().epoch(), self.last().epoch(), step).collect();

        epochs.into_par_iter().for_each_with(sender, |s, epoch| {
            // The `at` call will work because we only query within the start and end times of the trajectory
            let state = self.at(epoch).unwrap();
            if let Ok(this_eval) = event.eval(&state, almanac.clone()) {
                s.send((this_eval, state)).unwrap();
            }
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

    /// Identifies and pairs rising and falling edge events in a trajectory.
    ///
    /// This function processes a sequence of events in a trajectory and pairs each rising edge event with its subsequent falling edge event to form arcs.
    /// Each arc represents a complete cycle of an event rising above and then falling below a specified threshold.
    /// Use this to analyze a trajectory's behavior when understanding the complete cycle of an event (from rising to falling) is essential, e.g. ground station passes.
    ///
    /// # Arguments
    /// - `event`: A reference to an object implementing the `EventEvaluator<S>` trait, which is used to evaluate and classify events in the trajectory.
    ///
    /// # Returns
    /// - `Result<Vec<EventArc>, NyxError>`: On success, returns a vector of EventArc, where each struct contains a pair of `EventDetails` (one for the rising edge and one for the falling edge). Returns an error if any issues occur during the event evaluation process.
    ///
    /// # Logic
    /// - Sorts the events by their epoch to ensure chronological processing.
    /// - Iterates through the sorted events, identifying transitions from falling to rising edges and vice versa.
    /// - Pairs a rising edge with the subsequent falling edge to form an arc.
    /// - Handles edge cases where the trajectory starts or ends with a rising or falling edge.
    /// - Prints debug information for each event and arc.
    ///
    /// ## Note
    /// If no zero crossing happens in the trajectory, i.e. the there is "event is true" _and_ "event is false",
    /// then this function checks whether the event is true at the start and end of the trajectory. If so, it means
    /// that there is a single arc that spans the whole trajectory.
    ///
    pub fn find_arcs<E>(
        &self,
        event: &E,
        almanac: Arc<Almanac>,
    ) -> Result<Vec<EventArc<S>>, EventError>
    where
        E: EventEvaluator<S>,
    {
        let mut events = match self.find(event, almanac.clone()) {
            Ok(events) => events,
            Err(_) => {
                // We haven't found the start or end of an arc, i.e. no zero crossing on the event.
                // However, if the trajectory start and end are above the event value, then we found an arc.
                let first_eval = event.eval(self.first(), almanac.clone())?;
                let last_eval = event.eval(self.last(), almanac.clone())?;
                if first_eval > 0.0 && last_eval > 0.0 {
                    // No event crossing found, but from the start until the end of the trajectory, we're in the same arc
                    // because the evaluation of the event is above the zero crossing.
                    // Hence, there's a single arc, and it's from start until the end of the trajectory.
                    vec![
                        EventDetails::new(*self.first(), first_eval, event, self, almanac.clone())?,
                        EventDetails::new(*self.last(), last_eval, event, self, almanac.clone())?,
                    ]
                } else {
                    return Err(EventError::NotFound {
                        start: self.first().epoch(),
                        end: self.last().epoch(),
                        event: format!("{event}"),
                    });
                }
            }
        };
        events.sort_by_key(|event| event.state.epoch());

        // Now, let's pair the events.
        let mut arcs = Vec::new();

        if events.is_empty() {
            return Ok(arcs);
        }

        // If the first event isn't a rising edge, then we mark the start of the trajectory as a rising edge
        let mut prev_rise = if events[0].edge != EventEdge::Rising {
            let value = event.eval(self.first(), almanac.clone())?;
            Some(EventDetails::new(
                *self.first(),
                value,
                event,
                self,
                almanac.clone(),
            )?)
        } else {
            Some(events[0].clone())
        };

        let mut prev_fall = if events[0].edge == EventEdge::Falling {
            Some(events[0].clone())
        } else {
            None
        };

        for event in events {
            if event.edge == EventEdge::Rising {
                if prev_rise.is_none() && prev_fall.is_none() {
                    // This is a new rising edge
                    prev_rise = Some(event.clone());
                } else if prev_fall.is_some() {
                    // We've found a transition from a fall to a rise, so we can close this arc out.
                    if prev_rise.is_some() {
                        let arc = EventArc {
                            rise: prev_rise.clone().unwrap(),
                            fall: prev_fall.clone().unwrap(),
                        };
                        arcs.push(arc);
                    } else {
                        let arc = EventArc {
                            rise: event.clone(),
                            fall: prev_fall.clone().unwrap(),
                        };
                        arcs.push(arc);
                    }
                    prev_fall = None;
                    // We have a new rising edge since this is how we ended up here.
                    prev_rise = Some(event.clone());
                }
            } else if event.edge == EventEdge::Falling {
                prev_fall = Some(event.clone());
            }
        }

        // Add the final pass
        if prev_rise.is_some() {
            if prev_fall.is_some() {
                let arc = EventArc {
                    rise: prev_rise.clone().unwrap(),
                    fall: prev_fall.clone().unwrap(),
                };
                arcs.push(arc);
            } else {
                // Use the last trajectory as the end of the arc
                let value = event.eval(self.last(), almanac.clone())?;
                let fall = EventDetails::new(*self.last(), value, event, self, almanac.clone())?;
                let arc = EventArc {
                    rise: prev_rise.clone().unwrap(),
                    fall,
                };
                arcs.push(arc);
            }
        }

        Ok(arcs)
    }
}
