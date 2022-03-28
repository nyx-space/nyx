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

use super::error_ctrl::{ErrorCtrl, RSSCartesianStep};
use super::rayon::iter::ParallelBridge;
use super::rayon::prelude::ParallelIterator;
use super::{Dormand78, IntegrationDetails, RK, RK89};
use crate::dynamics::Dynamics;
use crate::errors::NyxError;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, OVector};
use crate::md::trajectory::spline::INTERPOLATION_SAMPLES;
use crate::md::trajectory::{interpolate, InterpState, Traj, TrajError};
use crate::md::EventEvaluator;
use crate::time::{Duration, Epoch, Unit};
use crate::State;
use std::collections::BTreeMap;
use std::f64;
use std::sync::mpsc::{channel, Sender};

/// A Propagator allows propagating a set of dynamics forward or backward in time.
/// It is an EventTracker, without any event tracking. It includes the options, the integrator
/// details of the previous step, and the set of coefficients used for the monomorphic instance.
#[derive(Clone, Debug)]
pub struct Propagator<'a, D: Dynamics, E: ErrorCtrl>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<usize, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>,
{
    pub dynamics: D, // Stores the dynamics used. *Must* use this to get the latest values
    pub opts: PropOpts<E>, // Stores the integration options (tolerance, min/max step, init step, etc.)
    order: u8,             // Order of the integrator
    stages: usize,         // Number of stages, i.e. how many times the derivatives will be called
    a_coeffs: &'a [f64],
    b_coeffs: &'a [f64],
}

/// The `Propagator` trait defines the functions of a propagator and of an event tracker.
impl<'a, D: Dynamics, E: ErrorCtrl> Propagator<'a, D, E>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<usize, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>,
{
    /// Each propagator must be initialized with `new` which stores propagator information.
    pub fn new<T: RK>(dynamics: D, opts: PropOpts<E>) -> Self {
        Self {
            dynamics,
            opts,
            stages: T::STAGES,
            order: T::ORDER,
            a_coeffs: T::A_COEFFS,
            b_coeffs: T::B_COEFFS,
        }
    }

    /// Set the tolerance for the propagator
    pub fn set_tolerance(&mut self, tol: f64) {
        self.opts.tolerance = tol;
    }

    /// Set the maximum step size for the propagator and sets the initial step to that value if currently greater
    pub fn set_max_step(&mut self, step: Duration) {
        self.opts.set_max_step(step);
    }

    pub fn set_min_step(&mut self, step: Duration) {
        self.opts.set_min_step(step);
    }

    /// An RK89 propagator (the default) with custom propagator options.
    pub fn rk89(dynamics: D, opts: PropOpts<E>) -> Self {
        Self::new::<RK89>(dynamics, opts)
    }

    /// A Dormand Prince 7-8 propagator with custom propagator options: it's about 20% faster than an RK98, and more stable in two body dynamics.
    /// WARNINGS: Dormand Prince may have issues with generating proper trajectories, leading to glitches in event finding.
    pub fn dp78(dynamics: D, opts: PropOpts<E>) -> Self {
        Self::new::<Dormand78>(dynamics, opts)
    }

    pub fn with(&'a self, state: D::StateType) -> PropInstance<'a, D, E> {
        // Pre-allocate the k used in the propagator
        let mut k = Vec::with_capacity(self.stages + 1);
        for _ in 0..self.stages {
            k.push(OVector::<f64, <D::StateType as State>::VecLength>::zeros());
        }
        PropInstance {
            state,
            prop: self,
            details: IntegrationDetails {
                step: self.opts.init_step,
                error: 0.0,
                attempts: 1,
            },
            step_size: self.opts.init_step,
            fixed_step: self.opts.fixed_step,
            k,
        }
    }
}

impl<'a, D: Dynamics> Propagator<'a, D, RSSCartesianStep>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<usize, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>,
{
    /// Default propagator is an RK89 with the default PropOpts.
    pub fn default(dynamics: D) -> Self {
        Self::new::<RK89>(dynamics, PropOpts::default())
    }

    /// A default Dormand Prince 78 propagator with the default PropOpts.
    /// Faster and more stable than an RK89 (`default`) but seems to cause issues for event finding.
    /// WARNINGS: Dormand Prince may have issues with generating proper trajectories, leading to glitches in event finding.
    pub fn default_dp78(dynamics: D) -> Self {
        Self::new::<Dormand78>(dynamics, PropOpts::default())
    }
}

/// A Propagator allows propagating a set of dynamics forward or backward in time.
/// It is an EventTracker, without any event tracking. It includes the options, the integrator
/// details of the previous step, and the set of coefficients used for the monomorphic instance.
#[derive(Debug)]
pub struct PropInstance<'a, D: Dynamics, E: ErrorCtrl>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<usize, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>,
{
    /// The state of this propagator instance
    pub state: D::StateType,
    /// The propagator setup (kind, stages, etc.)
    pub prop: &'a Propagator<'a, D, E>,
    /// Stores the details of the previous integration step
    pub details: IntegrationDetails,
    step_size: Duration, // Stores the adapted step for the _next_ call
    fixed_step: bool,
    // Allows us to do pre-allocation of the ki vectors
    k: Vec<OVector<f64, <D::StateType as State>::VecLength>>,
}

impl<'a, D: Dynamics, E: ErrorCtrl> PropInstance<'a, D, E>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<usize, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>,
{
    /// Allows setting the step size of the propagator
    pub fn set_step(&mut self, step_size: Duration, fixed: bool) {
        self.step_size = step_size;
        self.fixed_step = fixed;
    }

    #[allow(clippy::erasing_op)]
    fn for_duration_channel_option(
        &mut self,
        duration: Duration,
        maybe_tx_chan: Option<Sender<D::StateType>>,
    ) -> Result<D::StateType, NyxError> {
        if duration == 0 * Unit::Second {
            return Ok(self.state);
        }
        let stop_time = self.state.epoch() + duration;
        if duration > 2 * Unit::Minute || duration < -2 * Unit::Minute {
            // Prevent the print spam for orbit determination cases
            info!("Propagating for {} until {}", duration, stop_time);
        }
        // Call `finally` on the current state to set anything up
        self.state = self.prop.dynamics.finally(self.state)?;

        let backprop = duration < Unit::Nanosecond;
        if backprop {
            self.step_size = -self.step_size; // Invert the step size
        }
        loop {
            let dt = self.state.epoch();
            if (!backprop && dt + self.step_size > stop_time)
                || (backprop && dt + self.step_size <= stop_time)
            {
                if stop_time == dt {
                    // No propagation necessary
                    return Ok(self.state);
                }
                // Take one final step of exactly the needed duration until the stop time
                let prev_step_size = self.step_size;
                let prev_step_kind = self.fixed_step;
                self.set_step(stop_time - dt, true);

                self.single_step()?;

                // Publish to channel if provided
                if let Some(ref chan) = maybe_tx_chan {
                    if let Err(e) = chan.send(self.state) {
                        warn!("could not publish to channel: {}", e)
                    }
                }

                // Restore the step size for subsequent calls
                self.set_step(prev_step_size, prev_step_kind);
                if backprop {
                    self.step_size = -self.step_size; // Restore to a positive step size
                }
                return Ok(self.state);
            } else {
                self.single_step()?;
                // Publish to channel if provided
                if let Some(ref chan) = maybe_tx_chan {
                    if let Err(e) = chan.send(self.state) {
                        warn!("could not publish to channel: {}", e)
                    }
                }
            }
        }
    }

    /// This method propagates the provided Dynamics for the provided duration.
    pub fn for_duration(&mut self, duration: Duration) -> Result<D::StateType, NyxError> {
        self.for_duration_channel_option(duration, None)
    }

    /// This method propagates the provided Dynamics for the provided duration and publishes each state on the channel.
    pub fn for_duration_with_channel(
        &mut self,
        duration: Duration,
        tx_chan: Sender<D::StateType>,
    ) -> Result<D::StateType, NyxError> {
        self.for_duration_channel_option(duration, Some(tx_chan))
    }

    /// Propagates the provided Dynamics until the provided epoch. Returns the end state.
    pub fn until_epoch(&mut self, end_time: Epoch) -> Result<D::StateType, NyxError> {
        let duration: Duration = end_time - self.state.epoch();
        self.for_duration(duration)
    }

    /// Propagates the provided Dynamics until the provided epoch and publishes states on the provided channel. Returns the end state.
    pub fn until_epoch_with_channel(
        &mut self,
        end_time: Epoch,
        tx_chan: Sender<D::StateType>,
    ) -> Result<D::StateType, NyxError> {
        let duration: Duration = end_time - self.state.epoch();
        self.for_duration_with_channel(duration, tx_chan)
    }

    /// Propagates the provided Dynamics for the provided duration and generate the trajectory of these dynamics on its own thread.
    /// Returns the end state and the trajectory.
    /// Known bug #190: Cannot generate a valid trajectory when propagating backward
    #[allow(clippy::map_clone)]
    pub fn for_duration_with_traj(
        &mut self,
        duration: Duration,
    ) -> Result<(D::StateType, Traj<D::StateType>), NyxError>
    where
        <DefaultAllocator as Allocator<f64, <D::StateType as State>::VecLength>>::Buffer: Send,
        D::StateType: InterpState,
    {
        let start_state = self.state;
        let end_state;

        let rx = {
            // Channels that have the states in a bucket of the correct length
            let (tx_bucket, rx_bucket) = channel();

            let rx = {
                // Channels that have a single state for the propagator
                let (tx, rx) = channel();
                // Propagate the dynamics
                end_state = self.for_duration_with_channel(duration, tx)?;
                rx
            };

            /* *** */
            /* Map: bucket the states and send on a channel */
            /* *** */

            let items_per_segments = INTERPOLATION_SAMPLES;
            let mut window_states = Vec::with_capacity(2 * items_per_segments);
            // Push the initial state
            window_states.push(start_state);

            // Note that we're using the typical map+reduce pattern
            // Start receiving states on a blocking call
            while let Ok(state) = rx.recv() {
                window_states.push(state);
                if window_states.len() == 2 * items_per_segments {
                    // Publish the first items
                    let this_wdn = window_states[..items_per_segments]
                        .iter()
                        .map(|&x| x)
                        .collect::<Vec<D::StateType>>();

                    tx_bucket.send(this_wdn).map_err(|_| {
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
            // If there aren't enough states, set the propagator step size to make sure there is at least that many states
            if window_states.len() < items_per_segments {
                let step_size =
                    (end_state.epoch() - start_state.epoch()) / ((items_per_segments - 1) as f64);

                self.state = start_state;
                window_states.clear();
                self.set_step(step_size, true);
                let rx = {
                    // Channels that have a single state for the propagator
                    let (tx, rx) = channel();
                    // Propagate the dynamics
                    self.for_duration_with_channel(duration, tx)?;
                    rx
                };
                window_states.push(start_state);
                while let Ok(state) = rx.recv() {
                    window_states.push(state);
                }
            }
            // And interpolate the remaining states too, even if the buffer is not full!
            let mut start_idx = 0;
            loop {
                tx_bucket
                    .send(
                        window_states
                            [start_idx..(start_idx + items_per_segments).min(window_states.len())]
                            .iter()
                            .map(|&x| x)
                            .collect::<Vec<D::StateType>>(),
                    )
                    .map_err(|_| {
                        NyxError::from(TrajError::CreationError(
                            "could not send onto channel".to_string(),
                        ))
                    })?;
                if start_idx > 0 || window_states.len() < items_per_segments {
                    break;
                }
                start_idx = window_states.len() - items_per_segments;
                if start_idx == 0 {
                    // This means that the window states are exactly the correct size, break here
                    break;
                }
            }

            // Return the rx channel for these buckets
            rx_bucket
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
            traj.append_spline(spline)?;
        }

        Ok((end_state, traj))
    }

    /// Propagates the provided Dynamics until the provided epoch and generate the trajectory of these dynamics on its own thread.
    /// Returns the end state and the trajectory.
    /// Known bug #190: Cannot generate a valid trajectory when propagating backward
    pub fn until_epoch_with_traj(
        &mut self,
        end_time: Epoch,
    ) -> Result<(D::StateType, Traj<D::StateType>), NyxError>
    where
        <DefaultAllocator as Allocator<f64, <D::StateType as State>::VecLength>>::Buffer: Send,
        D::StateType: InterpState,
    {
        let duration: Duration = end_time - self.state.epoch();
        self.for_duration_with_traj(duration)
    }

    /// Propagate until a specific event is found once.
    /// Returns the state found and the trajectory until `max_duration`
    pub fn until_event<F: EventEvaluator<D::StateType>>(
        &mut self,
        max_duration: Duration,
        event: &F,
    ) -> Result<(D::StateType, Traj<D::StateType>), NyxError>
    where
        <DefaultAllocator as Allocator<f64, <D::StateType as State>::VecLength>>::Buffer: Send,
        D::StateType: InterpState,
    {
        self.until_nth_event(max_duration, event, 0)
    }

    /// Propagate until a specific event is found `trigger` times.
    /// Returns the state found and the trajectory until `max_duration`
    pub fn until_nth_event<F: EventEvaluator<D::StateType>>(
        &mut self,
        max_duration: Duration,
        event: &F,
        trigger: usize,
    ) -> Result<(D::StateType, Traj<D::StateType>), NyxError>
    where
        <DefaultAllocator as Allocator<f64, <D::StateType as State>::VecLength>>::Buffer: Send,
        D::StateType: InterpState,
    {
        info!("Searching for {}", event);

        let (_, traj) = self.for_duration_with_traj(max_duration)?;
        // Now, find the requested event
        let events = traj.find_all(event)?;
        match events.get(trigger) {
            Some(event_state) => Ok((*event_state, traj)),
            None => Err(NyxError::UnsufficientTriggers(trigger, events.len())),
        }
    }

    /// Take a single propagator step and emit the result on the TX channel (if enabled)
    pub fn single_step(&mut self) -> Result<(), NyxError> {
        let (t, state_vec) = self.derive()?;
        self.state.set(self.state.epoch() + t, &state_vec)?;
        self.state = self.prop.dynamics.finally(self.state)?;

        Ok(())
    }

    /// This method integrates whichever function is provided as `d_xdt`. Everything passed to this function is in **seconds**.
    ///
    /// This function returns the step sized used (as a Duration) and the new state as y_{n+1} = y_n + \frac{dy_n}{dt}.
    /// To get the integration details, check `self.latest_details`.
    fn derive(
        &mut self,
    ) -> Result<(Duration, OVector<f64, <D::StateType as State>::VecLength>), NyxError> {
        let state = &self.state.as_vector()?;
        let ctx = &self.state;
        // Reset the number of attempts used (we don't reset the error because it's set before it's read)
        self.details.attempts = 1;
        // Convert the step size to seconds -- it's mutable because we may change it below
        let mut step_size = self.step_size.in_seconds();
        loop {
            let ki = self.prop.dynamics.eom(0.0, state, ctx)?;
            self.k[0] = ki;
            let mut a_idx: usize = 0;
            for i in 0..(self.prop.stages - 1) {
                // Let's compute the c_i by summing the relevant items from the list of coefficients.
                // \sum_{j=1}^{i-1} a_ij  ∀ i ∈ [2, s]
                let mut ci: f64 = 0.0;
                // The wi stores the a_{s1} * k_1 + a_{s2} * k_2 + ... + a_{s, s-1} * k_{s-1} +
                let mut wi = OVector::<f64, <D::StateType as State>::VecLength>::from_element(0.0);
                for kj in &self.k[0..i + 1] {
                    let a_ij = self.prop.a_coeffs[a_idx];
                    ci += a_ij;
                    wi += a_ij * kj;
                    a_idx += 1;
                }

                let ki = self
                    .prop
                    .dynamics
                    .eom(ci * step_size, &(state + step_size * wi), ctx)?;
                self.k[i + 1] = ki;
            }
            // Compute the next state and the error
            let mut next_state = state.clone();
            // State error estimation from https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Adaptive_Runge%E2%80%93Kutta_methods
            // This is consistent with GMAT https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/propagator/RungeKutta.cpp#L537
            let mut error_est =
                OVector::<f64, <D::StateType as State>::VecLength>::from_element(0.0);
            for (i, ki) in self.k.iter().enumerate() {
                let b_i = self.prop.b_coeffs[i];
                if !self.fixed_step {
                    let b_i_star = self.prop.b_coeffs[i + self.prop.stages];
                    error_est += step_size * (b_i - b_i_star) * ki;
                }
                next_state += step_size * b_i * ki;
            }

            if self.fixed_step {
                // Using a fixed step, no adaptive step necessary
                self.details.step = self.step_size;
                return Ok(((self.details.step), next_state));
            } else {
                // Compute the error estimate.
                self.details.error = E::estimate(&error_est, &next_state, state);
                if self.details.error <= self.prop.opts.tolerance
                    || step_size <= self.prop.opts.min_step.in_seconds()
                    || self.details.attempts >= self.prop.opts.attempts
                {
                    if self.details.attempts >= self.prop.opts.attempts {
                        warn!(
                            "Could not further decrease step size: maximum number of attempts reached ({})",
                            self.details.attempts
                        );
                    }

                    self.details.step = step_size * Unit::Second;
                    if self.details.error < self.prop.opts.tolerance {
                        // Let's increase the step size for the next iteration.
                        // Error is less than tolerance, let's attempt to increase the step for the next iteration.
                        let proposed_step = 0.9
                            * step_size
                            * (self.prop.opts.tolerance / self.details.error)
                                .powf(1.0 / f64::from(self.prop.order));
                        step_size = if proposed_step > self.prop.opts.max_step.in_seconds() {
                            self.prop.opts.max_step.in_seconds()
                        } else {
                            proposed_step
                        };
                    }
                    // In all cases, let's update the step size to whatever was the adapted step size
                    self.step_size = step_size * Unit::Second;
                    return Ok((self.details.step, next_state));
                } else {
                    // Error is too high and we aren't using the smallest step, and we haven't hit the max number of attempts.
                    // So let's adapt the step size.
                    self.details.attempts += 1;
                    let proposed_step = 0.9
                        * step_size
                        * (self.prop.opts.tolerance / self.details.error)
                            .powf(1.0 / f64::from(self.prop.order - 1));
                    step_size = if proposed_step < self.prop.opts.min_step.in_seconds() {
                        self.prop.opts.min_step.in_seconds()
                    } else {
                        proposed_step
                    };
                    // Note that we don't set self.step_size, that will be updated right before we return
                }
            }
        }
    }

    /// Borrow the details of the latest integration step.
    pub fn latest_details(&self) -> &IntegrationDetails {
        &self.details
    }
}

/// PropOpts stores the integrator options, including the minimum and maximum step sizes, and the
/// max error size.
///
/// Note that different step sizes and max errors are only used for adaptive
/// methods. To use a fixed step integrator, initialize the options using `with_fixed_step`, and
/// use whichever adaptive step integrator is desired.  For example, initializing an RK45 with
/// fixed step options will lead to an RK4 being used instead of an RK45.
#[derive(Clone, Copy, Debug)]
pub struct PropOpts<E: ErrorCtrl> {
    pub init_step: Duration,
    pub min_step: Duration,
    pub max_step: Duration,
    pub tolerance: f64,
    pub attempts: u8,
    pub fixed_step: bool,
    pub _errctrl: E,
}

impl<E: ErrorCtrl> PropOpts<E> {
    /// `with_adaptive_step` initializes an `PropOpts` such that the integrator is used with an
    ///  adaptive step size. The number of attempts is currently fixed to 50 (as in GMAT).
    pub fn with_adaptive_step(
        min_step: Duration,
        max_step: Duration,
        tolerance: f64,
        errctrl: E,
    ) -> Self {
        PropOpts {
            init_step: max_step,
            min_step,
            max_step,
            tolerance,
            attempts: 50,
            fixed_step: false,
            _errctrl: errctrl,
        }
    }

    pub fn with_adaptive_step_s(min_step: f64, max_step: f64, tolerance: f64, errctrl: E) -> Self {
        Self::with_adaptive_step(
            min_step * Unit::Second,
            max_step * Unit::Second,
            tolerance,
            errctrl,
        )
    }

    /// Returns a string with the information about these options
    pub fn info(&self) -> String {
        format!(
            "[min_step: {:.e}, max_step: {:.e}, tol: {:.e}, attempts: {}]",
            self.min_step, self.max_step, self.tolerance, self.attempts,
        )
    }

    /// Set the maximum step size and sets the initial step to that value if currently greater
    pub fn set_max_step(&mut self, max_step: Duration) {
        if self.init_step > max_step {
            self.init_step = max_step;
        }
        self.max_step = max_step;
    }

    /// Set the minimum step size and sets the initial step to that value if currently smaller
    pub fn set_min_step(&mut self, min_step: Duration) {
        if self.init_step < min_step {
            self.init_step = min_step;
        }
        self.min_step = min_step;
    }
}

impl PropOpts<RSSCartesianStep> {
    /// `with_fixed_step` initializes an `PropOpts` such that the integrator is used with a fixed
    ///  step size.
    pub fn with_fixed_step(step: Duration) -> Self {
        PropOpts {
            init_step: step,
            min_step: step,
            max_step: step,
            tolerance: 0.0,
            fixed_step: true,
            attempts: 0,
            _errctrl: RSSCartesianStep {},
        }
    }

    pub fn with_fixed_step_s(step: f64) -> Self {
        Self::with_fixed_step(step * Unit::Second)
    }

    /// Returns the default options with a specific tolerance.
    #[allow(clippy::field_reassign_with_default)]
    pub fn with_tolerance(tolerance: f64) -> Self {
        let mut opts = Self::default();
        opts.tolerance = tolerance;
        opts
    }

    /// Creates a propagator with the provided max step, and sets the initial step to that value as well.
    #[allow(clippy::field_reassign_with_default)]
    pub fn with_max_step(max_step: Duration) -> Self {
        let mut opts = Self::default();
        opts.set_max_step(max_step);
        opts
    }
}

impl Default for PropOpts<RSSCartesianStep> {
    /// `default` returns the same default options as GMAT.
    fn default() -> PropOpts<RSSCartesianStep> {
        PropOpts {
            init_step: 60.0 * Unit::Second,
            min_step: 0.001 * Unit::Second,
            max_step: 2700.0 * Unit::Second,
            tolerance: 1e-12,
            attempts: 50,
            fixed_step: false,
            _errctrl: RSSCartesianStep {},
        }
    }
}

#[test]
fn test_options() {
    use super::error_ctrl::RSSStep;

    let opts = PropOpts::with_fixed_step_s(1e-1);
    assert_eq!(opts.min_step, 1e-1 * Unit::Second);
    assert_eq!(opts.max_step, 1e-1 * Unit::Second);
    assert!(opts.tolerance.abs() < f64::EPSILON);
    assert!(opts.fixed_step);

    let opts = PropOpts::with_adaptive_step_s(1e-2, 10.0, 1e-12, RSSStep {});
    assert_eq!(opts.min_step, 1e-2 * Unit::Second);
    assert_eq!(opts.max_step, 10.0 * Unit::Second);
    assert!((opts.tolerance - 1e-12).abs() < f64::EPSILON);
    assert!(!opts.fixed_step);

    let opts: PropOpts<RSSCartesianStep> = Default::default();
    assert_eq!(opts.init_step, 60.0 * Unit::Second);
    assert_eq!(opts.min_step, 0.001 * Unit::Second);
    assert_eq!(opts.max_step, 2700.0 * Unit::Second);
    assert!((opts.tolerance - 1e-12).abs() < f64::EPSILON);
    assert_eq!(opts.attempts, 50);
    assert!(!opts.fixed_step);

    let opts = PropOpts::with_max_step(1.0 * Unit::Second);
    assert_eq!(opts.init_step, 1.0 * Unit::Second);
    assert_eq!(opts.min_step, 0.001 * Unit::Second);
    assert_eq!(opts.max_step, 1.0 * Unit::Second);
    assert!((opts.tolerance - 1e-12).abs() < f64::EPSILON);
    assert_eq!(opts.attempts, 50);
    assert!(!opts.fixed_step);
}
