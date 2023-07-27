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

use super::error_ctrl::ErrorCtrl;
use super::{IntegrationDetails, Propagator};
use crate::dynamics::Dynamics;
use crate::errors::NyxError;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, OVector};
use crate::md::trajectory::{Interpolatable, Traj};
use crate::md::EventEvaluator;
use crate::time::{Duration, Epoch, Unit};
use crate::State;
use rayon::iter::ParallelBridge;
use rayon::prelude::ParallelIterator;
use std::f64;
use std::sync::mpsc::{channel, Sender};
#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant;

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
    pub(crate) step_size: Duration, // Stores the adapted step for the _next_ call
    pub(crate) fixed_step: bool,
    // Allows us to do pre-allocation of the ki vectors
    pub(crate) k: Vec<OVector<f64, <D::StateType as State>::VecLength>>,
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

        #[cfg(not(target_arch = "wasm32"))]
        let tick = Instant::now();
        let log_progress = duration.abs() >= 2 * Unit::Minute;

        if log_progress {
            // Prevent the print spam for orbit determination cases
            info!("Propagating for {} until {}", duration, stop_time);
        }
        // Call `finally` on the current state to set anything up
        self.state = self.prop.dynamics.finally(self.state)?;

        let backprop = duration.is_negative();
        if backprop {
            self.step_size = -self.step_size; // Invert the step size
        }
        loop {
            let epoch = self.state.epoch();
            if (!backprop && epoch + self.step_size > stop_time)
                || (backprop && epoch + self.step_size <= stop_time)
            {
                if stop_time == epoch {
                    // No propagation necessary
                    #[cfg(not(target_arch = "wasm32"))]
                    {
                        if log_progress {
                            let tock: Duration = tick.elapsed().into();
                            info!("Done in {}", tock);
                        }
                    }
                    return Ok(self.state);
                }
                // Take one final step of exactly the needed duration until the stop time
                let prev_step_size = self.step_size;
                let prev_step_kind = self.fixed_step;
                self.set_step(stop_time - epoch, true);

                self.single_step()?;

                // Publish to channel if provided
                if let Some(ref chan) = maybe_tx_chan {
                    if let Err(e) = chan.send(self.state) {
                        warn!("{} when sending on channel", e)
                    }
                }

                // Restore the step size for subsequent calls
                self.set_step(prev_step_size, prev_step_kind);

                if backprop {
                    self.step_size = -self.step_size; // Restore to a positive step size
                }

                #[cfg(not(target_arch = "wasm32"))]
                {
                    if log_progress {
                        let tock: Duration = tick.elapsed().into();
                        info!("Done in {}", tock);
                    }
                }

                return Ok(self.state);
            } else {
                self.single_step()?;
                // Publish to channel if provided
                if let Some(ref chan) = maybe_tx_chan {
                    if let Err(e) = chan.send(self.state) {
                        warn!("{} when sending on channel", e)
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
    #[allow(clippy::map_clone)]
    pub fn for_duration_with_traj(
        &mut self,
        duration: Duration,
    ) -> Result<(D::StateType, Traj<D::StateType>), NyxError>
    where
        <DefaultAllocator as Allocator<f64, <D::StateType as State>::VecLength>>::Buffer: Send,
        D::StateType: Interpolatable,
    {
        let end_state;
        let mut traj = Traj::new();
        let start_state = self.state;

        let rx = {
            // Channels that have a single state for the propagator
            let (tx, rx) = channel();
            // Propagate the dynamics
            // Note that the end state is also sent on the channel before the return of this function.
            end_state = self.for_duration_with_channel(duration, tx)?;
            rx
        };

        traj.states = rx.into_iter().par_bridge().collect();
        // Push the start state -- will be reordered in the finalize call.
        // For some reason, this must happen at the end -- can't figure out why.
        traj.states.push(start_state);

        traj.finalize();

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
        D::StateType: Interpolatable,
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
        D::StateType: Interpolatable,
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
        D::StateType: Interpolatable,
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
        let state_vec = &self.state.as_vector()?;
        let state_ctx = &self.state;
        // Reset the number of attempts used (we don't reset the error because it's set before it's read)
        self.details.attempts = 1;
        // Convert the step size to seconds -- it's mutable because we may change it below
        let mut step_size = self.step_size.to_seconds();
        loop {
            let ki = self.prop.dynamics.eom(0.0, state_vec, state_ctx)?;
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

                let ki = self.prop.dynamics.eom(
                    ci * step_size,
                    &(state_vec + step_size * wi),
                    state_ctx,
                )?;
                self.k[i + 1] = ki;
            }
            // Compute the next state and the error
            let mut next_state = state_vec.clone();
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
                self.details.error = E::estimate(&error_est, &next_state, state_vec);
                if self.details.error <= self.prop.opts.tolerance
                    || step_size <= self.prop.opts.min_step.to_seconds()
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
                        step_size = if proposed_step > self.prop.opts.max_step.to_seconds() {
                            self.prop.opts.max_step.to_seconds()
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
                    step_size = if proposed_step < self.prop.opts.min_step.to_seconds() {
                        self.prop.opts.min_step.to_seconds()
                    } else {
                        proposed_step
                    };
                    // Note that we don't set self.step_size, that will be updated right before we return
                }
            }
        }
    }

    /// Copy the details of the latest integration step.
    pub fn latest_details(&self) -> IntegrationDetails {
        self.details
    }
}
