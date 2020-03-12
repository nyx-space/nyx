extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, VectorN};
use super::error_ctrl::{ErrorCtrl, RSSStepPV};
use super::events::{ConvergenceError, EventTrackers, StopCondition};
use super::{IntegrationDetails, RK, RK89};
use dynamics::Dynamics;
use std::f64;
use std::f64::EPSILON;
use std::sync::mpsc::Sender;

/// A Propagator allows propagating a set of dynamics forward or backward in time.
/// It is an EventTracker, without any event tracking. It includes the options, the integrator
/// details of the previous step, and the set of coefficients used for the monomorphic instance.
#[derive(Debug)]
pub struct Propagator<'a, M: Dynamics, E: ErrorCtrl>
where
    DefaultAllocator: Allocator<f64, M::StateSize>,
{
    pub dynamics: &'a mut M, // Stores the dynamics used. *Must* use this to get the latest values
    // An output channel for all of the states computed by this propagator
    pub tx_chan: Option<&'a Sender<M::StateType>>,
    pub event_trackers: EventTrackers<M::StateType>,
    opts: PropOpts<E>, // Stores the integration options (tolerance, min/max step, init step, etc.)
    details: IntegrationDetails, // Stores the details of the previous integration step
    step_size: f64,    // Stores the adapted step for the _next_ call
    order: u8,         // Order of the integrator
    stages: usize,     // Number of stages, i.e. how many times the derivatives will be called
    a_coeffs: &'a [f64],
    b_coeffs: &'a [f64],
    fixed_step: bool,
}

/// The `Propagator` trait defines the functions of a propagator and of an event tracker.
impl<'a, M: Dynamics, E: ErrorCtrl> Propagator<'a, M, E>
where
    DefaultAllocator: Allocator<f64, M::StateSize>,
{
    /// Each propagator must be initialized with `new` which stores propagator information.
    pub fn new<T: RK>(dynamics: &'a mut M, opts: &PropOpts<E>) -> Self {
        Self {
            tx_chan: None,
            dynamics,
            event_trackers: EventTrackers::none(),
            opts: *opts,
            details: IntegrationDetails {
                step: 0.0,
                error: 0.0,
                attempts: 1,
            },
            step_size: opts.init_step,
            stages: T::stages(),
            order: T::order(),
            a_coeffs: T::a_coeffs(),
            b_coeffs: T::b_coeffs(),
            fixed_step: opts.fixed_step,
        }
    }

    /// Default propagator is an RK89.
    pub fn default(dynamics: &'a mut M, opts: &PropOpts<E>) -> Self {
        Self::new::<RK89>(dynamics, opts)
    }

    pub fn set_step(&mut self, step_size: f64, fixed: bool) {
        self.step_size = step_size;
        self.fixed_step = fixed;
    }

    /// Returns the time of the propagation
    ///
    /// WARNING: Do not use the dynamics to get the time, it will be the initial value!
    pub fn time(&self) -> f64 {
        self.dynamics.time()
    }

    /// Returns the state of the propagation
    ///
    /// WARNING: Do not use the dynamics to get the state, it will be the initial value!
    pub fn state_vector(&self) -> VectorN<f64, M::StateSize> {
        self.dynamics.state_vector()
    }

    /// This method propagates the provided Dynamics `dyn` for `elapsed_time` seconds. WARNING: This function has many caveats (please read detailed docs).
    ///
    /// ### IMPORTANT CAVEAT of `until_time_elapsed`
    /// - It is **assumed** that `self.dynamics.time()` returns a time in the same units as elapsed_time.
    pub fn until_time_elapsed(&mut self, elapsed_time: f64) -> M::StateType {
        let backprop = elapsed_time < 0.0;
        if backprop {
            self.step_size *= -1.0; // Invert the step size
        }
        let init_seconds = self.dynamics.time();
        let stop_time = init_seconds + elapsed_time;
        loop {
            let dt = self.dynamics.time();
            if (!backprop && dt + self.step_size > stop_time)
                || (backprop && dt + self.step_size <= stop_time)
            {
                if (stop_time - dt).abs() < f64::EPSILON {
                    // No propagation necessary
                    return self.dynamics.state();
                }
                // Take one final step of exactly the needed duration until the stop time
                let prev_step_size = self.step_size;
                let prev_step_kind = self.fixed_step;
                self.set_step(stop_time - dt, true);
                let (t, state) = self.derive(dt, &self.dynamics.state_vector());
                self.dynamics.set_state(t, &state);
                // Evaluate the event trackers
                self.event_trackers
                    .eval_and_save(dt, t, &self.dynamics.state());
                // Restore the step size for subsequent calls
                self.set_step(prev_step_size, prev_step_kind);
                if let Some(ref chan) = self.tx_chan {
                    if let Err(e) = chan.send(self.dynamics.state()) {
                        warn!("could not publish to channel: {}", e)
                    }
                }
                if backprop {
                    self.step_size *= -1.0; // Restore to a positive step size
                }
                return self.dynamics.state();
            } else {
                let (t, state) = self.derive(dt, &self.dynamics.state_vector());
                // We haven't passed the time based stopping condition.
                self.dynamics.set_state(t, &state);
                // Evaluate the event trackers
                self.event_trackers
                    .eval_and_save(dt, t, &self.dynamics.state());
                if let Some(ref chan) = self.tx_chan {
                    if let Err(e) = chan.send(self.dynamics.state()) {
                        warn!("could not publish to channel: {}", e)
                    }
                }
            }
        }
    }

    pub fn until_event(
        &mut self,
        condition: StopCondition<M::StateType>,
    ) -> Result<M::StateType, ConvergenceError> {
        // Store the initial time and state
        let init_time = self.dynamics.time();
        let init_state_vec = self.dynamics.state_vector();
        // Rewrite the event tracker
        if !self.event_trackers.events.is_empty() {
            warn!("Rewriting event tracker with the StopCondition");
        }
        self.event_trackers = EventTrackers::from_event(condition.event);
        self.until_time_elapsed(condition.max_prop_time);
        // Check if the event has been triggered
        if self.event_trackers.found_bounds[0].len() < condition.trigger {
            if condition.trigger == 1 {
                // Event was never triggered
                return Err(ConvergenceError::NeverTriggered);
            } else {
                // Event not triggered enough times
                return Err(ConvergenceError::UnsufficientTriggers(
                    condition.trigger,
                    self.event_trackers.found_bounds[0].len(),
                ));
            }
        }
        let (mut xa, mut xb) = self.event_trackers.found_bounds[0][condition.trigger - 1];
        // Reinitialize the dynamics and start the search
        self.dynamics.set_state(init_time, &init_state_vec);
        self.event_trackers.reset();
        // Compute the initial values of the condition at those bounds.
        self.until_time_elapsed(xa);
        let mut ya = self.event_trackers.events[0].eval(&self.dynamics.state());
        // And prop until next bound
        self.until_time_elapsed(xb - xa);
        let mut yb = self.event_trackers.events[0].eval(&self.dynamics.state());
        // The Brent solver, from the roots crate (sadly could not directly integrate it here)
        // Source: https://docs.rs/roots/0.0.5/src/roots/numerical/brent.rs.html#57-131

        // Helper lambdas, for f64s only
        let eps = condition.epsilon;
        let has_converged = |x1: f64, x2: f64| (x1 - x2).abs() < eps.abs();
        let arrange = |a: f64, ya: f64, b: f64, yb: f64| {
            if ya.abs() > yb.abs() {
                (a, ya, b, yb)
            } else {
                (b, yb, a, ya)
            }
        };

        let (mut c, mut yc, mut d) = (xa, ya, xa);
        let mut flag = true;
        let mut iter = 0;
        let closest_t;
        loop {
            if ya.abs() < condition.epsilon.abs() {
                closest_t = xa;
                break;
            }
            if yb.abs() < condition.epsilon.abs() {
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
            self.until_time_elapsed(s - self.dynamics.time());
            let ys = self.event_trackers.events[0].eval(&self.dynamics.state());
            d = c;
            c = xb;
            yc = yb;
            if ya * ys < 0.0 {
                // Root bracketed between a and s
                // Propagate until time xa
                self.until_time_elapsed(xa - self.dynamics.time());
                let ya_p = self.event_trackers.events[0].eval(&self.dynamics.state());
                match arrange(xa, ya_p, s, ys) {
                    (_a, _ya, _b, _yb) => {
                        xa = _a;
                        ya = _ya;
                        xb = _b;
                        yb = _yb;
                    }
                }
            } else {
                // Root bracketed between s and b
                // Propagate until time xb
                self.until_time_elapsed(xb - self.dynamics.time());
                let yb_p = self.event_trackers.events[0].eval(&self.dynamics.state());
                match arrange(s, ys, xb, yb_p) {
                    (_a, _ya, _b, _yb) => {
                        xa = _a;
                        ya = _ya;
                        xb = _b;
                        yb = _yb;
                    }
                }
            }
            iter += 1;
            if iter >= condition.max_iter {
                return Err(ConvergenceError::MaxIterReached(iter));
            }
        }

        // Now that we have the time at which the condition is matched, let's propagate until then
        self.until_time_elapsed(self.dynamics.time() - closest_t);

        Ok(self.dynamics.state())
    }

    /// This method integrates whichever function is provided as `d_xdt`.
    ///
    /// The `derive` method is monomorphic to increase speed. This function takes a time `t` and a current state `state`
    /// then derives the dynamics at that time (i.e. propagates for one time step). The `d_xdt` parameter is the derivative
    /// function which take a time t of type f64 and a reference to a state of type VectorN<f64, N>, and returns the
    /// result as VectorN<f64, N> of the derivative. The reference should preferrably only be borrowed.
    /// This function returns the next time (i.e. the previous time incremented by the timestep used) and
    /// the new state as y_{n+1} = y_n + \frac{dy_n}{dt}. To get the integration details, check `Self.latest_details`.
    /// Note: using VectorN<f64, N> instead of DVector implies that the function *must* always return a vector of the same
    /// size. This static allocation allows for high execution speeds.
    pub fn derive(
        &mut self,
        t: f64,
        state: &VectorN<f64, M::StateSize>,
    ) -> (f64, VectorN<f64, M::StateSize>) {
        // Reset the number of attempts used (we don't reset the error because it's set before it's read)
        self.details.attempts = 1;
        loop {
            let mut k = Vec::with_capacity(self.stages + 1); // Will store all the k_i.
            let ki = self.dynamics.eom(t, &state);
            k.push(ki);
            let mut a_idx: usize = 0;
            for _ in 0..(self.stages - 1) {
                // Let's compute the c_i by summing the relevant items from the list of coefficients.
                // \sum_{j=1}^{i-1} a_ij  ∀ i ∈ [2, s]
                let mut ci: f64 = 0.0;
                // The wi stores the a_{s1} * k_1 + a_{s2} * k_2 + ... + a_{s, s-1} * k_{s-1} +
                let mut wi = VectorN::<f64, M::StateSize>::from_element(0.0);
                for kj in &k {
                    let a_ij = self.a_coeffs[a_idx];
                    ci += a_ij;
                    wi += a_ij * kj;
                    a_idx += 1;
                }

                let ki = self
                    .dynamics
                    .eom(t + ci * self.step_size, &(state + self.step_size * wi));
                k.push(ki);
            }
            // Compute the next state and the error
            let mut next_state = state.clone();
            // State error estimation from https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods#Adaptive_Runge%E2%80%93Kutta_methods
            // This is consistent with GMAT https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/propagator/RungeKutta.cpp#L537
            let mut error_est = VectorN::<f64, M::StateSize>::from_element(0.0);
            for (i, ki) in k.iter().enumerate() {
                let b_i = self.b_coeffs[i];
                if !self.fixed_step {
                    let b_i_star = self.b_coeffs[i + self.stages];
                    error_est += self.step_size * (b_i - b_i_star) * ki;
                }
                next_state += self.step_size * b_i * ki;
            }

            if self.fixed_step {
                // Using a fixed step, no adaptive step necessary
                self.details.step = self.step_size;
                return ((t + self.details.step), next_state);
            } else {
                // Compute the error estimate.
                self.details.error = E::estimate(&error_est, &next_state.clone(), &state);
                if self.details.error <= self.opts.tolerance
                    || self.step_size <= self.opts.min_step
                    || self.details.attempts >= self.opts.attempts
                {
                    if self.details.attempts >= self.opts.attempts {
                        warn!(
                            "maximum number of attempts reached ({})",
                            self.details.attempts
                        );
                    }

                    self.details.step = self.step_size;
                    if self.details.error < self.opts.tolerance {
                        // Let's increase the step size for the next iteration.
                        // Error is less than tolerance, let's attempt to increase the step for the next iteration.
                        let proposed_step = 0.9
                            * self.step_size
                            * (self.opts.tolerance / self.details.error)
                                .powf(1.0 / f64::from(self.order));
                        self.step_size = if proposed_step > self.opts.max_step {
                            self.opts.max_step
                        } else {
                            proposed_step
                        };
                    }
                    return ((t + self.details.step), next_state);
                } else {
                    // Error is too high and we aren't using the smallest step, and we haven't hit the max number of attempts.
                    // So let's adapt the step size.
                    self.details.attempts += 1;
                    let proposed_step = 0.9
                        * self.step_size
                        * (self.opts.tolerance / self.details.error)
                            .powf(1.0 / f64::from(self.order - 1));
                    self.step_size = if proposed_step < self.opts.min_step {
                        self.opts.min_step
                    } else {
                        proposed_step
                    };
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
    init_step: f64,
    min_step: f64,
    max_step: f64,
    tolerance: f64,
    attempts: u8,
    fixed_step: bool,
    errctrl: E,
}

impl<E: ErrorCtrl> PropOpts<E> {
    /// `with_adaptive_step` initializes an `PropOpts` such that the integrator is used with an
    ///  adaptive step size. The number of attempts is currently fixed to 50 (as in GMAT).
    pub fn with_adaptive_step(
        min_step: f64,
        max_step: f64,
        tolerance: f64,
        errctrl: E,
    ) -> PropOpts<E> {
        PropOpts {
            init_step: max_step,
            min_step,
            max_step,
            tolerance,
            attempts: 50,
            fixed_step: false,
            errctrl,
        }
    }
}

impl PropOpts<RSSStepPV> {
    /// `with_fixed_step` initializes an `PropOpts` such that the integrator is used with a fixed
    ///  step size.
    pub fn with_fixed_step(step: f64) -> PropOpts<RSSStepPV> {
        PropOpts {
            init_step: step,
            min_step: step,
            max_step: step,
            tolerance: 0.0,
            fixed_step: true,
            attempts: 0,
            errctrl: RSSStepPV {},
        }
    }
}

impl Default for PropOpts<RSSStepPV> {
    /// `default` returns the same default options as GMAT.
    fn default() -> PropOpts<RSSStepPV> {
        PropOpts {
            init_step: 60.0,
            min_step: 0.001,
            max_step: 2700.0,
            tolerance: 1e-12,
            attempts: 50,
            fixed_step: false,
            errctrl: RSSStepPV {},
        }
    }
}

#[test]
fn test_options() {
    use std::f64::EPSILON;
    macro_rules! f64_eq {
        ($x:expr, $val:expr) => {
            assert!(($x - $val).abs() < EPSILON)
        };
    }

    use super::error_ctrl::RSSStep;

    let opts = PropOpts::with_fixed_step(1e-1);
    f64_eq!(opts.min_step, 1e-1);
    f64_eq!(opts.max_step, 1e-1);
    f64_eq!(opts.tolerance, 0.0);
    assert_eq!(opts.fixed_step, true);

    let opts = PropOpts::with_adaptive_step(1e-2, 10.0, 1e-12, RSSStep {});
    f64_eq!(opts.min_step, 1e-2);
    f64_eq!(opts.max_step, 10.0);
    f64_eq!(opts.tolerance, 1e-12);
    assert_eq!(opts.fixed_step, false);

    let opts: PropOpts<RSSStepPV> = Default::default();
    f64_eq!(opts.init_step, 60.0);
    f64_eq!(opts.min_step, 0.001);
    f64_eq!(opts.max_step, 2700.0);
    f64_eq!(opts.tolerance, 1e-12);
    assert_eq!(opts.attempts, 50);
    assert_eq!(opts.fixed_step, false);
}
