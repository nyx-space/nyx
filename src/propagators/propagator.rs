use super::error_ctrl::{ErrorCtrl, RSSStepPV};
use super::events::{EventTrackers, StopCondition};
use super::{IntegrationDetails, RK, RK89};
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, VectorN};
use crate::time::{Duration, Epoch, TimeUnit};
use crate::{State, TimeTagged};
use dynamics::Dynamics;
use errors::NyxError;
use std::f64;
use std::f64::EPSILON;
use std::sync::mpsc::Sender;

/// A Propagator allows propagating a set of dynamics forward or backward in time.
/// It is an EventTracker, without any event tracking. It includes the options, the integrator
/// details of the previous step, and the set of coefficients used for the monomorphic instance.
#[derive(Clone, Debug)]
pub struct Propagator<'a, D: Dynamics, E: ErrorCtrl>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::PropVecSize>,
{
    pub dynamics: &'a D, // Stores the dynamics used. *Must* use this to get the latest values
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
        + Allocator<f64, <D::StateType as State>::PropVecSize>,
{
    /// Each propagator must be initialized with `new` which stores propagator information.
    pub fn new<T: RK>(dynamics: &'a D, opts: PropOpts<E>) -> Self {
        Self {
            dynamics,
            opts,
            stages: T::stages(),
            order: T::order(),
            a_coeffs: T::a_coeffs(),
            b_coeffs: T::b_coeffs(),
        }
    }

    /// Set the tolerance for the propagator
    pub fn set_tolerance(&mut self, tol: f64) {
        self.opts.tolerance = tol;
    }

    /// Set the maximum step size for the propagator
    pub fn set_max_step(&mut self, step: Duration) {
        self.opts.max_step = step;
    }

    /// An RK89 propagator (the default) with custom propagator options.
    pub fn rk89(dynamics: &'a D, opts: PropOpts<E>) -> Self {
        Self::new::<RK89>(dynamics, opts)
    }

    pub fn with(&'a self, state: D::StateType) -> PropInstance<'a, D, E> {
        // let init_time = state.epoch();
        // let init_state_vec = dynamics.state_vector();
        // Pre-allocate the k used in the propagator
        let mut k = Vec::with_capacity(self.stages + 1);
        for _ in 0..self.stages {
            k.push(VectorN::<f64, <D::StateType as State>::PropVecSize>::zeros());
        }
        PropInstance {
            state,
            prop: &self,
            tx_chan: None,
            event_trackers: EventTrackers::none(),
            prevent_tx: false,
            details: IntegrationDetails {
                step: self.opts.init_step,
                error: 0.0,
                attempts: 1,
            },
            step_size: self.opts.init_step,
            fixed_step: self.opts.fixed_step,
            // init_time,
            k,
        }
    }
}

impl<'a, D: Dynamics> Propagator<'a, D, RSSStepPV>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::PropVecSize>,
{
    /// Default propagator is an RK89 with the default PropOpts.
    pub fn default(dynamics: &'a D) -> Self {
        Self::new::<RK89>(dynamics, PropOpts::default())
    }
}

/// A Propagator allows propagating a set of dynamics forward or backward in time.
/// It is an EventTracker, without any event tracking. It includes the options, the integrator
/// details of the previous step, and the set of coefficients used for the monomorphic instance.
#[derive(Debug)]
pub struct PropInstance<'a, D: Dynamics, E: ErrorCtrl>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::PropVecSize>,
{
    /// The state of this propagator instance
    pub state: D::StateType,
    /// The propagator setup (kind, stages, etc.)
    pub prop: &'a Propagator<'a, D, E>,
    /// An output channel for all of the states computed by this propagator instance
    pub tx_chan: Option<Sender<D::StateType>>,
    /// An event tracking instance
    pub event_trackers: EventTrackers<D::StateType>,
    /// Stores the details of the previous integration step
    pub details: IntegrationDetails,
    prevent_tx: bool, // Allows preventing publishing to channel even if channel is set
    step_size: Duration, // Stores the adapted step for the _next_ call
    fixed_step: bool,
    // init_time: Epoch,
    // init_state_vec: VectorN<f64, <D::StateType as State>::Size>,
    // Allows us to do pre-allocation of the ki vectors
    k: Vec<VectorN<f64, <D::StateType as State>::PropVecSize>>,
}

impl<'a, D: Dynamics, E: ErrorCtrl> PropInstance<'a, D, E>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::PropVecSize>
        + Allocator<f64, <D::StateType as State>::Size>,
{
    /// Allows setting the step size of the propagator
    pub fn set_step(&mut self, step_size: Duration, fixed: bool) {
        self.step_size = step_size;
        self.fixed_step = fixed;
    }

    /// Returns the time of the propagation
    ///
    /// WARNING: Do not use the dynamics to get the time, it will be the initial value!
    pub fn time(&self) -> f64 {
        self.state.epoch().as_tai_seconds()
    }

    /// Returns the state of the propagation
    ///
    /// WARNING: Do not use the dynamics to get the state, it will be the initial value!
    pub fn state_vector(&self) -> VectorN<f64, <D::StateType as State>::PropVecSize> {
        self.state.as_vector().unwrap()
    }

    /// A shortcut to dynamics.state()
    /// TODO: Rremove this
    pub fn state(&self) -> D::StateType {
        self.state
    }

    /// This method propagates the provided Dynamics `dyn` for `elapsed_time` seconds. WARNING: This function has many caveats (please read detailed docs).
    /// TODO: Rename this to `for` such that we have `prop.with(state).for(5 * TimeUnit::Hour)`.
    ///
    /// ### IMPORTANT CAVEAT of `until_time_elapsed`
    /// - It is **assumed** that `self.dynamics.time()` returns a time in the same units as elapsed_time.
    pub fn until_time_elapsed(&mut self, elapsed_time: Duration) -> Result<D::StateType, NyxError> {
        let backprop = elapsed_time < TimeUnit::Nanosecond;
        if backprop {
            self.step_size = -self.step_size; // Invert the step size
        }
        // let init_seconds = self.dynamics.time();
        // let stop_time = init_seconds + elapsed_time.in_seconds();
        let stop_time = self.state.epoch() + elapsed_time;
        // let stop_time_s = stop_time.as_tai_seconds();
        loop {
            let dt = self.state.epoch();
            // let dt_s = self.state.epoch().as_tai_seconds();
            if (!backprop && dt + self.step_size > stop_time)
                || (backprop && dt + self.step_size <= stop_time)
            /*if (!backprop && dt_s + self.step_size.in_seconds() > stop_time_s)
            || (backprop && dt_s + self.step_size.in_seconds() <= stop_time_s)*/
            {
                if stop_time == dt {
                    // No propagation necessary
                    return Ok(self.state);
                }
                // Take one final step of exactly the needed duration until the stop time
                let prev_step_size = self.step_size;
                let prev_step_kind = self.fixed_step;
                self.set_step(dbg!(stop_time - dt), true);
                // let state_vector = &self.state_vector();
                // let state = &self.state;
                // let (t, state_vec) = self.derive(dt.as_tai_seconds(), state_vector, state)?;
                let (t, state_vec) = self.derive(dt.as_tai_seconds())?;
                self.state.set(Epoch::from_tai_seconds(t), &state_vec)?;
                // Evaluate the event trackers
                self.event_trackers
                    .eval_and_save(dt, self.state.epoch(), &self.state);
                // Restore the step size for subsequent calls
                self.set_step(prev_step_size, prev_step_kind);
                if !self.prevent_tx {
                    if let Some(ref chan) = self.tx_chan {
                        if let Err(e) = chan.send(self.state) {
                            warn!("could not publish to channel: {}", e)
                        }
                    }
                }
                if backprop {
                    self.step_size = -self.step_size; // Restore to a positive step size
                }
                return Ok(self.state);
            } else {
                let (t, state_vec) = self.derive(dt.as_tai_seconds())?;
                // let (t, state_vec) =
                //     self.derive(dt.as_tai_seconds(), &self.state_vector(), &self.state)?;
                // We haven't passed the time based stopping condition.
                self.state.set(Epoch::from_tai_seconds(t), &state_vec)?;
                // Evaluate the event trackers
                self.event_trackers
                    .eval_and_save(dt, self.state.epoch(), &self.state);
                if !self.prevent_tx {
                    if let Some(ref chan) = self.tx_chan {
                        if let Err(e) = chan.send(self.state) {
                            warn!("could not publish to channel: {}", e)
                        }
                    }
                }
            }
        }
    }

    pub fn until_event(
        &mut self,
        condition: StopCondition<D::StateType>,
    ) -> Result<D::StateType, NyxError> {
        let init_state = self.state;
        self.prevent_tx = true;
        // Rewrite the event tracker
        if !self.event_trackers.events.is_empty() {
            warn!("Rewriting event tracker with the StopCondition");
        }

        self.event_trackers = EventTrackers::from_event(condition.event);
        self.until_time_elapsed(condition.max_prop_time)?;
        // Check if the event has been triggered
        if self.event_trackers.found_bounds[0].len() < condition.trigger {
            if condition.trigger == 1 {
                // Event was never triggered
                return Err(NyxError::ConditionNeverTriggered);
            } else {
                // Event not triggered enough times
                return Err(NyxError::UnsufficientTriggers(
                    condition.trigger,
                    self.event_trackers.found_bounds[0].len(),
                ));
            }
        }
        // TODO: Convert all of the epochs into TAI seconds, and do the math that way
        let (xa_e, xb_e) = self.event_trackers.found_bounds[0][condition.trigger - 1];
        let mut xa = xa_e.as_tai_seconds();
        let mut xb = xb_e.as_tai_seconds();
        // Reinitialize the dynamics and start the search
        self.state = init_state;
        self.event_trackers.reset();
        let initial_epoch = self.state.epoch();
        // Compute the initial values of the condition at those bounds.
        self.until_time_elapsed(xa_e - initial_epoch)?;
        let mut ya = self.event_trackers.events[0].eval(&self.state);
        // And prop until next bound
        self.until_time_elapsed((xb - xa) * TimeUnit::Second)?;
        let mut yb = self.event_trackers.events[0].eval(&self.state);
        // The Brent solver, from the roots crate (sadly could not directly integrate it here)
        // Source: https://docs.rs/roots/0.0.5/src/roots/numerical/brent.rs.html#57-131

        // Helper lambdas, for f64s only
        let eps = condition.epsilon;
        let has_converged = |x1: f64, x2: f64| (x1 - x2).abs() <= eps.in_seconds();
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
            if ya.abs() < condition.epsilon_eval.abs() {
                closest_t = xa;
                break;
            }
            if yb.abs() < condition.epsilon_eval.abs() {
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
            // self.until_time_elapsed(s - self.dynamics.time())?;
            self.until_time_elapsed(s * TimeUnit::Second)?;
            let ys = self.event_trackers.events[0].eval(&self.state);
            d = c;
            c = xb;
            yc = yb;
            if ya * ys < 0.0 {
                // Root bracketed between a and s
                // Propagate until time xa
                self.until_time_elapsed(xa * TimeUnit::Second)?;
                let ya_p = self.event_trackers.events[0].eval(&self.state);
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
                self.until_time_elapsed(xb * TimeUnit::Second)?;
                let yb_p = self.event_trackers.events[0].eval(&self.state);
                let (_a, _ya, _b, _yb) = arrange(s, ys, xb, yb_p);
                {
                    xa = _a;
                    ya = _ya;
                    xb = _b;
                    yb = _yb;
                }
            }
            iter += 1;
            if iter >= condition.max_iter {
                return Err(NyxError::MaxIterReached(iter));
            }
        }

        // Now that we have the time at which the condition is matched, let's propagate until then
        self.prevent_tx = false;
        let elapsed_time = -(closest_t * TimeUnit::Second);
        // self.until_time_elapsed(self.state.time() - closest_t)?;
        self.until_time_elapsed(elapsed_time)?;

        Ok(self.state)
    }

    /// This method integrates whichever function is provided as `d_xdt`. Everything passed to this function is in **seconds**.
    ///
    /// The `derive` method is monomorphic to increase speed. This function takes a time `t` and a current state `state`
    /// then derives the dynamics at that time (i.e. propagates for one time step). The `d_xdt` parameter is the derivative
    /// function which take a time t of type f64 and a reference to a state of type VectorN<f64, N>, and returns the
    /// result as VectorN<f64, N> of the derivative. The reference should preferrably only be borrowed.
    /// This function returns the next time (i.e. the previous time incremented by the timestep used) and
    /// the new state as y_{n+1} = y_n + \frac{dy_n}{dt}. To get the integration details, check `Self.latest_details`.
    /// Note: using VectorN<f64, N> instead of DVector implies that the function *must* always return a vector of the same
    /// size. This static allocation allows for high execution speeds.
    fn derive(
        &mut self,
        t: f64,
        // state: &VectorN<f64, <D::StateType as State>::PropVecSize>,
        // ctx: &D::StateType,
    ) -> Result<(f64, VectorN<f64, <D::StateType as State>::PropVecSize>), NyxError> {
        let state = &self.state_vector();
        let ctx = &self.state;
        // Reset the number of attempts used (we don't reset the error because it's set before it's read)
        self.details.attempts = 1;
        // Convert the step size to seconds;
        let step_size = self.step_size.in_seconds();
        loop {
            let ki = self.prop.dynamics.eom(0.0, state, ctx)?;
            self.k[0] = ki;
            let mut a_idx: usize = 0;
            for i in 0..(self.prop.stages - 1) {
                // Let's compute the c_i by summing the relevant items from the list of coefficients.
                // \sum_{j=1}^{i-1} a_ij  ∀ i ∈ [2, s]
                let mut ci: f64 = 0.0;
                // The wi stores the a_{s1} * k_1 + a_{s2} * k_2 + ... + a_{s, s-1} * k_{s-1} +
                let mut wi =
                    VectorN::<f64, <D::StateType as State>::PropVecSize>::from_element(0.0);
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
                VectorN::<f64, <D::StateType as State>::PropVecSize>::from_element(0.0);
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
                return Ok(((t + step_size), next_state));
            } else {
                // Compute the error estimate.
                self.details.error = E::estimate(&error_est, &next_state, &state);
                if self.details.error <= self.prop.opts.tolerance
                    || self.step_size <= self.prop.opts.min_step
                    || self.details.attempts >= self.prop.opts.attempts
                {
                    if self.details.attempts >= self.prop.opts.attempts {
                        warn!(
                            "maximum number of attempts reached ({})",
                            self.details.attempts
                        );
                    }

                    self.details.step = self.step_size;
                    if self.details.error < self.prop.opts.tolerance {
                        // Let's increase the step size for the next iteration.
                        // Error is less than tolerance, let's attempt to increase the step for the next iteration.
                        let proposed_step = 0.9
                            * self.step_size.in_seconds()
                            * (self.prop.opts.tolerance / self.details.error)
                                .powf(1.0 / f64::from(self.prop.order));
                        self.step_size = if proposed_step > self.prop.opts.max_step.in_seconds() {
                            self.prop.opts.max_step
                        } else {
                            proposed_step * TimeUnit::Second
                        };
                    }
                    return Ok(((t + self.details.step.in_seconds()), next_state));
                } else {
                    // Error is too high and we aren't using the smallest step, and we haven't hit the max number of attempts.
                    // So let's adapt the step size.
                    self.details.attempts += 1;
                    let proposed_step = 0.9
                        * self.step_size.in_seconds()
                        * (self.prop.opts.tolerance / self.details.error)
                            .powf(1.0 / f64::from(self.prop.order - 1));
                    assert!(!proposed_step.is_nan());
                    self.step_size = if proposed_step < self.prop.opts.min_step.in_seconds() {
                        self.prop.opts.min_step
                    } else {
                        proposed_step * TimeUnit::Second
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
    init_step: Duration,
    min_step: Duration,
    max_step: Duration,
    tolerance: f64,
    attempts: u8,
    fixed_step: bool,
    errctrl: E,
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
            errctrl,
        }
    }

    pub fn with_adaptive_step_s(min_step: f64, max_step: f64, tolerance: f64, errctrl: E) -> Self {
        Self::with_adaptive_step(
            min_step * TimeUnit::Second,
            max_step * TimeUnit::Second,
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
}

impl PropOpts<RSSStepPV> {
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
            errctrl: RSSStepPV {},
        }
    }

    pub fn with_fixed_step_s(step: f64) -> Self {
        Self::with_fixed_step(step * TimeUnit::Second)
    }

    /// Returns the default options with a specific tolerance.
    #[allow(clippy::field_reassign_with_default)]
    pub fn with_tolerance(tolerance: f64) -> Self {
        let mut opts = Self::default();
        opts.tolerance = tolerance;
        opts
    }
}

impl Default for PropOpts<RSSStepPV> {
    /// `default` returns the same default options as GMAT.
    fn default() -> PropOpts<RSSStepPV> {
        PropOpts {
            init_step: 60.0 * TimeUnit::Second,
            min_step: 0.001 * TimeUnit::Second,
            max_step: 2700.0 * TimeUnit::Second,
            tolerance: 1e-12,
            attempts: 50,
            fixed_step: false,
            errctrl: RSSStepPV {},
        }
    }
}

#[test]
fn test_options() {
    use super::error_ctrl::RSSStep;

    let opts = PropOpts::with_fixed_step_s(1e-1);
    assert_eq!(opts.min_step, 1e-1 * TimeUnit::Second);
    assert_eq!(opts.max_step, 1e-1 * TimeUnit::Second);
    assert!(opts.tolerance.abs() < std::f64::EPSILON);
    assert_eq!(opts.fixed_step, true);

    let opts = PropOpts::with_adaptive_step_s(1e-2, 10.0, 1e-12, RSSStep {});
    assert_eq!(opts.min_step, 1e-2 * TimeUnit::Second);
    assert_eq!(opts.max_step, 10.0 * TimeUnit::Second);
    assert!((opts.tolerance - 1e-12).abs() < std::f64::EPSILON);
    assert_eq!(opts.fixed_step, false);

    let opts: PropOpts<RSSStepPV> = Default::default();
    assert_eq!(opts.init_step, 60.0 * TimeUnit::Second);
    assert_eq!(opts.min_step, 0.001 * TimeUnit::Second);
    assert_eq!(opts.max_step, 2700.0 * TimeUnit::Second);
    assert!((opts.tolerance - 1e-12).abs() < std::f64::EPSILON);
    assert_eq!(opts.attempts, 50);
    assert_eq!(opts.fixed_step, false);
}
