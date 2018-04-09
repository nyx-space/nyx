extern crate nalgebra as na;

use std::f64;
use self::na::{DefaultAllocator, Dim, DimName, VectorN};
use self::na::allocator::Allocator;

// Re-Export
mod rk;
pub use self::rk::*;
mod dormand;
pub use self::dormand::*;
mod fehlberg;
pub use self::fehlberg::*;
mod verner;
pub use self::verner::*;

/// The `RK` trait defines a Runge Kutta integrator.
pub trait RK
where
    Self: Sized,
{
    /// Returns the order of this integrator (as u8 because there probably isn't an order greater than 255).
    /// The order is used for the adaptive step size only to compute the error between estimates.
    fn order() -> u8;

    /// Returns the stages of this integrator (as usize because it's used as indexing)
    fn stages() -> usize;

    /// Returns a pointer to a list of f64 corresponding to the A coefficients of the Butcher table for that RK.
    /// This module only supports *implicit* integrators, and as such, `Self.a_coeffs().len()` must be of
    /// size (order+1)*(order)/2.
    /// *Warning:* this RK trait supposes that the implementation is consistent, i.e. c_i = \sum_j a_{ij}.
    fn a_coeffs() -> &'static [f64];
    /// Returns a pointer to a list of f64 corresponding to the b_i and b^*_i coefficients of the
    /// Butcher table for that RK. `Self.a_coeffs().len()` must be of size (order+1)*2.
    fn b_coeffs() -> &'static [f64];
}

/// Store the details of the previous integration step of a given propagator. Access as `my_prop.clone().latest_details()`.
#[derive(Clone, Debug)]
pub struct IntegrationDetails {
    pub step: f64,
    pub error: f64,
    pub attempts: u8,
}

/// Includes the options, the integrator details of the previous step, and
/// the set of coefficients used for the monomorphic instance. **WARNING:** must be stored in a mutuable variable.
#[derive(Clone, Debug)]
pub struct Propagator<'a> {
    opts: Options,
    details: IntegrationDetails,
    order: u8,
    stages: usize,
    a_coeffs: &'a [f64],
    b_coeffs: &'a [f64],
}

/// The `Propagator` trait defines the functions of a propagator.
/// TODO: Add examples
impl<'a> Propagator<'a> {
    /// Each propagator must be initialized with `new` which stores propagator information.
    pub fn new<T: RK>(opts: &Options) -> Propagator<'a> {
        Propagator {
            opts: opts.clone(),
            details: IntegrationDetails {
                step: opts.max_step,
                error: 0.0,
                attempts: 1,
            },
            stages: T::stages(),
            order: T::order(),
            a_coeffs: T::a_coeffs(),
            b_coeffs: T::b_coeffs(),
        }
    }

    /// The `derive` method is monomorphic to increase speed. This function takes a time `t` and a current state `state`
    /// then derives the dynamics at that time (i.e. propagates for one time step). The `d_xdt` parameter is the derivative
    /// function which take a time t of type f64 and a reference to a state of type VectorN<f64, N>, and returns the
    /// result as VectorN<f64, N> of the derivative. The reference should preferrably only be borrowed.
    /// This function returns the next time (i.e. the previous time incremented by the timestep used) and
    /// the new state as y_{n+1} = y_n + \frac{dy_n}{dt}. To get the integration details, check `Self.latest_details`.
    /// Note: using VectorN<f64, N> instead of DVector implies that the function *must* always return a vector of the same
    /// size. This static allocation allows for high execution speeds.
    pub fn derive<F, N: Dim + DimName>(
        &mut self,
        t: f64,
        state: &VectorN<f64, N>,
        d_xdt: F,
    ) -> (f64, VectorN<f64, N>)
    where
        F: Fn(f64, &VectorN<f64, N>) -> VectorN<f64, N>,
        DefaultAllocator: Allocator<f64, N>,
    {
        loop {
            let mut k = Vec::with_capacity(self.stages + 1); // Will store all the k_i.
            let ki = d_xdt(t, state);
            k.push(ki);
            let mut a_idx: usize = 0;
            for _ in 0..(self.stages - 1) {
                // Let's compute the c_i by summing the relevant items from the list of coefficients.
                // \sum_{j=1}^{i-1} a_ij  ∀ i ∈ [2, s]
                let mut ci: f64 = 0.0;
                // The wi stores the a_{s1} * k_1 + a_{s2} * k_2 + ... + a_{s, s-1} * k_{s-1} +
                let mut wi = VectorN::from_element(0.0);
                for kj in &k {
                    let a_ij = self.a_coeffs[a_idx];
                    ci += a_ij;
                    wi += a_ij * kj;
                    a_idx += 1;
                }

                let ki = d_xdt(
                    t + ci * self.details.step,
                    &(state + self.details.step * wi),
                );
                k.push(ki);
            }
            // Compute the next state and the error
            let mut next_state = state.clone();
            let mut error_est = VectorN::from_element(0.0);
            for (i, ki) in k.iter().enumerate() {
                let b_i = self.b_coeffs[i];
                let b_i_star = self.b_coeffs[i + self.stages];
                error_est += b_i_star * ki;
                next_state += self.details.step * b_i * ki;
            }

            if !self.opts.fixed_step {
                for (i, error_est_i) in error_est.clone().iter().enumerate() {
                    let delta = next_state[(i, 0)] - state.clone()[(i, 0)];
                    let err = if delta > self.opts.tolerance {
                        // If greater than the relative tolerance, then we normalize it by the difference.
                        (error_est_i / delta).abs()
                    } else {
                        error_est_i.abs()
                    };
                    if i == 0 || err > self.details.error {
                        self.details.error = err;
                    }
                }
            } else {
                self.details.error = 0.0;
            }

            if self.opts.fixed_step
                || (self.details.error < self.opts.tolerance
                    || (self.details.step - self.opts.min_step).abs() <= f64::EPSILON)
                || self.details.attempts >= self.opts.attempts
            {
                // Using a fixed step, no adaptive step necessary, or
                // Error is within the desired tolerance, or it isn't but we've already reach the minimum step allowed
                return ((t + self.details.step), next_state);
            } else if !self.opts.fixed_step {
                self.details.attempts += 1;
                // Error is too high and using adaptive step size
                let proposed_step = 0.9 * self.details.step
                    * (self.opts.tolerance / self.details.error)
                        .powf(1.0 / (self.order as f64 - 1.0));
                self.details.step = if proposed_step < self.opts.min_step {
                    self.opts.min_step
                } else {
                    proposed_step
                };
            }
        }
    }
    pub fn latest_details(self) -> IntegrationDetails {
        self.details
    }
}

/// Options stores the integrator options, including the minimum and maximum step sizes, and the
/// max error size. Note that different step sizes and max errors are only used for adaptive
/// methods. To use a fixed step integrator, initialize the options using `with_fixed_step`, and
/// use whichever adaptive step integrator is desired.  For example, initializing an RK45 with
/// fixed step options will lead to an RK4 being used instead of an RK45.
#[derive(Clone, Debug)]
pub struct Options {
    min_step: f64,
    max_step: f64,
    tolerance: f64,
    fixed_step: bool,
    attempts: u8,
}

impl Options {
    /// `with_fixed_step` initializes an `Options` such that the integrator is used with a fixed
    ///  step size.
    pub fn with_fixed_step(step: f64) -> Options {
        Options {
            min_step: step,
            max_step: step,
            tolerance: 0.0,
            fixed_step: true,
            attempts: 0,
        }
    }

    /// `with_adaptive_step` initializes an `Options` such that the integrator is used with an
    ///  adaptive step size.
    pub fn with_adaptive_step(min_step: f64, max_step: f64, tolerance: f64) -> Options {
        Options {
            min_step: min_step,
            max_step: max_step,
            tolerance: tolerance,
            fixed_step: false,
            attempts: 50,
        }
    }
}

// TODO: export all RK methods here

#[test]
fn test_options() {
    let opts = Options::with_fixed_step(1e-1);
    assert_eq!(opts.min_step, 1e-1);
    assert_eq!(opts.max_step, 1e-1);
    assert_eq!(opts.tolerance, 0.0);
    assert_eq!(opts.fixed_step, true);

    let opts = Options::with_adaptive_step(1e-2, 10.0, 1e-12);
    assert_eq!(opts.min_step, 1e-2);
    assert_eq!(opts.max_step, 10.0);
    assert_eq!(opts.tolerance, 1e-12);
    assert_eq!(opts.fixed_step, false);
}
