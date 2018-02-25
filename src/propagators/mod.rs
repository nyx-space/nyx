extern crate nalgebra as na;

use self::na::{DefaultAllocator, Dim, DimName, VectorN};
use self::na::allocator::Allocator;

// Re-Export
mod rk;
pub use self::rk::*;
mod dormand;
pub use self::dormand::*;
mod ferhlberg;
pub use self::ferhlberg::*;

/// The `RK` trait defines a Runge Kutta
pub trait RK
where
    Self: Sized,
{
    // Returns the order of this integrator.
    fn order() -> usize;
    /// Returns a pointer to a list of f64 corresponding to the A coefficients of the Butcher table for that RK.
    /// This module only supports *implicit* integrators, and as such, `Self.a_coeffs().len()` must be of
    /// size (order+1)*(order)/2.
    /// *Warning:* this RK trait supposes that the implementation is consistent, i.e. c_i = \sum_j a_{ij}.
    fn a_coeffs() -> &'static [f64];
    /// Returns a pointer to a list of f64 corresponding to the b_i and b^*_i coefficients of the
    /// Butcher table for that RK. `Self.a_coeffs().len()` must be of size (order+1)*2.
    fn b_coeffs() -> &'static [f64];
}

#[derive(Clone, Debug)]
pub struct IntegrationDetails {
    pub step: f64,
    pub error: f64,
}

#[derive(Clone, Debug)]
pub struct Propagator<'a> {
    opts: Options,
    details: IntegrationDetails,
    order: usize,
    a_coeffs: &'a [f64],
    b_coeffs: &'a [f64],
}

/// The `Propagator` trait defines the functions of a propagator.
/// TODO: Add examples
impl<'a> Propagator<'a> {
    /// Each propagator must be initialized with `new` which stores propagator information.
    pub fn new<T: RK>(opts: Options) -> Propagator<'a> {
        Propagator {
            opts: opts.clone(),
            details: IntegrationDetails {
                step: opts.max_step,
                error: 0.0,
            },
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
        state: VectorN<f64, N>,
        d_xdt: F,
    ) -> (f64, VectorN<f64, N>)
    where
        F: Fn(f64, &VectorN<f64, N>) -> VectorN<f64, N>,
        DefaultAllocator: Allocator<f64, N>,
    {
        loop {
            let mut k = Vec::with_capacity(self.order + 1); // Will store all the k_i.
            let mut prev_end = 0;
            let ki = d_xdt(t, &state.clone());
            k.push(ki);
            let mut a_idx: usize = 0;
            for i in 0..self.order {
                // Let's compute the c_i by summing the relevant items from the list of coefficients.
                let mut ci: f64 = 0.0;
                for ak in prev_end..prev_end + i + 1 {
                    ci += self.a_coeffs[ak];
                }
                prev_end += i + 1;
                let mut wi = VectorN::from_element(0.0);
                for kj in &k {
                    let a_ij = self.a_coeffs[a_idx];
                    wi += a_ij * kj;
                    a_idx += 1;
                }

                let ki = d_xdt(
                    t + ci * self.details.step,
                    &(state.clone() + self.details.step * wi),
                );
                k.push(ki);
            }
            // Compute the next state and the error
            let mut next_state = state.clone();
            let mut next_state_star = state.clone();
            for i in 0..k.len() {
                let b_i = self.b_coeffs[i];
                let b_i_star = self.b_coeffs[i + self.order];
                next_state += self.details.step * b_i * &k[i];
                next_state_star += self.details.step * b_i_star * &k[i];
            }

            self.details.error = (next_state.clone() - next_state_star).norm();

            // TODO: Implement the adaptive step size as per https://en.wikipedia.org/wiki/Adaptive_stepsize
            // and check its efficiency compared to the current /=2 algo

            if self.opts.fixed_step
                || (!self.opts.fixed_step
                    && (self.details.error < self.opts.tolerance
                        || self.details.step == self.opts.min_step))
            {
                return ((t + self.details.step), next_state);
            } else if !self.opts.fixed_step {
                if self.details.step > self.opts.min_step {
                    // Let's compute the new step using the [WP](https://en.wikipedia.org/wiki/Adaptive_stepsize) algorithm.
                    // NOTE: We do a custom min and max implementation because std::cmp doesn't support f64s due to NaN and NNaN.
                    let eta = self.opts.tolerance / self.details.error;
                    self.details.step = 0.9 * self.details.step * if eta > 0.3 {
                        // max(self.opts.tolerance / self.details.error, 0.3)
                        if eta < 2.0 {
                            // min(..., 2.0)
                            eta
                        } else {
                            2.0
                        }
                    } else {
                        // min(0.3, 2.0)
                        0.3
                    };
                    // Check that our step isn't too small
                    if self.details.step < self.opts.min_step {
                        self.details.step = self.opts.min_step;
                    }
                }
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
        }
    }

    /// `with_adaptive_step` initializes an `Options` such that the integrator is used with an
    ///  adaptive step size. TODO: Add algorithms for step size computation (sigmoid, etc.)
    pub fn with_adaptive_step(min_step: f64, max_step: f64, tolerance: f64) -> Options {
        Options {
            min_step: min_step,
            max_step: max_step,
            tolerance: tolerance,
            fixed_step: false,
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
