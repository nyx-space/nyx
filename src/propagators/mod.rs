extern crate nalgebra as na;

use self::na::{DMatrix, DVector};

// Re-Export
mod classic_rk;
pub use self::classic_rk::*;

/// The `RK` trait defines a Runge Kutta
pub trait RK
where
    Self: Sized,
{
    /// Returns a Matrix corresponding to the A coefficients of the Butcher table for that RK.
    /// *Warning:* this RK trait supposes that the implementation is consistent, i.e. c_i = \sum_j a_{ij}.
    fn a_coeffs() -> DMatrix<f64>;
    /// Returns a Matrix corresponding to the b_i and b^*_i coefficients of the Butcher table for that RK
    fn b_coeffs() -> DMatrix<f64>;
    fn from_options(opts: Options) -> Self;
    fn options(self) -> Options;
}

#[derive(Clone)]
pub struct Propagator {
    ti: f64,
    state: DVector<f64>,
    opts: Options,
    latest_error: f64,
    a_coeffs: DMatrix<f64>,
    b_coeffs: DMatrix<f64>,
}

/// The `Propagator` trait defines the functions of a propagator.
/// TODO: Add examples
impl Propagator {
    /// Each propagator must be initialized with `new` in order to store the options needed for this propagation.
    pub fn new<T: RK>(t0: f64, state: DVector<f64>, method: T) -> Propagator {
        Propagator {
            ti: t0,
            state: state,
            opts: method.options(),
            latest_error: 0.0,
            a_coeffs: T::a_coeffs(),
            b_coeffs: T::b_coeffs(),
        }
    }

    /// The `prop` method is monomorphic to increase speed. The `d_xdt` parameter is the derivative function which take
    /// a time t of type f64 and a state of type DVector<f64>, and returns the result as DVector<f64> of the derivative.
    pub fn prop<F>(&mut self, d_xdt: F) -> DVector<f64>
    where
        F: Fn(f64, DVector<f64>) -> DVector<f64>,
    {
        let mut k = Vec::new(); // Will store all the k_i.
        let (order, _) = self.a_coeffs.shape();
        for i in 0..order {
            let ci: f64 = self.a_coeffs.row(i).iter().sum();
            let mut wi = DVector::from_element(self.state.shape().0 as usize, 0.0);
            for j in 0..k.len() {
                let a_ij = self.a_coeffs.slice((i, j), (1, 1)).0;
                wi += a_ij * &k[j];
            }
            let ki = d_xdt(
                self.ti + ci * self.opts.min_step,
                self.state.clone() + self.opts.min_step * wi,
            );
            k.push(ki);
        }
        // Compute the next state and the error
        let mut next_state = self.state.clone();
        let mut next_state_star = self.state.clone();
        for i in 0..k.len() {
            let b_ij = self.b_coeffs.slice((i, 0), (1, 1));
            let b_ij_star = self.b_coeffs.slice((i, 1), (1, 1));
            next_state += self.opts.min_step * b_ij * &k[i];
            next_state_star += self.opts.min_step * b_ij_star * &k[i];
        }
        // TODO: Adaptive step size
        self.latest_error = (next_state.clone() - next_state_star).norm();
        next_state
    }
    pub fn latest_error(self) -> f64 {
        self.latest_error
    }
    pub fn latest_state(self) -> DVector<f64> {
        self.state
    }
}

/// Options stores the integrator options, including the minimum and maximum step sizes, and the max error size.
/// Note that different step sizes and max errors are only used for adaptive methods. To use a fixed step
/// integrator, initialize the options using `with_fixed_step`, and use whichever adaptive step integrator is desired.
/// For example, initializing an RK45 with fixed step options will lead to an RK4 being used instead of an RK45.
#[derive(Clone, Debug)]
pub struct Options {
    min_step: f64,
    max_step: f64,
    tolerance: f64,
    fixed_step: bool,
    debug: bool,
}

impl Options {
    /// `with_fixed_step` initializes an `Options` such that the integrator is used with a fixed step size.
    pub fn with_fixed_step(step: f64) -> Options {
        Options {
            min_step: step,
            max_step: step,
            tolerance: 0.0,
            fixed_step: true,
            debug: false,
        }
    }

    /// `with_adaptive_step` initializes an `Options` such that the integrator is used with an adaptive step size.
    pub fn with_adaptive_step(min_step: f64, max_step: f64, tolerance: f64) -> Options {
        Options {
            min_step: min_step,
            max_step: max_step,
            tolerance: tolerance,
            fixed_step: false,
            debug: false,
        }
    }

    /// Calling enable_debug on a mutable `Options` will enable verbose debugging of the integrator.
    pub fn enable_debug(&mut self) {
        self.debug = true;
    }
}

// TODO: export all RK methods here

#[test]
fn test_options() {
    let mut opts = Options::with_fixed_step(1e-1);
    assert_eq!(opts.min_step, 1e-1);
    assert_eq!(opts.max_step, 1e-1);
    assert_eq!(opts.tolerance, 0.0);
    assert_eq!(opts.fixed_step, true);
    assert_eq!(opts.debug, false);
    opts.enable_debug();
    assert_eq!(opts.debug, true);

    let mut opts = Options::with_adaptive_step(1e-2, 10.0, 1e-12);
    assert_eq!(opts.min_step, 1e-2);
    assert_eq!(opts.max_step, 10.0);
    assert_eq!(opts.tolerance, 1e-12);
    assert_eq!(opts.fixed_step, false);
    assert_eq!(opts.debug, false);
    opts.enable_debug();
    assert_eq!(opts.debug, true);
}

#[test]
fn test_consistency() {
    // All the RKs should be consistent.
    println!("{}", RKF54::a_coeffs());
}
