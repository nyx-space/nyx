extern crate nalgebra as na;

use self::na::DVector;

/// The `Propagator` trait defines the functions of a propagator.
pub trait Propagator {
    /// Each propagator must be initialized with `new` in order to store the options needed for this propagation.
    fn new(options: Options) -> Self;
    /// The `prop` method is monomorphic to increase speed. The `d_xdt` parameter is the derivative function which take
    /// a time t of type f64 and a state of type DVector<f64>, and returns the result as DVector<f64> of the derivative.
    fn prop<F>(self, d_xdt: F) -> DVector<f64>
    where
        F: Fn(f64, DVector<f64>) -> DVector<f64>;
}

/// Options stores the integrator options, including the minimum and maximum step sizes, and the max error size.
/// Note that different step sizes and max errors are only used for adaptive methods. To use a fixed step
/// integrator, initialize the options using `with_fixed_step`, and use whichever adaptive step integrator is desired.
/// For example, initializing an RK45 with fixed step options will lead to an RK4 being used instead of an RK45.
#[derive(Debug)]
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
