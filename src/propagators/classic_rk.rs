extern crate nalgebra as na;

use self::na::DMatrix;
pub use super::{Options, RK};

pub struct RKF54 {
    opts: Options,
}

/// RKF54 is a [Runge Kutta Ferhlberg integrator](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method).
impl RK for RKF54 {
    fn a_coeffs() -> DMatrix<f64> {
        DMatrix::from_row_slice(
            5,
            5,
            &[
                1.0 / 4.0,
                0.0,
                0.0,
                0.0,
                0.0, // end row 1
                3.0 / 32.0,
                9.0 / 32.0,
                0.0,
                0.0,
                0.0, // end row 2
                1932.0 / 2197.0,
                -7200.0 / 2197.0,
                7296.0 / 2197.0,
                0.0,
                0.0, // end row 3
                439.0 / 216.0,
                -8.0,
                3680.0 / 513.0,
                -845.0 / 4104.0,
                0.0, // end row 4
                -8.0 / 27.0,
                2.0,
                -3544.0 / 2565.0,
                1859.0 / 4104.0,
                -11.0 / 40.0, // end row 5
            ],
        )
    }
    fn b_coeffs() -> DMatrix<f64> {
        DMatrix::from_row_slice(
            2,
            6,
            &[
                16.0 / 135.0,
                0.0,
                6656.0 / 12825.0,
                28561.0 / 56430.0,
                -9.0 / 50.0,
                2.0 / 55.0,
                25.0 / 216.0,
                0.0,
                1408.0 / 2565.0,
                2197.0 / 4104.0,
                -1.0 / 5.0,
                0.0,
            ],
        )
    }
    fn from_options(opts: Options) -> RKF54 {
        RKF54 { opts: opts }
    }
    fn options(self) -> Options {
        self.opts
    }
}
