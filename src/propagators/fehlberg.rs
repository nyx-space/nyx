use super::RK;

/// Fehlberg45 is a [Runge Kutta Fehlberg integrator](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method).
pub struct Fehlberg45 {}

impl RK for Fehlberg45 {
    fn order() -> u8 {
        5
    }
    fn stages() -> usize {
        6
    }
    fn a_coeffs() -> &'static [f64] {
        &[
            1.0 / 4.0,
            3.0 / 32.0,
            9.0 / 32.0,
            1932.0 / 2197.0,
            -7200.0 / 2197.0,
            7296.0 / 2197.0,
            439.0 / 216.0,
            -8.0,
            3680.0 / 513.0,
            -845.0 / 4104.0,
            -8.0 / 27.0,
            2.0,
            -3544.0 / 2565.0,
            1859.0 / 4104.0,
            -11.0 / 40.0,
        ]
    }
    fn b_coeffs() -> &'static [f64] {
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
        ]
    }
}

/// Fehlberg56 is a [Runge Kutta Fehlberg integrator](http://www.mymathlib.com/diffeq/embedded_runge_kutta/embedded_fehlberg_5_6.html).
/// NOTE: The adaptive step size is identical for all integrators, regardless of the information in this implementation. Moreover,
/// the coefficients used in this embedded implementation are slightly different than the WP definition when comparing the 54 order.
pub struct Fehlberg56 {}

impl RK for Fehlberg56 {
    fn order() -> u8 {
        6
    }
    fn stages() -> usize {
        8
    }
    fn a_coeffs() -> &'static [f64] {
        &[
            1.0 / 6.0,
            4.0 / 75.0,
            16.0 / 75.0,
            5.0 / 6.0,
            -8.0 / 3.0,
            5.0 / 2.0,
            -8.0 / 5.0,
            144.0 / 25.0,
            -4.0,
            16.0 / 25.0,
            361.0 / 320.0,
            -18.0 / 5.0,
            407.0 / 128.0,
            -11.0 / 80.0,
            55.0 / 128.0,
            -11.0 / 640.0,
            0.0,
            11.0 / 256.0,
            -11.0 / 160.0,
            11.0 / 256.0,
            0.0,
            93.0 / 640.0,
            -18.0 / 5.0,
            803.0 / 256.0,
            -11.0 / 160.0,
            99.0 / 256.0,
            0.0,
            1.0,
        ]
    }
    fn b_coeffs() -> &'static [f64] {
        &[
            7.0 / 1408.0,
            0.0,
            1125.0 / 2816.0,
            9.0 / 32.0,
            125.0 / 768.0,
            0.0,
            5.0 / 66.0,
            5.0 / 66.0,
            7.0 / 1408.0 - 31.0 / 384.0,
            0.0,
            0.0,
            0.0,
            0.0,
            -5.0 / 66.0,
            5.0 / 66.0,
            5.0 / 66.0,
        ]
    }
}
