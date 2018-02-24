extern crate nalgebra as na;

pub use super::RK;

pub struct RKDP54 {}

/// RKDP54 is a [Dormand-Prince integrator](https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method).
impl RK for RKDP54 {
    fn order() -> usize {
        5 as usize
    }

    fn a_coeffs() -> &'static [f64] {
        &[
            1.0 / 5.0,
            3.0 / 40.0,
            9.0 / 40.0,
            44.0 / 45.0,
            -56.0 / 15.0,
            32.0 / 9.0,
            19372.0 / 6561.0,
            -25360.0 / 2187.0,
            64448.0 / 6561.0,
            -212.0 / 729.0,
            9017.0 / 3168.0,
            -355.0 / 33.0,
            46732.0 / 5247.0,
            49.0 / 176.0,
            -5103.0 / 18656.0,
        ]
    }
    fn b_coeffs() -> &'static [f64] {
        &[
            35.0 / 384.0,
            0.0,
            500.0 / 1113.0,
            125.0 / 192.0,
            -2187.0 / 6784.0,
            11.0 / 84.0,
            5179.0 / 57600.0,
            0.0,
            7571.0 / 16695.0,
            393.0 / 640.0,
            -92097.0 / 339200.0,
            187.0 / 2100.0,
            1.0 / 40.0,
        ]
    }
}

pub struct RKCK54 {}

/// RKCK54 is a [Runge Kutta Cash Karp integrator](https://en.wikipedia.org/wiki/Cash%E2%80%93Karp_method).
impl RK for RKCK54 {
    fn order() -> usize {
        5 as usize
    }

    fn a_coeffs() -> &'static [f64] {
        &[
            1.0 / 5.0,
            3.0 / 40.0,
            9.0 / 40.0,
            3.0 / 10.0,
            -9.0 / 10.0,
            6.0 / 5.0,
            -11.0 / 54.0,
            5.0 / 2.0,
            -70.0 / 27.0,
            35.0 / 27.0,
            1631.0 / 55296.0,
            175.0 / 512.0,
            575.0 / 13824.0,
            44275.0 / 110592.0,
            253.0 / 4096.0,
        ]
    }
    fn b_coeffs() -> &'static [f64] {
        &[
            37.0 / 378.0,
            0.0,
            250.0 / 621.0,
            125.0 / 594.0,
            0.0,
            512.0 / 1771.0,
            2825.0 / 27648.0,
            0.0,
            18575.0 / 48384.0,
            13525.0 / 55296.0,
            277.0 / 14336.0,
            1.0 / 4.0,
        ]
    }
}

pub struct RKF54 {}

/// RKF54 is a [Runge Kutta Ferhlberg integrator](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method).
impl RK for RKF54 {
    fn order() -> usize {
        5 as usize
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

pub struct RK4Fixed {}

/// RKF54 is a [Runge Kutta Ferhlberg integrator](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method).
impl RK for RK4Fixed {
    fn order() -> usize {
        3 as usize
    }

    fn a_coeffs() -> &'static [f64] {
        &[0.5, 0.0, 0.5, 0.0, 0.0, 1.0]
    }
    fn b_coeffs() -> &'static [f64] {
        &[
            1.0 / 6.0,
            1.0 / 3.0,
            1.0 / 3.0,
            1.0 / 6.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ]
    }
}
