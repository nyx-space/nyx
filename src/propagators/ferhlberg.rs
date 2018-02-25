pub use super::RK;

pub struct Ferhlberg54 {}
pub struct Ferhlberg65 {}

/// Ferhlberg54 is a [Runge Kutta Ferhlberg integrator](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method).
impl RK for Ferhlberg54 {
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

/// Ferhlberg65 is a [Runge Kutta Ferhlberg integrator](http://www.mymathlib.com/diffeq/embedded_runge_kutta/embedded_fehlberg_5_6.html).
/// NOTE: The adaptive step size is identical for all integrators, regardless of the information in this implementation. Moreover,
/// the coefficients used in this embedded implementation are slightly different than the WP definition when comparing the 54 order.
impl RK for Ferhlberg65 {
    fn order() -> usize {
        6 as usize
    }

    //     k1 = f( x[i],y[i] ),                                                   //
    //     k2 = f( x[i]+h/6, y[i] + h*k1/6 ),                                     //
    //     k3 = f( x[i]+4h/15, y[i]+h/75*(4 k1 + 16 k2) ),                        //
    //     k4 = f( x[i]+2h/3, y[i]+h/6*(5 k1 - 16 k2 + 15 k3) ),                  //
    //     k5 = f( x[i]+4h/5, y[i]+h/25*(-40 k1 + 144 k2 - 100 k3 + 16 k4))       //
    //     k6 = f( x[i]+h, y[i]+h/640*(722 k1 - 2304 k2 + 2035 k3                 //
    //                                                       - 88 k4 + 275 k5) )  //
    //     k7 = f( x[i], y[i]+h/1280*( -22 k1 + 55 k3 - 88 k4 + 55 k5) )          //
    //     k8 = f( x[i]+h, y[i]+h/1280*( 186 k1 - 4608 k2 + 4015 k3 - 88 k4       //
    //                                             + 495 k5 + 1280 K7) )          //

    fn a_coeffs() -> &'static [f64] {
        &[
            1.0 / 6.0,
            4.0 / 75.0,
            16.0 / 75.0,
            5.0 / 6.0,
            -16.0 / 6.0,
            15.0 / 6.0,
            -40.0 / 25.0,
            144.0 / 25.0,
            -100.0 / 25.0,
            16.0 / 25.0,
            722.0 / 640.0,
            -2304.0 / 640.0,
            2035.0 / 640.0,
            -88.0 / 640.0,
            275.0 / 640.0,
            31.0 / 384.0,
            0.0,
            1125.0 / 2816.0,
            9.0 / 32.0,
            125.0 / 768.0,
            5.0 / 66.0,
            /*-22.0 / 1280.0,
            0.0,
            55.0 / 1280.0,
            -88.0 / 1280.0,
            55.0 / 1280.0,
            0.0,
            186.0 / 1280.0,
            -4608.0 / 1280.0,
            4015.0 / 1280.0,
            -88.0 / 1280.0,
            495.0 / 1280.0,
            0.0,
            1280.0 / 1280.0,*/
        ]
    }

    //        y[i+1] = y[i] +  h / 8448 * ( 682 * k1 + 3375 * k3 + 2376 * k4      //
    //                                                 + 1375 * k5 + 640 * k6 )   //
    //        err = - 5*h*( k1 + k6 - k7 - k8 ) / 66                              //

    fn b_coeffs() -> &'static [f64] {
        &[
            31.0 / 384.0,
            0.0,
            1125.0 / 2816.0,
            9.0 / 32.0,
            125.0 / 768.0,
            5.0 / 66.0,
            251.0 / 3072.0,
            0.0,
            8925.0 / 22528.0,
            75.0 / 256.0,
            925.0 / 6144.0,
            25.0 / 528.0,
            1.0 / 32.0, /*
            682.0 / 8448.0,
            0.0,
            3375.0 / 8448.0,
            2376.0 / 8448.0,
            1375.0 / 8448.0,
            640.0 / 8448.0,
            -5.0 / 66.0,
            0.0,
            0.0,
            0.0,
            0.0,
            5.0 / 66.0,
            5.0 / 66.0,*/
        ]
    }
}
