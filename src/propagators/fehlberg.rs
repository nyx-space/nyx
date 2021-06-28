/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

use super::RK;

/// `Fehlberg45` is a [Runge Kutta Fehlberg integrator](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method).
pub struct Fehlberg45 {}

impl RK for Fehlberg45 {
    const ORDER: u8 = 5;
    const STAGES: usize = 6;
    const A_COEFFS: &'static [f64] = &[
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
    ];
    const B_COEFFS: &'static [f64] = &[
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
    ];
}
