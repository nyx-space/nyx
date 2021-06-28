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

/// `Verner56` is an RK Verner integrator of order 5-6.
///
/// Coefficients taken from [here (PDF)](http://people.math.sfu.ca/~jverner/classify.1992.ps).
pub struct Verner56 {}

impl RK for Verner56 {
    const ORDER: u8 = 6;
    const STAGES: usize = 8;

    const A_COEFFS: &'static [f64] = &[
        1.0 / 6.0,
        4.0 / 75.0,
        16.0 / 75.0,
        5.0 / 6.0,
        -8.0 / 3.0,
        5.0 / 2.0,
        -165.0 / 64.0,
        55.0 / 6.0,
        -425.0 / 64.0,
        85.0 / 96.0,
        -8_263.0 / 15_000.0,
        124.0 / 75.0,
        -643.0 / 680.0,
        -81.0 / 250.0,
        2_484.0 / 10_625.0,
        3_501.0 / 1_720.0,
        -300.0 / 43.0,
        297_275.0 / 52_632.0,
        -319.0 / 2_322.0,
        24_068.0 / 84_065.0,
        3_850.0 / 26_703.0,
        12.0 / 5.0,
        -8.0,
        4_015.0 / 612.0,
        -11.0 / 36.0,
        88.0 / 255.0,
        0.0,
        0.0,
    ];

    const B_COEFFS: &'static [f64] = &[
        3.0 / 40.0,
        0.0,
        875.0 / 2_244.0,
        23.0 / 72.0,
        264.0 / 1_955.0,
        125.0 / 11_592.0,
        43.0 / 616.0,
        0.0,
        13.0 / 160.0,
        0.0,
        2_375.0 / 5_984.0,
        5.0 / 16.0,
        12.0 / 85.0,
        0.0,
        0.0,
        3.0 / 44.0,
    ];
}
