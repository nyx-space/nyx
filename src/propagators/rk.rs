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

pub use super::RK;

/// `CashKarp45` is a [Runge Kutta Cash Karp integrator](https://en.wikipedia.org/wiki/Cash%E2%80%93Karp_method).
pub struct CashKarp45 {}

impl RK for CashKarp45 {
    const ORDER: u8 = 5;
    const STAGES: usize = 6;
    const A_COEFFS: &'static [f64] = &[
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
        1_631.0 / 55_296.0,
        175.0 / 512.0,
        575.0 / 13_824.0,
        44_275.0 / 110_592.0,
        253.0 / 4_096.0,
    ];
    const B_COEFFS: &'static [f64] = &[
        37.0 / 378.0,
        0.0,
        250.0 / 621.0,
        125.0 / 594.0,
        0.0,
        512.0 / 1_771.0,
        2_825.0 / 27_648.0,
        0.0,
        18_575.0 / 48_384.0,
        13_525.0 / 55_296.0,
        277.0 / 14_336.0,
        1.0 / 4.0,
    ];
}

/// `RK4Fixed` is a fixed step RK4.
///
/// If initialized with an `PropOpts.with_adaptive_step`, the variable step will **not** be taken into consideration.
#[allow(clippy::upper_case_acronyms)]
pub struct RK4Fixed {}

impl RK for RK4Fixed {
    const ORDER: u8 = 4;
    const STAGES: usize = 4;
    const A_COEFFS: &'static [f64] = &[0.5, 0.0, 0.5, 0.0, 0.0, 1.0];
    const B_COEFFS: &'static [f64] = &[
        1.0 / 6.0,
        1.0 / 3.0,
        1.0 / 3.0,
        1.0 / 6.0,
        // NOTE: Duplicating the B coefficients for force the error to zero.
        1.0 / 6.0,
        1.0 / 3.0,
        1.0 / 3.0,
        1.0 / 6.0,
    ];
}

/// `RK2Fixed` is a fixed step RK4 (or midpoint method).
///
/// If initialized with an `PropOpts.with_adaptive_step`, the variable step will **not** be taken into consideration.
#[allow(clippy::upper_case_acronyms)]
pub struct RK2Fixed {}

impl RK for RK2Fixed {
    const ORDER: u8 = 2;
    const STAGES: usize = 2;
    const A_COEFFS: &'static [f64] = &[2.0 / 3.0];
    const B_COEFFS: &'static [f64] = &[
        1.0 / 4.0,
        3.0 / 4.0,
        // NOTE: Duplicating the B coefficients for force the error to zero.
        1.0 / 4.0,
        3.0 / 4.0,
    ];
}

const SQRT6: f64 = 2.449_489_742_783_178;

/// `RK89` is a Runge Kutta 8-9 integrator.
///
/// Coefficients taken from GMAT `src/base/propagator/RungeKutta89.cpp`.
#[allow(clippy::upper_case_acronyms)]
pub struct RK89 {}

impl RK for RK89 {
    const ORDER: u8 = 9;
    const STAGES: usize = 16;
    const A_COEFFS: &'static [f64] = &[
        1.0 / 12.0,
        1.0 / 27.0,
        2.0 / 27.0,
        1.0 / 24.0,
        0.0,
        1.0 / 8.0,
        (4.0 + 94.0 * SQRT6) / 375.0,
        0.0,
        (-94.0 - 84.0 * SQRT6) / 125.0,
        (328.0 + 208.0 * SQRT6) / 375.0,
        (9.0 - SQRT6) / 150.0,
        0.0,
        0.0,
        (312.0 + 32.0 * SQRT6) / 1425.0,
        (69.0 + 29.0 * SQRT6) / 570.0,
        (927.0 - 347.0 * SQRT6) / 1250.0,
        0.0,
        0.0,
        (-16248.0 + 7328.0 * SQRT6) / 9375.0,
        (-489.0 + 179.0 * SQRT6) / 3750.0,
        (14268.0 - 5798.0 * SQRT6) / 9375.0,
        2.0 / 27.0,
        0.0,
        0.0,
        0.0,
        0.0,
        (16.0 - SQRT6) / 54.0,
        (16.0 + SQRT6) / 54.0,
        19.0 / 256.0,
        0.0,
        0.0,
        0.0,
        0.0,
        (118.0 - 23.0 * SQRT6) / 512.0,
        (118.0 + 23.0 * SQRT6) / 512.0,
        -9.0 / 256.0,
        11.0 / 144.0,
        0.0,
        0.0,
        0.0,
        0.0,
        (266.0 - SQRT6) / 864.0,
        (266.0 + SQRT6) / 864.0,
        -1.0 / 16.0,
        -8.0 / 27.0,
        (5034.0 - 271.0 * SQRT6) / 61440.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        (7859.0 - 1626.0 * SQRT6) / 10240.0,
        (-2232.0 + 813.0 * SQRT6) / 20480.0,
        (-594.0 + 271.0 * SQRT6) / 960.0,
        (657.0 - 813.0 * SQRT6) / 5120.0,
        (5996.0 - 3794.0 * SQRT6) / 405.0,
        0.0,
        0.0,
        0.0,
        0.0,
        (-4342.0 - 338.0 * SQRT6) / 9.0,
        (154_922.0 - 40458.0 * SQRT6) / 135.0,
        (-4176.0 + 3794.0 * SQRT6) / 45.0,
        (-340_864.0 + 242_816.0 * SQRT6) / 405.0,
        (26304.0 - 15176.0 * SQRT6) / 45.0,
        -26624.0 / 81.0,
        (3793.0 + 2168.0 * SQRT6) / 103_680.0,
        0.0,
        0.0,
        0.0,
        0.0,
        (4042.0 + 2263.0 * SQRT6) / 13824.0,
        (-231_278.0 + 40717.0 * SQRT6) / 69120.0,
        (7947.0 - 2168.0 * SQRT6) / 11520.0,
        (1048.0 - 542.0 * SQRT6) / 405.0,
        (-1383.0 + 542.0 * SQRT6) / 720.0,
        2624.0 / 1053.0,
        3.0 / 1664.0,
        -137.0 / 1296.0,
        0.0,
        0.0,
        0.0,
        0.0,
        (5642.0 - 337.0 * SQRT6) / 864.0,
        (5642.0 + 337.0 * SQRT6) / 864.0,
        -299.0 / 48.0,
        184.0 / 81.0,
        -44.0 / 9.0,
        -5120.0 / 1053.0,
        -11.0 / 468.0,
        16.0 / 9.0,
        (33617.0 - 2168.0 * SQRT6) / 518_400.0,
        0.0,
        0.0,
        0.0,
        0.0,
        (-3846.0 + 31.0 * SQRT6) / 13824.0,
        (155_338.0 - 52807.0 * SQRT6) / 345_600.0,
        (-12537.0 + 2168.0 * SQRT6) / 57600.0,
        (92.0 + 542.0 * SQRT6) / 2025.0,
        (-1797.0 - 542.0 * SQRT6) / 3600.0,
        320.0 / 567.0,
        -1.0 / 1920.0,
        4.0 / 105.0,
        0.0,
        (-36487.0 - 30352.0 * SQRT6) / 279_600.0,
        0.0,
        0.0,
        0.0,
        0.0,
        (-29666.0 - 4499.0 * SQRT6) / 7456.0,
        (2_779_182.0 - 615_973.0 * SQRT6) / 186_400.0,
        (-94329.0 + 91056.0 * SQRT6) / 93200.0,
        (-232_192.0 + 121_408.0 * SQRT6) / 17475.0,
        (101_226.0 - 22764.0 * SQRT6) / 5825.0,
        -169_984.0 / 9087.0,
        -87.0 / 30290.0,
        492.0 / 1165.0,
        0.0,
        1260.0 / 233.0,
    ];
    const B_COEFFS: &'static [f64] = &[
        23.0 / 525.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        171.0 / 1400.0,
        86.0 / 525.0,
        93.0 / 280.0,
        -2048.0 / 6825.0,
        -3.0 / 18200.0,
        39.0 / 175.0,
        0.0,
        9.0 / 25.0,
        233.0 / 4200.0,
        // NOTE: The b_i_star are defined here as a subtraction since the coefficients come from GMAT,
        // and GMAT actually hard codes the b_i (as `cj[]`) and the errors.
        23.0 / 525.0 + 7.0 / 400.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        171.0 / 1400.0 - 63.0 / 200.0,
        86.0 / 525.0 + 14.0 / 25.0,
        93.0 / 280.0 - 21.0 / 20.0,
        -2048.0 / 6825.0 + 1024.0 / 975.0,
        -3.0 / 18200.0 + 21.0 / 36400.0,
        39.0 / 175.0 + 3.0 / 25.0,
        9.0 / 280.0,
        0.0,
        0.0,
    ];
}
