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

// Dumb factorial, avoid using this for large numbers
fn factorial(x: f64) -> f64 {
    if (x - 1.0).abs() < std::f64::EPSILON {
        1.0
    } else {
        x * factorial(x - 1.0)
    }
}

/// Taylor expansion of sin with x near zero, order 5
pub fn sin(x: f64) -> f64 {
    x - x.powi(3) / factorial(3.0) + x.powi(5) / factorial(5.0)
}

/// Taylor expansion of sin with x near zero, order 5
pub fn cos(x: f64) -> f64 {
    1.0 - x.powi(2) / 2.0 + x.powi(4) / factorial(4.0)
}

#[test]
fn sincos_taylor() {
    use std::f64::consts::{FRAC_PI_2, TAU};
    let increment = TAU * 0.05;
    let mut x = -FRAC_PI_2 + increment;
    loop {
        let sin_x = x.sin();
        let sin_x_te = sin(x);
        let cos_x = x.cos();
        let cos_x_te = cos(x);

        assert!(dbg!(sin_x - sin_x_te).abs() < 1e-2);
        assert!(dbg!(cos_x - cos_x_te).abs() < 1e-2);

        x += increment;
        if x > FRAC_PI_2 - increment {
            break;
        }
    }
}
