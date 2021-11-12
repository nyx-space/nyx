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

/*
 * SOURCE: bacon-sci, MIT licensed, Copyright (c) Wyatt Campbell.
 */

use crate::polyfit::polynomial::{multiply, Polynomial};
use crate::NyxError;

// mod spline;
// pub use spline::*;

/// Create a Lagrange interpolating polynomial.
///
/// Create an nth degree polynomial matching the n points (xs[i], ys[i])
/// using Neville's iterated method for Lagrange polynomials. The result will
/// match no derivatives.
///
/// # Examples
/// ```
/// use bacon_sci::interp::lagrange;
/// use bacon_sci::polynomial::Polynomial;
/// fn example() {
///     let xs: Vec<_> = (0..10).map(|i| i as f64).collect();
///     let ys: Vec<_> = xs.iter().map(|x| x.cos()).collect();
///     let poly = lagrange(&xs, &ys, 1e-6).unwrap();
///     for x in xs {
///         assert!((x.cos() - poly.evaluate(x)).abs() < 0.00001);
///     }
/// }
/// ```
// pub fn lagrange<const DEGREE: usize>(
//     xs: &[N],
//     ys: &[N],
//     tol: N::RealField,
// ) -> Result<Polynomial<DEGREE>, String> {
//     if xs.len() != ys.len() {
//         return Err("lagrange: slices have mismatched dimension".to_owned());
//     }

//     let mut qs = vec![Polynomial::with_tolerance(tol)?; xs.len() * xs.len()];
//     for (ind, y) in ys.iter().enumerate() {
//         qs[ind] = polynomial![*y];
//     }

//     for i in 1..xs.len() {
//         let mut poly_2 = polynomial![N::one(), -xs[i]];
//         poly_2.set_tolerance(tol)?;
//         for j in 1..=i {
//             let mut poly_1 = polynomial![N::one(), -xs[i - j]];
//             poly_1.set_tolerance(tol)?;
//             let idenom = N::one() / (xs[i] - xs[i - j]);
//             let numer =
//                 &poly_1 * &qs[i + xs.len() * (j - 1)] - &poly_2 * &qs[(i - 1) + xs.len() * (j - 1)];
//             qs[i + xs.len() * j] = numer * idenom;
//         }
//     }

//     for i in 0..=qs[xs.len() * xs.len() - 1].order() {
//         if qs[xs.len() * xs.len() - 1].get_coefficient(i).abs() < tol {
//             qs[xs.len() * xs.len() - 1].purge_coefficient(i);
//         }
//     }

//     qs[xs.len() * xs.len() - 1].purge_leading();
//     Ok(qs[xs.len() * xs.len() - 1].clone())
// }

pub fn hermite<const DEGREE: usize>(
    xs: &[f64],
    ys: &[f64],
    derivs: &[f64],
) -> Result<Polynomial<DEGREE>, NyxError> {
    if xs.len() != ys.len() {
        return Err(NyxError::InvalidInterpolationData(
            "Lengths of X and Y data differ".to_owned(),
        ));
    }
    if xs.len() != derivs.len() {
        return Err(NyxError::InvalidInterpolationData(
            "Lengths of X and its derivatives data differ".to_owned(),
        ));
    }

    if DEGREE < 2 * xs.len() - 1 {
        warn!(
            "Building Hermite interpolation of degree {} with {} samples, {} degree recommended",
            DEGREE,
            xs.len(),
            2 * xs.len() - 1
        );
    }

    let mut zs = vec![0.0; 2 * xs.len()];
    let mut qs = vec![0.0; 4 * xs.len() * xs.len()];

    for i in 0..xs.len() {
        zs[2 * i] = xs[i];
        zs[2 * i + 1] = xs[i];
        qs[2 * i] = ys[i];
        qs[2 * i + 1] = ys[i];
        qs[2 * i + 1 + (2 * xs.len())] = derivs[i];

        if i != 0 {
            qs[2 * i + (2 * xs.len())] = (qs[2 * i] - qs[2 * i - 1]) / (zs[2 * i] - zs[2 * i - 1]);
        }
    }

    for i in 2..2 * xs.len() {
        for j in 2..=i {
            qs[i + j * (2 * xs.len())] = (qs[i + (j - 1) * (2 * xs.len())]
                - qs[i - 1 + (j - 1) * (2 * xs.len())])
                / (zs[i] - zs[i - j]);
        }
    }

    let mut hermite = Polynomial::<DEGREE>::zeros();
    for i in (1..2 * xs.len()).rev() {
        hermite += qs[i + i * (2 * xs.len())];
        let new_poly = Polynomial::<2>::from_most_significant([1.0, -xs[(i - 1) / 2]]);
        hermite = multiply::<DEGREE, 2, DEGREE>(hermite, new_poly);
    }
    hermite += qs[0];

    Ok(hermite)
}

#[test]
fn hermite_sine_test() {
    let xs: Vec<_> = (0..8).map(|i| i as f64).collect();
    let ys: Vec<_> = xs.iter().map(|x| x.cos()).collect();
    let derivs: Vec<_> = xs.iter().map(|x| -x.sin()).collect();

    let tol = 1e-10;
    let poly = hermite::<16>(&xs, &ys, &derivs).unwrap();

    println!("{}", poly);

    let mut max_eval_err: f64 = 0.0;
    let mut max_deriv_err: f64 = 0.0;

    for x in xs {
        let (eval, deriv) = poly.eval_n_deriv(x);
        let eval_err = (eval - x.cos()).abs();
        assert!(eval_err < tol);
        max_eval_err = max_eval_err.max(eval_err);

        let deriv_err = (deriv - -x.sin()).abs();
        assert!(deriv_err < tol);
        max_deriv_err = max_eval_err.max(eval_err);
    }

    println!(
        "Max eval error: {:.e}\tMax deriv error: {:.e}\t",
        max_eval_err, max_deriv_err
    );
}
