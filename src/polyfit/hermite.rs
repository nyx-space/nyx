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

/// Builds a polynomial interpolation using the Hermite method
///
/// ```
/// use nyx_space::polyfit::hermite;
///
/// let xs: Vec<_> = (0..8).map(|i| i as f64).collect();
/// let ys: Vec<_> = xs.iter().map(|x| x.cos()).collect();
/// let derivs: Vec<_> = xs.iter().map(|x| -x.sin()).collect();
///
/// let tol = 1e-10;
/// let poly = hermite::<16>(&xs, &ys, &derivs).unwrap();
///
/// println!("{:x}", poly);
/// ```
pub fn hermite<const DEGREE: usize>(
    xs: &[f64],
    ys: &[f64],
    derivs: &[f64],
) -> Result<Polynomial<DEGREE>, NyxError> {
    if xs.len() == 0 {
        return Err(NyxError::InvalidInterpolationData(
            "No X data to interpolate".to_owned(),
        ));
    }
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

    if hermite.is_nan() {
        dbg!(xs, ys, derivs);
        return Err(NyxError::InvalidInterpolationData(format!(
            "Invalid interpolation {:x}",
            hermite
        )));
    }

    Ok(hermite)
}

#[test]
fn hermite_sine_test() {
    let xs: Vec<_> = (0..8).map(|i| i as f64).collect();
    let ys: Vec<_> = xs.iter().map(|x| x.cos()).collect();
    let derivs: Vec<_> = xs.iter().map(|x| -x.sin()).collect();

    let tol = 1e-10;
    let poly = hermite::<16>(&xs, &ys, &derivs).unwrap();

    println!("{:x}", poly);

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

#[test]
fn hermite_constant_test() {
    let xs: Vec<_> = (0..8).map(|i| i as f64).collect();
    let ys: Vec<_> = xs.iter().map(|_| 2.0159).collect();
    let derivs: Vec<_> = xs.iter().map(|_| 0.0).collect();

    let tol = 1e-10;
    let poly = hermite::<16>(&xs, &ys, &derivs).unwrap();

    println!("{:x}", poly);

    let mut max_eval_err: f64 = 0.0;
    let mut max_deriv_err: f64 = 0.0;

    for x in xs {
        let (eval, deriv) = poly.eval_n_deriv(x);
        let eval_err = (eval - 2.0159).abs();
        assert!(eval_err < tol);
        max_eval_err = max_eval_err.max(eval_err);

        let deriv_err = (deriv).abs();
        assert!(deriv_err < tol);
        max_deriv_err = max_eval_err.max(eval_err);
    }

    println!(
        "Max eval error: {:.e}\tMax deriv error: {:.e}\t",
        max_eval_err, max_deriv_err
    );
}

#[test]
fn hermite_ephem_spline_test() {
    let ts = [
        -1.0,
        -0.7142321608948587,
        -0.4284548929983568,
        -0.14272281352821248,
        0.1430009063036013,
        0.4286973024022658,
        0.714367019041751,
        1.0,
    ];
    let values = [
        -1200.6957374089038,
        -1649.3350718512218,
        -2088.1291193578113,
        -2514.3714789070427,
        -2925.5702772667646,
        -3319.240151300038,
        -3693.030156393982,
        -4044.695271513933,
    ];
    let values_dt = [
        -5.450221271198159,
        -5.3475633589540585,
        -5.212915678573803,
        -5.0471031201910135,
        -4.851091887968967,
        -4.626059429784994,
        -4.373345524123602,
        -4.094465775216765,
    ];

    let tol = 2e-7;
    let tol_deriv = 3e-6;
    let poly = hermite::<16>(&ts, &values, &values_dt).unwrap();

    println!("{:x}", poly);

    let mut max_eval_err: f64 = 0.0;
    let mut max_deriv_err: f64 = 0.0;

    for (i, t) in ts.iter().enumerate() {
        let (eval, deriv) = poly.eval_n_deriv(*t);
        let eval_err = (eval - values[i]).abs();
        assert!(dbg!(eval_err) < tol);
        max_eval_err = max_eval_err.max(eval_err);

        let deriv_err = (deriv - values_dt[i]).abs();
        assert!(dbg!(deriv_err) < tol_deriv);
        max_deriv_err = max_deriv_err.max(deriv_err);
    }

    println!(
        "Max eval error: {:.e}\tMax deriv error: {:.e}\t",
        max_eval_err, max_deriv_err
    );
}

#[test]
fn hermite_duplication_test() {
    let ts = [-1.0, 0.0, 1.0];
    let values = [239213.98224426163, 239342.1452415863, 239492.31122918683];
    let values_dt = [5.856883346456119, 1259.7108315572618, 737.5474327513627];

    let tol = 2e-16;
    let tol_deriv = 1e-11;
    let poly = hermite::<8>(&ts, &values, &values_dt).unwrap();

    println!("{:x}", poly);

    let mut max_eval_err: f64 = 0.0;
    let mut max_deriv_err: f64 = 0.0;

    for (i, t) in ts.iter().enumerate() {
        let (eval, deriv) = poly.eval_n_deriv(*t);
        let eval_err = (eval - values[i]).abs();
        assert!(dbg!(eval_err) < tol);
        max_eval_err = max_eval_err.max(eval_err);

        let deriv_err = (deriv - values_dt[i]).abs();
        assert!(dbg!(deriv_err) < tol_deriv);
        max_deriv_err = max_deriv_err.max(deriv_err);
    }

    println!(
        "Max eval error: {:.e}\tMax deriv error: {:.e}\t",
        max_eval_err, max_deriv_err
    );
}
