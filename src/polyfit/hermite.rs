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

use crate::linalg::{
    allocator::Allocator, Const, DefaultAllocator, DimMin, DimMinimum, DimSub, OMatrix, OVector,
    SVector,
};
use crate::polyfit::polynomial::{multiply, Polynomial};
use crate::NyxError;

/// Returns the pseudo-Vandermonde matrix of degree `deg` and sample points `xs`.
/// Fully statically allocated.
/// This is a translation from [numpy](https://github.com/numpy/numpy/blob/b235f9e701e14ed6f6f6dcba885f7986a833743f/numpy/polynomial/hermite.py#L1107)
pub(crate) fn hermvander<const VALS: usize, const DEGREE: usize>(
    xs: &[f64; VALS],
    // ) -> SMatrix<f64, { VALS }, { DEGREE }> {
) -> OMatrix<f64, Const<{ VALS }>, Const<{ DEGREE }>> {
    let mut v = OMatrix::<f64, Const<{ VALS }>, Const<{ DEGREE }>>::zeros();
    let x = SVector::<f64, { VALS }>::from_column_slice(xs);

    // Set the first row to one
    let mut one = SVector::<f64, { VALS }>::zeros();
    for i in 0..VALS {
        one[i] = 1.0;
    }
    v.set_column(0, &one);
    if DEGREE > 0 {
        let x2 = 2.0 * x;
        v.set_column(1, &x2);
        for i in 2..DEGREE {
            let col =
                v.column(i - 1).component_mul(&x2) - v.column(i - 2) * (2.0 * ((i as f64) - 1.0));
            // v[i] = (v[i-1]*x2 - v[i-2]*(2*(i - 1)))
            v.set_column(i, &col);
        }
    }
    v
}

pub fn hermfit<const VALS: usize, const DEGREE: usize>(
    xs: &[f64; VALS],
    ys: &[f64; VALS],
) -> Result<Polynomial<DEGREE>, NyxError>
where
    DefaultAllocator: Allocator<f64, Const<{ VALS }>, Const<{ DEGREE }>>
        + Allocator<f64, Const<{ VALS }>>
        + Allocator<f64, <Const<{ VALS }> as DimMin<Const<{ DEGREE }>>>::Output>
        + Allocator<f64, Const<{ VALS }>, <Const<{ VALS }> as DimMin<Const<{ DEGREE }>>>::Output>
        + Allocator<f64, <Const<{ VALS }> as DimMin<Const<{ DEGREE }>>>::Output, Const<{ DEGREE }>>
        + Allocator<f64, <Const<{ VALS }> as DimMin<Const<{ DEGREE }>>>::Output, Const<{ VALS }>>
        + Allocator<usize, <Const<{ VALS }> as DimMin<Const<{ DEGREE }>>>::Output, Const<{ DEGREE }>>
        + Allocator<
            f64,
            <<Const<{ VALS }> as DimMin<Const<{ DEGREE }>>>::Output as DimSub<Const<1>>>::Output,
        >,
    Const<{ VALS }>: DimMin<Const<{ DEGREE }>>,
    DimMinimum<Const<{ VALS }>, Const<{ DEGREE }>>: DimSub<Const<1>>,
{
    let vand = hermvander::<VALS, DEGREE>(&xs);
    let y = OVector::<f64, Const<{ VALS }>>::from_column_slice(ys);

    // Normalize the vand matrix
    // https://github.com/numpy/numpy/blob/b235f9e701e14ed6f6f6dcba885f7986a833743f/numpy/polynomial/polyutils.py#L652

    match vand.svd(true, true).solve(&y, 1e-10) {
        Ok(sol) => {
            let mut coeffs = [0.0; DEGREE];
            sol.iter().enumerate().for_each(|(i, x)| coeffs[i] = *x);
            // Build the Polynomial
            // Ok(Polynomial::from_most_significant(coeffs))
            Ok(Polynomial {
                coefficients: coeffs,
            })
        }
        Err(e) => Err(NyxError::CustomError(e.to_string())),
    }
}

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

#[test]
fn hermvander_numpy_test() {
    use crate::linalg::SMatrix;
    let rslt = hermvander::<3, 4>(&[-1.0, 0.0, 1.0]);
    let expect = SMatrix::<f64, 3, 4>::from_row_slice(&[
        1.0, -2.0, 2.0, 4.0, 1.0, 0.0, -2.0, 0.0, 1.0, 2.0, 2.0, -4.0,
    ]);
    assert!(
        (rslt - expect).norm() < 2e-16,
        "hermvander returned a result different than numpy test"
    );
}

#[test]
fn hermfit_numpy_test() {
    let xs = [-1., -0.75, -0.5, -0.25, 0., 0.25, 0.5, 0.75, 1.];
    let ys = [
        1.83541234,
        0.07808521,
        -1.06151568,
        -0.67510267,
        0.194533,
        0.12639794,
        -2.10571261,
        1.30597704,
        0.1143639,
    ];
    let poly = hermfit::<9, 8>(&xs, &ys).unwrap();
    println!("{:x}", poly);

    // Polynomial from numpy
    // `herm2poly(hermfit(x, y, 8))`
    let coeffs = [
        0.194533,
        3.61185323,
        -6.13353243,
        -37.53450715,
        -29.24982842,
        89.83425821,
        123.05798117,
        -56.77212851,
        -86.89426521,
    ];
    let npoly = Polynomial {
        coefficients: coeffs,
    };

    // Evaluate error
    for (i, x) in xs.iter().enumerate() {
        println!(
            "P({}) = {}\t\terr = {:e}\t\tP1(x) err = {:e}",
            x,
            poly.eval(*x),
            poly.eval(*x) - ys[i],
            npoly.eval(*x) - ys[i]
        );
    }
}
