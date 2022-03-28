/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

/// Stores a Hermite series
pub struct HermiteSeries<const N: usize> {
    coefficients: [f64; N],
}

impl<const N: usize> HermiteSeries<N> {
    /// Convert a Hermite series to a Polynomial
    pub fn to_polynomial(&self) -> Polynomial<N> {
        let mut rtn = Polynomial {
            coefficients: self.coefficients,
        };
        if N == 1 {
            // Do nothing more
            return rtn;
        } else if N == 2 {
            rtn.coefficients[1] *= 2.0;
        } else {
            let mut c0 = Polynomial::<N>::zeros();
            let mut c1 = Polynomial::<N>::zeros();
            c0.coefficients[0] = self.coefficients[self.coefficients.len() - 2];
            c1.coefficients[0] = self.coefficients[self.coefficients.len() - 1];

            for i in (2..self.coefficients.len()).rev() {
                let tmp = c0;
                let mut c_im2 = Polynomial::<N>::zeros();
                c_im2.coefficients[0] = self.coefficients[i - 2];
                c0 = c_im2 - c1 * (2 * (i - 1)) as f64;
                c1.shift_by_one();
                c1 = tmp + 2.0 * c1;
            }
            c1.shift_by_one();
            rtn = c0 + 2.0 * c1;
        }
        rtn
    }
}

/// Returns the pseudo-Vandermonde matrix of degree `DEGREE-1` and sample points `xs`.
/// Fully statically allocated. **WARNING** Degree must be set to 9 for a hermite of size 8!
/// This is a translation from [numpy](https://github.com/numpy/numpy/blob/b235f9e701e14ed6f6f6dcba885f7986a833743f/numpy/polynomial/hermite.py#L1107)
fn hermvander<const VALS: usize, const DEGREE: usize>(
    xs: &[f64],
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
            v.set_column(i, &col);
        }
    }
    v
}

pub fn hermfit<const VALS: usize, const DEGREE: usize>(
    xs: &[f64],
    ys: &[f64],
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
        > + Allocator<f64, Const<{ DEGREE }>, <Const<{ VALS }> as DimMin<Const<{ DEGREE }>>>::Output>
        + Allocator<
            f64,
            <<Const<{ DEGREE }> as DimMin<Const<{ VALS }>>>::Output as DimSub<Const<1>>>::Output,
        > + Allocator<f64, <Const<{ DEGREE }> as DimMin<Const<{ VALS }>>>::Output, Const<{ VALS }>>
        + Allocator<f64, <Const<{ DEGREE }> as DimMin<Const<{ VALS }>>>::Output>
        + Allocator<f64, Const<{ DEGREE }>, <Const<{ DEGREE }> as DimMin<Const<{ VALS }>>>::Output>
        + Allocator<(f64, usize), <Const<{ VALS }> as DimMin<Const<{ DEGREE }>>>::Output>
        + Allocator<(usize, usize), <Const<{ VALS }> as DimMin<Const<{ DEGREE }>>>::Output>,

    Const<{ VALS }>: DimMin<Const<{ DEGREE }>>,
    Const<{ DEGREE }>: DimMin<Const<{ VALS }>>,
    DimMinimum<Const<{ VALS }>, Const<{ DEGREE }>>: DimSub<Const<1>>,
    DimMinimum<Const<{ DEGREE }>, Const<{ VALS }>>: DimSub<Const<1>>,
{
    let vand = hermvander::<VALS, DEGREE>(xs);
    let y = OVector::<f64, Const<{ VALS }>>::from_column_slice(ys);

    // Normalize the vand matrix
    let mut lhs = vand.transpose();
    let mut scl = lhs.component_mul(&lhs).column_sum();
    for c in scl.iter_mut() {
        if c.abs() < 2e-16 {
            *c = 1.0;
        }
        *c = c.sqrt();
    }

    for mut col in lhs.column_iter_mut() {
        col.component_div_assign(&scl);
    }

    match lhs.transpose().svd(true, true).solve(&y, 2e-16) {
        Ok(mut sol) => {
            for mut col in sol.column_iter_mut() {
                col.component_div_assign(&scl);
            }

            let mut coeffs = [0.0; DEGREE];
            sol.iter().enumerate().for_each(|(i, x)| coeffs[i] = *x);
            // Build the Hermite Series and convert to a polynomial
            let h = HermiteSeries {
                coefficients: coeffs,
            };
            Ok(h.to_polynomial())
        }
        Err(e) => Err(NyxError::InvalidInterpolationData(e.to_string())),
    }
}

pub fn hermite<const DEGREE: usize>(
    xs: &[f64],
    ys: &[f64],
    derivs: &[f64],
) -> Result<Polynomial<DEGREE>, NyxError> {
    if xs.is_empty() {
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
    let vand = hermvander::<3, 5>(&[-1.0, 0.0, 1.0]);
    let expect = SMatrix::<f64, 3, 5>::from_row_slice(&[
        1.0, -2.0, 2.0, 4.0, -20.0, 1.0, 0.0, -2.0, 0.0, 12.0, 1.0, 2.0, 2.0, -4.0, -20.0,
    ]);
    assert!(
        (vand - expect).norm() < 2e-16,
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
    let poly = hermfit::<9, 9>(&xs, &ys).unwrap();
    println!("{:x}", poly);

    // Polynomial from numpy
    // `herm2poly(hermfit(x, y, 8))`
    let npoly = Polynomial {
        coefficients: [
            0.1945330000000354,
            3.61185323000015,
            -6.133532429999718,
            -37.53450715000004,
            -29.24982842000058,
            89.83425820999997,
            123.0579811700001,
            -56.77212851,
            -86.89426521,
        ],
    };

    println!("{:x}\n", npoly);

    // Evaluate error
    for (i, x) in xs.iter().enumerate() {
        assert!((poly.eval(*x) - ys[i]).abs() < 1e-11, "hermfit failed");
        assert!(
            (npoly.eval(*x) - ys[i]).abs() < 1e-7,
            "check numpy data, it does not fit!"
        );
        println!(
            "P({}) = {}\t\terr = {:e}\t\tP1(x) err = {:e}",
            x,
            poly.eval(*x),
            poly.eval(*x) - ys[i],
            npoly.eval(*x) - ys[i]
        );
    }
}

#[test]
fn herm2poly() {
    let series = HermiteSeries {
        coefficients: [
            -364.319505276875,
            -230.472812950625,
            -817.857413263125,
            -134.8289486859375,
            -229.266493323125,
            -15.82103409828125,
            -17.08533955890625,
            -0.443532253984375,
            -0.3394307234765625,
        ],
    };
    let expected = Polynomial {
        coefficients: [
            0.1945330000000354,
            3.61185323000015,
            -6.133532429999718,
            -37.53450715000004,
            -29.24982842000058,
            89.83425820999997,
            123.0579811700001,
            -56.77212851,
            -86.89426521,
        ],
    };
    let poly = series.to_polynomial();
    println!("{}", poly);
    println!("{}", expected);
    let delta = poly - expected;
    println!("DELTA = {}", delta);
    for c in delta.coefficients {
        assert!(c.abs() < 1e-10);
    }
}
