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

/* NOTE: This code is effectively a clone of bacon-sci, MIT License, by Wyatt Campbell. */

use std::f64::EPSILON;
use std::fmt;
use std::ops;

/// Polynomial is a statically allocated polynomial.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct Polynomial<const SIZE: usize> {
    /// Coefficients are orders by their power, e.g. index 0 is to the power 0, 1 is linear, 2 is quadratic, etc.
    pub coefficients: [f64; SIZE],
}

impl<const SIZE: usize> Polynomial<SIZE> {
    pub fn from_most_significant(mut coeffs: [f64; SIZE]) -> Self {
        coeffs.reverse();
        Self {
            coefficients: coeffs,
        }
    }

    /// Get the order of the polynomial
    pub const fn order(&self) -> usize {
        SIZE - 1
    }

    /// Evaluate the polynomial at the provided position
    pub fn eval(&self, x: f64) -> f64 {
        self.eval_n_deriv(x).0
    }

    /// Evaluate the derivative at the provided position
    pub fn deriv(&self, x: f64) -> f64 {
        self.eval_n_deriv(x).1
    }

    /// Evaluate the polynomial and its derivative at the provided position
    pub fn eval_n_deriv(&self, x: f64) -> (f64, f64) {
        if SIZE == 1 {
            return (self.coefficients[0], 0.0);
        }

        // Start with biggest coefficients
        let mut acc_eval = *self.coefficients.last().unwrap();
        let mut acc_deriv = *self.coefficients.last().unwrap();
        // For every coefficient except the constant and largest
        for val in self.coefficients.iter().skip(1).rev().skip(1) {
            acc_eval = acc_eval * x + *val;
            acc_deriv = acc_deriv * x + acc_eval;
        }
        // Do the constant for the polynomial evaluation
        acc_eval = x * acc_eval + self.coefficients[0];

        (acc_eval, acc_deriv)
    }

    /// Initializes a Polynomial with only zeros
    pub fn zeros() -> Self {
        Self {
            coefficients: [0.0; SIZE],
        }
    }

    /// Set the i-th power of this polynomial to zero (e.g. if i=0, set the x^0 coefficient to zero, i.e. the constant part goes to zero)
    pub fn zero_power(&mut self, i: usize) {
        if i < SIZE {
            self.coefficients[i] = 0.0;
        }
    }

    /// Set all of the coefficients below this tolerance to zero
    pub fn zero_below_tolerance(&mut self, tol: f64) {
        for i in 0..=self.order() {
            if self.coefficients[i].abs() < tol {
                self.zero_power(i);
            }
        }
    }

    /// Returns true if any of the coefficients are NaN
    pub fn is_nan(&self) -> bool {
        for c in self.coefficients {
            if c.is_nan() {
                return true;
            }
        }
        false
    }

    /// Shifts all of the coefficients by one degree, dropping the largest degree.
    /// For example:
    /// P(x) = 10x^3 -6.13353243x^2 + 3.61185323x + 0.194533 .. becomes ...
    /// P(x) = -6.13353243x^3 + 3.61185323x^2 + 0.194533x
    pub(crate) fn shift_by_one(&mut self) {
        let prev_coeff = self.coefficients;
        self.coefficients[1..((prev_coeff.len() - 1) + 1)]
            .clone_from_slice(&prev_coeff[..(prev_coeff.len() - 1)]);
        self.coefficients[0] = 0.0;
    }

    fn fmt_with_var(&self, f: &mut fmt::Formatter, var: String) -> fmt::Result {
        write!(f, "P({}) = ", var)?;
        let mut data = Vec::with_capacity(SIZE);

        for (i, c) in self.coefficients.iter().enumerate().rev() {
            if c.abs() <= EPSILON {
                continue;
            }

            let mut d;
            if c.abs() > 100.0 || c.abs() < 0.01 {
                // Use scientific notation
                if c > &0.0 {
                    d = format!("+{:e}", c);
                } else {
                    d = format!("{:e}", c);
                }
            } else if c > &0.0 {
                d = format!("+{}", c);
            } else {
                d = format!("{}", c);
            }
            // Add the power
            let p = i;
            match p {
                0 => {} // Show nothing for zero
                1 => d = format!("{}{}", d, var),
                _ => d = format!("{}{}^{}", d, var, p),
            }
            data.push(d);
        }
        write!(f, "{}", data.join(" "))
    }
}

/// In-place multiplication of a polynomial with an f64
impl<const SIZE: usize> ops::Mul<f64> for Polynomial<SIZE> {
    type Output = Polynomial<SIZE>;

    fn mul(mut self, rhs: f64) -> Self::Output {
        for val in &mut self.coefficients {
            *val *= rhs;
        }
        self
    }
}

/// Clone current polynomial and then multiply it with an f64
impl<const SIZE: usize> ops::Mul<f64> for &Polynomial<SIZE> {
    type Output = Polynomial<SIZE>;

    fn mul(self, rhs: f64) -> Self::Output {
        *self * rhs
    }
}

/// In-place multiplication of a polynomial with an f64
impl<const SIZE: usize> ops::Mul<Polynomial<SIZE>> for f64 {
    type Output = Polynomial<SIZE>;

    fn mul(self, rhs: Polynomial<SIZE>) -> Self::Output {
        let mut me = rhs;
        for val in &mut me.coefficients {
            *val *= self;
        }
        me
    }
}

impl<const SIZE: usize> ops::AddAssign<f64> for Polynomial<SIZE> {
    fn add_assign(&mut self, rhs: f64) {
        self.coefficients[0] += rhs;
    }
}

impl<const SIZE: usize> fmt::Display for Polynomial<SIZE> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.fmt_with_var(f, "t".to_string())
    }
}

impl<const SIZE: usize> fmt::LowerHex for Polynomial<SIZE> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        self.fmt_with_var(f, "x".to_string())
    }
}

pub(crate) fn add<const S1: usize, const S2: usize>(
    p1: Polynomial<S1>,
    p2: Polynomial<S2>,
) -> Polynomial<S1> {
    if S1 < S2 {
        panic!();
    }
    let mut rtn = Polynomial::zeros();
    for (i, c1) in p1.coefficients.iter().enumerate() {
        rtn.coefficients[i] = match p2.coefficients.get(i) {
            Some(c2) => c1 + c2,
            None => *c1,
        };
    }
    rtn
}

impl<const S1: usize, const S2: usize> ops::Add<Polynomial<S1>> for Polynomial<S2> {
    type Output = Polynomial<S1>;
    /// Add Self and Other, _IF_ S2 >= S1 (else panic!)
    fn add(self, other: Polynomial<S1>) -> Self::Output {
        add(other, self)
    }
}

/// Subtracts p1 from p2 (p3 = p1 - p2)
pub(crate) fn sub<const S1: usize, const S2: usize>(
    p1: Polynomial<S1>,
    p2: Polynomial<S2>,
) -> Polynomial<S1> {
    if S1 < S2 {
        panic!();
    }
    let mut rtn = Polynomial::zeros();
    for (i, c1) in p1.coefficients.iter().enumerate() {
        rtn.coefficients[i] = match p2.coefficients.get(i) {
            Some(c2) => c1 - c2,
            None => *c1,
        };
    }
    rtn
}

impl<const S1: usize, const S2: usize> ops::Sub<Polynomial<S2>> for Polynomial<S1> {
    type Output = Polynomial<S1>;
    fn sub(self, other: Polynomial<S2>) -> Self::Output {
        sub(self, other)
    }
}

/// Multiply two polynomials. First parameter is the size of the first polynomial, second is the size of the second, and third is the sum of both minus one.
/// Implementation is naive and has a complexity of O(n*m) where n and m are the sizes of the polynomials.
pub(crate) fn multiply<const S1: usize, const S2: usize, const S3: usize>(
    p1: Polynomial<S1>,
    p2: Polynomial<S2>,
) -> Polynomial<S3> {
    let mut rslt = Polynomial::<S3>::zeros();
    for (exponent, val) in p2.coefficients.iter().enumerate() {
        if (*val).abs() < std::f64::EPSILON {
            // Skip any zeros to allow multiplying large polynomials with themselves.
            continue;
        }
        let if_was_scalar = *val * p1;
        for (pos, ival) in if_was_scalar.coefficients.iter().enumerate() {
            if (*ival).abs() < std::f64::EPSILON {
                // Skip any zeros to allow multiplying large polynomials with themselves.
                continue;
            }
            rslt.coefficients[pos + exponent] += *ival;
        }
    }

    rslt
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum CommonPolynomial {
    Constant(f64),
    /// Linear(a, b) <=> f(x) = ax + b (order is FLIPPED from Polynomial<N> structure)
    Linear(f64, f64),
    /// Quadratic(a, b, c) <=> f(x) = ax^2 + bx + c (order is FLIPPED from Polynomial<N> structure)
    Quadratic(f64, f64, f64),
}

impl CommonPolynomial {
    pub fn eval(&self, x: f64) -> f64 {
        match *self {
            Self::Constant(a) => Polynomial::<1> { coefficients: [a] }.eval(x),
            Self::Linear(a, b) => Polynomial::<2> {
                coefficients: [b, a],
            }
            .eval(x),
            Self::Quadratic(a, b, c) => Polynomial::<3> {
                coefficients: [c, b, a],
            }
            .eval(x),
        }
    }

    pub fn deriv(&self, x: f64) -> f64 {
        match *self {
            Self::Constant(a) => Polynomial::<1> { coefficients: [a] }.deriv(x),
            Self::Linear(a, b) => Polynomial::<2> {
                coefficients: [b, a],
            }
            .deriv(x),
            Self::Quadratic(a, b, c) => Polynomial::<3> {
                coefficients: [c, b, a],
            }
            .deriv(x),
        }
    }
}

impl fmt::Display for CommonPolynomial {
    /// Prints the polynomial with the least significant coefficients first
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Self::Constant(a) => write!(f, "{}", Polynomial::<1> { coefficients: [a] }),
            Self::Linear(a, b) => write!(
                f,
                "{}",
                Polynomial::<2> {
                    coefficients: [b, a],
                }
            ),
            Self::Quadratic(a, b, c) => write!(
                f,
                "{}",
                Polynomial::<3> {
                    coefficients: [c, b, a],
                }
            ),
        }
    }
}

#[test]
fn poly_constant() {
    let c = CommonPolynomial::Constant(10.0);
    for i in -100..=100 {
        assert!(
            (c.eval(i as f64) - 10.0).abs() < 2e-16,
            "Constant polynomial returned wrong value"
        );
    }
}

#[test]
fn poly_linear() {
    let c = CommonPolynomial::Linear(2.0, 10.0);
    for i in -100..=100 {
        let x = i as f64;
        let expect = 2.0 * x + 10.0;
        assert!(
            (c.eval(x) - expect).abs() < 2e-16,
            "Constant polynomial returned wrong value"
        );
    }
}

#[test]
fn poly_quadratic() {
    let p = Polynomial {
        coefficients: [101.0, -2.0, 3.0],
    };
    let p2 = 2.0 * p;
    let c = CommonPolynomial::Quadratic(3.0, -2.0, 101.0);
    for i in -100..=100 {
        let x = i as f64;
        let expect = 3.0 * x.powi(2) - 2.0 * x + 101.0;
        let expect_deriv = 6.0 * x - 2.0;
        assert!(
            (c.eval(x) - expect).abs() < 2e-16,
            "Polynomial returned wrong value"
        );
        assert!(
            (p.deriv(x) - expect_deriv).abs() < 2e-16,
            "Polynomial derivative returned wrong value"
        );

        assert!(
            (p.eval(x) - expect).abs() < 2e-16,
            "Polynomial returned wrong value"
        );
        assert!(
            (p2.eval(x) - 2.0 * expect).abs() < 2e-16,
            "Polynomial returned wrong value"
        );
    }
}

#[test]
fn poly_print() {
    let p = Polynomial {
        coefficients: [101.0, -2.0, 3.0],
    };
    println!("{}", p);
    assert_eq!(
        format!("{}", p),
        format!("{}", CommonPolynomial::Quadratic(3.0, -2.0, 101.0))
    );
}

#[test]
fn poly_add() {
    let p1 = Polynomial {
        coefficients: [4.0, -2.0, 3.0],
    };
    let p2 = Polynomial {
        coefficients: [0.0, -5.0, 0.0, 2.0],
    };
    //      P(x) = (3x^2 - 2x + 4) + (2x^3 - 5x)
    // <=>  P(x) = 2x^3 + 3x^2 -7x + 4
    let p_expected = Polynomial {
        coefficients: [4.0, -7.0, 3.0, 2.0],
    };

    // let p3 = add::<4, 3>(p2, p1);
    let p3 = p1 + p2;
    println!("p3 = {:x}\npe = {:x}", p3, p_expected);
    assert_eq!(p3, p_expected);
    // Check this is correct
    for i in -100..=100 {
        let x = i as f64;
        let expect = p1.eval(x) + p2.eval(x);
        assert!(
            (p3.eval(x) - expect).abs() < 2e-16,
            "Constant polynomial returned wrong value"
        );
    }
}

#[test]
fn poly_sub() {
    let p2 = Polynomial {
        coefficients: [4.0, -2.0, 3.0],
    };
    let p1 = Polynomial {
        coefficients: [0.0, -5.0, 0.0, 2.0],
    };
    //      P(x) = (3x^2 - 2x + 4) + (2x^3 - 5x)
    // <=>  P(x) = 2x^3 + 3x^2 -7x + 4
    let p_expected = Polynomial {
        coefficients: [-4.0, -3.0, -3.0, 2.0],
    };

    let p3 = p1 - p2;
    println!("p3 = {:x}\npe = {:x}", p3, p_expected);
    assert_eq!(p3, p_expected);
    // Check this is correct
    for i in -100..=100 {
        let x = i as f64;
        let expect = p1.eval(x) - p2.eval(x);
        assert!(
            (p3.eval(x) - expect).abs() < 2e-16,
            "Constant polynomial returned wrong value"
        );
    }
}

#[test]
fn poly_multiply() {
    let p1 = Polynomial {
        coefficients: [4.0, -2.0, 3.0],
    };
    let p2 = Polynomial {
        coefficients: [0.0, -5.0, 0.0, 2.0],
    };
    //      P(x) = (3x^2 - 2x + 4) * (2x^3 - 5x)
    // <=>  P(x) = (3x^2 - 2x + 4) * (2x^3) + (- 5x) * (3x^2 - 2x + 4)
    // <=>  P(x) = (6x^5 - 4x^4 + 8x^3) + (-15x^3 + 10x^2 -20x)
    // <=>  P(x) = 6x^5 - 4x^4 -7x^3 + 10x^2 -20x
    let p_expected = Polynomial {
        coefficients: [0.0, -20.0, 10.0, -7.0, -4.0, 6.0],
    };

    let p3 = multiply::<3, 4, 6>(p1, p2);
    println!("p3 = {:x}\npe = {:x}", p3, p_expected);
    assert_eq!(p3, p_expected);
    // Check this is correct
    for i in -100..=100 {
        let x = i as f64;
        let expect = p1.eval(x) * p2.eval(x);
        assert!(
            (p3.eval(x) - expect).abs() < 2e-16,
            "Constant polynomial returned wrong value"
        );
    }
}

#[test]
fn poly_shift_mulx() {
    let mut p1 = Polynomial {
        coefficients: [0.194533, 3.61185323, -6.13353243, 10.0],
    };

    let pe = Polynomial {
        coefficients: [0.0, 0.194533, 3.61185323, -6.13353243],
    };

    println!("p1 = {:x}", p1);
    p1.shift_by_one();
    println!("p1 = {:x}\npe = {:x}", p1, pe);
    assert_eq!(p1, pe);
}
