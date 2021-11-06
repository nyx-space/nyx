use std::fmt;

/// Polynomial is a statically allocated polynomial.
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

    /// Evaluate the polynomial at the provided position
    pub fn eval(&self, x: f64) -> f64 {
        let mut acc = *self.coefficients.last().unwrap();
        for val in self.coefficients.iter().rev().skip(1) {
            acc *= x;
            acc += *val;
        }

        acc
    }
}

impl<const SIZE: usize> fmt::Display for Polynomial<SIZE> {
    /// Prints the polynomial with the least significant coefficients first
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        Ok(for (i, c) in self.coefficients.iter().enumerate() {
            if c.abs() > 100.0 {
                // Use scientific notation
                write!(f, "{:e} ", c)?;
            } else {
                write!(f, "{}", c)?;
            }
            // Add the power
            let p = self.coefficients.len() - i - 1;
            match p {
                0 => {} // Do nothing for zero
                1 => write!(f, "x")?,
                _ => write!(f, "x^{}", p)?,
            }
            if i < self.coefficients.len() - 1 {
                if *c > 0.0 {
                    write!(f, " ")?;
                } else {
                    write!(f, " + ")?;
                }
            }
        })
    }
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
        match self {
            &Self::Constant(a) => Polynomial::<1> { coefficients: [a] }.eval(x),
            &Self::Linear(a, b) => Polynomial::<2> {
                coefficients: [b, a],
            }
            .eval(x),
            &Self::Quadratic(a, b, c) => Polynomial::<3> {
                coefficients: [c, b, a],
            }
            .eval(x),
        }
    }
}

impl fmt::Display for CommonPolynomial {
    /// Prints the polynomial with the least significant coefficients first
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            &Self::Constant(a) => write!(f, "{}", Polynomial::<1> { coefficients: [a] }),
            &Self::Linear(a, b) => write!(
                f,
                "{}",
                Polynomial::<2> {
                    coefficients: [b, a],
                }
            ),
            &Self::Quadratic(a, b, c) => write!(
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
    let c = CommonPolynomial::Quadratic(3.0, -2.0, 101.0);
    for i in -100..=100 {
        let x = i as f64;
        let expect = 3.0 * x.powi(2) - 2.0 * x + 101.0;
        assert!(
            (c.eval(x) - expect).abs() < 2e-16,
            "Constant polynomial returned wrong value"
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
