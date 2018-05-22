extern crate nalgebra as na;
use self::na::{Matrix3, Vector3};
use std::f64;

/// Returns the tilde matrix from the provided Vector3.
pub fn tilde_matrix(v: &Vector3<f64>) -> Matrix3<f64> {
    Matrix3::new(
        0.0,
        -v[(2, 0)],
        v[(1, 0)],
        v[(2, 0)],
        0.0,
        -v[(0, 0)],
        -v[(1, 0)],
        v[(0, 0)],
        0.0,
    )
}

/// Returns whether the provided square matrix (3x3) is diagonal
pub fn is_diagonal(m: &Matrix3<f64>) -> bool {
    for i in 1..2 {
        for j in 0..i {
            if (i == j && (m[(i, j)] - m[(0, 0)]) > f64::EPSILON)
                || (i != j && (m[(i, j)].abs() > f64::EPSILON || (m[(i, j)] - m[(j, i)]).abs() > f64::EPSILON))
            {
                return false;
            }
        }
    }
    true
}

/// Returns the provided angle bounded between 0.0 and 360.0
pub fn between_0_360(angle: f64) -> f64 {
    let mut bounded = angle;
    while bounded > 360.0 {
        bounded -= 360.0;
    }
    while bounded < 0.0 {
        bounded += 360.0;
    }
    bounded
}

/// Returns the provided angle bounded between -180.0 and +180.0
pub fn between_pm_180(angle: f64) -> f64 {
    let mut bounded = angle;
    while bounded > 180.0 {
        bounded -= 360.0;
    }
    while bounded < -180.0 {
        bounded += 360.0;
    }
    bounded
}

pub fn factorial(num: f64) -> f64 {
    if num <= f64::EPSILON || (num - 1.0).abs() <= f64::EPSILON {
        1.0
    } else {
        num * factorial(num - 1.0)
    }
}

pub fn kronecker(a: f64, b: f64) -> f64 {
    if (a - b).abs() <= f64::EPSILON {
        1.0
    } else {
        0.0
    }
}
