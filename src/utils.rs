extern crate nalgebra as na;
use self::na::{Matrix3, Vector3, Vector6, U3};
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
            if i == j && (m[(i, j)] - m[(0, 0)]) > f64::EPSILON
                || i != j
                    && (m[(i, j)].abs() > f64::EPSILON
                        || (m[(i, j)] - m[(j, i)]).abs() > f64::EPSILON)
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

/// Returns a rotation about the X axis. The angle must be provided in radians.
pub fn r1(angle: f64) -> Matrix3<f64> {
    let (s, c) = angle.sin_cos();
    Matrix3::new(1.0, 0.0, 0.0, 0.0, c, s, 0.0, -s, c)
}

/// Returns a rotation about the Y axis. The angle must be provided in radians.
pub fn r2(angle: f64) -> Matrix3<f64> {
    let (s, c) = angle.sin_cos();
    Matrix3::new(c, 0.0, -s, 0.0, 1.0, 0.0, s, 0.0, c)
}

/// Returns a rotation about the Z axis. The angle must be provided in radians.
pub fn r3(angle: f64) -> Matrix3<f64> {
    let (s, c) = angle.sin_cos();
    Matrix3::new(c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0)
}

/// Computes the RSS state errors in position and in velocity
pub fn rss_state_errors(prop_err: &Vector6<f64>, cur_state: &Vector6<f64>) -> (f64, f64) {
    let err_radius = (prop_err.fixed_rows::<U3>(0) - cur_state.fixed_rows::<U3>(0)).norm();

    let err_velocity = (prop_err.fixed_rows::<U3>(3) - cur_state.fixed_rows::<U3>(3)).norm();

    (err_radius, err_velocity)
}

#[test]
fn test_tilde_matrix() {
    let vec = Vector3::new(1.0, 2.0, 3.0);
    let rslt = Matrix3::new(0.0, -3.0, 2.0, 3.0, 0.0, -1.0, -2.0, 1.0, 0.0);
    assert_eq!(tilde_matrix(&vec), rslt);
}

#[test]
fn test_diagonality() {
    assert_eq!(
        is_diagonal(&Matrix3::new(10.0, 0.0, 0.0, 1.0, 5.0, 0.0, 0.0, 0.0, 2.0)),
        false,
        "lower triangular"
    );

    assert_eq!(
        is_diagonal(&Matrix3::new(10.0, 1.0, 0.0, 1.0, 5.0, 0.0, 0.0, 0.0, 2.0)),
        false,
        "symmetric but not diag"
    );

    assert_eq!(
        is_diagonal(&Matrix3::new(10.0, 1.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 2.0)),
        false,
        "upper triangular"
    );

    assert_eq!(
        is_diagonal(&Matrix3::new(10.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 2.0)),
        true,
        "diagonal"
    );
}
