/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use crate::cosmic::Orbit;
use crate::linalg::{
    allocator::Allocator, DefaultAllocator, DimName, Matrix3, Matrix6, OVector, Vector3, Vector6,
};
use nalgebra::Complex;

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

/// Checks if the provided 3x3 matrix is diagonal.
///
/// A square matrix is considered diagonal if all its off-diagonal elements are zero.
/// This function verifies this property for a given 3x3 matrix.
/// It checks each off-diagonal element of the matrix and returns `false` if any of them
/// is not approximately zero, considering a tolerance defined by `f64::EPSILON`.
///
/// # Arguments
///
/// * `m` - A 3x3 matrix of `f64` elements to be checked.
///
/// # Returns
///
/// * `bool` - Returns `true` if the matrix is diagonal, `false` otherwise.
///
/// # Example
///
/// ```
/// let m = Matrix3::new(1.0, 0.0, 0.0,
///                      0.0, 2.0, 0.0,
///                      0.0, 0.0, 3.0);
/// assert_eq!(is_diagonal(&m), true);
/// ```
///
/// # Note
///
/// This function uses `f64::EPSILON` as the tolerance for checking if an element is approximately zero.
/// This means that elements with absolute value less than `f64::EPSILON` are considered zero.
pub fn is_diagonal(m: &Matrix3<f64>) -> bool {
    for i in 0..3 {
        for j in 0..3 {
            if i != j && m[(i, j)].abs() > f64::EPSILON {
                return false;
            }
        }
    }
    true
}

/// Returns whether this matrix represent a stable linear system looking at its eigen values.
/// NOTE: This code is not super pretty on purpose to make it clear that we're covering all of the cases
#[allow(clippy::if_same_then_else)]
#[allow(clippy::branches_sharing_code)]
pub fn are_eigenvalues_stable<N: DimName>(eigenvalues: OVector<Complex<f64>, N>) -> bool
where
    DefaultAllocator: Allocator<Complex<f64>, N>,
{
    // Source: https://eng.libretexts.org/Bookshelves/Industrial_and_Systems_Engineering/Book%3A_Chemical_Process_Dynamics_and_Controls_(Woolf)/10%3A_Dynamical_Systems_Analysis/10.04%3A_Using_eigenvalues_and_eigenvectors_to_find_stability_and_solve_ODEs#Summary_of_Eigenvalue_Graphs
    for ev in &eigenvalues {
        if ev.im.abs() > 0.0 {
            // There is an imaginary part to this eigenvalue
            if ev.re > 0.0 {
                // At least one of the EVs has a positive real part and a non-zero imaginary part
                // This is sufficient condition for unstability
                return false;
            } else {
                // Option 1: The real part is zero, the system is oscilliatory
                // Option 2: The real part is negative, so the systems tends toward stability
                continue;
            }
        } else {
            // There is no imaginary part
            if ev.re > 0.0 {
                // Real value is positive, so the system is unstable
                return false;
            } else {
                // Option 1: The real value is negative, so the system is stable
                // Option 2: The real value is zero, so the system is invariant
                continue;
            }
        }
    }

    true
}

/// Returns the provided angle bounded between 0.0 and 360.0.
///
/// This function takes an angle (in degrees) and normalizes it to the range [0, 360).
/// If the angle is negative, it will be converted to a positive angle in the equivalent position.
/// For example, an angle of -90 degrees will be converted to 270 degrees.
///
/// # Arguments
///
/// * `angle` - An angle in degrees.
///
pub fn between_0_360(angle: f64) -> f64 {
    let mut bounded = angle % 360.0;
    if bounded < 0.0 {
        bounded += 360.0;
    }
    bounded
}

/// Returns the provided angle bounded between -180.0 and +180.0
pub fn between_pm_180(angle: f64) -> f64 {
    between_pm_x(angle, 180.0)
}

/// Returns the provided angle bounded between -x and +x.
///
/// This function takes an angle (in degrees) and normalizes it to the range [-x, x).
/// If the angle is outside this range, it will be converted to an equivalent angle within this range.
/// For example, if x is 180, an angle of 270 degrees will be converted to -90 degrees.
///
/// # Arguments
///
/// * `angle` - An angle in degrees.
/// * `x` - The boundary for the angle normalization.
pub fn between_pm_x(angle: f64, x: f64) -> f64 {
    let mut bounded = angle % (2.0 * x);
    if bounded > x {
        bounded -= 2.0 * x;
    }
    if bounded < -x {
        bounded += 2.0 * x;
    }
    bounded
}

/// The Kronecker delta function
pub fn kronecker(a: f64, b: f64) -> f64 {
    if (a - b).abs() <= f64::EPSILON {
        1.0
    } else {
        0.0
    }
}

/// Returns a rotation about the X axis. The angle must be provided in radians.
/// WARNING: this is a COORDINATE SYSTEM rotation by x radians; this matrix, when applied to a vector, rotates the vector by -x  radians, not x radians.
/// Applying the matrix to a vector yields the vector's representation relative to the rotated coordinate system.
/// Source: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/eul2xf_c.html
pub fn r1(angle: f64) -> Matrix3<f64> {
    let (s, c) = angle.sin_cos();
    Matrix3::new(1.0, 0.0, 0.0, 0.0, c, s, 0.0, -s, c)
}

/// Returns a rotation about the Y axis. The angle must be provided in radians.
/// WARNING: this is a COORDINATE SYSTEM rotation by x radians; this matrix, when applied to a vector, rotates the vector by -x  radians, not x radians.
/// Applying the matrix to a vector yields the vector's representation relative to the rotated coordinate system.
/// Source: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/eul2xf_c.html
pub fn r2(angle: f64) -> Matrix3<f64> {
    let (s, c) = angle.sin_cos();
    Matrix3::new(c, 0.0, -s, 0.0, 1.0, 0.0, s, 0.0, c)
}

/// Returns a rotation about the Z axis. The angle must be provided in radians.
/// WARNING: this is a COORDINATE SYSTEM rotation by x radians; this matrix, when applied to a vector, rotates the vector by -x  radians, not x radians.
/// Applying the matrix to a vector yields the vector's representation relative to the rotated coordinate system.
/// Source: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/eul2xf_c.html
pub fn r3(angle: f64) -> Matrix3<f64> {
    let (s, c) = angle.sin_cos();
    Matrix3::new(c, s, 0.0, -s, c, 0.0, 0.0, 0.0, 1.0)
}

/// Rotate a vector about a given axis
pub fn rotv(v: &Vector3<f64>, axis: &Vector3<f64>, theta: f64) -> Vector3<f64> {
    let k_hat = axis / axis.norm();
    v * theta.cos() + k_hat.cross(v) * theta.sin() + k_hat.dot(v) * k_hat * (1.0 - theta.cos())
}

/// Returns the components of vector a orthogonal to b
pub fn perpv(a: &Vector3<f64>, b: &Vector3<f64>) -> Vector3<f64> {
    let big_a = a[0].abs().max(a[1].abs().max(a[2].abs()));
    let big_b = b[0].abs().max(b[1].abs().max(b[2].abs()));
    if big_a < f64::EPSILON {
        Vector3::zeros()
    } else if big_b < f64::EPSILON {
        *a
    } else {
        let a_scl = a / big_a;
        let b_scl = b / big_b;
        let v = projv(&a_scl, &b_scl);
        big_a * (a_scl - v)
    }
}

/// Returns the projection of a onto b
pub fn projv(a: &Vector3<f64>, b: &Vector3<f64>) -> Vector3<f64> {
    b * a.dot(b) / b.dot(b)
}

/// Computes the RSS state errors in two provided vectors
pub fn rss_errors<N: DimName>(prop_err: &OVector<f64, N>, cur_state: &OVector<f64, N>) -> f64
where
    DefaultAllocator: Allocator<f64, N>,
{
    let mut v = 0.0;
    for i in 0..N::dim() {
        v += (prop_err[i] - cur_state[i]).powi(2);
    }
    v.sqrt()
}

/// Returns the RSS orbit errors in kilometers and kilometers per second
pub fn rss_orbit_errors(prop_err: &Orbit, cur_state: &Orbit) -> (f64, f64) {
    (
        rss_errors(&prop_err.radius(), &cur_state.radius()),
        rss_errors(&prop_err.velocity(), &cur_state.velocity()),
    )
}

/// Computes the RSS state errors in position and in velocity of two orbit vectors [P V]
pub fn rss_orbit_vec_errors(prop_err: &Vector6<f64>, cur_state: &Vector6<f64>) -> (f64, f64) {
    let err_radius = (prop_err.fixed_rows::<3>(0) - cur_state.fixed_rows::<3>(0)).norm();

    let err_velocity = (prop_err.fixed_rows::<3>(3) - cur_state.fixed_rows::<3>(3)).norm();

    (err_radius, err_velocity)
}

// Normalize between -1.0 and 1.0
pub fn normalize(x: f64, min_x: f64, max_x: f64) -> f64 {
    2.0 * (x - min_x) / (max_x - min_x) - 1.0
}

// Denormalize between -1.0 and 1.0
pub fn denormalize(xp: f64, min_x: f64, max_x: f64) -> f64 {
    (max_x - min_x) * (xp + 1.0) / 2.0 + min_x
}

// Source: https://stackoverflow.com/questions/38406793/why-is-capitalizing-the-first-letter-of-a-string-so-convoluted-in-rust
pub fn capitalize(s: &str) -> String {
    let mut c = s.chars();
    match c.next() {
        None => String::new(),
        Some(f) => f.to_uppercase().collect::<String>() + c.as_str(),
    }
}

/// Builds a 6x6 DCM from the current, previous, and post DCMs, assuming that the previous and post DCMs are exactly one second before and one second after the current DCM.
pub(crate) fn dcm_finite_differencing(
    dcm_pre: Matrix3<f64>,
    dcm_cur: Matrix3<f64>,
    dcm_post: Matrix3<f64>,
) -> Matrix6<f64> {
    let drdt = 0.5 * dcm_post - 0.5 * dcm_pre;

    dcm_assemble(dcm_cur, drdt)
}

pub(crate) fn dcm_assemble(r: Matrix3<f64>, drdt: Matrix3<f64>) -> Matrix6<f64> {
    let mut full_dcm = Matrix6::zeros();
    for i in 0..6 {
        for j in 0..6 {
            if (i < 3 && j < 3) || (i >= 3 && j >= 3) {
                full_dcm[(i, j)] = r[(i % 3, j % 3)];
            } else if i >= 3 && j < 3 {
                full_dcm[(i, j)] = drdt[(i - 3, j)];
            }
        }
    }

    full_dcm
}

#[macro_export]
macro_rules! pseudo_inverse {
    ($mat:expr) => {{
        use $crate::NyxError;
        let (rows, cols) = $mat.shape();
        if rows < cols {
            match ($mat * $mat.transpose()).try_inverse() {
                Some(m1_inv) => Ok($mat.transpose() * m1_inv),
                None => Err(NyxError::SingularJacobian),
            }
        } else {
            match ($mat.transpose() * $mat).try_inverse() {
                Some(m2_inv) => Ok(m2_inv * $mat.transpose()),
                None => Err(NyxError::SingularJacobian),
            }
        }
    }};
}

/// Returns the order of mangitude of the provided value
/// ```
/// use nyx_space::utils::mag_order;
/// assert_eq!(mag_order(1000.0), 3);
/// assert_eq!(mag_order(-5000.0), 3);
/// assert_eq!(mag_order(-0.0005), -4);
/// ```
pub fn mag_order(value: f64) -> i32 {
    value.abs().log10().floor() as i32
}

/// Returns the unit vector of the moved input vector
pub fn unitize(v: Vector3<f64>) -> Vector3<f64> {
    if v.norm() < f64::EPSILON {
        v
    } else {
        v / v.norm()
    }
}

/// Converts the input vector V from Cartesian coordinates to spherical coordinates
/// Returns ρ, θ, φ where the range ρ is in the units of the input vector and the angles are in radians
pub fn cartesian_to_spherical(v: &Vector3<f64>) -> (f64, f64, f64) {
    if v.norm() < f64::EPSILON {
        (0.0, 0.0, 0.0)
    } else {
        let range_ρ = v.norm();
        let θ = v.y.atan2(v.x);
        let φ = (v.z / range_ρ).acos();
        (range_ρ, θ, φ)
    }
}

/// Converts the input vector V from Cartesian coordinates to spherical coordinates
/// Returns ρ, θ, φ where the range ρ is in the units of the input vector and the angles are in radians
pub fn spherical_to_cartesian(range_ρ: f64, θ: f64, φ: f64) -> Vector3<f64> {
    if range_ρ < f64::EPSILON {
        // Treat a negative range as a zero vector
        Vector3::zeros()
    } else {
        let x = range_ρ * φ.sin() * θ.cos();
        let y = range_ρ * φ.sin() * θ.sin();
        let z = range_ρ * φ.cos();
        Vector3::new(x, y, z)
    }
}

#[test]
fn test_tilde_matrix() {
    let vec = Vector3::new(1.0, 2.0, 3.0);
    let rslt = Matrix3::new(0.0, -3.0, 2.0, 3.0, 0.0, -1.0, -2.0, 1.0, 0.0);
    assert_eq!(tilde_matrix(&vec), rslt);
}

#[rustfmt::skip]
#[test]
fn test_diagonality() {
    assert!(!is_diagonal(&Matrix3::new(10.0, 0.0, 0.0,
                                       1.0, 5.0, 0.0,
                                       0.0, 0.0, 2.0)),
        "lower triangular"
    );
    assert!(!is_diagonal(&Matrix3::new(10.0, 1.0, 0.0,
                                       1.0, 5.0, 0.0,
                                       0.0, 0.0, 2.0)),
        "symmetric but not diag"
    );
    assert!(!is_diagonal(&Matrix3::new(10.0, 1.0, 0.0,
                                       0.0, 5.0, 0.0,
                                       0.0, 0.0, 2.0)),
        "upper triangular"
    );
    assert!(is_diagonal(&Matrix3::new(10.0, 0.0, 0.0,
                                       0.0, 0.0, 0.0,
                                       0.0, 0.0, 2.0)),
        "diagonal with zero diagonal element"
    );
    assert!(is_diagonal(&Matrix3::new(10.0, 0.0, 0.0,
                                      0.0, 5.0, 0.0,
                                      0.0, 0.0, 2.0)),
        "diagonal"
    );
}

#[test]
fn test_perpv() {
    assert_eq!(
        perpv(&Vector3::new(6.0, 6.0, 6.0), &Vector3::new(2.0, 0.0, 0.0)),
        Vector3::new(0.0, 6.0, 6.0)
    );
    assert_eq!(
        perpv(&Vector3::new(6.0, 6.0, 6.0), &Vector3::new(-3.0, 0.0, 0.0)),
        Vector3::new(0.0, 6.0, 6.0)
    );
    assert_eq!(
        perpv(&Vector3::new(6.0, 6.0, 0.0), &Vector3::new(0.0, 7.0, 0.0)),
        Vector3::new(6.0, 0.0, 0.0)
    );
    assert_eq!(
        perpv(&Vector3::new(6.0, 0.0, 0.0), &Vector3::new(0.0, 0.0, 9.0)),
        Vector3::new(6.0, 0.0, 0.0)
    );
}

#[test]
fn test_projv() {
    assert_eq!(
        projv(&Vector3::new(6.0, 6.0, 6.0), &Vector3::new(2.0, 0.0, 0.0)),
        Vector3::new(6.0, 0.0, 0.0)
    );
    assert_eq!(
        projv(&Vector3::new(6.0, 6.0, 6.0), &Vector3::new(-3.0, 0.0, 0.0)),
        Vector3::new(6.0, 0.0, 0.0)
    );
    assert_eq!(
        projv(&Vector3::new(6.0, 6.0, 0.0), &Vector3::new(0.0, 7.0, 0.0)),
        Vector3::new(0.0, 6.0, 0.0)
    );
    assert_eq!(
        projv(&Vector3::new(6.0, 0.0, 0.0), &Vector3::new(0.0, 0.0, 9.0)),
        Vector3::new(0.0, 0.0, 0.0)
    );
}

#[test]
fn test_angle_bounds() {
    assert!((between_pm_180(181.0) - -179.0).abs() < f64::EPSILON);
    assert!((between_0_360(-179.0) - 181.0).abs() < f64::EPSILON);
}

#[test]
fn test_positive_angle() {
    assert_eq!(between_0_360(450.0), 90.0);
    assert_eq!(between_pm_x(270.0, 180.0), -90.0);
}

#[test]
fn test_negative_angle() {
    assert_eq!(between_0_360(-90.0), 270.0);
    assert_eq!(between_pm_x(-270.0, 180.0), 90.0);
}

#[test]
fn test_angle_in_range() {
    assert_eq!(between_0_360(180.0), 180.0);
    assert_eq!(between_pm_x(90.0, 180.0), 90.0);
}

#[test]
fn test_zero_angle() {
    assert_eq!(between_0_360(0.0), 0.0);
    assert_eq!(between_pm_x(0.0, 180.0), 0.0);
}

#[test]
fn test_full_circle_angle() {
    assert_eq!(between_0_360(360.0), 0.0);
    assert_eq!(between_pm_x(360.0, 180.0), 0.0);
}

#[test]
fn test_pseudo_inv() {
    use crate::linalg::{DMatrix, SMatrix};
    let mut mat = DMatrix::from_element(1, 3, 0.0);
    mat[(0, 0)] = -1407.273208782421;
    mat[(0, 1)] = -2146.3100013104886;
    mat[(0, 2)] = 84.05022886527551;

    println!("{}", pseudo_inverse!(&mat).unwrap());

    let mut mat = SMatrix::<f64, 1, 3>::zeros();
    mat[(0, 0)] = -1407.273208782421;
    mat[(0, 1)] = -2146.3100013104886;
    mat[(0, 2)] = 84.05022886527551;

    println!("{}", pseudo_inverse!(&mat).unwrap());

    let mut mat = SMatrix::<f64, 3, 1>::zeros();
    mat[(0, 0)] = -1407.273208782421;
    mat[(1, 0)] = -2146.3100013104886;
    mat[(2, 0)] = 84.05022886527551;

    println!("{}", pseudo_inverse!(&mat).unwrap());

    // Compare a pseudo inverse with a true inverse
    let mat = SMatrix::<f64, 2, 2>::new(3.0, 4.0, -2.0, 1.0);
    println!("{}", mat.try_inverse().unwrap());

    println!("{}", pseudo_inverse!(&mat).unwrap());
}

#[test]
fn spherical() {
    for v in &[
        Vector3::<f64>::x(),
        Vector3::<f64>::y(),
        Vector3::<f64>::z(),
        Vector3::<f64>::zeros(),
        Vector3::<f64>::new(159.1, 561.2, 756.3),
    ] {
        let (range_ρ, θ, φ) = cartesian_to_spherical(v);
        let v_prime = spherical_to_cartesian(range_ρ, θ, φ);

        assert!(rss_errors(v, &v_prime) < 1e-12, "{} != {}", v, &v_prime);
    }
}
