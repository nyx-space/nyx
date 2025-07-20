/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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

use core::f64::consts::PI;

use crate::errors::LambertError;

use super::{LambertSolution, TransferKind, Vector3, LAMBERT_EPSILON, MAX_ITERATIONS};

/// Solve the Lambert boundary problem using Izzo's method.
///
/// This is an implementation of D. Izzo's method for solving Lambert's problem, as described in "Revisiting Lambert’s problem".
/// The code was adapted from the Python version available in jorgepiloto's [lamberthub](https://github.com/jorgepiloto/lamberthub/blob/main/src/lamberthub/universal_solvers/izzo.py),
/// which is released under the GPL v3 license, compatible with Nyx's AGPL v3 license.
///
/// Given the initial and final radii, a time of flight, and a gravitational parameters, it returns the needed initial and final velocities.
///
/// # Arguments
///
/// * `r_init` - The initial radius vector.
/// * `r_final` - The final radius vector.
/// * `tof_s` - The time of flight in seconds.
/// * `mu_km3_s2` - The gravitational parameter in km^3/s^2.
/// * `kind` - The kind of transfer (auto, short way, long way, or number of revolutions).
///
/// # Returns
///
/// `Result<LambertSolution, NyxError>` - The solution to the Lambert problem or an error if the problem could not be solved.
pub fn izzo(
    r_init: Vector3<f64>,
    r_final: Vector3<f64>,
    tof_s: f64,
    mu_km3_s2: f64,
    kind: TransferKind,
) -> Result<LambertSolution, LambertError> {
    let ri_norm = r_init.norm();
    let rf_norm = r_final.norm();

    let chord = r_final - r_init;
    let c_norm = chord.norm();

    // Semi parameter
    let s = (ri_norm + rf_norm + c_norm) * 0.5;

    // Versors
    let i_r1 = r_init / ri_norm;
    let i_r2 = r_final / rf_norm;

    let mut i_h = i_r1.cross(&i_r2);
    i_h /= i_h.norm(); // Ensure normalization

    // Geometry of the problem
    let lambda = 1.0 - c_norm / s;
    let mut m_lambda = lambda.sqrt();

    let retrograde = matches!(kind, TransferKind::LongWay);

    let (mut i_t1, mut i_t2) = if i_h.z < 0.0 || retrograde {
        // Transfer angle greater than 180 degrees
        m_lambda = -m_lambda;
        (i_r1.cross(&i_h), i_r2.cross(&i_h))
    } else {
        (i_h.cross(&i_r1), i_h.cross(&i_r2))
    };
    // Ensure unit vector
    i_t1 /= i_t1.norm();
    i_t2 /= i_t2.norm();

    // Always assume prograde for now
    let t = (2.0 * mu_km3_s2 / s.powi(3)).sqrt() * tof_s;

    // Find then filter solutions.
    let (x, y) = find_xy(
        m_lambda,
        t,
        match kind {
            TransferKind::NRevs(revs) => revs.into(),
            _ => 0,
        },
        MAX_ITERATIONS,
        LAMBERT_EPSILON,
        LAMBERT_EPSILON,
        true,
    )?;
    // Reconstruct
    let gamma = (mu_km3_s2 * s / 2.0).sqrt();
    let rho = (ri_norm - rf_norm) / c_norm;
    let sigma = (1.0 - rho.powi(2)).sqrt();

    // Compute the radial and tangential components at initial and final
    // position vectors
    let (v_r1, v_r2, v_t1, v_t2) = reconstruct(x, y, ri_norm, rf_norm, m_lambda, gamma, rho, sigma);

    // Solve for the initial and final velocity
    let v_init = v_r1 * (r_init / ri_norm) + v_t1 * i_t1;
    let v_final = v_r2 * (r_final / rf_norm) + v_t2 * i_t2;

    Ok(LambertSolution {
        v_init,
        v_final,
        phi: 0.0,
    })
}

/// Reconstructs solution velocity vectors.
/// Returns a tuple of (V_r1, V_r2, V_t1, V_t2).
pub fn reconstruct(
    x: f64,
    y: f64,
    r1: f64,
    r2: f64,
    ll: f64,
    gamma: f64,
    rho: f64,
    sigma: f64,
) -> (f64, f64, f64, f64) {
    let v_r1 = gamma * ((ll * y - x) - rho * (ll * y + x)) / r1;
    let v_r2 = -gamma * ((ll * y - x) + rho * (ll * y + x)) / r2;
    let v_t1 = gamma * sigma * (y + ll * x) / r1;
    let v_t2 = gamma * sigma * (y + ll * x) / r2;

    (v_r1, v_r2, v_t1, v_t2)
}

/// Computes x and y for a given number of revolutions.
///
/// This function orchestrates the solving process, from refining the number
/// of possible revolutions to running the Householder root-finding method.
/// It returns a Result containing a tuple of (x, y) on success.
pub fn find_xy(
    ll: f64,
    t: f64,
    m: u32,
    maxiter: usize,
    atol: f64,
    rtol: f64,
    low_path: bool,
) -> Result<(f64, f64), LambertError> {
    // For abs(ll) == 1 the derivative is not continuous.
    // An assertion is used for documenting an unrecoverable precondition.
    assert!(ll.abs() < 1.0, "|ll| must be less than 1");

    let mut m_max = (t / PI).floor() as u32;
    let t_00 = ll.acos() + ll * (1.0 - ll.powi(2)).sqrt(); // T_xM

    // Refine maximum number of revolutions if necessary
    if t < t_00 + (m_max as f64) * PI && m_max > 0 {
        let (_, t_min) = compute_t_min(ll, m_max, maxiter, atol, rtol)?;
        if t < t_min {
            m_max -= 1;
        }
    }

    // Check if a feasible solution exists for the given number of revolutions
    if m > m_max {
        return Err(LambertError::MultiRevNotFeasible { m, m_max });
    }

    // Initial guess
    let x_0 = initial_guess(t, ll, m, low_path);

    // Start Householder iterations from x_0 to find x.
    // The '?' will propagate any error from the solver.
    let x = householder(x_0, t, ll, m, atol, rtol, maxiter)?;
    let y = compute_y(x, ll);

    Ok((x, y))
}

fn compute_y(x: f64, ll: f64) -> f64 {
    (1.0 - ll.powi(2) * (1.0 - x.powi(2))).sqrt()
}

/// Computes psi.
/// "The auxiliary angle psi is computed using Eq.(17) by the appropriate
/// inverse function"
fn compute_psi(x: f64, y: f64, ll: f64) -> f64 {
    if x > 1.0 {
        // Hyperbolic motion
        ((y - x * ll) * (x.powi(2) - 1.0).sqrt()).asinh()
    } else if x < 1.0 {
        // Catches the original `-1 <= x < 1`
        // Elliptic motion
        (x * y + ll * (1.0 - x.powi(2))).acos()
    } else {
        // Parabolic motion (x = 1.0)
        0.0
    }
}

/// Time of flight equation.
fn tof_equation(x: f64, t0: f64, ll: f64, m: u32) -> f64 {
    let y = compute_y(x, ll);
    tof_equation_y(x, y, t0, ll, m)
}

/// Time of flight equation with externally computed y.
fn tof_equation_y(x: f64, y: f64, t0: f64, ll: f64, m: u32) -> f64 {
    let t_ = if m == 0 && x.powi(2) > 0.6 && x.powi(2) < 1.4 {
        let eta = y - ll * x;
        let s_1 = (1.0 - ll - x * eta) * 0.5;
        let q = 4.0 / 3.0 * hyp2f1b(s_1);
        (eta.powi(3) * q + 4.0 * ll * eta) * 0.5
    } else {
        let psi = compute_psi(x, y, ll);
        let m_float = m as f64;
        let den = 1.0 - x.powi(2);

        // In Rust, division by zero on floats results in `inf` or `NaN`,
        // which matches the behavior of `np.divide` in this context.
        ((psi + m_float * PI) / den.abs().sqrt() - x + ll * y) / den
    };

    t_ - t0
}

/// Derivative of the time of flight equation.
fn tof_equation_p(x: f64, y: f64, t: f64, ll: f64) -> f64 {
    (3.0 * t * x - 2.0 + 2.0 * ll.powi(3) * x / y) / (1.0 - x.powi(2))
}

/// Second derivative of the time of flight equation.
fn tof_equation_p2(x: f64, y: f64, t: f64, dt: f64, ll: f64) -> f64 {
    (3.0 * t + 5.0 * x * dt + 2.0 * (1.0 - ll.powi(2)) * ll.powi(3) / y.powi(3)) / (1.0 - x.powi(2))
}

/// Third derivative of the time of flight equation.
fn tof_equation_p3(x: f64, y: f64, _t: f64, dt: f64, ddt: f64, ll: f64) -> f64 {
    (7.0 * x * ddt + 8.0 * dt - 6.0 * (1.0 - ll.powi(2)) * ll.powi(5) * x / y.powi(5))
        / (1.0 - x.powi(2))
}

/// Computes minimum T.
/// Returns a tuple of `(x_T_min, T_min)`.
fn compute_t_min(
    ll: f64,
    m: u32,
    maxiter: usize,
    atol: f64,
    rtol: f64,
) -> Result<(f64, f64), LambertError> {
    // Use an epsilon for floating point comparison
    if (ll - 1.0).abs() < 1e-9 {
        let x_t_min = 0.0;
        let t_min = tof_equation(x_t_min, 0.0, ll, m);
        Ok((x_t_min, t_min))
    } else if m == 0 {
        let x_t_min = f64::INFINITY;
        let t_min = 0.0;
        Ok((x_t_min, t_min))
    } else {
        // Set x_i > 0 to avoid problems at ll = -1
        let x_i = 0.1;
        let t_i = tof_equation(x_i, 0.0, ll, m);
        let x_t_min = halley(x_i, t_i, ll, atol, rtol, maxiter)?;
        let t_min = tof_equation(x_t_min, 0.0, ll, m);
        Ok((x_t_min, t_min))
    }
}

/// Calculates the initial guess for the Lambert solver.
///
/// This function is a Rust translation of the Python `_initial_guess` method,
/// using f64 for floating-point calculations.
///
/// # Arguments
///
/// * `t_of_f` - The time of flight.
/// * `ll` - The lambda parameter, related to the geometry of the transfer.
/// * `m` - The number of complete revolutions (0 for the short path).
/// * `low_path` - A boolean indicating whether to select the low-path solution for multi-revolution cases.
///
/// # Returns
///
/// An `f64` value representing the initial guess `x_0`.
///
fn initial_guess(tof_s: f64, ll: f64, m: u32, low_path: bool) -> f64 {
    if m == 0 {
        // --- Single revolution case ---
        let t_0 = ll.acos() + ll * (1.0 - ll.powi(2)).sqrt(); // Eq. 19 (for m=0)
        let t_1 = 2.0 * (1.0 - ll.powi(3)) / 3.0; // Eq. 21

        // The if/else block is an expression that returns a value directly.
        if tof_s >= t_0 {
            (t_0 / tof_s).powf(2.0 / 3.0) - 1.0
        } else if tof_s < t_1 {
            5.0 / 2.0 * t_1 / tof_s * (t_1 - tof_s) / (1.0 - ll.powi(5)) + 1.0
        } else {
            // Condition for T_1 <= T < T_0
            // Uses the corrected initial guess formula.
            (2.0f64.ln() * (tof_s / t_0).ln() / (t_1 / t_0).ln()).exp() - 1.0
        }
    } else {
        // --- Multiple revolution case ---
        let m_float = m as f64;

        // The Python code calculates a common term twice for both x_0l and x_0r.
        // We can pre-calculate these terms.
        let term_l = ((m_float * PI + PI) / (8.0 * tof_s)).powf(2.0 / 3.0);
        let x_0l = (term_l - 1.0) / (term_l + 1.0);

        let term_r = ((8.0 * tof_s) / (m_float * PI)).powf(2.0 / 3.0);
        let x_0r = (term_r - 1.0) / (term_r + 1.0);

        // Filter out the solution using .max() and .min() methods
        if low_path {
            x_0l.max(x_0r)
        } else {
            x_0l.min(x_0r)
        }
    }
}

/// Hypergeometric function 2F1(3, 1, 5/2, x), see [Battin].
pub fn hyp2f1b(x: f64) -> f64 {
    if x >= 1.0 {
        return f64::INFINITY;
    }

    let mut res = 1.0;
    let mut term = 1.0;
    let mut ii = 0.0_f64; // Use float for loop counter to avoid casting
    loop {
        term *= (3.0 + ii) * (1.0 + ii) / (2.5 + ii) * x / (ii + 1.0);
        let res_old = res;
        res += term;

        // Convergence is reached when the result no longer changes.
        if res == res_old {
            return res;
        }
        ii += 1.0;
    }
}

// --- Solvers and High-Level Functions ---

/// Finds a minimum of time of flight equation using Halley's method.
/// Returns a Result, which is Ok(value) on success or Err(message) on failure.
pub fn halley(
    mut p0: f64,
    t0: f64,
    ll: f64,
    atol: f64,
    rtol: f64,
    maxiter: usize,
) -> Result<f64, LambertError> {
    for _ in 1..=maxiter {
        let y = compute_y(p0, ll);
        let fder = tof_equation_p(p0, y, t0, ll);
        let fder2 = tof_equation_p2(p0, y, t0, fder, ll);

        // Avoid division by zero
        if fder2.abs() < 1e-14 {
            return Err(LambertError::TargetsTooClose);
        }

        let fder3 = tof_equation_p3(p0, y, t0, fder, fder2, ll);

        // Halley step (cubic)
        let p = p0 - 2.0 * fder * fder2 / (2.0 * fder2.powi(2) - fder * fder3);

        if (p - p0).abs() < rtol * p0.abs() + atol {
            return Ok(p);
        }
        p0 = p;
    }

    Err(LambertError::SolverMaxIter { maxiter })
}

/// Finds a zero of time of flight equation using Householder's method.
/// Returns a Result, which is Ok(value) on success or Err(message) on failure.
pub fn householder(
    mut p0: f64,
    t0: f64,
    ll: f64,
    m: u32,
    atol: f64,
    rtol: f64,
    maxiter: usize,
) -> Result<f64, LambertError> {
    for _ in 1..=maxiter {
        let y = compute_y(p0, ll);
        let fval = tof_equation_y(p0, y, t0, ll, m);
        let t = fval + t0;
        let fder = tof_equation_p(p0, y, t, ll);
        let fder2 = tof_equation_p2(p0, y, t, fder, ll);
        let fder3 = tof_equation_p3(p0, y, t, fder, fder2, ll);

        // Householder step (quartic)
        let num = fder.powi(2) - fval * fder2 / 2.0;
        let den = fder * (fder.powi(2) - fval * fder2) + fder3 * fval.powi(2) / 6.0;

        if den.abs() < 1e-14 {
            return Err(LambertError::TargetsTooClose);
        }

        let p = p0 - fval * (num / den);

        if (p - p0).abs() < rtol * p0.abs() + atol {
            return Ok(p);
        }
        p0 = p;
    }
    Err(LambertError::SolverMaxIter { maxiter })
}

#[test]
fn test_lambert_izzo_shortway() {
    // Test case from Vallado, Example 7-1, p. 462
    let ri = Vector3::new(15945.34, 0.0, 0.0);
    let rf = Vector3::new(12214.83899, 10249.46731, 0.0);
    let tof_s = 76.0 * 60.0;
    let mu_km3_s2 = 3.98600433e5;

    let exp_vi = Vector3::new(2.058913, 2.915965, 0.0);
    let exp_vf = Vector3::new(-3.451565, 0.910315, 0.0);

    let sol = izzo(ri, rf, tof_s, mu_km3_s2, TransferKind::ShortWay).unwrap();

    println!("{sol:?}\t{exp_vi}\t{exp_vf}");

    assert!((sol.v_init - exp_vi).norm() < 1e-5);
    assert!((sol.v_final - exp_vf).norm() < 1e-5);
}

#[test]
fn test_lambert_izzo_longway() {
    // Test case from Vallado, Example 7-1, p. 462
    let ri = Vector3::new(15945.34, 0.0, 0.0);
    let rf = Vector3::new(12214.83899, 10249.46731, 0.0);
    let tof_s = 76.0 * 60.0;
    let mu_km3_s2 = 3.98600433e5;

    let exp_vi = Vector3::new(-3.811158, -2.003854, 0.0);
    let exp_vf = Vector3::new(4.207569, 0.914724, 0.0);

    let sol = izzo(ri, rf, tof_s, mu_km3_s2, TransferKind::LongWay).unwrap();

    assert!((sol.v_init - exp_vi).norm() < 1e-5);
    assert!((sol.v_final - exp_vf).norm() < 1e-5);
}
