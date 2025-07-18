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

use crate::errors::NyxError;
use crate::linalg::Vector3;
use std::f64::consts::PI;

const TAU: f64 = 2.0 * PI;
const LAMBERT_EPSILON: f64 = 1e-4; // General epsilon
const LAMBERT_EPSILON_TIME: f64 = 1e-4; // Time epsilon
const LAMBERT_EPSILON_RAD: f64 = (5e-5 / 180.0) * PI; // 0.00005 degrees
/// Maximum number of iterations allowed in the Lambert problem solver.
/// This is a safety measure to prevent infinite loops in case a solution cannot be found.
const MAX_ITERATIONS: usize = 1000;

/// Define the transfer kind for a Lambert
pub enum TransferKind {
    Auto,
    ShortWay,
    LongWay,
    NRevs(u8),
}

impl TransferKind {
    /// Calculate the direction multiplier based on the transfer kind.
    ///
    /// # Arguments
    ///
    /// * `r_final` - The final radius vector.
    /// * `r_init` - The initial radius vector.
    ///
    /// # Returns
    ///
    /// * `Result<f64, NyxError>` - The direction multiplier or an error if the transfer kind is not supported.
    fn direction_of_motion(
        self,
        r_final: &Vector3<f64>,
        r_init: &Vector3<f64>,
    ) -> Result<f64, NyxError> {
        match self {
            TransferKind::Auto => {
                let mut dnu = r_final[1].atan2(r_final[0]) - r_init[1].atan2(r_final[1]);
                if dnu > TAU {
                    dnu -= TAU;
                } else if dnu < 0.0 {
                    dnu += TAU;
                }

                if dnu > std::f64::consts::PI {
                    Ok(-1.0)
                } else {
                    Ok(1.0)
                }
            }
            TransferKind::ShortWay => Ok(1.0),
            TransferKind::LongWay => Ok(-1.0),
            _ => Err(NyxError::LambertMultiRevNotSupported),
        }
    }
}

#[derive(Debug)]
pub struct LambertSolution {
    pub v_init: Vector3<f64>,
    pub v_final: Vector3<f64>,
    pub phi: f64,
}

/// Solve the Lambert boundary problem using Gooding's method.
///
/// This is an implementation of R. H. Gooding's method for solving Lambert's problem, as described in "A procedure for the solution of Lambert's orbital boundary-value problem".
/// Given the initial and final radii, a time of flight, and a gravitational parameters, it returns the needed initial and final velocities
/// along with φ which is the square of the difference in eccentric anomaly. Note that the direction of motion
/// is computed directly in this function to simplify the generation of Pork chop plots.
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
pub fn gooding(
    r_init: Vector3<f64>,
    r_final: Vector3<f64>,
    tof_s: f64,
    mu_km3_s2: f64,
    kind: TransferKind,
) -> Result<LambertSolution, NyxError> {
    let r_init_norm = r_init.norm();
    let r_final_norm = r_final.norm();
    let r_norm_product = r_init_norm * r_final_norm;
    let cos_dnu = r_init.dot(&r_final) / r_norm_product;

    let dm = kind.direction_of_motion(&r_final, &r_init)?;

    let nu_init = r_init[1].atan2(r_init[0]);
    let nu_final = r_final[1].atan2(r_final[0]);

    let a = dm * (r_norm_product * (1.0 + cos_dnu)).sqrt();

    if nu_final - nu_init < LAMBERT_EPSILON_RAD && a.abs() < LAMBERT_EPSILON {
        return Err(NyxError::TargetsTooClose);
    }

    let mut phi_upper = 4.0 * PI.powi(2);
    let mut phi_lower = -4.0 * PI.powi(2);
    let mut phi = 0.0;

    let mut c2: f64 = 1.0 / 2.0;
    let mut c3: f64 = 1.0 / 6.0;
    let mut iter: usize = 0;
    let mut cur_tof: f64 = 0.0;
    let mut y = 0.0;

    while (cur_tof - tof_s).abs() > LAMBERT_EPSILON_TIME {
        if iter > MAX_ITERATIONS {
            return Err(NyxError::MaxIterReached {
                msg: format!("Lambert solver failed after {MAX_ITERATIONS} iterations"),
            });
        }
        iter += 1;

        y = r_init_norm + r_final_norm + a * (phi * c3 - 1.0) / c2.sqrt();
        if a > 0.0 && y < 0.0 {
            for _ in 0..500 {
                phi += 0.1;
                y = r_init_norm + r_final_norm + a * (phi * c3 - 1.0) / c2.sqrt();
                if y >= 0.0 {
                    break;
                }
            }
            if y < 0.0 {
                return Err(NyxError::LambertNotReasonablePhi);
            }
        }

        let chi = (y / c2).sqrt();
        cur_tof = (chi.powi(3) * c3 + a * y.sqrt()) / mu_km3_s2.sqrt();

        if cur_tof < tof_s {
            phi_lower = phi;
        } else {
            phi_upper = phi;
        }

        phi = (phi_upper + phi_lower) / 2.0;

        if phi > LAMBERT_EPSILON {
            let sqrt_phi = phi.sqrt();
            let (s_sphi, c_sphi) = sqrt_phi.sin_cos();
            c2 = (1.0 - c_sphi) / phi;
            c3 = (sqrt_phi - s_sphi) / phi.powi(3).sqrt();
        } else if phi < -LAMBERT_EPSILON {
            let sqrt_phi = (-phi).sqrt();
            c2 = (1.0 - sqrt_phi.cosh()) / phi;
            c3 = (sqrt_phi.sinh() - sqrt_phi) / (-phi).powi(3).sqrt();
        } else {
            c2 = 0.5;
            c3 = 1.0 / 6.0;
        }
    }

    let f = 1.0 - y / r_init_norm;
    let g_dot = 1.0 - y / r_final_norm;
    let g = a * (y / mu_km3_s2).sqrt();

    Ok(LambertSolution {
        v_init: (r_final - f * r_init) / g,
        v_final: (1.0 / g) * (g_dot * r_final - r_init),
        phi,
    })
}

#[test]
fn test_lambert_vallado_shortway() {
    let ri = Vector3::new(15945.34, 0.0, 0.0);
    let rf = Vector3::new(12214.83899, 10249.46731, 0.0);
    let tof_s = 76.0 * 60.0;
    let mu_km3_s2 = 3.98600433e5;

    let exp_vi = Vector3::new(2.058913, 2.915965, 0.0);
    let exp_vf = Vector3::new(-3.451565, 0.910315, 0.0);

    let sol = gooding(ri, rf, tof_s, mu_km3_s2, TransferKind::ShortWay).unwrap();

    assert!((sol.v_init - exp_vi).norm() < 1e-6);
    assert!((sol.v_final - exp_vf).norm() < 1e-6);
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

#[test]
fn test_lambert_izzo_multi_rev() {
    // Test case from Izzo 2015 paper, Table 1, test case 10
    let r1 = Vector3::new(1.0, 0.0, 0.0);
    let r2 = Vector3::new(1.0, 1.0, 0.0);
    let tof = 2.5;
    let mu = 1.0;

    let exp_v1 = Vector3::new(0.345809286638515, 1.13280141942031, 0.0);
    let exp_v2 = Vector3::new(-0.898583134103195, 0.509849202315269, 0.0);

    let sol = izzo(r1, r2, tof, mu, TransferKind::NRevs(1)).unwrap();

    assert!((sol.v_init - exp_v1).norm() < 1e-5);
    assert!((sol.v_final - exp_v2).norm() < 1e-5);
}

#[test]
fn test_lambert_izzo_retrograde() {
    // Test case from Izzo 2015 paper, Table 1, test case 14
    let r1 = Vector3::new(1.0, 0.0, 0.0);
    let r2 = Vector3::new(1.0, 0.1, 0.0);
    let tof = 0.1;
    let mu = 1.0;

    let exp_v1 = Vector3::new(0.099512313361485, -1.98048259275143, 0.0);
    let exp_v2 = Vector3::new(-0.099512313361485, -1.98048259275143, 0.0);

    let sol = izzo(r1, r2, tof, mu, TransferKind::Auto).unwrap();

    assert!((sol.v_init - exp_v1).norm() < 1e-5);
    assert!((sol.v_final - exp_v2).norm() < 1e-5);
}

/// Solve the Lambert boundary problem using Izzo's method.
///
/// This is an implementation of D. Izzo's method for solving Lambert's problem, as described in "Revisiting Lambert’s problem".
/// The code was adapted from the Python version available in ESA's lamberthub, which is released under the MIT license.
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
) -> Result<LambertSolution, NyxError> {
    let r1_norm = r_init.norm();
    let r2_norm = r_final.norm();

    let cos_dnu = r_init.dot(&r_final) / (r1_norm * r2_norm);

    let dm = kind.direction_of_motion(&r_final, &r_init)?;

    let c = (r1_norm.powi(2) + r2_norm.powi(2) - 2.0 * r1_norm * r2_norm * cos_dnu).sqrt();
    let s = (r1_norm + r2_norm + c) / 2.0;
    let T = (mu_km3_s2 / (8.0 * s.powi(3))).sqrt() * tof_s;

    let lambda = (1.0 - c / s).sqrt();

    let x;
    let y;

    if dm < 0.0 {
        // Long way
        let x0 = (0.5 * PI).tan();
        let T0 = (1.0 / 3.0) * (1.0 - lambda.powi(3));
        let T1 = (2.0 / 3.0) * (1.0 - lambda);

        if T >= T0 {
            x = (T0 / T).powf(2.0 / 3.0) - 1.0;
        } else if T <= T1 {
            x = 5.0 / 2.0 * T1 / T * (T1 - T) / (1.0 - lambda.powi(5)) + 1.0;
        } else {
            x = (T0 / T).log2() * (x0 - (T0 / T1).powf(2.0 / 3.0) + 1.0)
                + (T0 / T1).powf(2.0 / 3.0)
                - 1.0;
        }
    } else {
        // Short way
        let x0 = (0.5 * PI).tan();
        let T0 = (1.0 - lambda.powi(3)) / 3.0;
        let T1 = (2.0 / 3.0) * (1.0 - lambda);

        if T >= T0 {
            x = (T0 / T).powf(2.0 / 3.0) - 1.0;
        } else if T <= T1 {
            x = 5.0 / 2.0 * T1 / T * (T1 - T) / (1.0 - lambda.powi(5)) + 1.0;
        } else {
            x = (T0 / T).log2() * (x0 - (T0 / T1).powf(2.0 / 3.0) + 1.0)
                + (T0 / T1).powf(2.0 / 3.0)
                - 1.0;
        }
    }

    let mut x_val = x;
    for _ in 0..MAX_ITERATIONS {
        let tof_x = tof(x_val, lambda, dm);

        if ((tof_x - T) / T).abs() < LAMBERT_EPSILON_TIME {
            break;
        }

        let dtdx = dtof_dx(x_val, lambda, dm);
        x_val -= (tof_x - T) / dtdx;
    }

    y = ((1.0 - lambda.powi(2) * (1.0 - x_val.powi(2)))).sqrt();

    let vr1 = 1.0 / (1.0 - lambda) * (y - x_val * lambda);
    let vt1 = (2.0 * (1.0 + x_val * y * lambda + x_val * y / lambda) / (1.0 + x_val)).sqrt();
    let vt2 = c / r2_norm * vt1;
    let vr2 = (vt1 - vt2) / (c / r1_norm);

    let v_init = Vector3::new(vr1, vt1, 0.0);
    let v_final = Vector3::new(vr2, vt2, 0.0);

    Ok(LambertSolution {
        v_init,
        v_final,
        phi: 0.0,
    })
}

fn tof(x: f64, lambda: f64, dm: f64) -> f64 {
    let a = 1.0 / (1.0 - x.powi(2));
    if a > 0.0 {
        // Elliptic
        let alpha = 2.0 * a.acos();
        let beta = dm * 2.0 * (a * (1.0 - lambda.powi(2))).acos();
        (a.powf(1.5) * (alpha - (alpha).sin()) - (a.powf(1.5) * (beta - (beta).sin())))
    } else {
        // Hyperbolic
        let alpha = 2.0 * (-a).acosh();
        let beta = dm * 2.0 * (-a * (1.0 - lambda.powi(2))).acosh();
        (a.powf(1.5) * ((alpha).sinh() - alpha) - (a.powf(1.5) * ((beta).sinh() - beta)))
    }
}

fn dtof_dx(x: f64, lambda: f64, dm: f64) -> f64 {
    let a = 1.0 / (1.0 - x.powi(2));
    let y = (1.0 - lambda.powi(2) * (1.0 - x.powi(2))).sqrt();

    if a > 0.0 {
        // Elliptic
        3.0 * tof(x, lambda, dm) * a / (1.0 - x.powi(2))
            - 2.0 * x * a.powf(1.5)
            + dm * 2.0 * lambda * a.powf(1.5) / y
    } else {
        // Hyperbolic
        3.0 * tof(x, lambda, dm) * a / (1.0 - x.powi(2))
            - 2.0 * x * a.powf(1.5)
            + dm * 2.0 * lambda * a.powf(1.5) / y
    }
}

#[test]
fn test_lambert_vallado_lonway() {
    let ri = Vector3::new(15945.34, 0.0, 0.0);
    let rf = Vector3::new(12214.83899, 10249.46731, 0.0);
    let tof_s = 76.0 * 60.0;
    let mu_km3_s2 = 3.98600433e5;

    let exp_vi = Vector3::new(-3.811158, -2.003854, 0.0);
    let exp_vf = Vector3::new(4.207569, 0.914724, 0.0);

    let sol = gooding(ri, rf, tof_s, mu_km3_s2, TransferKind::LongWay).unwrap();

    assert!((sol.v_init - exp_vi).norm() < 1e-6);
    assert!((sol.v_final - exp_vf).norm() < 1e-6);
}
