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

/// Solve the Lambert boundary problem using a standard secant method.
///
/// Given the initial and final radii, a time of flight, and a gravitational parameters, it returns the needed initial and final velocities
/// along with φ which is the square of the difference in eccentric anomaly. Note that the direction of motion
/// is computed directly in this function to simplify the generation of Pork chop plots.
///
/// # Arguments
///
/// * `r_init` - The initial radius vector.
/// * `r_final` - The final radius vector.
/// * `tof` - The time of flight.
/// * `gm` - The gravitational parameter.
/// * `kind` - The kind of transfer (auto, short way, long way, or number of revolutions).
///
/// # Returns
///
/// `Result<LambertSolution, NyxError>` - The solution to the Lambert problem or an error if the problem could not be solved.
pub fn standard(
    r_init: Vector3<f64>,
    r_final: Vector3<f64>,
    tof: f64,
    gm: f64,
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

    while (cur_tof - tof).abs() > LAMBERT_EPSILON_TIME {
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
        cur_tof = (chi.powi(3) * c3 + a * y.sqrt()) / gm.sqrt();

        if cur_tof < tof {
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
    let g = a * (y / gm).sqrt();

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
    let gm = 3.98600433e5;

    let exp_vi = Vector3::new(2.058913, 2.915965, 0.0);
    let exp_vf = Vector3::new(-3.451565, 0.910315, 0.0);

    let sol = standard(ri, rf, tof_s, gm, TransferKind::ShortWay).unwrap();

    assert!((sol.v_init - exp_vi).norm() < 1e-6);
    assert!((sol.v_final - exp_vf).norm() < 1e-6);
}

#[test]
fn test_lambert_vallado_lonway() {
    let ri = Vector3::new(15945.34, 0.0, 0.0);
    let rf = Vector3::new(12214.83899, 10249.46731, 0.0);
    let tof_s = 76.0 * 60.0;
    let gm = 3.98600433e5;

    let exp_vi = Vector3::new(-3.811158, -2.003854, 0.0);
    let exp_vf = Vector3::new(4.207569, 0.914724, 0.0);

    let sol = standard(ri, rf, tof_s, gm, TransferKind::LongWay).unwrap();

    assert!((sol.v_init - exp_vi).norm() < 1e-6);
    assert!((sol.v_final - exp_vf).norm() < 1e-6);
}
