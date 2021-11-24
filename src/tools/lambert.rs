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

use crate::errors::NyxError;
use crate::linalg::Vector3;
use std::f64::consts::PI;

const TAU: f64 = 2.0 * PI;
const LAMBERT_EPSILON: f64 = 1e-4; // General epsilon
const LAMBERT_EPSILON_TIME: f64 = 1e-4; // Time epsilon
const LAMBERT_EPSILON_RAD: f64 = (5e-5 / 180.0) * PI; // 0.00005 degrees

/// Define the transfer kind for a Lambert
pub enum TransferKind {
    Auto,
    ShortWay,
    LongWay,
    NRevs(u8),
}

#[derive(Debug)]
pub struct LambertSolution {
    pub v_init: Vector3<f64>,
    pub v_final: Vector3<f64>,
    pub phi: f64,
}

/// Solves the Lambert boundary problem using a standard secant method.
/// Given the initial and final radii, a time of flight, and a gravitational parameters, it returns the needed initial and final velocities
/// along with φ which is the square of the difference in eccentric anomaly. Note that the direction of motion
/// is computed directly in this function to simplify the generation of Pork chop plots.
pub fn standard(
    r_init: Vector3<f64>,
    r_final: Vector3<f64>,
    tof: f64,
    gm: f64,
    kind: TransferKind,
) -> Result<LambertSolution, NyxError> {
    let r_init_norm = r_init.norm();
    let r_final_norm = r_final.norm();

    let cos_dnu = r_init.dot(&r_final) / (r_init_norm * r_final_norm);

    let dm = match kind {
        TransferKind::Auto => {
            let mut dnu = r_final[1].atan2(r_final[0]) - r_init[1].atan2(r_final[1]);
            if dnu > TAU {
                dnu -= TAU;
            } else if dnu < 0.0 {
                dnu += TAU;
            }

            if dnu > std::f64::consts::PI {
                -1.0
            } else {
                1.0
            }
        }
        TransferKind::ShortWay => 1.0,
        TransferKind::LongWay => -1.0,
        _ => return Err(NyxError::LambertMultiRevNotSupported),
    };

    // Compute the direction of motion
    let nu_init = r_init[1].atan2(r_init[0]);
    let nu_final = r_final[1].atan2(r_final[0]);

    let a = dm * (r_init_norm * r_final_norm * (1.0 + cos_dnu)).sqrt();

    if nu_final - nu_init < LAMBERT_EPSILON_RAD && a.abs() < LAMBERT_EPSILON {
        return Err(NyxError::TargetsTooClose);
    }

    // Define the search space (note that we do not support multirevs in this algorithm)
    let mut phi_upper = 4.0 * PI.powi(2);
    let mut phi_lower = -4.0 * PI.powi(2);
    let mut phi = 0.0; // ??!?

    // Initial guesses for c2 and c3
    let mut c2: f64 = 1.0 / 2.0;
    let mut c3: f64 = 1.0 / 6.0;
    let mut iter: usize = 0;
    let mut cur_tof: f64 = 0.0;
    let mut y = 0.0;

    while (cur_tof - tof).abs() > LAMBERT_EPSILON_TIME {
        if iter > 1000 {
            return Err(NyxError::MaxIterReached(format!(
                "Lambert solver failed after {} iterations",
                1000
            )));
        }
        iter += 1;

        y = r_init_norm + r_final_norm + a * (phi * c3 - 1.0) / c2.sqrt();
        if a > 0.0 && y < 0.0 {
            // Try to increase phi
            for _ in 0..500 {
                phi += 0.1;
                // Recompute y
                y = r_init_norm + r_final_norm + a * (phi * c3 - 1.0) / c2.sqrt();
                if y >= 0.0 {
                    break;
                }
            }
            if y < 0.0 {
                // If y is still negative, then our attempts have failed.
                return Err(NyxError::LambertNotReasonablePhi);
            }
        }

        let chi = (y / c2).sqrt();
        // Compute the current time of flight
        cur_tof = (chi.powi(3) * c3 + a * y.sqrt()) / gm.sqrt();
        // Update the next TOF we should use
        if cur_tof < tof {
            phi_lower = phi;
        } else {
            phi_upper = phi;
        }

        // Compute the next phi
        phi = (phi_upper + phi_lower) / 2.0;

        // Update c2 and c3
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
            // Reset c2 and c3 and try again
            c2 = 0.5;
            c3 = 1.0 / 6.0;
        }
    }

    // Time of flight matches!

    let f = 1.0 - y / r_init_norm;
    let g_dot = 1.0 - y / r_final_norm;
    let g = a * (y / gm).sqrt();

    // Compute velocities
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
