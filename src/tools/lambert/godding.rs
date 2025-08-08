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

use super::{
    LambertError, LambertInput, LambertSolution, TransferKind, LAMBERT_EPSILON,
    LAMBERT_EPSILON_RAD, LAMBERT_EPSILON_TIME, MAX_ITERATIONS,
};

use core::f64::consts::PI;

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
pub fn gooding(input: LambertInput, kind: TransferKind) -> Result<LambertSolution, LambertError> {
    let r_init = input.initial_state.radius_km;
    let r_final = input.final_state.radius_km;
    let tof_s = (input.final_state.epoch - input.initial_state.epoch).to_seconds();
    let mu_km3_s2 = input.mu_km2_s3();

    let r_init_norm = r_init.norm();
    let r_final_norm = r_final.norm();
    let r_norm_product = r_init_norm * r_final_norm;
    let cos_dnu = r_init.dot(&r_final) / r_norm_product;

    let dm = kind.direction_of_motion(&r_final, &r_init)?;

    let nu_init = r_init[1].atan2(r_init[0]);
    let nu_final = r_final[1].atan2(r_final[0]);

    let a = dm * (r_norm_product * (1.0 + cos_dnu)).sqrt();

    if nu_final - nu_init < LAMBERT_EPSILON_RAD && a.abs() < LAMBERT_EPSILON {
        return Err(LambertError::TargetsTooClose);
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
            return Err(LambertError::SolverMaxIter {
                maxiter: MAX_ITERATIONS,
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
                return Err(LambertError::NotReasonablePhi);
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
        v_init_km_s: (r_final - f * r_init) / g,
        v_final_km_s: (1.0 / g) * (g_dot * r_final - r_init),
        phi_rad: phi,
        input,
    })
}

#[cfg(test)]
mod ut_lambert_gooding {

    use super::{gooding, TransferKind};
    use crate::{linalg::Vector3, tools::lambert::LambertInput};
    use anise::{frames::Frame, prelude::Orbit};
    use hifitime::{Epoch, Unit};
    #[test]
    fn test_lambert_vallado_shortway() {
        let frame = Frame {
            ephemeris_id: 301,
            orientation_id: 0,
            mu_km3_s2: Some(3.98600433e5),
            shape: None,
        };
        let initial_state = Orbit {
            radius_km: Vector3::new(15945.34, 0.0, 0.0),
            velocity_km_s: Vector3::zeros(),
            epoch: Epoch::from_gregorian_utc_at_midnight(2025, 1, 1),
            frame,
        };
        let final_state = Orbit {
            radius_km: Vector3::new(12214.83899, 10249.46731, 0.0),
            velocity_km_s: Vector3::zeros(),
            epoch: Epoch::from_gregorian_utc_at_midnight(2025, 1, 1) + Unit::Minute * 76.0,
            frame,
        };

        let input = LambertInput::from_planetary_states(initial_state, final_state).unwrap();

        let exp_vi = Vector3::new(2.058913, 2.915965, 0.0);
        let exp_vf = Vector3::new(-3.451565, 0.910315, 0.0);

        let sol = gooding(input, TransferKind::ShortWay).unwrap();

        assert!((sol.v_init_km_s - exp_vi).norm() < 1e-6);
        assert!((sol.v_final_km_s - exp_vf).norm() < 1e-6);

        println!("{sol:?}");
        println!("{}\n{}", sol.transfer_orbit(), sol.arrival_orbit());
        println!(
            "v_inf+ = {} km/s -- c3 = {} km^2/s^2",
            sol.v_inf_outgoing_km_s().norm(),
            sol.c3_km2_s2()
        );
        println!("v_inf- = {} km/s", sol.v_inf_incoming_km_s().norm());
    }

    #[test]
    fn test_lambert_vallado_longway() {
        let frame = Frame {
            ephemeris_id: 301,
            orientation_id: 0,
            mu_km3_s2: Some(3.98600433e5),
            shape: None,
        };
        let initial_state = Orbit {
            radius_km: Vector3::new(15945.34, 0.0, 0.0),
            velocity_km_s: Vector3::zeros(),
            epoch: Epoch::from_gregorian_utc_at_midnight(2025, 1, 1),
            frame,
        };
        let final_state = Orbit {
            radius_km: Vector3::new(12214.83899, 10249.46731, 0.0),
            velocity_km_s: Vector3::zeros(),
            epoch: Epoch::from_gregorian_utc_at_midnight(2025, 1, 1) + Unit::Minute * 76.0,
            frame,
        };

        let input = LambertInput::from_planetary_states(initial_state, final_state).unwrap();
        let exp_vi = Vector3::new(-3.811158, -2.003854, 0.0);
        let exp_vf = Vector3::new(4.207569, 0.914724, 0.0);

        let sol = gooding(input, TransferKind::LongWay).unwrap();

        assert!((sol.v_init_km_s - exp_vi).norm() < 1e-6);
        assert!((sol.v_final_km_s - exp_vf).norm() < 1e-6);
    }
}
