/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::{Distribution, Normal, Rng, Uniform};
use crate::linalg::Vector3;
use crate::NyxError;

/// Returns a unit vector from a normal distribution.
/// Implements the Sphere Point Picking method: https://mathworld.wolfram.com/SpherePointPicking.html
pub fn unit_vector_from_seed<R: Rng>(rng: &mut R) -> Vector3<f64> {
    let distr = Uniform::new_inclusive(0.0, 1.0);
    let u = distr.sample(rng);
    let v = distr.sample(rng);
    let theta = std::f64::consts::TAU * u;
    let phi = (2.0 * v - 1.0).acos();
    Vector3::new(theta.cos() * phi.sin(), theta.sin() * phi.sin(), phi.cos())
}

/// Apply a Normal pointing error to the provided delta-v vector (must be in km/s).
/// Requires a 3-sigma error percentage (e.g. 0.05 is a 5% pointing error) and a pseudo random number generator
/// Returns the delta-v vector with the applying pointing error
pub fn dv_pointing_error<R: Rng>(
    cur_pointing: &Vector3<f64>,
    dv: Vector3<f64>,
    error_prct3s: f64,
    rng: &mut R,
) -> Result<Vector3<f64>, NyxError> {
    if !(0.0..1.0).contains(&error_prct3s) {
        return Err(NyxError::MonteCarlo(format!(
            "Pointing error percentage must be between 0 and 1, got {}",
            error_prct3s
        )));
    }

    let dv_mag = dv.norm();
    if dv_mag < std::f64::EPSILON {
        return Err(NyxError::MonteCarlo(format!(
            "Delta-v vector is nil, cannot apply a pointing error: {}",
            dv
        )));
    }

    let dv_hat = dv / dv_mag;
    let cur_pointing_mag = cur_pointing.norm();
    let cur_angle = (cur_pointing.dot(&dv_hat) / cur_pointing_mag).acos();

    let new_angle = Normal::new(cur_angle, error_prct3s / 3.0)
        .unwrap()
        .sample(rng);
    // Return initial delta-v with pointing error
    Ok(dv_hat * new_angle.cos() * dv_mag)
}

/// Returns the provided Delta V vector with a magnitude and pointing error
pub fn dv_execution_error<R: Rng>(
    cur_pointing: &Vector3<f64>,
    dv: Vector3<f64>,
    pointing_3s: f64,
    mag_3s: f64,
    rng: &mut R,
) -> Result<Vector3<f64>, NyxError> {
    let dv = dv_pointing_error(cur_pointing, dv, pointing_3s, rng)?;

    let new_mag = Normal::new(dv.norm(), mag_3s / 3.0).unwrap().sample(rng);

    Ok(new_mag * (dv / dv.norm()))
}

#[test]
fn test_dv_mag_fixed() {
    use super::thread_rng;
    use crate::cosmic::{Cosm, Orbit};
    use crate::time::Epoch;
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let orbit = Orbit::cartesian(
        -2436.45,
        -2436.45,
        6891.037,
        5.088_611,
        -5.088_611,
        0.0,
        Epoch::from_gregorian_tai_at_noon(2021, 3, 24),
        eme2k,
    );

    let dv_mag_distr = Normal::new(5e-3, 5e-4).unwrap();

    for _ in 0..=1000 {
        let dv_mag = dv_mag_distr.sample(&mut thread_rng());
        let dv_point = unit_vector_from_seed(&mut thread_rng());
        let dv = dv_point * dv_mag;
        let dv_w_err = dv_pointing_error(&orbit.velocity(), dv, 0.1, &mut thread_rng()).unwrap();
        assert!(
            (dv_w_err.norm() - dv_mag) < std::f64::EPSILON,
            "{:.1e}",
            (dv_w_err.norm() - dv_mag)
        );
    }
}
