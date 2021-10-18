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

use super::hyperdual::linalg::norm;
use super::hyperdual::{Float, Hyperdual};
use super::{Frame, Orbit, OrbitDual, OrbitPartial};
use crate::linalg::{Matrix2, Matrix3, Vector2, Vector3};
use crate::md::targeter::Objective;
use crate::md::StateParameter;
use crate::time::{Duration, Epoch, TimeUnit};
use crate::utils::between_pm_180;
use crate::NyxError;

use std::convert::From;
use std::fmt;

/// Stores a B-Plane
#[derive(Copy, Clone, Debug)]
pub struct BPlane {
    /// The $B_T$ component, in kilometers
    pub b_t: OrbitPartial,
    /// The $B_R$ component, in kilometers
    pub b_r: OrbitPartial,
    /// The Linearized Time of Flight
    pub ltof_s: OrbitPartial,
    /// The B-Plane rotation matrix
    pub str_dcm: Matrix3<f64>,
    /// The frame in which this B Plane was computed
    pub frame: Frame,
    /// The time of computation
    pub epoch: Epoch,
}

impl BPlane {
    /// Returns a newly define B-Plane if the orbit is hyperbolic and already in Dual form
    pub fn from_dual(orbit: OrbitDual) -> Result<Self, NyxError> {
        if orbit.ecc().real() <= 1.0 {
            Err(NyxError::NotHyperbolic(
                "Orbit is not hyperbolic. Convert to target object first".to_string(),
            ))
        } else {
            let one = Hyperdual::from(1.0);
            let zero = Hyperdual::from(0.0);

            let e_hat = orbit.evec() / orbit.ecc().dual;
            let h_hat = orbit.hvec() / orbit.hmag().dual;
            let n_hat = h_hat.cross(&e_hat);

            // The reals implementation (which was initially validated) was:
            // let s = e_hat / orbit.ecc() + (1.0 - (1.0 / orbit.ecc()).powi(2)).sqrt() * n_hat;
            // let s_hat = s / s.norm();

            let incoming_asymptote_fact = (one - (one / orbit.ecc().dual).powi(2)).sqrt();

            let s = Vector3::new(
                e_hat[0] / orbit.ecc().dual + incoming_asymptote_fact * n_hat[0],
                e_hat[1] / orbit.ecc().dual + incoming_asymptote_fact * n_hat[1],
                e_hat[2] / orbit.ecc().dual + incoming_asymptote_fact * n_hat[2],
            );

            let s_hat = s / norm(&s); // Just to make sure to renormalize everything

            // The reals implementation (which was initially validated) was:
            // let b_vec = orbit.semi_minor_axis()
            //     * ((1.0 - (1.0 / orbit.ecc()).powi(2)).sqrt() * e_hat
            //         - (1.0 / orbit.ecc() * n_hat));
            let b_vec = Vector3::new(
                orbit.semi_minor_axis().dual
                    * (incoming_asymptote_fact * e_hat[0] - ((one / orbit.ecc().dual) * n_hat[0])),
                orbit.semi_minor_axis().dual
                    * (incoming_asymptote_fact * e_hat[1] - ((one / orbit.ecc().dual) * n_hat[1])),
                orbit.semi_minor_axis().dual
                    * (incoming_asymptote_fact * e_hat[2] - ((one / orbit.ecc().dual) * n_hat[2])),
            );

            let t = s_hat.cross(&Vector3::new(zero, zero, one));
            let t_hat = t / norm(&t);
            let r_hat = s_hat.cross(&t_hat);

            // Build the rotation matrix from inertial to B Plane
            let str_rot = Matrix3::new(
                s_hat[0].real(),
                s_hat[1].real(),
                s_hat[2].real(),
                t_hat[0].real(),
                t_hat[1].real(),
                t_hat[2].real(),
                r_hat[0].real(),
                r_hat[1].real(),
                r_hat[2].real(),
            );

            Ok(BPlane {
                b_r: OrbitPartial {
                    dual: b_vec.dot(&r_hat),
                    param: StateParameter::BdotR,
                },
                b_t: OrbitPartial {
                    dual: b_vec.dot(&t_hat),
                    param: StateParameter::BdotT,
                },
                ltof_s: OrbitPartial {
                    dual: b_vec.dot(&s_hat) / orbit.vmag().dual,
                    param: StateParameter::BLTOF,
                },
                str_dcm: str_rot,
                frame: orbit.frame,
                epoch: orbit.dt,
            })
        }
    }

    /// Returns a newly defined B-Plane if the orbit is hyperbolic.
    pub fn new(orbit: Orbit) -> Result<Self, NyxError> {
        // Convert to OrbitDual so we can target it
        Self::from_dual(OrbitDual::from(orbit))
    }

    pub fn b_dot_t(&self) -> f64 {
        self.b_t.real()
    }

    pub fn b_dot_r(&self) -> f64 {
        self.b_r.real()
    }

    pub fn ltof(&self) -> Duration {
        self.ltof_s.real() * TimeUnit::Second
    }

    /// Returns the B plane angle in degrees between -180 and 180
    pub fn angle(&self) -> f64 {
        between_pm_180(self.b_dot_r().atan2(self.b_dot_t()).to_degrees())
    }

    /// Returns the B plane vector magnitude, in kilometers
    pub fn mag(&self) -> f64 {
        (self.b_dot_t().powi(2) + self.b_dot_r().powi(2)).sqrt()
    }

    /// Returns the DCM to convert to the B Plane from the inertial frame
    pub fn inertial_to_bplane(&self) -> Matrix3<f64> {
        self.str_dcm
    }

    /// Returns the Jacobian of the B plane (BT, BR, LTOF) with respect to the velocity
    pub fn jacobian(&self) -> Matrix3<f64> {
        Matrix3::new(
            self.b_r.wtr_vx(),
            self.b_r.wtr_vy(),
            self.b_r.wtr_vz(),
            self.b_t.wtr_vx(),
            self.b_t.wtr_vy(),
            self.b_t.wtr_vz(),
            self.ltof_s.wtr_vx(),
            self.ltof_s.wtr_vy(),
            self.ltof_s.wtr_vz(),
        )
    }

    /// Returns the Jacobian of the B plane (BT, BR) with respect to two of the velocity components
    pub fn jacobian2(&self, invariant: StateParameter) -> Result<Matrix2<f64>, NyxError> {
        match invariant {
            StateParameter::VX => Ok(Matrix2::new(
                self.b_t.wtr_vy(),
                self.b_t.wtr_vz(),
                self.b_r.wtr_vy(),
                self.b_r.wtr_vz(),
            )),
            StateParameter::VY => Ok(Matrix2::new(
                self.b_t.wtr_vx(),
                self.b_t.wtr_vz(),
                self.b_r.wtr_vx(),
                self.b_r.wtr_vz(),
            )),
            StateParameter::VZ => Ok(Matrix2::new(
                self.b_t.wtr_vx(),
                self.b_t.wtr_vy(),
                self.b_r.wtr_vx(),
                self.b_r.wtr_vy(),
            )),
            _ => Err(NyxError::CustomError(
                "B Plane jacobian invariant must be either VX, VY or VZ".to_string(),
            )),
        }
    }
}

impl fmt::Display for BPlane {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "[{}] {} B-Plane: B∙R = {:.3} km\tB∙T = {:.3} km\tAngle = {:.3} deg",
            self.frame,
            self.epoch,
            self.b_dot_r(),
            self.b_dot_t(),
            self.angle()
        )
    }
}

#[derive(Copy, Clone, Debug)]
pub struct BPlaneTarget {
    /// The $B_T$ component, in kilometers
    pub b_t_km: f64,
    /// The $B_R$ component, in kilometers
    pub b_r_km: f64,
    /// The Linearized Time of Flight, in seconds
    pub ltof_s: f64,

    /// The tolerance on the $B_T$ component, in kilometers
    pub tol_b_t_km: f64,
    /// The tolerance on the $B_R$ component, in kilometers
    pub tol_b_r_km: f64,
    /// The tolerance on the Linearized Time of Flight, in seconds
    pub tol_ltof_s: f64,
}

impl BPlaneTarget {
    /// Initializes a new B Plane target with only the targets and the default tolerances.
    /// Default tolerances are 1 millimeter in positions and 1 second in LTOF
    pub fn from_targets(b_r_km: f64, b_t_km: f64, ltof: Duration) -> Self {
        let tol_ltof: Duration = 6.0 * TimeUnit::Hour;
        Self {
            b_t_km,
            b_r_km,
            ltof_s: ltof.in_seconds(),
            tol_b_t_km: 1e-6,
            tol_b_r_km: 1e-6,
            tol_ltof_s: tol_ltof.in_seconds(),
        }
    }

    /// Initializes a new B Plane target with only the B Plane targets (not LTOF constraint) and the default tolerances.
    /// Default tolerances are 1 millimeter in positions. Here, the LTOF tolerance is set to 100 days.
    pub fn from_bt_br(b_t_km: f64, b_r_km: f64) -> Self {
        let ltof_tol: Duration = 100 * TimeUnit::Day;
        Self {
            b_t_km,
            b_r_km,
            ltof_s: 0.0,
            tol_b_t_km: 1e-6,
            tol_b_r_km: 1e-6,
            tol_ltof_s: ltof_tol.in_seconds(),
        }
    }

    pub fn ltof_target_set(&self) -> bool {
        self.ltof_s.abs() > 1e-10
    }

    pub fn to_objectives(self) -> Vec<Objective> {
        self.to_objectives_with_tolerance(1.0)
    }

    pub fn to_objectives_with_tolerance(self, tol_km: f64) -> Vec<Objective> {
        let mut objs = vec![
            Objective::within_tolerance(StateParameter::BdotR, self.b_r_km, tol_km),
            Objective::within_tolerance(StateParameter::BdotT, self.b_t_km, tol_km),
        ];

        if self.ltof_s.abs() > std::f64::EPSILON {
            objs.push(Objective::within_tolerance(
                StateParameter::BLTOF,
                self.ltof_s,
                self.tol_ltof_s * 1e5,
            ));
        }

        objs
    }
}

impl fmt::Display for BPlaneTarget {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "B-Plane target: B∙R = {:.3} km (+/- {:.1} m)\tB∙T = {:.3} km (+/- {:.1} m)",
            self.b_r_km,
            self.tol_b_r_km * 1e-3,
            self.b_t_km,
            self.tol_b_t_km * 1e-3,
        )
    }
}

/// Returns the Delta V (in km/s) needed to achieve the B Plane specified by B dot R and B dot T.
/// If no LTOF target is set, this method will fix VX, VY and VZ successively and use the minimum of those as a seed for the LTOF variation finding.
/// If the 3x3 search is worse than any of the 2x2s, then a 2x2 will be returned.
/// This uses the hyperdual formulation of the Jacobian and will also vary the linearize time of flight (LTOF).
pub fn try_achieve_b_plane(
    orbit: Orbit,
    target: BPlaneTarget,
) -> Result<(Vector3<f64>, BPlane), NyxError> {
    let mut total_dv = Vector3::zeros();
    let mut attempt_no = 0;
    let max_iter = 10;

    let mut real_orbit = orbit;
    let mut prev_b_plane_err = std::f64::INFINITY;

    if !target.ltof_target_set() {
        // If no LTOF is targeted, we'll solve this with a least squared approach.
        loop {
            if attempt_no > max_iter {
                return Err(NyxError::MaxIterReached(format!(
                    "Error norm of {} km after {} iterations",
                    prev_b_plane_err, max_iter
                )));
            }

            // Build current B Plane
            let b_plane = BPlane::new(real_orbit)?;

            // Check convergence
            let br_err = target.b_r_km - b_plane.b_dot_r();
            let bt_err = target.b_t_km - b_plane.b_dot_t();

            if br_err.abs() < target.tol_b_r_km && bt_err.abs() < target.tol_b_t_km {
                return Ok((total_dv, b_plane));
            }

            // Build the error vector
            let b_plane_err = Vector2::new(br_err, bt_err);

            if b_plane_err.norm() >= prev_b_plane_err {
                // If the error is not going down, we'll raise an error
                return Err(NyxError::CorrectionIneffective(
                    format!("Delta-V correction is ineffective at reducing the B-Plane error:\nprev err norm: {:.3} km\tcur err norm: {:.3} km", prev_b_plane_err, b_plane_err.norm())
                ));
            }
            prev_b_plane_err = b_plane_err.norm();

            // Grab the first two rows of the Jacobian (discard the rest).
            let full_jac = b_plane.jacobian();
            let jac = full_jac.fixed_rows::<2>(0);
            // Solve the Least Squares / compute the delta-v
            let dv = jac.transpose() * (jac * jac.transpose()).try_inverse().unwrap() * b_plane_err;

            total_dv[0] += dv[0];
            total_dv[1] += dv[1];
            total_dv[2] += dv[2];

            // Rebuild a new orbit
            real_orbit.vx += dv[0];
            real_orbit.vy += dv[1];
            real_orbit.vz += dv[2];

            attempt_no += 1;
        }
    } else {
        // The LTOF targeting seems to break often, but it's still implemented
        loop {
            if attempt_no > max_iter {
                return Err(NyxError::MaxIterReached(format!(
                    "Error norm of {} km after {} iterations",
                    prev_b_plane_err, max_iter
                )));
            }

            // Build current B Plane
            let b_plane = BPlane::new(real_orbit)?;

            // Check convergence
            let br_err = target.b_r_km - b_plane.b_dot_r();
            let bt_err = target.b_t_km - b_plane.b_dot_t();
            let ltof_err = target.ltof_s - b_plane.ltof_s.real();

            if br_err.abs() < target.tol_b_r_km
                && bt_err.abs() < target.tol_b_t_km
                && ltof_err.abs() < target.tol_ltof_s
            {
                return Ok((total_dv, b_plane));
            }

            // Build the error vector
            let b_plane_err = Vector3::new(bt_err, br_err, ltof_err);

            if b_plane_err.norm() >= prev_b_plane_err {
                return Err(NyxError::CorrectionIneffective(
                    format!("LTOF enabled correction is failing. Try to not set an LTOF target. Delta-V correction is ineffective are reducing the B-Plane error:\nprev err norm: {:.3} km\tcur err norm: {:.3} km", prev_b_plane_err, b_plane_err.norm()),
                ));
            }
            prev_b_plane_err = b_plane_err.norm();

            // Compute the delta-v
            let dv = b_plane.jacobian() * b_plane_err;

            total_dv[0] += dv[0];
            total_dv[1] += dv[1];
            total_dv[2] += dv[2];

            // Rebuild a new orbit
            real_orbit.vx += dv[0];
            real_orbit.vy += dv[1];
            real_orbit.vz += dv[2];

            attempt_no += 1;
        }
    }
}
