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
use super::{Frame, Orbit, OrbitDual};
use crate::dimensions::{Matrix3, Vector3, U7};
use crate::time::{Duration, Epoch, TimeUnit};
use crate::utils::between_0_360;
use crate::NyxError;

use std::convert::From;
use std::fmt;

/// Stores a B-Plane
#[derive(Copy, Clone, Debug)]
pub struct BPlane {
    /// The $B_T$ component, in kilometers
    pub b_t: Hyperdual<f64, U7>,
    /// The $B_R$ component, in kilometers
    pub b_r: Hyperdual<f64, U7>,
    /// The Linearized Time of Flight
    pub ltof_s: Hyperdual<f64, U7>,
    /// The B-Plane rotation matrix
    pub str_dcm: Matrix3<f64>,
    /// The frame in which this B Plane was computed
    pub frame: Frame,
    /// The time of computation
    pub epoch: Epoch,
}

impl BPlane {
    /// Returns a newly defined B-Plane if the orbit is hyperbolic.
    pub fn new(orbit_real: Orbit) -> Result<Self, NyxError> {
        if orbit_real.ecc() <= 1.0 {
            Err(NyxError::NotHyperbolic(
                "Orbit is not hyperbolic. Convert to target object first".to_string(),
            ))
        } else {
            // Convert to OrbitDual so we can target it
            let orbit = OrbitDual::from(orbit_real);

            let one = Hyperdual::from(1.0);
            let zero = Hyperdual::from(0.0);

            let e_hat = orbit.evec() / orbit.ecc();
            let h_hat = orbit.hvec() / orbit.hmag();
            let n_hat = h_hat.cross(&e_hat);

            // The reals implementation (which was initially validation) was:
            // let s = e_hat / orbit.ecc() + (one - (one / orbit.ecc()).powi(2)).sqrt() * n_hat;
            // let s_hat = s / s.norm(); // Just to make sure to renormalize everything

            let s = Vector3::new(
                e_hat[0] / orbit.ecc() + (one - (one / orbit.ecc()).powi(2)).sqrt() * n_hat[0],
                e_hat[1] / orbit.ecc() + (one - (one / orbit.ecc()).powi(2)).sqrt() * n_hat[1],
                e_hat[2] / orbit.ecc() + (one - (one / orbit.ecc()).powi(2)).sqrt() * n_hat[2],
            );
            let s_hat = s / norm(&s); // Just to make sure to renormalize everything

            // The reals implementation (which was initially validation) was:
            // let b_vec = orbit.semi_minor_axis()
            //     * ((one - (one / orbit.ecc()).powi(2)).sqrt() * e_hat
            //         - (one / orbit.ecc() * n_hat));
            let b_vec = Vector3::new(
                orbit.semi_minor_axis()
                    * ((one - (one / orbit.ecc()).powi(2)).sqrt() * e_hat[0]
                        - (one / orbit.ecc() * n_hat[0])),
                orbit.semi_minor_axis()
                    * ((one - (one / orbit.ecc()).powi(2)).sqrt() * e_hat[1]
                        - (one / orbit.ecc() * n_hat[1])),
                orbit.semi_minor_axis()
                    * ((one - (one / orbit.ecc()).powi(2)).sqrt() * e_hat[2]
                        - (one / orbit.ecc() * n_hat[2])),
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

            // Compute the LTOF in seconds
            let f = (one
                + (orbit.vmag().powi(2) / Hyperdual::from(orbit.frame.gm()))
                    * (orbit.semi_parameter()
                        / (one + orbit.ecc() * orbit.ta().to_radians().cos())))
            .acosh();

            let ltof = (Hyperdual::from(orbit.frame.gm()) / orbit.vmag().powi(3)) * (f.sinh() - f);

            Ok(BPlane {
                b_r: b_vec.dot(&r_hat),
                b_t: b_vec.dot(&t_hat),
                ltof_s: ltof,
                str_dcm: str_rot,
                frame: orbit.frame,
                epoch: orbit.dt,
            })
        }
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

    /// Returns the B plane angle in degrees between 0 and 360
    pub fn angle(&self) -> f64 {
        between_0_360(self.b_dot_r().atan2(self.b_dot_t()).to_degrees())
    }

    /// Returns the B plane vector magnitude, in kilometers
    pub fn mag(&self) -> f64 {
        (self.b_dot_t().powi(2) + self.b_dot_r().powi(2)).sqrt()
    }

    /// Returns the DCM to convert to the B Plane from the inertial frame
    pub fn inertial_to_bplane(&self) -> Matrix3<f64> {
        self.str_dcm
    }
}

impl fmt::Display for BPlane {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "[{}] {} B-Plane: B∙R = {:.3} km\tB∙T = {:.3} km",
            self.frame,
            self.epoch,
            self.b_dot_r(),
            self.b_dot_t()
        )
    }
}
