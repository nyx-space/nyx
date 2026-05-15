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

use anise::constants::frames::{IAU_EARTH_FRAME, MOON_J2000, SUN_J2000};
use anise::errors::OrientationSnafu;
use anise::prelude::Almanac;
use serde::{Deserialize, Serialize};
use serde_dhall::StaticType;
use snafu::ResultExt;

use crate::cosmic::{AstroPhysicsSnafu, Epoch, Frame, Orbit};
use crate::dynamics::{AccelModel, DynamicsAlmanacSnafu, DynamicsAstroSnafu, DynamicsError};
use crate::linalg::{Matrix3, Vector3, Vector4, U7};
use hyperdual::linalg::norm;
use hyperdual::{hyperspace_from_vector, OHyperdual};
use std::fmt;
use std::sync::Arc;

/// `SolidTides` implements the solid tide acceleration model.
/// It accounts for the crust deformation due to the Moon and the Sun.
/// Formulas are based on IERS 2010 Conventions.
#[derive(Clone, Debug, Serialize, Deserialize, StaticType)]
pub struct SolidTides {
    /// Frame in which the gravity field is defined (e.g. IAU Earth)
    pub compute_frame: Frame,
    /// Mean equatorial radius of the central body (km)
    pub eq_radius_km: f64,
    /// Gravitational parameter of the central body (km^3/s^2)
    pub mu_km3_s2: f64,
    /// Degree-2 Love number
    pub k2: f64,
    /// Degree-3 Love number
    pub k3: f64,
    /// Frame for querying Moon position
    pub moon_frame: Frame,
    /// Frame for querying Sun position
    pub sun_frame: Frame,
}

impl Default for SolidTides {
    fn default() -> Self {
        Self {
            compute_frame: IAU_EARTH_FRAME.into(),
            eq_radius_km: 6378.1363, // Default Earth radius
            mu_km3_s2: 398600.4415,  // Default Earth GM
            k2: 0.30190,
            k3: 0.093,
            moon_frame: MOON_J2000.into(),
            sun_frame: SUN_J2000.into(),
        }
    }
}

impl SolidTides {
    /// Create a new SolidTides model for Earth using default values from IERS 2010.
    pub fn new(_almanac: Arc<Almanac>) -> Result<Arc<Self>, DynamicsError> {
        let compute_frame = Frame::from(IAU_EARTH_FRAME);
        let eq_radius_km = compute_frame
            .mean_equatorial_radius_km()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;
        let mu_km3_s2 = compute_frame
            .mu_km3_s2()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;

        Ok(Arc::new(Self {
            compute_frame,
            eq_radius_km,
            mu_km3_s2,
            ..Default::default()
        }))
    }

    /// Internal helper to compute tidal delta coefficients
    fn compute_deltas(
        &self,
        epoch: Epoch,
        almanac: Arc<Almanac>,
    ) -> Result<([[f64; 4]; 4], [[f64; 4]; 4]), DynamicsError> {
        let moon_in_ecef = almanac
            .transform(self.moon_frame, self.compute_frame, epoch, None)
            .context(DynamicsAlmanacSnafu {
                action: "Moon position in ECEF",
            })?;
        let sun_in_ecef = almanac
            .transform(self.sun_frame, self.compute_frame, epoch, None)
            .context(DynamicsAlmanacSnafu {
                action: "Sun position in ECEF",
            })?;

        let moon_gm = self
            .moon_frame
            .mu_km3_s2()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;
        let sun_gm = self
            .sun_frame
            .mu_km3_s2()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;

        let mut delta_c = [[0.0f64; 4]; 4];
        let mut delta_s = [[0.0f64; 4]; 4];

        for (body_pos, body_gm) in [
            (moon_in_ecef.radius_km, moon_gm),
            (sun_in_ecef.radius_km, sun_gm),
        ] {
            let r_body = body_pos.norm();
            let s_body = body_pos.x / r_body;
            let t_body = body_pos.y / r_body;
            let u_body = body_pos.z / r_body;

            let sin_phi = u_body;
            let cos_phi = (1.0 - sin_phi.powi(2)).max(0.0).sqrt();
            let cos_lambda = if cos_phi > 1e-12 {
                s_body / cos_phi
            } else {
                1.0
            };
            let sin_lambda = if cos_phi > 1e-12 {
                t_body / cos_phi
            } else {
                0.0
            };

            // Fully normalized Associated Legendre Polynomials P_nm(sin_phi) for n=2,3
            let p20 = 0.5 * (3.0 * sin_phi.powi(2) - 1.0) * 5.0f64.sqrt();
            let p21 = 3.0 * sin_phi * cos_phi * (5.0 / 3.0f64).sqrt();
            let p22 = 3.0 * cos_phi.powi(2) * (5.0 / 12.0f64).sqrt();

            let p30 = 0.5 * (5.0 * sin_phi.powi(3) - 3.0 * sin_phi) * 7.0f64.sqrt();
            let p31 = 1.5 * (5.0 * sin_phi.powi(2) - 1.0) * cos_phi * (7.0 / 6.0f64).sqrt();
            let p32 = 15.0 * sin_phi * cos_phi.powi(2) * (7.0 / 60.0f64).sqrt();
            let p33 = 15.0 * cos_phi.powi(3) * (7.0 / 360.0f64).sqrt();

            let gm_ratio = body_gm / self.mu_km3_s2;
            let r_ratio = self.eq_radius_km / r_body;

            for n in 2..=3 {
                let kn = if n == 2 { self.k2 } else { self.k3 };
                let common = kn / (2.0 * n as f64 + 1.0) * gm_ratio * r_ratio.powi(n as i32 + 1);

                for m in 0..=n {
                    let p_nm = match (n, m) {
                        (2, 0) => p20,
                        (2, 1) => p21,
                        (2, 2) => p22,
                        (3, 0) => p30,
                        (3, 1) => p31,
                        (3, 2) => p32,
                        (3, 3) => p33,
                        _ => 0.0,
                    };

                    let (cos_ml, sin_ml) = match m {
                        0 => (1.0, 0.0),
                        1 => (cos_lambda, sin_lambda),
                        2 => (
                            cos_lambda.powi(2) - sin_lambda.powi(2),
                            2.0 * sin_lambda * cos_lambda,
                        ),
                        3 => (
                            cos_lambda * (cos_lambda.powi(2) - 3.0 * sin_lambda.powi(2)),
                            sin_lambda * (3.0 * cos_lambda.powi(2) - sin_lambda.powi(2)),
                        ),
                        _ => (0.0, 0.0),
                    };

                    delta_c[n][m] += common * p_nm * cos_ml;
                    delta_s[n][m] += common * p_nm * sin_ml;
                }
            }
        }
        Ok((delta_c, delta_s))
    }
}

impl AccelModel for SolidTides {
    fn eom(&self, osc: &Orbit, almanac: Arc<Almanac>) -> Result<Vector3<f64>, DynamicsError> {
        let (delta_c, delta_s) = self.compute_deltas(osc.epoch, almanac.clone())?;

        let sc_in_ecef = almanac
            .transform_to(*osc, self.compute_frame, None)
            .context(DynamicsAlmanacSnafu {
                action: "transforming sc state to ECEF",
            })?;

        let r_sc = sc_in_ecef.rmag_km();
        let s_sc = sc_in_ecef.radius_km.x / r_sc;
        let t_sc = sc_in_ecef.radius_km.y / r_sc;
        let u_sc = sc_in_ecef.radius_km.z / r_sc;

        // Associated Legendre polynomials a_nm (scaled as in sph_harmonics.rs)
        let mut a_nm = [[0.0f64; 6]; 6];
        a_nm[0][0] = 1.0;
        for n in 1..=4 {
            a_nm[n][n] = (1.0 + 1.0 / (2.0 * n as f64)).sqrt() * a_nm[n - 1][n - 1];
        }
        a_nm[1][0] = u_sc * 3.0f64.sqrt();
        for n in 1..=4 {
            a_nm[n + 1][n] = (2.0 * n as f64 + 3.0).sqrt() * u_sc * a_nm[n][n];
        }

        let b_nm = |n: usize, m: usize| {
            (((2.0 * n as f64 + 1.0) * (2.0 * n as f64 - 1.0)) / ((n as f64 + m as f64) * (n as f64 - m as f64))).sqrt()
        };
        let c_nm = |n: usize, m: usize| {
            (((2.0 * n as f64 + 1.0) * (n as f64 + m as f64 - 1.0) * (n as f64 - m as f64 - 1.0))
                / ((n as f64 - m as f64) * (n as f64 + m as f64) * (2.0 * n as f64 - 3.0)))
                .sqrt()
        };

        for m in 0..=3 {
            for n in (m + 2)..=4 {
                a_nm[n][m] = u_sc * b_nm(n, m) * a_nm[n - 1][m] - c_nm(n, m) * a_nm[n - 2][m];
            }
        }

        let mut r_m = [0.0f64; 4];
        let mut i_m = [0.0f64; 4];
        r_m[0] = 1.0;
        i_m[0] = 0.0;
        for m in 1..=3 {
            r_m[m] = s_sc * r_m[m - 1] - t_sc * i_m[m - 1];
            i_m[m] = s_sc * i_m[m - 1] + t_sc * r_m[m - 1];
        }

        let rho = self.eq_radius_km / r_sc;
        let mut rho_np1 = self.mu_km3_s2 / r_sc * rho;
        let mut accel4 = Vector4::zeros();

        let vr01 = |n: usize, m: usize| {
            let mut val = ((n as f64 - m as f64) * (n as f64 + m as f64 + 1.0)).sqrt();
            if m == 0 {
                val /= 2.0f64.sqrt();
            }
            val
        };
        let vr11 = |n: usize, m: usize| {
            let mut val = (((2.0 * n as f64 + 1.0)
                * (n as f64 + m as f64 + 2.0)
                * (n as f64 + m as f64 + 1.0))
                / (2.0 * n as f64 + 3.0))
                .sqrt();
            if m == 0 {
                val /= 2.0f64.sqrt();
            }
            val
        };

        let sqrt2 = 2.0f64.sqrt();

        for n in 1..=3 {
            rho_np1 *= rho;
            if n < 2 {
                continue;
            } // only degree 2 and 3

            let mut sum = Vector4::zeros();
            for m in 0..=n {
                let c_val = delta_c[n][m];
                let s_val = delta_s[n][m];

                let d_ = (c_val * r_m[m] + s_val * i_m[m]) * sqrt2;
                let e_ = if m == 0 {
                    0.0
                } else {
                    (c_val * r_m[m - 1] + s_val * i_m[m - 1]) * sqrt2
                };
                let f_ = if m == 0 {
                    0.0
                } else {
                    (s_val * r_m[m - 1] - c_val * i_m[m - 1]) * sqrt2
                };

                sum.x += (m as f64) * a_nm[n][m] * e_;
                sum.y += (m as f64) * a_nm[n][m] * f_;
                sum.z += vr01(n, m) * a_nm[n][m + 1] * d_;
                sum.w -= vr11(n, m) * a_nm[n + 1][m + 1] * d_;
            }
            accel4 += (rho_np1 / self.eq_radius_km) * sum;
        }

        let accel_ecef = Vector3::new(
            accel4.x + accel4.w * s_sc,
            accel4.y + accel4.w * t_sc,
            accel4.z + accel4.w * u_sc,
        );

        let dcm = almanac
            .rotate(self.compute_frame, osc.frame, osc.epoch)
            .context(OrientationSnafu {
                action: "rotating accel back to integration frame",
            })
            .context(DynamicsAlmanacSnafu {
                action: "rotating accel back to integration frame",
            })?
            .rot_mat;

        Ok(dcm * accel_ecef)
    }

    fn gradient(
        &self,
        osc: &Orbit,
        almanac: Arc<Almanac>,
    ) -> Result<(Vector3<f64>, Matrix3<f64>), DynamicsError> {
        let (delta_c, delta_s) = self.compute_deltas(osc.epoch, almanac.clone())?;

        let sc_in_ecef = almanac
            .transform_to(*osc, self.compute_frame, None)
            .context(DynamicsAlmanacSnafu {
                action: "transforming sc state to ECEF",
            })?;

        let radius: Vector3<OHyperdual<f64, U7>> = hyperspace_from_vector(&sc_in_ecef.radius_km);
        let r_sc = norm(&radius);
        let s_sc = radius[0] / r_sc;
        let t_sc = radius[1] / r_sc;
        let u_sc = radius[2] / r_sc;

        // Legendre polynomials recursion in Hyperdual
        let mut a_nm = [[OHyperdual::<f64, U7>::from(0.0); 6]; 6];
        a_nm[0][0] = OHyperdual::from(1.0);
        for n in 1..=4 {
            a_nm[n][n] = OHyperdual::from((1.0 + 1.0 / (2.0 * n as f64)).sqrt()) * a_nm[n - 1][n - 1];
        }
        a_nm[1][0] = u_sc * OHyperdual::from(3.0f64.sqrt());
        for n in 1..=4 {
            a_nm[n + 1][n] = OHyperdual::from((2.0 * n as f64 + 3.0).sqrt()) * u_sc * a_nm[n][n];
        }

        let b_nm = |n: usize, m: usize| {
            (((2.0 * n as f64 + 1.0) * (2.0 * n as f64 - 1.0)) / ((n as f64 + m as f64) * (n as f64 - m as f64))).sqrt()
        };
        let c_nm = |n: usize, m: usize| {
            (((2.0 * n as f64 + 1.0) * (n as f64 + m as f64 - 1.0) * (n as f64 - m as f64 - 1.0))
                / ((n as f64 - m as f64) * (n as f64 + m as f64) * (2.0 * n as f64 - 3.0)))
                .sqrt()
        };

        for m in 0..=3 {
            for n in (m + 2)..=4 {
                a_nm[n][m] = u_sc * OHyperdual::from(b_nm(n, m)) * a_nm[n - 1][m]
                    - OHyperdual::from(c_nm(n, m)) * a_nm[n - 2][m];
            }
        }

        let mut r_m = [OHyperdual::<f64, U7>::from(0.0); 4];
        let mut i_m = [OHyperdual::<f64, U7>::from(0.0); 4];
        r_m[0] = OHyperdual::from(1.0);
        i_m[0] = OHyperdual::from(0.0);
        for m in 1..=3 {
            r_m[m] = s_sc * r_m[m - 1] - t_sc * i_m[m - 1];
            i_m[m] = s_sc * i_m[m - 1] + t_sc * r_m[m - 1];
        }

        let eq_radius = OHyperdual::<f64, U7>::from(self.eq_radius_km);
        let rho = eq_radius / r_sc;
        let mut rho_np1 = OHyperdual::<f64, U7>::from(self.mu_km3_s2) / r_sc * rho;

        let mut a0 = OHyperdual::<f64, U7>::from(0.0);
        let mut a1 = OHyperdual::<f64, U7>::from(0.0);
        let mut a2 = OHyperdual::<f64, U7>::from(0.0);
        let mut a3 = OHyperdual::<f64, U7>::from(0.0);

        let vr01 = |n: usize, m: usize| {
            let mut val = ((n as f64 - m as f64) * (n as f64 + m as f64 + 1.0)).sqrt();
            if m == 0 {
                val /= 2.0f64.sqrt();
            }
            val
        };
        let vr11 = |n: usize, m: usize| {
            let mut val = (((2.0 * n as f64 + 1.0)
                * (n as f64 + m as f64 + 2.0)
                * (n as f64 + m as f64 + 1.0))
                / (2.0 * n as f64 + 3.0))
                .sqrt();
            if m == 0 {
                val /= 2.0f64.sqrt();
            }
            val
        };
        let sqrt2 = OHyperdual::<f64, U7>::from(2.0f64.sqrt());

        for n in 1..=3 {
            rho_np1 *= rho;
            if n < 2 {
                continue;
            }

            let mut sum0 = OHyperdual::from(0.0);
            let mut sum1 = OHyperdual::from(0.0);
            let mut sum2 = OHyperdual::from(0.0);
            let mut sum3 = OHyperdual::from(0.0);

            for m in 0..=n {
                let c_val = OHyperdual::from(delta_c[n][m]);
                let s_val = OHyperdual::from(delta_s[n][m]);

                let d_ = (c_val * r_m[m] + s_val * i_m[m]) * sqrt2;
                let e_ = if m == 0 {
                    OHyperdual::from(0.0)
                } else {
                    (c_val * r_m[m - 1] + s_val * i_m[m - 1]) * sqrt2
                };
                let f_ = if m == 0 {
                    OHyperdual::from(0.0)
                } else {
                    (s_val * r_m[m - 1] - c_val * i_m[m - 1]) * sqrt2
                };

                sum0 += OHyperdual::from(m as f64) * a_nm[n][m] * e_;
                sum1 += OHyperdual::from(m as f64) * a_nm[n][m] * f_;
                sum2 += OHyperdual::from(vr01(n, m)) * a_nm[n][m + 1] * d_;
                sum3 += OHyperdual::from(vr11(n, m)) * a_nm[n + 1][m + 1] * d_;
            }
            let rr = rho_np1 / eq_radius;
            a0 += rr * sum0;
            a1 += rr * sum1;
            a2 += rr * sum2;
            a3 -= rr * sum3;
        }

        let accel_local = Vector3::new(a0 + a3 * s_sc, a1 + a3 * t_sc, a2 + a3 * u_sc);

        let dcm = almanac
            .rotate(self.compute_frame, osc.frame, osc.epoch)
            .context(OrientationSnafu {
                action: "rotating accel back to integration frame",
            })
            .context(DynamicsAlmanacSnafu {
                action: "rotating accel back to integration frame",
            })?
            .rot_mat;

        let dx = dcm
            * Vector3::new(
                accel_local[0].real(),
                accel_local[1].real(),
                accel_local[2].real(),
            );

        let mut grad_local = Matrix3::zeros();
        for i in 0..3 {
            for j in 1..4 {
                grad_local[(i, j - 1)] += accel_local[i][j];
            }
        }
        let grad = dcm * grad_local * dcm.transpose();
        Ok((dx, grad))
    }
}

impl fmt::Display for SolidTides {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Solid tides for {}", self.compute_frame)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use anise::constants::frames::EARTH_J2000;
    use crate::cosmic::Orbit;
    use std::str::FromStr;

    #[test]
    #[ignore = "Requires LFS kernels to be loaded"]
    fn test_solid_tides_earth() {
        let mut almanac = Almanac::default();
        // Load kernels
        almanac = almanac.load("data/01_planetary/de440s.bsp").unwrap();
        almanac = almanac.load("data/01_planetary/pck08.pca").unwrap();
        let almanac = Arc::new(almanac);

        let epoch = Epoch::from_str("2024-01-01T12:00:00 UTC").unwrap();
        let frame = Frame::from(EARTH_J2000);
        let sc_orbit = Orbit::cartesian(
            7000.0, 0.0, 0.0,
            0.0, 7.5, 0.0,
            epoch,
            frame
        );

        let tides = SolidTides::new(almanac.clone()).unwrap();
        let acc = tides.eom(&sc_orbit, almanac.clone()).unwrap();

        println!("Solid tides acceleration: {:?}", acc);
        // Typical solid tide acceleration for LEO is around 1e-9 to 1e-7 km/s^2
        assert!(acc.norm() > 0.0);
        assert!(acc.norm() < 1e-6);

        let (acc_grad, grad) = tides.gradient(&sc_orbit, almanac.clone()).unwrap();
        assert!((acc - acc_grad).norm() < 1e-12);
        assert!(grad.norm() > 0.0);
    }
}
