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

use anise::constants::frames::SUN_J2000;
use anise::errors::OrientationSnafu;
use anise::prelude::Almanac;
use serde::{Deserialize, Serialize};
use serde_dhall::StaticType;
use snafu::ResultExt;
use typed_builder::TypedBuilder;

use crate::cosmic::{AstroPhysicsSnafu, Epoch, Frame, Orbit};
use crate::dynamics::{
    AccelModel, DynamicsAlmanacSnafu, DynamicsAstroSnafu, DynamicsError, DynamicsPlanetarySnafu,
};
use crate::linalg::{Matrix3, Vector3, Vector4, U7};
use hyperdual::linalg::norm;
use hyperdual::{hyperspace_from_vector, OHyperdual};
use std::fmt;
use std::sync::Arc;

/// `SolidTides` implements the solid tide acceleration model.
/// It accounts for the crust deformation due to the Moon and the Sun.
/// Formulas are based on IERS 2010 Conventions.
#[derive(Clone, Debug, Serialize, Deserialize, StaticType, TypedBuilder)]
pub struct SolidTides {
    /// The body-fixed frame of the central body being deformed.
    pub frame: Frame,
    /// 2nd degree Love number
    pub k2: f64,
    /// 3rd degree Love number
    pub k3: f64,
    /// The collection of celestial bodies raising the tide.
    pub perturbers: Vec<TidalPerturber>,
}

#[derive(Clone, Debug, Serialize, Deserialize, StaticType, TypedBuilder)]
pub struct TidalPerturber {
    /// The frame used to resolve the state of the perturber relative to central_frame.
    pub frame: Frame,
    /// Optimization flag: true only if (R_eq / r_j)^4 is large enough to warrant k3.
    /// Set to True for the Earth system
    pub compute_degree_3: bool,
}

impl TidalPerturber {
    fn compute_pert(
        &self,
        epoch: Epoch,
        tidal_model: &SolidTides,
        almanac: &Almanac,
        delta_c: &mut [[f64; 4]; 4],
        delta_s: &mut [[f64; 4]; 4],
    ) -> Result<(), DynamicsError> {
        let radius_km = almanac
            .transform(self.frame, tidal_model.frame, epoch, None)
            .context(DynamicsAlmanacSnafu {
                action: "Moon position in ECEF",
            })?
            .radius_km;

        let r_body = radius_km.norm();
        let s_body = radius_km.x / r_body;
        let t_body = radius_km.y / r_body;
        let u_body = radius_km.z / r_body;

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

        let primary_mu_km3_s2 = tidal_model
            .frame
            .mu_km3_s2()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;

        let primary_eq_radius_km = tidal_model
            .frame
            .mean_equatorial_radius_km()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;

        let secondary_mu_km3_s2 = self
            .frame
            .mu_km3_s2()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;

        let gm_ratio = secondary_mu_km3_s2 / primary_mu_km3_s2;
        let r_ratio = primary_eq_radius_km / r_body;

        let m = if self.compute_degree_3 { 3 } else { 2 };

        for n in 2..=m {
            let kn = if n == 2 {
                tidal_model.k2
            } else {
                tidal_model.k3
            };
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

        Ok(())
    }
}

impl SolidTides {
    /// Initializes solid tides with the Moon and the Sun, where the k3 is only computed for the Moon.
    /// Sets the k2 Love number to 0.3019 and the k3 Love number to 0.093
    pub fn earth_moon_system(
        mut earth_frame: Frame,
        mut moon_frame: Frame,
        almanac: Arc<Almanac>,
    ) -> Result<Arc<Self>, DynamicsError> {
        let mut sun_j2k = almanac
            .frame_info(SUN_J2000)
            .context(DynamicsPlanetarySnafu {
                action: "fetching sun frame",
            })?;

        // Repeat for the Earth and Moon
        for frame in [&mut earth_frame, &mut moon_frame, &mut sun_j2k] {
            if frame.mu_km3_s2.is_none() {
                *frame = almanac.frame_info(*frame).context(DynamicsPlanetarySnafu {
                    action: "fetching sun frame",
                })?;
            }

            // Ensure the gravitational parameter is set.
            frame
                .mu_km3_s2()
                .context(AstroPhysicsSnafu)
                .context(DynamicsAstroSnafu)?;

            // Ensure the equatorial radius is set.
            frame
                .mean_equatorial_radius_km()
                .context(AstroPhysicsSnafu)
                .context(DynamicsAstroSnafu)?;
        }

        let me = Self::builder()
            .k2(0.3019)
            .k3(0.093)
            .frame(earth_frame)
            .perturbers(vec![
                TidalPerturber::builder()
                    .frame(moon_frame)
                    .compute_degree_3(true)
                    .build(),
                TidalPerturber::builder()
                    .frame(sun_j2k)
                    .compute_degree_3(true)
                    .build(),
            ])
            .build();

        Ok(Arc::new(me))
    }

    /// Internal helper to compute tidal delta coefficients
    fn accumulate_deltas(
        &self,
        epoch: Epoch,
        almanac: Arc<Almanac>,
    ) -> Result<([[f64; 4]; 4], [[f64; 4]; 4]), DynamicsError> {
        let mut delta_c = [[0.0f64; 4]; 4];
        let mut delta_s = [[0.0f64; 4]; 4];

        for pert in &self.perturbers {
            pert.compute_pert(epoch, self, &almanac, &mut delta_c, &mut delta_s)?;
        }

        Ok((delta_c, delta_s))
    }
}

#[allow(clippy::needless_range_loop)]
impl AccelModel for SolidTides {
    fn eom(&self, osc: &Orbit, almanac: Arc<Almanac>) -> Result<Vector3<f64>, DynamicsError> {
        let (delta_c, delta_s) = self.accumulate_deltas(osc.epoch, almanac.clone())?;

        // Convert the osculating orbit to the correct frame (needed for multiple harmonic fields)
        let state = almanac
            .transform_to(*osc, self.frame, None)
            .context(DynamicsAlmanacSnafu {
                action: "transforming into solid tides frame",
            })?;

        // Using the GMAT notation, with extra character for ease of highlight
        let r_ = state.rmag_km();
        let s_ = state.radius_km.x / r_;
        let t_ = state.radius_km.y / r_;
        let u_ = state.radius_km.z / r_;

        // Associated Legendre polynomials a_nm (scaled as in sph_harmonics.rs)
        let mut a_nm = [[0.0f64; 6]; 6];
        a_nm[0][0] = 1.0;
        for n in 1..=4 {
            a_nm[n][n] = (1.0 + 1.0 / (2.0 * n as f64)).sqrt() * a_nm[n - 1][n - 1];
        }
        a_nm[1][0] = u_ * 3.0f64.sqrt();
        for n in 1..=4 {
            a_nm[n + 1][n] = (2.0 * n as f64 + 3.0).sqrt() * u_ * a_nm[n][n];
        }

        let b_nm = |n: usize, m: usize| {
            (((2.0 * n as f64 + 1.0) * (2.0 * n as f64 - 1.0))
                / ((n as f64 + m as f64) * (n as f64 - m as f64)))
                .sqrt()
        };
        let c_nm = |n: usize, m: usize| {
            (((2.0 * n as f64 + 1.0) * (n as f64 + m as f64 - 1.0) * (n as f64 - m as f64 - 1.0))
                / ((n as f64 - m as f64) * (n as f64 + m as f64) * (2.0 * n as f64 - 3.0)))
                .sqrt()
        };

        for m in 0..=3 {
            for n in (m + 2)..=4 {
                a_nm[n][m] = u_ * b_nm(n, m) * a_nm[n - 1][m] - c_nm(n, m) * a_nm[n - 2][m];
            }
        }

        let mut r_m = [0.0f64; 4];
        let mut i_m = [0.0f64; 4];
        r_m[0] = 1.0;
        i_m[0] = 0.0;
        for m in 1..=3 {
            r_m[m] = s_ * r_m[m - 1] - t_ * i_m[m - 1];
            i_m[m] = s_ * i_m[m - 1] + t_ * r_m[m - 1];
        }

        let eq_radius_km = self
            .frame
            .mean_equatorial_radius_km()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;

        let mu_km3_s2 = self
            .frame
            .mu_km3_s2()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;

        let rho = eq_radius_km / r_;

        let mut rho_np1 = mu_km3_s2 / r_ * rho;
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
            accel4 += (rho_np1 / eq_radius_km) * sum;
        }

        let accel_ecef = Vector3::new(
            accel4.x + accel4.w * s_,
            accel4.y + accel4.w * t_,
            accel4.z + accel4.w * u_,
        );

        let dcm = almanac
            .rotate(self.frame, osc.frame, osc.epoch)
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
        let (delta_c, delta_s) = self.accumulate_deltas(osc.epoch, almanac.clone())?;

        // Convert the osculating orbit to the correct frame (needed for multiple harmonic fields)
        let state = almanac
            .transform_to(*osc, self.frame, None)
            .context(DynamicsAlmanacSnafu {
                action: "transforming into gravity field frame",
            })?;

        let radius: Vector3<OHyperdual<f64, U7>> = hyperspace_from_vector(&state.radius_km);

        // Using the GMAT notation, with extra character for ease of highlight
        let r_ = norm(&radius);
        let s_ = radius[0] / r_;
        let t_ = radius[1] / r_;
        let u_ = radius[2] / r_;

        // Legendre polynomials recursion in Hyperdual
        let mut a_nm = [[OHyperdual::<f64, U7>::from(0.0); 6]; 6];
        a_nm[0][0] = OHyperdual::from(1.0);
        for n in 1..=4 {
            a_nm[n][n] =
                OHyperdual::from((1.0 + 1.0 / (2.0 * n as f64)).sqrt()) * a_nm[n - 1][n - 1];
        }
        a_nm[1][0] = u_ * OHyperdual::from(3.0f64.sqrt());
        for n in 1..=4 {
            a_nm[n + 1][n] = OHyperdual::from((2.0 * n as f64 + 3.0).sqrt()) * u_ * a_nm[n][n];
        }

        let b_nm = |n: usize, m: usize| {
            (((2.0 * n as f64 + 1.0) * (2.0 * n as f64 - 1.0))
                / ((n as f64 + m as f64) * (n as f64 - m as f64)))
                .sqrt()
        };
        let c_nm = |n: usize, m: usize| {
            (((2.0 * n as f64 + 1.0) * (n as f64 + m as f64 - 1.0) * (n as f64 - m as f64 - 1.0))
                / ((n as f64 - m as f64) * (n as f64 + m as f64) * (2.0 * n as f64 - 3.0)))
                .sqrt()
        };

        for m in 0..=3 {
            for n in (m + 2)..=4 {
                a_nm[n][m] = u_ * OHyperdual::from(b_nm(n, m)) * a_nm[n - 1][m]
                    - OHyperdual::from(c_nm(n, m)) * a_nm[n - 2][m];
            }
        }

        let mut r_m = [OHyperdual::<f64, U7>::from(0.0); 4];
        let mut i_m = [OHyperdual::<f64, U7>::from(0.0); 4];
        r_m[0] = OHyperdual::from(1.0);
        i_m[0] = OHyperdual::from(0.0);
        for m in 1..=3 {
            r_m[m] = s_ * r_m[m - 1] - t_ * i_m[m - 1];
            i_m[m] = s_ * i_m[m - 1] + t_ * r_m[m - 1];
        }

        let real_eq_radius_km = self
            .frame
            .mean_equatorial_radius_km()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;

        let real_mu_km3_s2 = self
            .frame
            .mu_km3_s2()
            .context(AstroPhysicsSnafu)
            .context(DynamicsAstroSnafu)?;

        let eq_radius = OHyperdual::<f64, U7>::from(real_eq_radius_km);
        let rho = eq_radius / r_;
        let mut rho_np1 = OHyperdual::<f64, U7>::from(real_mu_km3_s2) / r_ * rho;

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

        let accel_local = Vector3::new(a0 + a3 * s_, a1 + a3 * t_, a2 + a3 * u_);

        let dcm = almanac
            .rotate(self.frame, osc.frame, osc.epoch)
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
        write!(f, "Solid tides for {}", self.frame)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cosmic::Orbit;
    use anise::constants::frames::{EARTH_J2000, IAU_EARTH_FRAME, IAU_MOON_FRAME};
    use std::str::FromStr;

    #[test]
    fn test_solid_tides_earth() {
        let mut almanac = Almanac::default();
        // Load kernels
        almanac = almanac.load("data/01_planetary/de440s.bsp").unwrap();
        almanac = almanac.load("data/01_planetary/pck08.pca").unwrap();
        let almanac = Arc::new(almanac);

        let epoch = Epoch::from_str("2024-01-01T12:00:00 UTC").unwrap();
        let frame = Frame::from(EARTH_J2000);
        let sc_orbit = Orbit::cartesian(7000.0, 0.0, 0.0, 0.0, 7.5, 0.0, epoch, frame);

        let tides = SolidTides::earth_moon_system(IAU_EARTH_FRAME, IAU_MOON_FRAME, almanac.clone())
            .expect("could not init solid tides");
        let acc = tides.eom(&sc_orbit, almanac.clone()).unwrap();

        println!("Solid tides acceleration: {:?}", acc);
        // Typical solid tide acceleration for LEO is around 1e-9 to 1e-7 km/s^2
        assert!(acc.norm() > 0.0);
        assert!(acc.norm() < 1e-6);

        let (acc_grad, grad) = tides.gradient(&sc_orbit, almanac).unwrap();
        assert!((acc - acc_grad).norm() < 1e-12);
        assert!(grad.norm() > 0.0);
    }
}
