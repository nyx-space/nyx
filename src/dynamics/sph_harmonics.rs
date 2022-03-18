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

use super::hyperdual::linalg::norm;
use super::hyperdual::{hyperspace_from_vector, Float, OHyperdual};
use crate::cosmic::{Cosm, Frame, Orbit};
use crate::dynamics::AccelModel;
use crate::errors::NyxError;
use crate::io::gravity::GravityPotentialStor;
use crate::linalg::{DMatrix, Matrix3, Vector3, U7};
use std::cmp::min;
use std::fmt;
use std::sync::Arc;

#[derive(Clone)]
pub struct Harmonics<S>
where
    S: GravityPotentialStor,
{
    cosm: Arc<Cosm>,
    compute_frame: Frame,
    stor: S,
    a_nm: DMatrix<f64>,
    b_nm: DMatrix<f64>,
    c_nm: DMatrix<f64>,
    vr01: DMatrix<f64>,
    vr11: DMatrix<f64>,
    a_nm_h: DMatrix<OHyperdual<f64, U7>>,
    b_nm_h: DMatrix<OHyperdual<f64, U7>>,
    c_nm_h: DMatrix<OHyperdual<f64, U7>>,
    vr01_h: DMatrix<OHyperdual<f64, U7>>,
    vr11_h: DMatrix<OHyperdual<f64, U7>>,
}

impl<S> Harmonics<S>
where
    S: GravityPotentialStor,
{
    /// Create a new Harmonics dynamical model from the provided gravity potential storage instance.
    pub fn from_stor(compute_frame: Frame, stor: S, cosm: Arc<Cosm>) -> Arc<Self> {
        assert!(
            compute_frame.is_geoid(),
            "harmonics only work around geoids"
        );
        let degree_np2 = stor.max_degree_n() + 2;
        let mut a_nm = DMatrix::from_element(degree_np2 + 1, degree_np2 + 1, 0.0);
        let mut b_nm = DMatrix::from_element(degree_np2, degree_np2, 0.0);
        let mut c_nm = DMatrix::from_element(degree_np2, degree_np2, 0.0);
        let mut vr01 = DMatrix::from_element(degree_np2, degree_np2, 0.0);
        let mut vr11 = DMatrix::from_element(degree_np2, degree_np2, 0.0);

        // Initialize the diagonal elements (not a function of the input)
        a_nm[(0, 0)] = 1.0;
        for n in 1..=degree_np2 {
            let nf64 = n as f64;
            // Diagonal element
            a_nm[(n, n)] = (1.0 + 1.0 / (2.0 * nf64)).sqrt() * a_nm[(n - 1, n - 1)];
        }

        // Pre-compute the B_nm, C_nm, vr01 and vr11 storages
        for n in 0..degree_np2 {
            for m in 0..degree_np2 {
                let nf64 = n as f64;
                let mf64 = m as f64;
                // Compute c_nm, which is B_nm/B_(n-1,m) in Jones' dissertation
                c_nm[(n, m)] = (((2.0 * nf64 + 1.0) * (nf64 + mf64 - 1.0) * (nf64 - mf64 - 1.0))
                    / ((nf64 - mf64) * (nf64 + mf64) * (2.0 * nf64 - 3.0)))
                    .sqrt();

                b_nm[(n, m)] = (((2.0 * nf64 + 1.0) * (2.0 * nf64 - 1.0))
                    / ((nf64 + mf64) * (nf64 - mf64)))
                    .sqrt();

                vr01[(n, m)] = ((nf64 - mf64) * (nf64 + mf64 + 1.0)).sqrt();
                vr11[(n, m)] = (((2.0 * nf64 + 1.0) * (nf64 + mf64 + 2.0) * (nf64 + mf64 + 1.0))
                    / (2.0 * nf64 + 3.0))
                    .sqrt();

                if m == 0 {
                    vr01[(n, m)] /= 2.0_f64.sqrt();
                    vr11[(n, m)] /= 2.0_f64.sqrt();
                }
            }
        }

        // Repeat for the hyperdual part in case we need to super the partials
        let mut a_nm_h =
            DMatrix::from_element(degree_np2 + 1, degree_np2 + 1, OHyperdual::from(0.0));
        let mut b_nm_h = DMatrix::from_element(degree_np2, degree_np2, OHyperdual::from(0.0));
        let mut c_nm_h = DMatrix::from_element(degree_np2, degree_np2, OHyperdual::from(0.0));
        let mut vr01_h = DMatrix::from_element(degree_np2, degree_np2, OHyperdual::from(0.0));
        let mut vr11_h = DMatrix::from_element(degree_np2, degree_np2, OHyperdual::from(0.0));

        // initialize the diagonal elements (not a function of the input)
        a_nm_h[(0, 0)] = OHyperdual::from(1.0);
        for n in 1..=degree_np2 {
            // Diagonal element
            a_nm_h[(n, n)] = OHyperdual::from(a_nm[(n, n)]);
        }

        // Pre-compute the B_nm, C_nm, vr01 and vr11 storages
        for n in 0..degree_np2 {
            for m in 0..degree_np2 {
                vr01_h[(n, m)] = OHyperdual::from(vr01[(n, m)]);
                vr11_h[(n, m)] = OHyperdual::from(vr11[(n, m)]);
                b_nm_h[(n, m)] = OHyperdual::from(b_nm[(n, m)]);
                c_nm_h[(n, m)] = OHyperdual::from(c_nm[(n, m)]);
            }
        }

        Arc::new(Self {
            cosm,
            compute_frame,
            stor,
            a_nm,
            b_nm,
            c_nm,
            vr01,
            vr11,
            a_nm_h,
            b_nm_h,
            c_nm_h,
            vr01_h,
            vr11_h,
        })
    }
}

impl<S: GravityPotentialStor + Send> fmt::Display for Harmonics<S> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{} gravity field {}x{} (order x degree)",
            self.compute_frame,
            self.stor.max_order_m(),
            self.stor.max_degree_n(),
        )
    }
}

impl<S: GravityPotentialStor + Send> AccelModel for Harmonics<S> {
    fn eom(&self, osc: &Orbit) -> Result<Vector3<f64>, NyxError> {
        // Convert the osculating orbit to the correct frame (needed for multiple harmonic fields)
        let state = self.cosm.frame_chg(osc, self.compute_frame);

        // Using the GMAT notation, with extra character for ease of highlight
        let r_ = state.rmag();
        let s_ = state.x / r_;
        let t_ = state.y / r_;
        let u_ = state.z / r_;
        let max_degree = self.stor.max_degree_n() as usize; // In GMAT, the degree is NN
        let max_order = self.stor.max_order_m() as usize; // In GMAT, the order is MM

        // Create the associated Legendre polynomials. Note that we add three items as per GMAT (this may be useful for the STM)
        let mut a_nm = self.a_nm.clone();

        // Initialize the diagonal elements (not a function of the input)
        a_nm[(1, 0)] = u_ * 3.0f64.sqrt();
        for n in 1..=max_degree + 1 {
            let nf64 = n as f64;
            // Off diagonal
            a_nm[(n + 1, n)] = (2.0 * nf64 + 3.0).sqrt() * u_ * a_nm[(n, n)];
        }

        for m in 0..=max_order + 1 {
            for n in (m + 2)..=max_degree + 1 {
                let hm_idx = (n, m);
                a_nm[hm_idx] = u_ * self.b_nm[hm_idx] * a_nm[(n - 1, m)]
                    - self.c_nm[hm_idx] * a_nm[(n - 2, m)];
            }
        }

        // Generate r_m and i_m
        let mut r_m = Vec::with_capacity(min(max_degree, max_order) + 1);
        let mut i_m = Vec::with_capacity(min(max_degree, max_order) + 1);

        r_m.push(1.0);
        i_m.push(0.0);

        for m in 1..=min(max_degree, max_order) {
            r_m.push(s_ * r_m[m - 1] - t_ * i_m[m - 1]);
            i_m.push(s_ * i_m[m - 1] + t_ * r_m[m - 1]);
        }

        let rho = self.compute_frame.equatorial_radius() / r_;
        let mut rho_np1 = self.compute_frame.gm() / r_ * rho;
        let mut a0 = 0.0;
        let mut a1 = 0.0;
        let mut a2 = 0.0;
        let mut a3 = 0.0;

        for n in 1..max_degree {
            let mut sum0 = 0.0;
            let mut sum1 = 0.0;
            let mut sum2 = 0.0;
            let mut sum3 = 0.0;
            rho_np1 *= rho;

            for m in 0..=min(n, max_order) {
                let (c_val, s_val) = self.stor.cs_nm(n, m);
                let d_ = (c_val * r_m[m] + s_val * i_m[m]) * 2.0.sqrt();
                let e_ = if m == 0 {
                    0.0
                } else {
                    (c_val * r_m[m - 1] + s_val * i_m[m - 1]) * 2.0.sqrt()
                };
                let f_ = if m == 0 {
                    0.0
                } else {
                    (s_val * r_m[m - 1] - c_val * i_m[m - 1]) * 2.0.sqrt()
                };

                sum0 += (m as f64) * a_nm[(n, m)] * e_;
                sum1 += (m as f64) * a_nm[(n, m)] * f_;
                sum2 += self.vr01[(n, m)] * a_nm[(n, m + 1)] * d_;
                sum3 += self.vr11[(n, m)] * a_nm[(n + 1, m + 1)] * d_;
            }
            let rr = rho_np1 / self.compute_frame.equatorial_radius();
            a0 += rr * sum0;
            a1 += rr * sum1;
            a2 += rr * sum2;
            a3 -= rr * sum3;
        }
        let accel = Vector3::new(a0 + a3 * s_, a1 + a3 * t_, a2 + a3 * u_);
        // Rotate this acceleration vector back into the integration frame (no center change needed, it's just a vector)
        // As discussed with Sai, if the Earth was spinning faster, would the acceleration due to the harmonics be any different?
        // No. Therefore, we do not need to account for the transport theorem here.
        let dcm = self
            .cosm
            .try_position_dcm_from_to(&self.compute_frame, &osc.frame, osc.dt)?;
        Ok(dcm * accel)
    }

    fn dual_eom(&self, osc: &Orbit) -> Result<(Vector3<f64>, Matrix3<f64>), NyxError> {
        // Convert the osculating orbit to the correct frame (needed for multiple harmonic fields)
        let state = self.cosm.frame_chg(osc, self.compute_frame);

        let radius: Vector3<OHyperdual<f64, U7>> = hyperspace_from_vector(&state.radius());

        // Using the GMAT notation, with extra character for ease of highlight
        let r_ = norm(&radius);
        let s_ = radius[0] / r_;
        let t_ = radius[1] / r_;
        let u_ = radius[2] / r_;
        let max_degree = self.stor.max_degree_n() as usize; // In GMAT, the order is NN
        let max_order = self.stor.max_order_m() as usize; // In GMAT, the order is MM

        // Create the associated Legendre polynomials. Note that we add three items as per GMAT (this may be useful for the STM)
        let mut a_nm = self.a_nm_h.clone();

        // Initialize the diagonal elements (not a function of the input)
        a_nm[(1, 0)] = u_ * 3.0f64.sqrt();
        for n in 1..=max_degree + 1 {
            let nf64 = n as f64;
            // Off diagonal
            a_nm[(n + 1, n)] = OHyperdual::from((2.0 * nf64 + 3.0).sqrt()) * u_ * a_nm[(n, n)];
        }

        for m in 0..=max_order + 1 {
            for n in (m + 2)..=max_degree + 1 {
                let hm_idx = (n, m);
                a_nm[hm_idx] = u_ * self.b_nm_h[hm_idx] * a_nm[(n - 1, m)]
                    - self.c_nm_h[hm_idx] * a_nm[(n - 2, m)];
            }
        }

        // Generate r_m and i_m
        let mut r_m = Vec::with_capacity(min(max_degree, max_order) + 1);
        let mut i_m = Vec::with_capacity(min(max_degree, max_order) + 1);

        r_m.push(OHyperdual::<f64, U7>::from(1.0));
        i_m.push(OHyperdual::<f64, U7>::from(0.0));

        for m in 1..=min(max_degree, max_order) {
            r_m.push(s_ * r_m[m - 1] - t_ * i_m[m - 1]);
            i_m.push(s_ * i_m[m - 1] + t_ * r_m[m - 1]);
        }

        let eq_radius = OHyperdual::<f64, U7>::from(self.compute_frame.equatorial_radius());
        let rho = eq_radius / r_;
        let mut rho_np1 = OHyperdual::<f64, U7>::from(self.compute_frame.gm()) / r_ * rho;

        let mut a0 = OHyperdual::<f64, U7>::from(0.0);
        let mut a1 = OHyperdual::<f64, U7>::from(0.0);
        let mut a2 = OHyperdual::<f64, U7>::from(0.0);
        let mut a3 = OHyperdual::<f64, U7>::from(0.0);
        let sqrt2 = OHyperdual::<f64, U7>::from(2.0.sqrt());

        for n in 1..max_degree {
            let mut sum0 = OHyperdual::from(0.0);
            let mut sum1 = OHyperdual::from(0.0);
            let mut sum2 = OHyperdual::from(0.0);
            let mut sum3 = OHyperdual::from(0.0);
            rho_np1 *= rho;

            for m in 0..=min(n, max_order) {
                let (c_valf64, s_valf64) = self.stor.cs_nm(n, m);
                let c_val = OHyperdual::<f64, U7>::from(c_valf64);
                let s_val = OHyperdual::<f64, U7>::from(s_valf64);

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

                sum0 += OHyperdual::from(m as f64) * a_nm[(n, m)] * e_;
                sum1 += OHyperdual::from(m as f64) * a_nm[(n, m)] * f_;
                sum2 += self.vr01_h[(n, m)] * a_nm[(n, m + 1)] * d_;
                sum3 += self.vr11_h[(n, m)] * a_nm[(n + 1, m + 1)] * d_;
            }
            let rr = rho_np1 / eq_radius;
            a0 += rr * sum0;
            a1 += rr * sum1;
            a2 += rr * sum2;
            a3 -= rr * sum3;
        }

        let dcm = self
            .cosm
            .try_position_dcm_from_to(&self.compute_frame, &osc.frame, osc.dt)?;

        // Convert DCM to OHyperdual DCMs
        let mut dcm_d = Matrix3::<OHyperdual<f64, U7>>::zeros();
        for i in 0..3 {
            for j in 0..3 {
                dcm_d[(i, j)] = OHyperdual::from_fn(|k| {
                    if k == 0 {
                        dcm[(i, j)]
                    } else if i + 1 == k {
                        1.0
                    } else {
                        0.0
                    }
                })
            }
        }

        let accel = dcm_d * Vector3::new(a0 + a3 * s_, a1 + a3 * t_, a2 + a3 * u_);
        // Extract data
        let mut dx = Vector3::zeros();
        let mut grad = Matrix3::zeros();
        for i in 0..3 {
            dx[i] += accel[i].real();
            // NOTE: Although the hyperdual state is of size 7, we're only setting the values up to 3 (Matrix3)
            for j in 1..4 {
                grad[(i, j - 1)] += accel[i][j];
            }
        }
        Ok((dx, grad))
    }
}
