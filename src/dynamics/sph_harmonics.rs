use super::hyperdual::linalg::norm;
use super::hyperdual::{Float, Hyperdual};
use crate::celestia::{Cosm, Frame, Orbit};
use crate::dimensions::{DMatrix, Matrix3, Vector3, U7};
use crate::dynamics::AccelModel;
use crate::errors::NyxError;
use crate::io::gravity::GravityPotentialStor;
use crate::TimeTagged;
use std::cmp::min;
use std::sync::Arc;

#[derive(Clone)]
pub struct Harmonics<'a, S>
where
    S: GravityPotentialStor,
{
    cosm: &'a Cosm,
    compute_frame: Frame,
    stor: S,
    a_nm: DMatrix<f64>,
    b_nm: DMatrix<f64>,
    c_nm: DMatrix<f64>,
    vr01: DMatrix<f64>,
    vr11: DMatrix<f64>,
    a_nm_h: DMatrix<Hyperdual<f64, U7>>,
    b_nm_h: DMatrix<Hyperdual<f64, U7>>,
    c_nm_h: DMatrix<Hyperdual<f64, U7>>,
    vr01_h: DMatrix<Hyperdual<f64, U7>>,
    vr11_h: DMatrix<Hyperdual<f64, U7>>,
}

impl<'a, S> Harmonics<'a, S>
where
    S: GravityPotentialStor,
{
    /// Create a new Harmonics dynamical model from the provided gravity potential storage instance.
    pub fn from_stor(compute_frame: Frame, stor: S, cosm: &'a Cosm) -> Arc<Self> {
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
        a_nm[(1, 1)] = 3.0f64.sqrt(); // This is the same formulation as below for n=1
        for n in 2..=degree_np2 {
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
            DMatrix::from_element(degree_np2 + 1, degree_np2 + 1, Hyperdual::from_real(0.0));
        let mut b_nm_h = DMatrix::from_element(degree_np2, degree_np2, Hyperdual::from_real(0.0));
        let mut c_nm_h = DMatrix::from_element(degree_np2, degree_np2, Hyperdual::from_real(0.0));
        let mut vr01_h = DMatrix::from_element(degree_np2, degree_np2, Hyperdual::from_real(0.0));
        let mut vr11_h = DMatrix::from_element(degree_np2, degree_np2, Hyperdual::from_real(0.0));

        // initialize the diagonal elements (not a function of the input)
        a_nm_h[(0, 0)] = Hyperdual::from_real(1.0);
        a_nm_h[(1, 1)] = Hyperdual::from_real(3.0f64.sqrt()); // This is the same formulation as below for n=1
        for n in 2..=degree_np2 {
            let nf64 = n as f64;
            // Diagonal element
            a_nm_h[(n, n)] =
                Hyperdual::from_real((1.0 + 1.0 / (2.0 * nf64)).sqrt()) * a_nm_h[(n - 1, n - 1)];
        }

        // Pre-compute the B_nm, C_nm, vr01 and vr11 storages
        for n in 0..degree_np2 {
            for m in 0..degree_np2 {
                let nf64 = n as f64;
                let mf64 = m as f64;
                // Compute c_nm, which is B_nm/B_(n-1,m) in Jones' dissertation
                c_nm_h[(n, m)] = Hyperdual::from_real(
                    (((2.0 * nf64 + 1.0) * (nf64 + mf64 - 1.0) * (nf64 - mf64 - 1.0))
                        / ((nf64 - mf64) * (nf64 + mf64) * (2.0 * nf64 - 3.0)))
                        .sqrt(),
                );

                b_nm_h[(n, m)] = Hyperdual::from_real(
                    (((2.0 * nf64 + 1.0) * (2.0 * nf64 - 1.0)) / ((nf64 + mf64) * (nf64 - mf64)))
                        .sqrt(),
                );

                vr01_h[(n, m)] = Hyperdual::from_real(((nf64 - mf64) * (nf64 + mf64 + 1.0)).sqrt());
                vr11_h[(n, m)] = Hyperdual::from_real(
                    (((2.0 * nf64 + 1.0) * (nf64 + mf64 + 2.0) * (nf64 + mf64 + 1.0))
                        / (2.0 * nf64 + 3.0))
                        .sqrt(),
                );

                if m == 0 {
                    vr01_h[(n, m)] /= 2.0_f64.sqrt();
                    vr11_h[(n, m)] /= 2.0_f64.sqrt();
                }
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

impl<'a, S: GravityPotentialStor + Send> AccelModel for Harmonics<'a, S> {
    fn eom(&self, osc: &Orbit) -> Result<Vector3<f64>, NyxError> {
        // Get the DCM to convert from the integration state to the computation frame of the harmonics
        let dcm = self
            .cosm
            .try_frame_chg_dcm_from_to(&osc.frame, &self.compute_frame, osc.dt)?;
        // Convert to the computation frame
        let mut state = *osc;
        state.apply_dcm(dcm);

        // Using the GMAT notation, with extra character for ease of highlight
        let r_ = state.rmag();
        let s_ = state.x / r_;
        let t_ = state.y / r_;
        let u_ = state.z / r_;
        let max_degree = self.stor.max_degree_n() as usize; // In GMAT, the order is NN
        let max_order = self.stor.max_order_m() as usize; // In GMAT, the order is MM

        // Create the associated Legendre polynomials. Note that we add three items as per GMAT (this may be useful for the STM)
        let mut a_nm = self.a_nm.clone();

        // initialize the diagonal elements (not a function of the input)
        a_nm[(1, 0)] = u_ * 3.0f64.sqrt();
        for n in 1..=max_degree + 1 {
            let nf64 = n as f64;
            // Off diagonal
            a_nm[(n + 1, n)] = (2.0 * nf64 + 3.0).sqrt() * u_ * a_nm[(n, n)];
        }

        for m in 0..=max_order + 1 {
            for n in (m + 2)..=max_degree + 1 {
                let hm_idx = (n, m);
                a_nm[(n, m)] = u_ * self.b_nm[hm_idx] * a_nm[(n - 1, m)]
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
        let mut a0 = 0.0;
        let mut a1 = 0.0;
        let mut a2 = 0.0;
        let mut a3 = 0.0;

        for n in 1..max_degree {
            let mut sum0 = 0.0;
            let mut sum1 = 0.0;
            let mut sum2 = 0.0;
            let mut sum3 = 0.0;

            for m in 0..=min(n, max_order) {
                let (c_val, s_val) = self.stor.cs_nm(n, m);
                let d_ = c_val * r_m[m] + s_val * i_m[m];
                let e_ = if m == 0 {
                    0.0
                } else {
                    c_val * r_m[m - 1] + s_val * i_m[m - 1]
                };
                let f_ = if m == 0 {
                    0.0
                } else {
                    s_val * r_m[m - 1] - c_val * i_m[m - 1]
                };

                sum0 += (m as f64) * a_nm[(n, m)] * e_;
                sum1 += (m as f64) * a_nm[(n, m)] * f_;
                sum2 += self.vr01[(n, m)] * a_nm[(n, m + 1)] * d_;
                sum3 += self.vr11[(n, m)] * a_nm[(n + 1, m + 1)] * d_;
            }
            let rr = rho.powi(n as i32 + 1);
            a0 += rr * sum0;
            a1 += rr * sum1;
            a2 += rr * sum2;
            a3 += rr * sum3;
        }
        let mu_fact = self.compute_frame.gm() / (self.compute_frame.equatorial_radius() * r_);
        a0 *= mu_fact;
        a1 *= mu_fact;
        a2 *= mu_fact;
        a3 *= -mu_fact;
        let accel = Vector3::new(a0 + a3 * s_, a1 + a3 * t_, a2 + a3 * u_);
        // Convert back to integration frame
        Ok(dcm.transpose() * accel)
    }

    fn dual_eom(
        &self,
        radius: &Vector3<Hyperdual<f64, U7>>,
        ctx: &Orbit,
    ) -> Result<(Vector3<f64>, Matrix3<f64>), NyxError> {
        // Get the DCM to convert from the integration state to the computation frame of the harmonics
        let dcm =
            self.cosm
                .try_frame_chg_dcm_from_to(&ctx.frame, &self.compute_frame, ctx.epoch())?;
        // Convert DCM to Hyperdual DCMs
        let mut dcm_d = Matrix3::<Hyperdual<f64, U7>>::zeros();
        for i in 0..3 {
            for j in 0..3 {
                dcm_d[(i, j)] = Hyperdual::from_fn(|k| {
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

        // Convert to the computation frame
        let radius = dcm_d * radius;

        // Using the GMAT notation, with extra character for ease of highlight
        let r_ = norm(&radius);
        let s_ = radius[0] / r_;
        let t_ = radius[1] / r_;
        let u_ = radius[2] / r_;
        let max_degree = self.stor.max_degree_n() as usize; // In GMAT, the order is NN
        let max_order = self.stor.max_order_m() as usize; // In GMAT, the order is MM

        // Create the associated Legendre polynomials. Note that we add three items as per GMAT (this may be useful for the STM)
        let mut a_nm = self.a_nm_h.clone();

        // initialize the diagonal elements (not a function of the input)
        a_nm[(1, 0)] = u_ * Hyperdual::<f64, U7>::from_real(3.0f64.sqrt());
        for n in 1..=max_degree + 1 {
            let nf64 = n as f64;
            // Off diagonal
            a_nm[(n + 1, n)] =
                Hyperdual::<f64, U7>::from_real((2.0 * nf64 + 3.0).sqrt()) * u_ * a_nm[(n, n)];
        }

        for m in 0..=max_order + 1 {
            for n in (m + 2)..=max_degree + 1 {
                let hm_idx = (n, m);
                a_nm[(n, m)] = u_ * self.b_nm_h[hm_idx] * a_nm[(n - 1, m)]
                    - self.c_nm_h[hm_idx] * a_nm[(n - 2, m)];
            }
        }

        // Generate r_m and i_m
        let mut r_m = Vec::with_capacity(min(max_degree, max_order) + 1);
        let mut i_m = Vec::with_capacity(min(max_degree, max_order) + 1);

        r_m.push(Hyperdual::<f64, U7>::from_real(1.0));
        i_m.push(Hyperdual::<f64, U7>::from_real(0.0));

        for m in 1..=min(max_degree, max_order) {
            r_m.push(s_ * r_m[m - 1] - t_ * i_m[m - 1]);
            i_m.push(s_ * i_m[m - 1] + t_ * r_m[m - 1]);
        }

        let eq_radius = Hyperdual::<f64, U7>::from_real(self.compute_frame.equatorial_radius());
        let rho = eq_radius / r_;
        let mut a1 = Hyperdual::<f64, U7>::from_real(0.0);
        let mut a2 = Hyperdual::<f64, U7>::from_real(0.0);
        let mut a3 = Hyperdual::<f64, U7>::from_real(0.0);
        let mut a4 = Hyperdual::<f64, U7>::from_real(0.0);

        for n in 1..max_degree {
            let mut sum1 = Hyperdual::<f64, U7>::from_real(0.0);
            let mut sum2 = Hyperdual::<f64, U7>::from_real(0.0);
            let mut sum3 = Hyperdual::<f64, U7>::from_real(0.0);
            let mut sum4 = Hyperdual::<f64, U7>::from_real(0.0);

            for m in 0..=min(n, max_order) {
                let (c_valf64, s_valf64) = self.stor.cs_nm(n, m);
                let c_val = Hyperdual::<f64, U7>::from_real(c_valf64);
                let s_val = Hyperdual::<f64, U7>::from_real(s_valf64);
                let d_ = c_val * r_m[m] + s_val * i_m[m];
                let e_ = if m == 0 {
                    Hyperdual::<f64, U7>::from_real(0.0)
                } else {
                    c_val * r_m[m - 1] + s_val * i_m[m - 1]
                };
                let f_ = if m == 0 {
                    Hyperdual::<f64, U7>::from_real(0.0)
                } else {
                    s_val * r_m[m - 1] - c_val * i_m[m - 1]
                };

                sum1 += Hyperdual::<f64, U7>::from_real(m as f64) * a_nm[(n, m)] * e_;
                sum2 += Hyperdual::<f64, U7>::from_real(m as f64) * a_nm[(n, m)] * f_;
                sum3 += self.vr01_h[(n, m)] * a_nm[(n, m + 1)] * d_;
                sum4 += self.vr11_h[(n, m)] * a_nm[(n + 1, m + 1)] * d_;
            }
            let rr = rho.powi(n as i32 + 1);
            a1 += rr * sum1;
            a2 += rr * sum2;
            a3 += rr * sum3;
            a4 += rr * sum4;
        }
        let mu_fact = Hyperdual::<f64, U7>::from_real(self.compute_frame.gm()) / (eq_radius * r_);
        a1 *= mu_fact;
        a2 *= mu_fact;
        a3 *= mu_fact;
        a4 *= -mu_fact;
        // Convert back to integration frame
        let accel = dcm_d.transpose() * Vector3::new(a1 + a4 * s_, a2 + a4 * t_, a3 + a4 * u_);
        // Extract data
        let mut fx = Vector3::zeros();
        let mut grad = Matrix3::zeros();
        for i in 0..3 {
            fx[i] += accel[i][0];
            // NOTE: Although the hyperdual state is of size 7, we're only setting the values up to 3 (Matrix3)
            for j in 1..4 {
                grad[(i, j - 1)] += accel[i][j];
            }
        }
        Ok((fx, grad))
    }
}
