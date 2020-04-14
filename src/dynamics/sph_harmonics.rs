use crate::celestia::{Cosm, Frame, State};
use crate::dimensions::{DMatrix, Vector3};
use crate::dynamics::AccelModel;
use crate::io::gravity::GravityPotentialStor;
use std::cmp::min;

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
}

impl<'a, S> Harmonics<'a, S>
where
    S: GravityPotentialStor,
{
    /// Create a new Harmonics dynamical model from the provided gravity potential storage instance.
    pub fn from_stor(compute_frame: Frame, stor: S, cosm: &'a Cosm) -> Self {
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

        // initialize the diagonal elements (not a function of the input)
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

        Harmonics {
            cosm,
            compute_frame,
            stor,
            a_nm,
            b_nm,
            c_nm,
            vr01,
            vr11,
        }
    }
}

impl<'a, S: GravityPotentialStor> AccelModel for Harmonics<'a, S> {
    fn eom(&self, osc: &State) -> Vector3<f64> {
        // Get the DCM to convert from the integration state to the computation frame of the harmonics
        let dcm = self
            .cosm
            .try_frame_chg_dcm_from_to(&osc.frame, &self.compute_frame, osc.dt)
            .unwrap();
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

        let rho = self.compute_frame.equatorial_radius() / r_;
        let mut a1 = 0.0;
        let mut a2 = 0.0;
        let mut a3 = 0.0;
        let mut a4 = 0.0;

        for n in 1..max_degree {
            let mut sum1 = 0.0;
            let mut sum2 = 0.0;
            let mut sum3 = 0.0;
            let mut sum4 = 0.0;

            for m in 0..=min(n, max_order) {
                let (c_val, s_val) = self.stor.cs_nm(n, m);
                let d_ = c_val * r_m(m as u16, s_, t_) + s_val * i_m(m as u16, s_, t_);
                let e_ = if m == 0 {
                    0.0
                } else {
                    c_val * r_m(m as u16 - 1, s_, t_) + s_val * i_m(m as u16 - 1, s_, t_)
                };
                let f_ = if m == 0 {
                    0.0
                } else {
                    s_val * r_m(m as u16 - 1, s_, t_) - c_val * i_m(m as u16 - 1, s_, t_)
                };

                sum1 += (m as f64) * a_nm[(n, m)] * e_;
                sum2 += (m as f64) * a_nm[(n, m)] * f_;
                sum3 += self.vr01[(n, m)] * a_nm[(n, m + 1)] * d_;
                sum4 += self.vr11[(n, m)] * a_nm[(n + 1, m + 1)] * d_;
            }
            let rr = rho.powi(n as i32 + 1);
            a1 += rr * sum1;
            a2 += rr * sum2;
            a3 += rr * sum3;
            a4 += rr * sum4;
        }
        let mu_fact = self.compute_frame.gm() / (self.compute_frame.equatorial_radius() * r_);
        a1 *= mu_fact;
        a2 *= mu_fact;
        a3 *= mu_fact;
        a4 *= -mu_fact;
        let accel = Vector3::new(a1 + a4 * s_, a2 + a4 * t_, a3 + a4 * u_);
        // Convert back to integration frame
        dcm.transpose() * accel
    }
}

fn r_m(m: u16, s: f64, t: f64) -> f64 {
    if m == 0 {
        1.0
    } else {
        s * r_m(m - 1, s, t) - t * i_m(m - 1, s, t)
    }
}

fn i_m(m: u16, s: f64, t: f64) -> f64 {
    if m == 0 {
        0.0
    } else {
        s * i_m(m - 1, s, t) + t * r_m(m - 1, s, t)
    }
}
