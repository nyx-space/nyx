use super::Dynamics;
use super::na::{U3, U6, Vector6, VectorN};
use celestia::CelestialBody;
use io::gravity::GravityPotentialStor;
use std::cmp::min;
use utils::{factorial, kronecker};

/// `TwoBody` exposes the equations of motion for a simple two body propagation.
#[derive(Clone, Copy)]
pub struct Harmonics<S>
where
    S: GravityPotentialStor,
{
    neg_mu: f64,
    body_radius: f64,
    stor: S,
}

impl<S> Harmonics<S>
where
    S: GravityPotentialStor,
{
    /// Create a new Harmonics dynamical model from the provided gravity potential storage instance.
    pub fn from_stor<B: CelestialBody>(stor: S) -> Harmonics<S> {
        Harmonics {
            neg_mu: -B::gm(),
            body_radius: B::eq_radius(),
            stor: stor,
        }
    }
}

impl<S: GravityPotentialStor> Dynamics for Harmonics<S> {
    type StateSize = U6;

    /// NOTE: No state is associated with Harmonics, always return zero time
    fn time(&self) -> f64 {
        0.0
    }

    /// NOTE: No state is associated with Harmonics, always return zero
    fn state(&self) -> VectorN<f64, Self::StateSize> {
        Vector6::zeros()
    }

    /// NOTE: Nothing happens in this `set_state` since there is no state of spherical harmonics.
    fn set_state(&mut self, _new_t: f64, _new_state: &VectorN<f64, Self::StateSize>) {}

    /// This provides a **DELTA** of the state, which must be added to the result of the TwoBody propagator being used.
    /// However, the provided `state` must be the position and velocity.
    fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        // NOTE: All this code is a conversion from GMAT's CalculateField1
        let radius = state.fixed_rows::<U3>(0).into_owned();
        // Using the GMAT notation, with extra character for ease of highlight
        let r_ = radius.norm();
        let s_ = radius[(0, 0)] / radius.norm();
        let t_ = radius[(1, 0)] / radius.norm();
        let u_ = radius[(2, 0)] / radius.norm();
        let max_degree = self.stor.max_degree() as usize; // In GMAT, the order is NN
        let max_order = self.stor.max_order() as usize; // In GMAT, the order is MM

        let matrix_size = max_degree + 2;

        // Create the associated Legendre polynomials. Note that we add three items as per GMAT (this may be useful for the STM)
        let mut a_matrix: Vec<Vec<f64>> = (0..matrix_size).map(|_| Vec::with_capacity(matrix_size)).collect();
        let mut vr01: Vec<Vec<f64>> = (0..matrix_size).map(|_| Vec::with_capacity(matrix_size)).collect();
        let mut vr11: Vec<Vec<f64>> = (0..matrix_size).map(|_| Vec::with_capacity(matrix_size)).collect();

        // NOTE: In Jones 2009, the re are referred to as p_m and im as q_m.
        let mut re = Vec::with_capacity(matrix_size);
        let mut im = Vec::with_capacity(matrix_size);

        // Now that we have requested that capacity, let's set everything to zero so we can populate as a_matrix[n][m] (instead of pushing values).
        for n in 0..=max_degree + 1 {
            for _m in 0..=(n + 1) {
                a_matrix[n].push(0.0);
                vr01[n].push(0.0);
                vr11[n].push(0.0);
            }
        }

        for _m in 0..=max_order + 1 {
            re.push(0.0);
            im.push(0.0);
        }

        let b_mn = |n: f64, m: f64| (((2.0 * n + 1.0) * (2.0 * n - 1.0)) / ((n + m) * (n - m))).sqrt();

        a_matrix[0][0] = 1.0;
        a_matrix[1][0] = u_ * 3.0f64.sqrt();
        a_matrix[1][1] = 3.0f64.sqrt();

        for nu16 in 0..=max_degree {
            let n = nu16 as f64;
            for mu16 in 0..=nu16 {
                let m = mu16 as f64;
                vr01[nu16][mu16] = ((n - m) * (n + m + 1.0)).sqrt();
                vr11[nu16][mu16] = (((2.0 * n + 1.0) * (n + m + 2.0) * (n + m + 1.0)) / (2.0 * n + 3.0)).sqrt();
                if mu16 == 0 {
                    vr01[nu16][mu16] /= 2.0f64.sqrt();
                    vr11[nu16][mu16] /= 2.0f64.sqrt();
                }
                if nu16 >= 1 && mu16 == nu16 {
                    a_matrix[nu16][nu16] = ((2.0 * n + 1.0) / (2.0 * n)).sqrt() * a_matrix[nu16 - 1][nu16 - 1];
                } else if nu16 >= 2 {
                    a_matrix[nu16][mu16] =
                        u_ * b_mn(n, m) * a_matrix[nu16 - 1][mu16] - (b_mn(n, m) / b_mn(n - 1.0, m)) * a_matrix[nu16 - 2][mu16];
                }
            }
        }

        // apply column-fill recursion formula (Table 2, Row I, Ref.[1])
        for m in 0..=max_order + 1 {
            re[m] = if m == 0 { 1.0 } else { s_ * re[m - 1] - t_ * im[m - 1] };
            im[m] = if m == 0 { 0.0 } else { s_ * im[m - 1] + t_ * re[m - 1] };
        }

        // Now do summation ------------------------------------------------
        // initialize recursion
        let rho = self.body_radius / r_;
        let mut rho_np1 = -self.neg_mu / r_ * rho;
        // NOTE: There currently is no rho_np2 because that's only used when computing the STM
        let common_fact = -self.neg_mu / (self.body_radius * r_); // TODO: I'm taking the negative of the negative, remove this
        let mut a1 = 0.0;
        let mut a2 = 0.0;
        let mut a3 = 0.0;
        let mut a4 = 0.0;
        let sqrt2 = 2.0f64.sqrt();

        for n in 1..=max_degree {
            rho_np1 *= rho;
            let mut sum1 = 0.0;
            let mut sum2 = 0.0;
            let mut sum3 = 0.0;
            let mut sum4 = 0.0;
            for m in 0..=min(n, max_order) {
                let (c_val, s_val) = self.stor.cs_nm(n as u16, m as u16);
                // Pines Equation 27 (Part of)
                let d_ = (c_val * re[m] + s_val * im[m]) * sqrt2;
                let e_ = if m == 0 {
                    0.0
                } else {
                    (c_val * re[m - 1] + s_val * im[m - 1]) * sqrt2
                };
                let f_ = if m == 0 {
                    0.0
                } else {
                    (s_val * re[m - 1] - c_val * im[m - 1]) * sqrt2
                };

                sum1 += (m as f64) * a_matrix[n][m] * e_;
                sum2 += (m as f64) * a_matrix[n][m] * f_;
                sum3 += vr01[n][m] * a_matrix[n][m + 1] * d_;
                sum4 += vr11[n][m] * a_matrix[n + 1][m + 1] * d_;
            }
            let rr = rho_np1 / self.body_radius;
            a1 += rr * sum1;
            a2 += rr * sum2;
            a3 += rr * sum3;
            a4 -= rr * sum4;
        }

        Vector6::new(0.0, 0.0, 0.0, a1 + a4 * s_, a2 + a4 * t_, a3 + a4 * u_)
    }
}
