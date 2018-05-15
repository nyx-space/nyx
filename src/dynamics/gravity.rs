// use super::hifitime::SECONDS_PER_DAY;
use super::Dynamics;
use super::na::{U1, U3, U6, Vector6, VectorN};
use celestia::CelestialBody;
use io::gravity::GravityPotentialStor;

/// `TwoBody` exposes the equations of motion for a simple two body propagation.
#[derive(Clone, Copy)]
pub struct Harmonics<S>
where
    S: GravityPotentialStor,
{
    mu: f64,
    radius: f64,
    stor: S,
}

impl<S> Harmonics<S>
where
    S: GravityPotentialStor,
{
    /// Create a new Harmonics dynamical model from the provided gravity potential storage instance.
    pub fn from_stor<B: CelestialBody>(stor: S) -> Harmonics<S> {
        Harmonics {
            mu: B::gm(),
            radius: B::eq_radius(),
            stor: stor,
        }
    }
}

impl<S: GravityPotentialStor> Dynamics for Harmonics<S> {
    type StateSize = U6;

    fn time(&self) -> f64 {
        0.0
    }

    fn state(&self) -> VectorN<f64, Self::StateSize> {
        // No state is associated with Harmonics, always return zero
        Vector6::zeros()
    }

    /// WARNING: The *new_t* parameter is considered a TIME INCREASE from the initial julian days.
    /// This is likely a bad assumption (breaks back prop unless new_t is negative) -- I will need to clarify that in the docs.
    fn set_state(&mut self, new_t: f64, _new_state: &VectorN<f64, Self::StateSize>) {
        // self.time += new_t / SECONDS_PER_DAY;
        unimplemented!();
    }

    /// WARNING: This provides a DELTA of the state, which must be added to the result of the TwoBody propagator being used.
    /// However, the provided `state` must be the position and velocity.
    fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        // NOTE: All this code is a conversion from GMAT's CalculateField1
        let radius = state.fixed_slice::<U3, U1>(0, 0);
        // Using the GMAT notation, with extra character for ease of highlight
        let r_ = radius.norm();
        let s_ = radius[(0, 0)] / radius.norm();
        let t_ = radius[(1, 0)] / radius.norm();
        let u_ = radius[(2, 0)] / radius.norm();
        let max_order = self.stor.max_order() as usize;
        let max_degree = self.stor.max_degree() as usize;
        // Create the associated Legendre polynomials. Note that we add three items as per GMAT (this may be useful for the STM)
        let mut a_matrix: Vec<Vec<f64>> = (0..max_order + 3)
            .map(|_| Vec::with_capacity(max_degree + 3))
            .collect();
        // Now that we have requested that capacity, let's set everything to zero so we can populate with [n][m].
        for n in 0..max_order + 3 {
            for _m in 0..max_degree + 3 {
                a_matrix[n].push(0.0);
            }
        }

        // initialize the diagonal elements (not a function of the input)
        a_matrix[0][0] = 1.0; // Temp value for this first initialization
        for n in 1..max_degree {
            let nf64 = n as f64;
            a_matrix[n][n] = ((2.0 * nf64 + 1.0).sqrt() / (2.0 * nf64)) * a_matrix[n - 1][n - 1];
        }

        // generate the off-diagonal elements
        a_matrix[1][0] = u_ * 3.0f64.sqrt();
        for n in 1..max_degree + 1 {
            a_matrix[n + 1][n] = u_ * ((2 * n + 3) as f64).sqrt() * a_matrix[n][n];
        }

        let mut re = Vec::with_capacity(max_degree);
        let mut im = Vec::with_capacity(max_degree);

        for m in 0..max_degree + 1 {
            for n in (m + 2)..max_order + 1 {
                // XX: This looks like a hack to avoid J<2 because they don't exist (but I do not store them)
                // normalization factors
                let n1 = ((((2 * n + 1) * (2 * n - 1)) / ((n - m) * (n + m))) as f64).sqrt();
                let n2 = ((((2 * n + 1) * (n - m - 1) * (n + m - 1))
                    / ((2 * n - 3) * (n + m) * (n - m))) as f64)
                    .sqrt();
                a_matrix[n][m] = u_ * n1 * a_matrix[n - 1][m] - n2 * a_matrix[n - 2][m];
            }
            // real part of (s + i*t)^m -- note, we can't borrow mutability, so we can't use clever Rust notation here
            if m == 0 {
                re.push(1.0);
            } else {
                let re_m = re[(m - 1)];
                re.push(s_ * re_m - t_ * im[(m - 1)]);
            }

            if m == 0 {
                im.push(0.0);
            } else {
                let im_m = im[(m - 1)];
                im.push(s_ * im_m + t_ * re[(m - 1)]);
            }; // imaginary part of (s + i*t)^m
        }

        // Now do summation ------------------------------------------------
        // initialize recursion
        let rho = self.radius / r_;
        // NOTE: There currently is no rho_np2 because that's only used when computing the STM
        let mut rho_np1 = -self.mu / r_ * rho; // rho(0) ,Ref[3], Eq 26 , factor = mu for gravity
        let mut a1 = 0.0;
        let mut a2 = 0.0;
        let mut a3 = 0.0;
        let mut a4 = 0.0;
        let sqrt2 = 2.0f64.sqrt();
        for n in 1..max_order {
            rho_np1 *= rho;
            let mut sum1 = 0.0;
            let mut sum2 = 0.0;
            let mut sum3 = 0.0;
            let mut sum4 = 0.0;

            let mut m = 0;
            while m <= n && m <= max_degree {
                let (c_val, s_val) = self.stor.cs_nm(n as u16, m as u16);
                // let c_val = self.Cnm(self.julian_days, n, m);
                // let s_val = self.Snm(self.julian_days, n, m);
                // Pines Equation 27 (Part of)
                let d_ = (c_val * re[m] + s_val * im[m]) * sqrt2;
                let e_ = if m == 0 {
                    0.0
                } else {
                    (c_val * re[(m - 1)] + s_val * im[(m - 1)]) * sqrt2
                };
                let f_ = if m == 0 {
                    0.0
                } else {
                    (s_val * re[(m - 1)] - c_val * im[(m - 1)]) * sqrt2
                };
                // Correct for normalization
                let mut vr01 = (((n - m) * (n + m + 1)) as f64).sqrt();
                let mut vr11 = (((2 * n + 1) * (n + m + 2) * (n + m + 1)) as f64
                    / ((2 * n + 3) as f64))
                    .sqrt();
                if m == 0 {
                    vr01 /= (2.0f64).sqrt();
                    vr11 /= (2.0f64).sqrt();
                }

                let avv01 = vr01 * a_matrix[n][m + 1];
                let avv11 = vr11 * a_matrix[n + 1][m + 1];
                // Pines Equation 30 and 30b (Part of)
                sum1 += (m as f64) * a_matrix[n][m] * e_;
                sum2 += (m as f64) * a_matrix[n][m] * f_;
                sum3 += avv01 * d_;
                sum4 += avv11 * d_;
                m += 1;
            }
            // Pines Equation 30 and 30b (Part of)
            let rr = rho_np1 / self.radius;
            a1 += rr * sum1;
            a2 += rr * sum2;
            a3 += rr * sum3;
            a4 -= rr * sum4;
        }

        // Pines Equation 31
        Vector6::new(0.0, 0.0, 0.0, a1 + a4 * s_, a2 + a4 * t_, a3 + a4 * u_)
    }
}
