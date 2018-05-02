// use super::hifitime::SECONDS_PER_DAY;
use super::Dynamics;
use super::na::{U1, U3, U6, Vector6, VectorN};
use celestia::{CelestialBody, EARTH};
use std::collections::HashMap;

/// `TwoBody` exposes the equations of motion for a simple two body propagation.
pub struct Harmonics {
    julian_days: f64, // XXX: This is not actually used in GMAT
    mu: f64,
    radius: f64,
    order: u16,
    degree: u16,
    // data is (degree, order) -> (C_nm, S_nm)
    data: HashMap<(u16, u16), (f64, f64)>,
}

impl Harmonics {
    /// Initialize `Harmonics` as J<sub>2</sub> only using the JGM3 model (available in GMAT)
    ///
    /// Use the embedded Earth parameter. If others are needed, low from `from_gunzip`.
    /// *WARNING:* This is an EARTH gravity model, and _should not_ be used around any other body.
    pub fn J2_JGM3(start_julian_days: f64) -> Harmonics {
        let mut data = HashMap::new();
        data.insert((2, 0), (-4.84165374886470e-04, 0.0));
        Harmonics {
            julian_days: start_julian_days,
            mu: EARTH::gm(),
            radius: EARTH::eq_radius(),
            order: 2,
            degree: 0,
            data: data,
        }
    }

    /// Initialize `Harmonics` as J<sub>2</sub> only using the EGM2008 model (from the GRACE mission, best model as of 2018)
    ///
    /// *WARNING:* This is an EARTH gravity model, and _should not_ be used around any other body.
    pub fn J2_EGM2008(start_julian_days: f64) -> Harmonics {
        let mut data = HashMap::new();
        data.insert((2, 0), (-0.484165143790815e-03, 0.0));
        Harmonics {
            julian_days: start_julian_days,
            mu: EARTH::gm(),
            radius: EARTH::eq_radius(),
            order: 2,
            degree: 0,
            data: data,
        }
    }

    /// Initialize `Harmonics` from the file path (must be a gunzipped file)
    ///
    /// Gravity models provided by `nyx`: TODO: add github links and examples
    /// + EMG2008 to 2190 for Earth (tide free)
    /// + Moon to 1500 (from SHADR file)
    /// + Mars to 120 (from SHADR file)
    /// + Venus to 150 (from SHADR file)
    pub fn from_gunzip(
        start_julian_days: f64,
        degree: u16,
        order: u16,
        filepath: String,
    ) -> Harmonics {
        unimplemented!();
    }
}

impl Dynamics for Harmonics {
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
        // Create the associated Legendre polynomials using a DMatrix (for runtime flexibility).
        let mut A_ = [[0f64; self.order + 1]; self.order];
        let mut Re = Vec::with_capacity(self.degree as usize);
        let mut Im = Vec::with_capacity(self.degree as usize);

        // Populate the off-diagonal elements
        A_[(1, 0)] = u_ * (3.0f64).sqrt();
        for i in 1..self.order {
            A_[(i + 1, i)] = u_ * (2.0 * f64::from(i) + 3.0).sqrt() * A_[(i, i)];
        }

        for m in 0..self.degree {
            for n in (m + 2)..self.order {
                // normalization factors
                let N1 = f64::from(((2 * n + 1) * (2 * n - 1)) / ((n - m) * (n + m))).sqrt();
                let N2 = f64::from(
                    ((2 * n + 1) * (n - m - 1) * (n + m - 1)) / ((2 * n - 3) * (n + m) * (n - m)),
                ).sqrt();
                A_[(n, m)] = u_ * N1 * A_[(n - 1, m)] - N2 * A_[(n - 2, m)];
            }
            Re.push(if m == 0 {
                1.0
            } else {
                s_ * Re[(m - 1) as usize] - t_ * Im[(m - 1) as usize]
            }); // real part of (s + i*t)^m

            Im.push(if m == 0 {
                0.0
            } else {
                s_ * Im[(m - 1) as usize] + t_ * Re[(m - 1) as usize]
            }); // imaginary part of (s + i*t)^m
        }

        // Now do summation ------------------------------------------------
        // initialize recursion
        let rho = self.radius / r_;
        let mut rho_np1 = -self.mu / r_ * rho; // rho(0) ,Ref[3], Eq 26 , factor = mu for gravity
        let mut rho_np2 = rho_np1 * rho;
        let mut a1 = 0.0;
        let mut a2 = 0.0;
        let mut a3 = 0.0;
        let mut a4 = 0.0;
        let mut sqrt2 = 2.0f64.sqrt();
        for n in 1..self.order {
            rho_np1 *= rho;
            rho_np2 *= rho;
            let mut sum1 = 0.0;
            let mut sum2 = 0.0;
            let mut sum3 = 0.0;
            let mut sum4 = 0.0;

            let mut m = 0;
            while m <= n && m <= self.degree {
                let &(Cval, Sval) = self.data.get(&(n, m)).unwrap(); // We expect the data to be there
                                                                     // let Cval = self.Cnm(self.julian_days, n, m);
                                                                     // let Sval = self.Snm(self.julian_days, n, m);
                                                                     // Pines Equation 27 (Part of)
                let D = (Cval * Re[m as usize] + Sval * Im[m as usize]) * sqrt2;
                let E = if m == 0 {
                    0.0
                } else {
                    (Cval * Re[(m - 1) as usize] + Sval * Im[(m - 1) as usize]) * sqrt2
                };
                let F = if m == 0 {
                    0.0
                } else {
                    (Sval * Re[(m - 1) as usize] - Cval * Im[(m - 1) as usize]) * sqrt2
                };
                // Correct for normalization
                let mut VR01 = f64::from((n - m) * (n + m + 1)).sqrt();
                let mut VR11 = (f64::from((2 * n + 1) * (n + m + 2) * (n + m + 1))
                    / (2.0 * f64::from(n) + 3.0))
                    .sqrt();
                if m == 0 {
                    VR01 /= (2.0f64).sqrt();
                    VR11 /= (2.0f64).sqrt();
                }

                let Avv01 = VR01 * A_[(n, m + 1)];
                let Avv11 = VR11 * A_[(n + 1, m + 1)];
                // Pines Equation 30 and 30b (Part of)
                sum1 += m * A_[(n, m)] * E;
                sum2 += m * A_[(n, m)] * F;
                sum3 += Avv01 * D;
                sum4 += Avv11 * D;
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
