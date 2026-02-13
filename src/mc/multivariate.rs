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

use super::{DispersedState, StateDispersion};
use crate::cosmic::AstroPhysicsSnafu;
use crate::errors::StateError;
use crate::md::prelude::BPlane;
use crate::md::{AstroSnafu, StateParameter};
use crate::{pseudo_inverse, NyxError, Spacecraft, State};
use anise::analysis::prelude::OrbitalElement;
use anise::astro::orbit_gradient::OrbitGrad;
use nalgebra::{DMatrix, DVector, SMatrix, SVector};
use rand_distr::{Distribution, Normal};
use snafu::ResultExt;
use std::error::Error;

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// A multivariate spacecraft state generator for Monte Carlo analyses. Ensures that the covariance is properly applied on all provided state variables.
#[derive(Clone, Debug)]
#[cfg_attr(feature = "python", pyclass)]
pub struct MvnSpacecraft {
    /// The template state
    pub template: Spacecraft,
    pub dispersions: Vec<StateDispersion>,
    /// The mean of the multivariate normal distribution
    pub mean: SVector<f64, 9>,
    /// The dot product \sqrt{\vec s} \cdot \vec v, where S is the singular values and V the V matrix from the SVD decomp of the covariance of multivariate normal distribution
    pub sqrt_s_v: SMatrix<f64, 9, 9>,
    /// The standard normal distribution used to seed the multivariate normal distribution
    pub std_norm_distr: Normal<f64>,
}

impl MvnSpacecraft {
    /// Creates a new mulivariate state generator from a mean and covariance on the set of state parameters.
    /// The covariance must be positive semi definite.
    ///
    /// # Algorithm
    /// This function will build the rotation matrix to rotate the requested dispersions into the Spacecraft state space using [OrbitDual].
    /// If there are any dispersions on the Cr and Cd, then these are dispersed independently (because they are iid).
    pub fn new(
        template: Spacecraft,
        dispersions: Vec<StateDispersion>,
    ) -> Result<Self, Box<dyn Error>> {
        let mut cov = SMatrix::<f64, 9, 9>::zeros();
        let mut mean = SVector::<f64, 9>::zeros();

        let orbit_dual = OrbitGrad::from(template.orbit);
        let mut b_plane = None;
        for obj in &dispersions {
            if obj.param.is_b_plane() {
                b_plane = Some(
                    BPlane::from_dual(orbit_dual)
                        .context(AstroSnafu)
                        .map_err(Box::new)?,
                );
                break;
            }
        }

        let num_orbital = dispersions
            .iter()
            .filter(|disp| disp.param.is_orbital())
            .count();

        if num_orbital > 0 {
            // Build the rotation matrix from the orbital dispersions to the Cartesian state.
            let mut jac = DMatrix::from_element(num_orbital, 6, 0.0);
            let mut covar = DMatrix::from_element(num_orbital, num_orbital, 0.0);
            let mut means = DVector::from_element(num_orbital, 0.0);
            let orbit_dual = OrbitGrad::from(template.orbit);

            for (rno, disp) in dispersions
                .iter()
                .filter(|d| d.param.is_orbital())
                .enumerate()
            {
                let partial = if let StateParameter::Element(oe) = disp.param {
                    orbit_dual.partial_for(oe).context(AstroPhysicsSnafu)?
                } else if disp.param.is_b_plane() {
                    match disp.param {
                        StateParameter::BdotR() => b_plane.unwrap().b_r_km,
                        StateParameter::BdotT() => b_plane.unwrap().b_t_km,
                        StateParameter::BLTOF() => b_plane.unwrap().ltof_s,
                        _ => unreachable!(),
                    }
                } else {
                    unreachable!()
                };

                println!("{partial}");

                for (cno, val) in [
                    partial.wrt_x(),
                    partial.wrt_y(),
                    partial.wrt_z(),
                    partial.wrt_vx(),
                    partial.wrt_vy(),
                    partial.wrt_vz(),
                ]
                .iter()
                .copied()
                .enumerate()
                {
                    jac[(rno, cno)] = val;
                }
                covar[(rno, rno)] = disp.std_dev.unwrap_or(0.0).powi(2);
                means[rno] = disp.mean.unwrap_or(0.0);
            }

            // Now that we have the Jacobian that rotates from the Cartesian elements to the dispersions parameters,
            // let's compute the inverse of this Jacobian to rotate from the dispersion params into the Cartesian elements.
            let jac_inv = pseudo_inverse!(&jac)?;

            // Rotate the orbital covariance back into the Cartesian state space, making this a 6x6.
            let orbit_cov = &jac_inv * &covar * jac_inv.transpose();
            // Rotate the means into the Cartesian space
            let cartesian_mean = jac_inv * means;
            for ii in 0..6 {
                for jj in 0..6 {
                    cov[(ii, jj)] = orbit_cov[(ii, jj)];
                }
                mean[ii] = cartesian_mean[ii];
            }
        };

        if dispersions.len() > num_orbital {
            for disp in &dispersions {
                if disp.param.is_orbital() {
                    continue;
                } else {
                    match disp.param {
                        StateParameter::Cr() => {
                            cov[(7, 7)] = disp.mean.unwrap_or(0.0).powi(2);
                        }
                        StateParameter::Cd() => {
                            cov[(8, 8)] = disp.mean.unwrap_or(0.0).powi(2);
                        }
                        StateParameter::DryMass() | StateParameter::PropMass() => {
                            cov[(9, 9)] = disp.mean.unwrap_or(0.0).powi(2);
                        }
                        _ => return Err(Box::new(StateError::ReadOnly { param: disp.param })),
                    }
                }
            }
        }

        // At this point, the cov matrix is a 9x9 with all dispersions transformed into the Cartesian state space.
        let svd = cov.svd_unordered(false, true);
        if svd.v_t.is_none() {
            return Err(Box::new(NyxError::CovarianceMatrixNotPsd));
        }

        let sqrt_s = svd.singular_values.map(|x| x.sqrt());
        let mut sqrt_s_v_t = svd.v_t.unwrap().transpose();

        for (i, mut col) in sqrt_s_v_t.column_iter_mut().enumerate() {
            col *= sqrt_s[i];
        }

        Ok(Self {
            template,
            dispersions,
            mean,
            sqrt_s_v: sqrt_s_v_t,
            std_norm_distr: Normal::new(0.0, 1.0).unwrap(),
        })
    }

    /// Same as `new` but with a zero mean
    pub fn zero_mean(
        template: Spacecraft,
        mut dispersions: Vec<StateDispersion>,
    ) -> Result<Self, Box<dyn Error>> {
        for disp in &mut dispersions {
            disp.mean = Some(0.0);
        }

        Self::new(template, dispersions)
    }

    /// Initializes a new multivariate distribution using the state data in the spacecraft state space.
    pub fn from_spacecraft_cov(
        template: Spacecraft,
        cov: SMatrix<f64, 9, 9>,
        mean: SVector<f64, 9>,
    ) -> Result<Self, Box<dyn Error>> {
        // Check that covariance is PSD by ensuring that all the eigenvalues are positive or nil
        match cov.eigenvalues() {
            None => return Err(Box::new(NyxError::CovarianceMatrixNotPsd)),
            Some(evals) => {
                for eigenval in &evals {
                    if *eigenval < 0.0 {
                        return Err(Box::new(NyxError::CovarianceMatrixNotPsd));
                    }
                }
            }
        };

        let svd = cov.svd_unordered(false, true);
        if svd.v_t.is_none() {
            return Err(Box::new(NyxError::CovarianceMatrixNotPsd));
        }

        let s = svd.singular_values;
        // Item by item multiplication
        let mut sqrt_s_v = svd.v_t.unwrap().transpose();
        for (i, mut col) in sqrt_s_v.column_iter_mut().enumerate() {
            col *= s[i].sqrt();
        }

        let dispersions = vec![
            StateDispersion::builder()
                .param(StateParameter::Element(OrbitalElement::X))
                .std_dev(cov[(0, 0)])
                .build(),
            StateDispersion::builder()
                .param(StateParameter::Element(OrbitalElement::Y))
                .std_dev(cov[(1, 1)])
                .build(),
            StateDispersion::builder()
                .param(StateParameter::Element(OrbitalElement::Z))
                .std_dev(cov[(2, 2)])
                .build(),
            StateDispersion::builder()
                .param(StateParameter::Element(OrbitalElement::VX))
                .std_dev(cov[(3, 3)])
                .build(),
            StateDispersion::builder()
                .param(StateParameter::Element(OrbitalElement::VY))
                .std_dev(cov[(4, 4)])
                .build(),
            StateDispersion::builder()
                .param(StateParameter::Element(OrbitalElement::VZ))
                .std_dev(cov[(5, 5)])
                .build(),
            StateDispersion::builder()
                .param(StateParameter::Cr())
                .std_dev(cov[(6, 6)])
                .build(),
            StateDispersion::builder()
                .param(StateParameter::Cd())
                .std_dev(cov[(7, 7)])
                .build(),
            StateDispersion::builder()
                .param(StateParameter::PropMass())
                .std_dev(cov[(8, 8)])
                .build(),
        ];

        Ok(Self {
            template,
            dispersions,
            mean,
            sqrt_s_v,
            std_norm_distr: Normal::new(0.0, 1.0).unwrap(),
        })
    }
}

impl Distribution<DispersedState<Spacecraft>> for MvnSpacecraft {
    fn sample<R: rand::Rng + ?Sized>(&self, rng: &mut R) -> DispersedState<Spacecraft> {
        // Generate the vector representing the state
        let x_rng = SVector::<f64, 9>::from_fn(|_, _| self.std_norm_distr.sample(rng));
        let x = self.sqrt_s_v * x_rng + self.mean;
        let mut state = self.template;

        // Set the new state data
        for (coord, val) in x.iter().copied().enumerate() {
            if coord < 3 {
                state.orbit.radius_km[coord] += val;
            } else if coord < 6 {
                state.orbit.velocity_km_s[coord % 3] += val;
            } else if coord == 6 {
                state.srp.coeff_reflectivity += val;
            } else if coord == 7 {
                state.drag.coeff_drag += val;
            } else if coord == 8 {
                state.mass.prop_mass_kg += val;
            }
        }

        let mut actual_dispersions = Vec::new();
        for disp in &self.dispersions {
            // Compute the delta
            let delta = self.template.value(disp.param).unwrap() - state.value(disp.param).unwrap();
            actual_dispersions.push((disp.param, delta));
        }

        DispersedState {
            state,
            actual_dispersions,
        }
    }
}

#[cfg(test)]
mod multivariate_ut {
    use super::*;
    use crate::time::Epoch;
    use crate::Spacecraft;
    use crate::GMAT_EARTH_GM;
    use anise::constants::frames::EARTH_J2000;
    use anise::prelude::Orbit;
    use rand::Rng;
    use statrs;
    use statrs::distribution::ContinuousCDF;

    /// Returns the Mahalanobis distance of a random variable x from a multivariate normal distribution
    /// with mean mu and covariance matrix S.
    fn mahalanobis_distance<const T: usize>(
        x: &SVector<f64, T>,
        mu: &SVector<f64, T>,
        cov_inv: &SMatrix<f64, T, T>,
    ) -> f64 {
        let delta = x - mu;
        (delta.transpose() * cov_inv * delta)[(0, 0)]
    }

    #[test]
    fn mahalanobis_distance_test() {
        // Simple 2D case with diagonal covariance
        let x = SVector::<f64, 2>::new(1.0, 2.0);
        let mu = SVector::<f64, 2>::new(0.0, 0.0);
        let cov = SMatrix::<f64, 2, 2>::from_diagonal(&SVector::<f64, 2>::new(2.0, 3.0));
        let cov_inv = cov.pseudo_inverse(1e-12).unwrap();

        // Expected distance: (1-0)^2/2 + (2-0)^2/3 = 1/2 + 4/3 = 0.5 + 1.333... = 1.8333...
        let md = mahalanobis_distance(&x, &mu, &cov_inv);
        assert!((md - 1.8333333333333333).abs() < 1e-9);
    }

    #[test]
    fn test_mvn_generator() {
        // Generate a random covariance matrix.
        let mut rng = rand::rng();
        let cov = SMatrix::<f64, 6, 6>::from_fn(|_, _| rng.random());
        let cov = cov * cov.transpose();
        let mut cov_resized = SMatrix::<f64, 9, 9>::zeros();
        cov_resized.fixed_view_mut::<6, 6>(0, 0).copy_from(&cov);

        let eme2k = EARTH_J2000.with_mu_km3_s2(GMAT_EARTH_GM);

        let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
        let state = Spacecraft::builder()
            .orbit(Orbit::keplerian(
                8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k,
            ))
            .build();

        let mvn = MvnSpacecraft::from_spacecraft_cov(state, cov_resized, SVector::zeros()).unwrap();

        // Generate N samples from this distribution
        let n = 1000;
        let samples = mvn.sample_iter(&mut rng).take(n);

        // Compute the Mahalanobis distance for each sample
        let cov_inv = cov_resized.pseudo_inverse(1e-12).unwrap();
        let md: Vec<f64> = samples
            .map(|sample| {
                let mut v = SVector::<f64, 9>::zeros();
                for i in 0..6 {
                    v[i] = sample.state.orbit.to_cartesian_pos_vel()[i]
                        - state.orbit.to_cartesian_pos_vel()[i];
                }
                mahalanobis_distance(&v, &SVector::zeros(), &cov_inv)
            })
            .collect();

        // Perform a chi-squared test on the Mahalanobis distances.
        // The squared Mahalanobis distance of a sample from the mean of a multivariate normal distribution
        // with k degrees of freedom follows a chi-squared distribution with k degrees of freedom.
        // We will test if the 95th percentile of the Mahalanobis distances is within the 95th percentile
        // of the chi-squared distribution.
        let mut md = md;
        md.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let p95_md = md[(0.95 * n as f64) as usize];

        let chi_squared = statrs::distribution::ChiSquared::new(6.0).unwrap();
        let p95_chi_squared = chi_squared.inverse_cdf(0.95);

        assert!((p95_md - p95_chi_squared).abs() / p95_chi_squared < 0.2);
    }

    #[test]
    fn disperse_r_mag() {
        use anise::constants::frames::EARTH_J2000;
        use anise::prelude::Orbit;

        use crate::time::Epoch;
        use rand_pcg::Pcg64Mcg;

        // Ensure that this worked: a 3 sigma deviation around 1 km means we shouldn't have 99.7% of samples within those bounds.
        // Create a reproducible fast seed
        let seed = 0;

        let eme2k = EARTH_J2000.with_mu_km3_s2(GMAT_EARTH_GM);

        let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
        let state = Spacecraft::builder()
            .orbit(Orbit::keplerian(
                8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k,
            ))
            .build();

        // Check that we can modify the radius magnitude
        let std_dev = 1.0;
        let generator = MvnSpacecraft::new(
            state,
            vec![StateDispersion::builder()
                .param(StateParameter::Element(OrbitalElement::Rmag))
                .std_dev(std_dev)
                .build()],
        )
        .unwrap();

        let rng = Pcg64Mcg::new(seed);
        let init_rmag = state.orbit.rmag_km();
        let cnt_too_far: u16 = generator
            .sample_iter(rng)
            .take(1000)
            .map(|dispersed_state| {
                if (init_rmag - dispersed_state.state.orbit.rmag_km()).abs() >= 3.0 * std_dev {
                    1
                } else {
                    0
                }
            })
            .sum::<u16>();

        // We specified a seed so we know exactly what to expect and we've reset the seed to 0.
        assert_eq!(
            cnt_too_far,
            6, // Mathematically, this should be 3!
            "Should have about 3% of samples being more than 3 sigma away, got {cnt_too_far}"
        );
    }

    #[test]
    fn disperse_full_cartesian() {
        use anise::constants::frames::EARTH_J2000;
        use anise::prelude::Orbit;

        use crate::time::Epoch;
        use crate::Spacecraft;
        use crate::GMAT_EARTH_GM;

        use rand_pcg::Pcg64Mcg;

        let eme2k = EARTH_J2000.with_mu_km3_s2(GMAT_EARTH_GM);

        let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
        let state = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);

        let std_dev = [10.0, 10.0, 10.0, 0.2, 0.2, 0.2, 0.0, 0.0, 0.0];

        let generator = MvnSpacecraft::new(
            Spacecraft {
                orbit: state,
                ..Default::default()
            },
            vec![
                StateDispersion::builder()
                    .param(StateParameter::Element(OrbitalElement::X))
                    .std_dev(10.0)
                    .build(),
                StateDispersion::builder()
                    .param(StateParameter::Element(OrbitalElement::Y))
                    .std_dev(10.0)
                    .build(),
                StateDispersion::builder()
                    .param(StateParameter::Element(OrbitalElement::Z))
                    .std_dev(10.0)
                    .build(),
                StateDispersion::builder()
                    .param(StateParameter::Element(OrbitalElement::VX))
                    .std_dev(0.2)
                    .build(),
                StateDispersion::builder()
                    .param(StateParameter::Element(OrbitalElement::VY))
                    .std_dev(0.2)
                    .build(),
                StateDispersion::builder()
                    .param(StateParameter::Element(OrbitalElement::VZ))
                    .std_dev(0.2)
                    .build(),
            ],
        )
        .unwrap();

        // Ensure that this worked: a 3 sigma deviation around 1 km means we shouldn't have 99.7% of samples within those bounds.
        // Create a reproducible fast seed
        let seed = 0;
        let rng = Pcg64Mcg::new(seed);

        let cnt_too_far: u16 = generator
            .sample_iter(rng)
            .take(1000)
            .map(|dispersed_state| {
                let mut cnt = 0;
                for (idx, val_std_dev) in std_dev.iter().take(6).enumerate() {
                    let cur_val = dispersed_state.state.to_vector()[idx];
                    let nom_val = state.to_cartesian_pos_vel()[idx];
                    if (cur_val - nom_val).abs() > *val_std_dev {
                        cnt += 1;
                    }
                }
                cnt
            })
            .sum::<u16>();

        // We specified a seed so we know exactly what to expect
        assert_eq!(
            cnt_too_far / 6,
            312,
            "Should have about 3% of samples being more than 3 sigma away, got {cnt_too_far}"
        );
    }

    #[test]
    fn disperse_raan_only() {
        use anise::constants::frames::EARTH_J2000;
        use anise::prelude::Orbit;

        use crate::time::Epoch;
        use rand_pcg::Pcg64Mcg;

        let eme2k = EARTH_J2000.with_mu_km3_s2(GMAT_EARTH_GM);

        let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
        let state = Spacecraft::builder()
            .orbit(Orbit::keplerian(
                8_100.0, 1e-6, 12.85, 356.614, 14.19, 199.887_7, dt, eme2k,
            ))
            .build();

        let angle_sigma_deg = 0.2;

        assert!(StateParameter::Element(OrbitalElement::RAAN).is_orbital());

        let generator = MvnSpacecraft::new(
            state,
            vec![StateDispersion::zero_mean(
                StateParameter::Element(OrbitalElement::RAAN),
                angle_sigma_deg,
            )],
        )
        .unwrap();

        // Ensure that this worked: a 3 sigma deviation around 1 km means we shouldn't have 99.7% of samples within those bounds.
        // Create a reproducible fast seed
        let seed = 0;
        let rng = Pcg64Mcg::new(seed);
        let n = 1000;
        let samples = generator.sample_iter(rng).take(n);

        let cov = SMatrix::<f64, 1, 1>::from_diagonal(&SVector::<f64, 1>::from_vec(vec![
            angle_sigma_deg.powi(2),
        ]));
        let cov_inv = cov.pseudo_inverse(1e-12).unwrap();

        let md: Vec<f64> = samples
            .map(|dispersed_state| {
                // For all other orbital parameters, make sure that we have not changed things dramatically.
                for param in [
                    StateParameter::Element(OrbitalElement::SemiMajorAxis),
                    StateParameter::Element(OrbitalElement::Inclination),
                ] {
                    let orig = state.value(param).unwrap();
                    let new = dispersed_state.state.value(param).unwrap();
                    let prct_change = 100.0 * (orig - new).abs() / orig;
                    assert!(
                        prct_change < 5.0,
                        "{param} changed by {prct_change:.3} % ({orig:.3e} -> {new:.3e})"
                    );
                }

                let delta =
                    SVector::<f64, 1>::from_vec(vec![dispersed_state.actual_dispersions[0].1]);
                mahalanobis_distance(&delta, &SVector::zeros(), &cov_inv)
            })
            .collect();

        let mut md = md;
        md.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let p95_md = md[(0.95 * n as f64) as usize];

        let chi_squared = statrs::distribution::ChiSquared::new(1.0).unwrap();
        let p95_chi_squared = chi_squared.inverse_cdf(0.95);

        assert!((p95_md - p95_chi_squared).abs() / p95_chi_squared < 0.2);
    }

    #[test]
    fn disperse_keplerian() {
        use anise::constants::frames::EARTH_J2000;
        use anise::prelude::Orbit;

        use crate::time::Epoch;
        use rand_pcg::Pcg64Mcg;

        let eme2k = EARTH_J2000.with_mu_km3_s2(GMAT_EARTH_GM);

        let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
        let state = Spacecraft::builder()
            .orbit(Orbit::keplerian(
                8_100.0, 1e-6, 12.85, 356.614, 14.19, 199.887_7, dt, eme2k,
            ))
            .build();

        let sma_sigma_km = 10.0;
        let inc_sigma_deg = 0.15;
        let angle_sigma_deg = 0.02;

        let generator = MvnSpacecraft::new(
            state,
            vec![
                StateDispersion::zero_mean(
                    StateParameter::Element(OrbitalElement::SemiMajorAxis),
                    sma_sigma_km,
                ),
                StateDispersion::zero_mean(
                    StateParameter::Element(OrbitalElement::Inclination),
                    inc_sigma_deg,
                ),
                StateDispersion::zero_mean(
                    StateParameter::Element(OrbitalElement::RAAN),
                    angle_sigma_deg,
                ),
                StateDispersion::zero_mean(
                    StateParameter::Element(OrbitalElement::AoP),
                    angle_sigma_deg,
                ),
            ],
        )
        .unwrap();

        // The generator computes the equivalent Cartesian covariance. Let's retrieve it.
        let expected_cart_cov = &generator.sqrt_s_v * generator.sqrt_s_v.transpose();

        // Create a reproducible fast seed, and generate the samples in the Cartesian space
        let seed = 0;
        let rng = Pcg64Mcg::new(seed);
        let n = 2000; // Increase samples for better stats
        let samples: Vec<SVector<f64, 6>> = generator
            .sample_iter(rng)
            .take(n)
            .map(|s| s.state.orbit.to_cartesian_pos_vel())
            .collect();

        let nominal_cart = state.orbit.to_cartesian_pos_vel();

        // Check that the mean of the Cartesian states is close to the nominal
        let mut mean_cart = SVector::<f64, 6>::zeros();
        for sample in &samples {
            mean_cart += sample;
        }
        mean_cart /= n as f64;

        let mean_diff = mean_cart - nominal_cart;
        // We expect some deviation because of the non-linearities.
        assert!(mean_diff.norm() < 1.0);

        // Check that the covariance of the Cartesian states is close to the expected one
        let mut sample_cov = SMatrix::<f64, 6, 6>::zeros();
        for sample in &samples {
            let disp_vec = sample - &mean_cart;
            sample_cov += &disp_vec * disp_vec.transpose();
        }
        sample_cov /= (n - 1) as f64;

        let expected_cov_6x6 = expected_cart_cov.fixed_view::<6, 6>(0, 0);

        let diff = sample_cov - expected_cov_6x6;
        assert!(diff.norm() / expected_cov_6x6.norm() < 0.2);
    }
}
