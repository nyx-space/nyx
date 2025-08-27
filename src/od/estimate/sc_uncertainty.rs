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

use core::fmt;

use anise::errors::MathError;
use anise::{astro::PhysicsResult, errors::PhysicsError};
use nalgebra::{SMatrix, SVector};
use rand::SeedableRng;
use rand_distr::Distribution;
use rand_pcg::Pcg64Mcg;
use typed_builder::TypedBuilder;

use crate::mc::MvnSpacecraft;
use crate::{dynamics::guidance::LocalFrame, Spacecraft};

use super::KfEstimate;

#[derive(Clone, Copy, Debug, TypedBuilder)]
/// Builds a spacecraft uncertainty in different local frames, dispersing any of the parameters of the spacecraft state.
///
/// # Usage
/// Use the `TypeBuilder` trait, e.g `SpacecraftUncertainty::builder().nominal(spacecraft).frame(LocalFrame::RIC).x_km(0.5).y_km(0.5).z_km(0.5).build()`
/// to build an uncertainty on position in the RIC frame of 500 meters on R, I, and C, and zero on all other parameters (velocity components, Cr, Cd, mass).
pub struct SpacecraftUncertainty {
    pub nominal: Spacecraft,
    #[builder(default, setter(strip_option))]
    pub frame: Option<LocalFrame>,
    #[builder(default = 0.5)]
    pub x_km: f64,
    #[builder(default = 0.5)]
    pub y_km: f64,
    #[builder(default = 0.5)]
    pub z_km: f64,
    #[builder(default = 50e-5)]
    pub vx_km_s: f64,
    #[builder(default = 50e-5)]
    pub vy_km_s: f64,
    #[builder(default = 50e-5)]
    pub vz_km_s: f64,
    #[builder(default)]
    pub coeff_reflectivity: f64,
    #[builder(default)]
    pub coeff_drag: f64,
    #[builder(default)]
    pub mass_kg: f64,
}

impl SpacecraftUncertainty {
    /// Builds a Kalman filter estimate for a spacecraft state, ready to ingest into an OD Process.
    ///
    /// Note: this function will rotate from the provided local frame into the inertial frame with the same central body.
    pub fn to_estimate(&self) -> PhysicsResult<KfEstimate<Spacecraft>> {
        if self.x_km < 0.0
            || self.y_km < 0.0
            || self.z_km < 0.0
            || self.vx_km_s < 0.0
            || self.vy_km_s < 0.0
            || self.vz_km_s < 0.0
            || self.coeff_drag < 0.0
            || self.coeff_reflectivity < 0.0
            || self.mass_kg < 0.0
        {
            return Err(PhysicsError::AppliedMath {
                source: MathError::DomainError {
                    value: -0.0,
                    msg: "uncertainties must be positive ",
                },
            });
        }

        // Build the orbit state vector as provided.
        let orbit_vec = SVector::<f64, 6>::new(
            self.x_km,
            self.y_km,
            self.z_km,
            self.vx_km_s,
            self.vy_km_s,
            self.vz_km_s,
        );

        // Rotate into the correct frame.
        let dcm_local2inertial = match self.frame {
            None => LocalFrame::Inertial.dcm_to_inertial(self.nominal.orbit)?,
            Some(frame) => frame.dcm_to_inertial(self.nominal.orbit)?,
        };

        let mut init_covar =
            SMatrix::<f64, 9, 9>::from_diagonal(&SVector::<f64, 9>::from_iterator([
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                self.coeff_reflectivity.powi(2),
                self.coeff_drag.powi(2),
                self.mass_kg.powi(2),
            ]));

        let other_cov = SMatrix::<f64, 6, 6>::from_diagonal(&SVector::<f64, 6>::from_iterator([
            orbit_vec[0].powi(2),
            orbit_vec[1].powi(2),
            orbit_vec[2].powi(2),
            orbit_vec[3].powi(2),
            orbit_vec[4].powi(2),
            orbit_vec[5].powi(2),
        ]));

        let rot_covar =
            dcm_local2inertial.state_dcm().transpose() * other_cov * dcm_local2inertial.state_dcm();

        for i in 0..6 {
            for j in 0..6 {
                init_covar[(i, j)] = rot_covar[(i, j)];
            }
        }

        Ok(KfEstimate::from_covar(self.nominal, init_covar))
    }

    /// Returns an estimation whose nominal state is dispersed.
    /// Known bug in MultiVariate dispersion: https://github.com/nyx-space/nyx/issues/339
    pub fn to_estimate_randomized(
        &self,
        seed: Option<u128>,
    ) -> PhysicsResult<KfEstimate<Spacecraft>> {
        let mut estimate = self.to_estimate()?;
        let mvn =
            MvnSpacecraft::from_spacecraft_cov(self.nominal, estimate.covar, SVector::zeros())
                .unwrap();

        // Setup the RNG
        let mut rng = match seed {
            Some(seed) => Pcg64Mcg::new(seed),
            None => Pcg64Mcg::from_os_rng(),
        };

        // Generate the states, forcing the borrow as specified in the `sample_iter` docs.
        let disp_state = mvn.sample(&mut rng);

        estimate.nominal_state = disp_state.state;

        Ok(estimate)
    }
}

impl fmt::Display for SpacecraftUncertainty {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let frame = match self.frame {
            None => format!("{}", self.nominal.orbit.frame),
            Some(frame) => match frame {
                LocalFrame::Inertial => format!("{}", self.nominal.orbit.frame),
                _ => format!("{frame:?}"),
            },
        };
        writeln!(f, "{}", self.nominal)?;
        writeln!(
            f,
            "{frame}  σ_x = {} km  σ_y = {} km  σ_z = {} km",
            self.x_km, self.y_km, self.z_km
        )?;
        writeln!(
            f,
            "{frame}  σ_vx = {} km/s  σ_vy = {} km/s  σ_vz = {} km/s",
            self.vx_km_s, self.vy_km_s, self.vz_km_s
        )?;
        writeln!(
            f,
            "σ_cr = {}  σ_cd = {}  σ_mass = {} kg",
            self.coeff_reflectivity, self.coeff_drag, self.mass_kg
        )
    }
}

#[cfg(test)]
mod ut_sc_uncertainty {

    use super::{Spacecraft, SpacecraftUncertainty};
    use crate::dynamics::guidance::LocalFrame;
    use crate::md::StateParameter;
    use crate::GMAT_EARTH_GM;
    use anise::constants::frames::EME2000;
    use anise::prelude::{Epoch, Orbit};

    use rstest::*;
    #[fixture]
    fn spacecraft() -> Spacecraft {
        let eme2k = EME2000.with_mu_km3_s2(GMAT_EARTH_GM);

        Spacecraft::builder()
            .orbit(Orbit::keplerian(
                7000.0,
                0.01,
                28.5,
                15.0,
                55.0,
                123.0,
                Epoch::from_gregorian_utc_hms(2024, 2, 29, 1, 2, 3),
                eme2k,
            ))
            .build()
    }

    #[rstest]
    fn test_inertial_frame(spacecraft: Spacecraft) {
        let uncertainty = SpacecraftUncertainty::builder()
            .nominal(spacecraft)
            .x_km(0.5)
            .y_km(0.5)
            .z_km(0.5)
            .vx_km_s(0.5e-3)
            .vy_km_s(0.5e-3)
            .vz_km_s(0.5e-3)
            .build();

        assert!((uncertainty.x_km - 0.5).abs() < f64::EPSILON);
        assert!((uncertainty.y_km - 0.5).abs() < f64::EPSILON);
        assert!((uncertainty.z_km - 0.5).abs() < f64::EPSILON);
        assert!((uncertainty.vx_km_s - 0.5e-3).abs() < f64::EPSILON);
        assert!((uncertainty.vy_km_s - 0.5e-3).abs() < f64::EPSILON);
        assert!((uncertainty.vz_km_s - 0.5e-3).abs() < f64::EPSILON);

        println!("{uncertainty}");

        let estimate = uncertainty.to_estimate().unwrap();

        // Ensure that the covariance is a diagonal.
        for i in 0..6 {
            for j in 0..6 {
                if i == j {
                    if i < 3 {
                        assert!((estimate.covar[(i, j)] - 0.5_f64.powi(2)).abs() < f64::EPSILON);
                    } else {
                        assert!((estimate.covar[(i, j)] - 0.5e-3_f64.powi(2)).abs() < f64::EPSILON);
                    }
                } else {
                    assert_eq!(estimate.covar[(i, j)], 0.0);
                }
            }
        }
    }

    #[rstest]
    fn test_ric_frame(spacecraft: Spacecraft) {
        let uncertainty = SpacecraftUncertainty::builder()
            .nominal(spacecraft)
            .frame(LocalFrame::RIC)
            .x_km(0.5)
            .y_km(0.5)
            .z_km(0.5)
            .vx_km_s(0.5e-3)
            .vy_km_s(0.5e-3)
            .vz_km_s(0.5e-3)
            .build();

        assert!((uncertainty.x_km - 0.5).abs() < f64::EPSILON);
        assert!((uncertainty.y_km - 0.5).abs() < f64::EPSILON);
        assert!((uncertainty.z_km - 0.5).abs() < f64::EPSILON);
        assert!((uncertainty.vx_km_s - 0.5e-3).abs() < f64::EPSILON);
        assert!((uncertainty.vy_km_s - 0.5e-3).abs() < f64::EPSILON);
        assert!((uncertainty.vz_km_s - 0.5e-3).abs() < f64::EPSILON);

        println!("{uncertainty}");

        let estimate = uncertainty.to_estimate().unwrap();

        println!("{estimate}");

        let sma_sigma_km = estimate.sigma_for(StateParameter::SMA).unwrap();
        assert!((sma_sigma_km - 1.3528088887049263).abs() < f64::EPSILON);

        let covar_keplerian = estimate.keplerian_covar();
        for i in 1..=3 {
            let val = covar_keplerian[(i, i)].sqrt();
            assert!(val.abs() < 1e-1);
        }
        for i in 4..6 {
            let val = covar_keplerian[(i, i)].sqrt();
            assert!(val.abs() < 0.8);
        }

        // Rotate back into the RIC frame.
        let dcm_ric2inertial = estimate
            .nominal_state
            .orbit
            .dcm_from_ric_to_inertial()
            .unwrap()
            .state_dcm();

        // Build the matrix view of the orbit part of the covariance.
        let orbit_cov = estimate.covar.fixed_view::<6, 6>(0, 0);

        // Rotate back into the RIC frame
        let ric_covar = dcm_ric2inertial * orbit_cov * dcm_ric2inertial.transpose();

        println!("{ric_covar:.9}");
        for i in 0..6 {
            for j in 0..6 {
                if i == j {
                    // We don't abs the data in the sqrt to ensure that the diag is positive.
                    if i < 3 {
                        assert!((ric_covar[(i, j)].sqrt() - 0.5).abs() < 1e-9);
                    } else {
                        assert!((ric_covar[(i, j)].sqrt() - 0.5e-3).abs() < 1e-9);
                    }
                }
            }
        }
    }
}
