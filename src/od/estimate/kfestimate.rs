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

use super::{Estimate, State};
use crate::cosmic::AstroError;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, Matrix, OMatrix, OVector};
use crate::mc::{MultivariateNormal, StateDispersion};
use crate::md::prelude::OrbitDual;
use crate::md::StateParameter;
use crate::Spacecraft;
use na::SMatrix;
use nalgebra::Const;
use rand::SeedableRng;
use rand_distr::Distribution;
use rand_pcg::Pcg64Mcg;
use std::cmp::PartialEq;
use std::error::Error;
use std::fmt;
use std::ops::Mul;

/// Kalman filter Estimate
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct KfEstimate<T: State>
where
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<<T as State>::Size>
        + Allocator<<T as State>::VecLength>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    /// The estimated state
    pub nominal_state: T,
    /// The state deviation
    pub state_deviation: OVector<f64, <T as State>::Size>,
    /// The Covariance of this estimate
    pub covar: OMatrix<f64, <T as State>::Size, <T as State>::Size>,
    /// The predicted covariance of this estimate
    pub covar_bar: OMatrix<f64, <T as State>::Size, <T as State>::Size>,
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    pub predicted: bool,
    /// The STM used to compute this Estimate
    pub stm: OMatrix<f64, <T as State>::Size, <T as State>::Size>,
}

impl<T: State> KfEstimate<T>
where
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<<T as State>::Size>
        + Allocator<<T as State>::VecLength>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    /// Initializes a new filter estimate from the nominal state (not dispersed) and the full covariance
    pub fn from_covar(
        nominal_state: T,
        covar: OMatrix<f64, <T as State>::Size, <T as State>::Size>,
    ) -> Self {
        Self {
            nominal_state,
            state_deviation: OVector::<f64, <T as State>::Size>::zeros(),
            covar,
            covar_bar: covar,
            predicted: true,
            stm: OMatrix::<f64, <T as State>::Size, <T as State>::Size>::identity(),
        }
    }

    /// Initializes a new filter estimate from the nominal state (not dispersed) and the diagonal of the covariance
    pub fn from_diag(nominal_state: T, diag: OVector<f64, <T as State>::Size>) -> Self {
        let covar = Matrix::from_diagonal(&diag);
        Self {
            nominal_state,
            state_deviation: OVector::<f64, <T as State>::Size>::zeros(),
            covar,
            covar_bar: covar,
            predicted: true,
            stm: OMatrix::<f64, <T as State>::Size, <T as State>::Size>::identity(),
        }
    }
}

impl KfEstimate<Spacecraft> {
    /// Generates an initial Kalman filter state estimate dispersed from the nominal state using the provided standard deviation parameters.
    ///
    /// The resulting estimate will have a diagonal covariance matrix constructed from the variances of each parameter.
    /// *Limitation:* This method may not work correctly for all Keplerian orbital elements, refer to
    /// <https://github.com/nyx-space/nyx/issues/339> for details.
    pub fn disperse_from_diag(
        nominal_state: Spacecraft,
        dispersions: Vec<StateDispersion>,
        seed: Option<u128>,
    ) -> Result<Self, Box<dyn Error>> {
        let generator = MultivariateNormal::new(nominal_state, dispersions)?;

        let mut rng = match seed {
            Some(seed) => Pcg64Mcg::new(seed),
            None => Pcg64Mcg::from_entropy(),
        };
        let dispersed_state = generator.sample(&mut rng);

        // Compute the difference between both states
        let delta_orbit = (nominal_state.orbit - dispersed_state.state.orbit).unwrap();

        // Build the covariance as three times the absolute value of the error, squared.

        let diag_data = [
            (3.0 * delta_orbit.radius_km.x.abs()).powi(2),
            (3.0 * delta_orbit.radius_km.y.abs()).powi(2),
            (3.0 * delta_orbit.radius_km.z.abs()).powi(2),
            (3.0 * delta_orbit.velocity_km_s.x.abs()).powi(2),
            (3.0 * delta_orbit.velocity_km_s.y.abs()).powi(2),
            (3.0 * delta_orbit.velocity_km_s.z.abs()).powi(2),
            (3.0 * (nominal_state.srp.cr - dispersed_state.state.srp.cr).abs()).powi(2),
            (3.0 * (nominal_state.drag.cd - dispersed_state.state.drag.cd).abs()).powi(2),
            (3.0 * (nominal_state.fuel_mass_kg - dispersed_state.state.fuel_mass_kg).abs()).powi(2),
        ];

        let diag = OVector::<f64, Const<9>>::from_iterator(diag_data);

        // Build the covar from the diagonal
        let covar = Matrix::from_diagonal(&diag);

        Ok(Self {
            nominal_state: dispersed_state.state,
            state_deviation: OVector::<f64, Const<9>>::zeros(),
            covar,
            covar_bar: covar,
            predicted: true,
            stm: OMatrix::<f64, Const<9>, Const<9>>::identity(),
        })
    }

    /// Builds a multivariate random variable from this estimate's nominal state and covariance, zero mean.
    pub fn to_random_variable(&self) -> Result<MultivariateNormal, Box<dyn Error>> {
        MultivariateNormal::from_spacecraft_cov(
            self.nominal_state,
            self.covar,
            self.state_deviation,
        )
    }

    /// Returns the 1-sigma uncertainty for a given parameter, in that parameter's unit
    ///
    /// This method uses the [OrbitDual] structure to compute the estimate in the hyperdual space
    /// and rotate the nominal covariance into that space.
    pub fn sigma_for(&self, param: StateParameter) -> Result<f64, AstroError> {
        // Build the rotation matrix using Orbit Dual.
        let mut rotmat = SMatrix::<f64, 1, 6>::zeros();
        let orbit_dual = OrbitDual::from(self.nominal_state.orbit);

        let xf_partial = orbit_dual.partial_for(param)?;
        for (cno, val) in [
            xf_partial.wtr_x(),
            xf_partial.wtr_y(),
            xf_partial.wtr_z(),
            xf_partial.wtr_vx(),
            xf_partial.wtr_vy(),
            xf_partial.wtr_vz(),
        ]
        .iter()
        .copied()
        .enumerate()
        {
            rotmat[(0, cno)] = val;
        }

        Ok((rotmat * self.covar.fixed_view::<6, 6>(0, 0) * rotmat.transpose())[(0, 0)].sqrt())
    }

    /// Returns the 6x6 covariance (i.e. square of the sigma/uncertainty) of the SMA, ECC, INC, RAAN, AOP, and True Anomaly.
    pub fn keplerian_covar(&self) -> SMatrix<f64, 6, 6> {
        // Build the rotation matrix using Orbit Dual.
        let mut rotmat = SMatrix::<f64, 6, 6>::zeros();
        let orbit_dual = OrbitDual::from(self.nominal_state.orbit);
        for (pno, param) in [
            StateParameter::SMA,
            StateParameter::Eccentricity,
            StateParameter::Inclination,
            StateParameter::RAAN,
            StateParameter::AoP,
            StateParameter::TrueAnomaly,
        ]
        .iter()
        .copied()
        .enumerate()
        {
            let xf_partial = orbit_dual.partial_for(param).unwrap();
            for (cno, val) in [
                xf_partial.wtr_x(),
                xf_partial.wtr_y(),
                xf_partial.wtr_z(),
                xf_partial.wtr_vx(),
                xf_partial.wtr_vy(),
                xf_partial.wtr_vz(),
            ]
            .iter()
            .copied()
            .enumerate()
            {
                rotmat[(pno, cno)] = val;
            }
        }

        rotmat * self.covar.fixed_view::<6, 6>(0, 0) * rotmat.transpose()
    }
}

impl<T: State> Estimate<T> for KfEstimate<T>
where
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<<T as State>::Size>
        + Allocator<<T as State>::VecLength>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    fn zeros(nominal_state: T) -> Self {
        Self {
            nominal_state,
            state_deviation: OVector::<f64, <T as State>::Size>::zeros(),
            covar: OMatrix::<f64, <T as State>::Size, <T as State>::Size>::zeros(),
            covar_bar: OMatrix::<f64, <T as State>::Size, <T as State>::Size>::zeros(),
            predicted: true,
            stm: OMatrix::<f64, <T as State>::Size, <T as State>::Size>::identity(),
        }
    }

    fn nominal_state(&self) -> T {
        self.nominal_state
    }

    fn state_deviation(&self) -> OVector<f64, <T as State>::Size> {
        self.state_deviation
    }

    fn covar(&self) -> OMatrix<f64, <T as State>::Size, <T as State>::Size> {
        self.covar
    }

    fn predicted_covar(&self) -> OMatrix<f64, <T as State>::Size, <T as State>::Size> {
        self.covar_bar
    }

    fn predicted(&self) -> bool {
        self.predicted
    }
    fn stm(&self) -> &OMatrix<f64, <T as State>::Size, <T as State>::Size> {
        &self.stm
    }
    fn set_state_deviation(&mut self, new_state: OVector<f64, <T as State>::Size>) {
        self.state_deviation = new_state;
    }
    fn set_covar(&mut self, new_covar: OMatrix<f64, <T as State>::Size, <T as State>::Size>) {
        self.covar = new_covar;
    }
}

impl<T: State> fmt::Display for KfEstimate<T>
where
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<<T as State>::Size>
        + Allocator<<T as State>::VecLength>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let dim = <T as State>::Size::dim();
        let word = if self.predicted {
            "Prediction"
        } else {
            "Estimate"
        };
        let mut fmt_cov = Vec::with_capacity(dim);
        for i in 0..dim {
            let unit = if i < 3 {
                "km"
            } else if i < 6 {
                "km/s"
            } else {
                ""
            };
            fmt_cov.push(format!("{:.6} {unit}", &self.covar[(i, i)]));
        }
        write!(
            f,
            "=== {} @ {} -- within 3 sigma: {} ===\nstate {}\nsigmas [{}]\n",
            word,
            &self.epoch(),
            self.within_3sigma(),
            &self.state(),
            fmt_cov.join(", ")
        )
    }
}

impl<T: State> fmt::LowerExp for KfEstimate<T>
where
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<<T as State>::Size>
        + Allocator<<T as State>::VecLength>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "=== PREDICTED: {} ===\nEstState {:e} Covariance {:e}\n=====================",
            &self.predicted, &self.state_deviation, &self.covar
        )
    }
}

impl<T: State> Mul<f64> for KfEstimate<T>
where
    DefaultAllocator: Allocator<<T as State>::Size>
        + Allocator<<T as State>::Size, <T as State>::Size>
        + Allocator<<T as State>::Size>
        + Allocator<<T as State>::VecLength>,
    <DefaultAllocator as Allocator<<T as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<T as State>::Size, <T as State>::Size>>::Buffer<f64>: Copy,
{
    type Output = Self;

    fn mul(mut self, rhs: f64) -> Self::Output {
        self.covar *= rhs;
        self
    }
}

#[cfg(test)]
mod ut_kfest {
    use crate::{
        mc::StateDispersion, md::StateParameter, od::estimate::KfEstimate, Spacecraft,
        GMAT_EARTH_GM,
    };
    use anise::{constants::frames::EARTH_J2000, prelude::Orbit};
    use hifitime::Epoch;

    #[test]
    fn test_estimate_from_disp() {
        let eme2k = EARTH_J2000.with_mu_km3_s2(GMAT_EARTH_GM);
        let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
        let initial_state = Spacecraft::builder()
            .orbit(Orbit::keplerian(
                22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k,
            ))
            .build();

        let initial_estimate = KfEstimate::disperse_from_diag(
            initial_state,
            vec![
                StateDispersion::builder()
                    .param(StateParameter::SMA)
                    .std_dev(1.1)
                    .build(),
                StateDispersion::zero_mean(StateParameter::Inclination, 0.2),
                StateDispersion::zero_mean(StateParameter::RAAN, 0.2),
                StateDispersion::zero_mean(StateParameter::AoP, 0.2),
            ],
            Some(0),
        )
        .unwrap();

        let initial_state_dev = initial_estimate.nominal_state;

        let (init_rss_pos_km, init_rss_vel_km_s, _) =
            initial_state.rss(&initial_state_dev).unwrap();

        let delta = (initial_state.orbit - initial_state_dev.orbit).unwrap();

        println!("Truth initial state:\n{initial_state}\n{initial_state:x}");
        println!("Filter initial state:\n{initial_state_dev}\n{initial_state_dev:x}");
        println!(
            "Initial state dev:\t{init_rss_pos_km:.6} km\t{init_rss_vel_km_s:.6} km/s\n{delta}",
        );
        println!("covariance: {:.6}", initial_estimate.covar);

        // Check that the error is in the square root of the covariance
        assert!(delta.radius_km.x < initial_estimate.covar[(0, 0)].sqrt());
        assert!(delta.radius_km.y < initial_estimate.covar[(1, 1)].sqrt());
        assert!(delta.radius_km.z < initial_estimate.covar[(2, 2)].sqrt());
        assert!(delta.velocity_km_s.x < initial_estimate.covar[(3, 3)].sqrt());
        assert!(delta.velocity_km_s.y < initial_estimate.covar[(4, 4)].sqrt());
        assert!(delta.velocity_km_s.z < initial_estimate.covar[(5, 5)].sqrt());
    }
}
