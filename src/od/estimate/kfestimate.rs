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
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, Matrix, OMatrix, OVector};
use crate::mc::GaussianGenerator;
use crate::md::StateParameter;
use crate::Spacecraft;
use nalgebra::Const;
use rand::SeedableRng;
use rand_distr::Distribution;
use rand_pcg::Pcg64Mcg;
use std::cmp::PartialEq;
use std::fmt;

/// Kalman filter Estimate
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct KfEstimate<T: State>
where
    DefaultAllocator: Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<usize, <T as State>::Size>
        + Allocator<f64, <T as State>::VecLength>
        + Allocator<usize, <T as State>::Size, <T as State>::Size>,
    <DefaultAllocator as Allocator<f64, <T as State>::Size>>::Buffer: Copy,
    <DefaultAllocator as Allocator<f64, <T as State>::Size, <T as State>::Size>>::Buffer: Copy,
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
    DefaultAllocator: Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<usize, <T as State>::Size>
        + Allocator<f64, <T as State>::VecLength>
        + Allocator<usize, <T as State>::Size, <T as State>::Size>,
    <DefaultAllocator as Allocator<f64, <T as State>::Size>>::Buffer: Copy,
    <DefaultAllocator as Allocator<f64, <T as State>::Size, <T as State>::Size>>::Buffer: Copy,
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
    /// *Limitation:* This method incorrectly assumes all parameters are statistically independent.
    pub fn disperse_from_diag(
        nominal_state: Spacecraft,
        params: &[(StateParameter, f64)],
        seed: Option<u128>,
    ) -> Self {
        // Build a generator.
        let gen = GaussianGenerator::from_3std_devs(nominal_state, params).unwrap();

        let mut rng = match seed {
            Some(seed) => Pcg64Mcg::new(seed),
            None => Pcg64Mcg::from_entropy(),
        };
        let dispersed_state = gen.sample(&mut rng);

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

        Self {
            nominal_state: dispersed_state.state,
            state_deviation: OVector::<f64, Const<9>>::zeros(),
            covar,
            covar_bar: covar,
            predicted: true,
            stm: OMatrix::<f64, Const<9>, Const<9>>::identity(),
        }
    }
}

impl<T: State> Estimate<T> for KfEstimate<T>
where
    DefaultAllocator: Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<usize, <T as State>::Size>
        + Allocator<f64, <T as State>::VecLength>
        + Allocator<usize, <T as State>::Size, <T as State>::Size>,
    <DefaultAllocator as Allocator<f64, <T as State>::Size>>::Buffer: Copy,
    <DefaultAllocator as Allocator<f64, <T as State>::Size, <T as State>::Size>>::Buffer: Copy,
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
    DefaultAllocator: Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<usize, <T as State>::Size>
        + Allocator<f64, <T as State>::VecLength>
        + Allocator<usize, <T as State>::Size, <T as State>::Size>,
    <DefaultAllocator as Allocator<f64, <T as State>::Size>>::Buffer: Copy,
    <DefaultAllocator as Allocator<f64, <T as State>::Size, <T as State>::Size>>::Buffer: Copy,
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
            fmt_cov.push(format!("{:e}", &self.covar[(i, i)]));
        }
        write!(
            f,
            "=== {} @ {} -- within 3 sigma: {} ===\nstate {}\nsigmas [{}]\n",
            word,
            &self.epoch(),
            self.within_3sigma(),
            &self.state(),
            fmt_cov.join(",")
        )
    }
}

impl<T: State> fmt::LowerExp for KfEstimate<T>
where
    DefaultAllocator: Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<usize, <T as State>::Size>
        + Allocator<f64, <T as State>::VecLength>
        + Allocator<usize, <T as State>::Size, <T as State>::Size>,
    <DefaultAllocator as Allocator<f64, <T as State>::Size>>::Buffer: Copy,
    <DefaultAllocator as Allocator<f64, <T as State>::Size, <T as State>::Size>>::Buffer: Copy,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "=== PREDICTED: {} ===\nEstState {:e} Covariance {:e}\n=====================",
            &self.predicted, &self.state_deviation, &self.covar
        )
    }
}

#[cfg(test)]
mod ut_kfest {
    use crate::{md::StateParameter, od::estimate::KfEstimate, Spacecraft, GMAT_EARTH_GM};
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
            &[
                (StateParameter::SMA, 1.1),
                (StateParameter::Inclination, 0.0025),
                (StateParameter::RAAN, 0.022),
                (StateParameter::AoP, 0.02),
            ],
            Some(0),
        );

        let initial_state_dev = initial_estimate.nominal_state;

        let (init_rss_pos_km, init_rss_vel_km_s, _) =
            initial_state.rss(&initial_state_dev).unwrap();

        let delta = (initial_state.orbit - initial_state_dev.orbit).unwrap();

        println!("Truth initial state:\n{initial_state}\n{initial_state:x}");
        println!("Filter initial state:\n{initial_state_dev}\n{initial_state_dev:x}");
        println!(
            "Initial state dev:\t{init_rss_pos_km:.3} km\t{init_rss_vel_km_s:.3} km/s\n{delta}",
        );
        println!("covariance: {}", initial_estimate.covar);

        // Check that the error is in the square root of the covariance
        assert!(delta.radius_km.x < initial_estimate.covar[(0, 0)].sqrt());
        assert!(delta.radius_km.y < initial_estimate.covar[(1, 1)].sqrt());
        assert!(delta.radius_km.z < initial_estimate.covar[(2, 2)].sqrt());
        assert!(delta.velocity_km_s.x < initial_estimate.covar[(3, 3)].sqrt());
        assert!(delta.velocity_km_s.y < initial_estimate.covar[(4, 4)].sqrt());
        assert!(delta.velocity_km_s.z < initial_estimate.covar[(5, 5)].sqrt());
    }
}
