/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::msr::StdMeasurement;
use super::{CovarFormat, EpochFormat, Measurement};
use super::{EstimateFrom, State};
use crate::cosmic::Orbit;
use crate::hifitime::Epoch;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, Matrix, OMatrix, OVector, Vector6, U6};
use crate::mc::GaussianGenerator;
use crate::md::StateParameter;
use crate::Spacecraft;
use nalgebra::Matrix2x6;
use rand::SeedableRng;
use rand_distr::Distribution;
use rand_pcg::Pcg64Mcg;
use serde::ser::SerializeSeq;
use serde::{Serialize, Serializer};
use std::cmp::PartialEq;
use std::fmt;

/// Stores an Estimate, as the result of a `time_update` or `measurement_update`.
pub trait Estimate<T: State>
where
    Self: Clone + PartialEq + Sized + fmt::Display,
    DefaultAllocator: Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<f64, <T as State>::VecLength>,
{
    /// An empty estimate. This is useful if wanting to store an estimate outside the scope of a filtering loop.
    fn zeros(state: T) -> Self;
    /// Epoch of this Estimate
    fn epoch(&self) -> Epoch {
        self.state().epoch()
    }
    // Sets the epoch
    fn set_epoch(&mut self, dt: Epoch) {
        self.state().set_epoch(dt);
    }
    /// The estimated state
    fn state(&self) -> T {
        self.nominal_state().add(self.state_deviation())
    }
    /// The state deviation as computed by the filter.
    fn state_deviation(&self) -> OVector<f64, <T as State>::Size>;
    /// The nominal state as reported by the filter dynamics
    fn nominal_state(&self) -> T;
    /// The Covariance of this estimate. Will return the predicted covariance if this is a time update/prediction.
    fn covar(&self) -> OMatrix<f64, <T as State>::Size, <T as State>::Size>;
    /// The predicted covariance of this estimate from the time update
    fn predicted_covar(&self) -> OMatrix<f64, <T as State>::Size, <T as State>::Size>;
    /// Sets the state deviation.
    fn set_state_deviation(&mut self, new_state: OVector<f64, <T as State>::Size>);
    /// Sets the Covariance of this estimate
    fn set_covar(&mut self, new_covar: OMatrix<f64, <T as State>::Size, <T as State>::Size>);
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    fn predicted(&self) -> bool;
    /// The STM used to compute this Estimate
    fn stm(&self) -> &OMatrix<f64, <T as State>::Size, <T as State>::Size>;
    /// The Epoch format upon serialization
    fn epoch_fmt(&self) -> EpochFormat;
    /// The covariance format upon serialization
    fn covar_fmt(&self) -> CovarFormat;
    /// Returns whether this estimate is within some bound
    /// The 68-95-99.7 rule is a good way to assess whether the filter is operating normally
    fn within_sigma(&self, sigma: f64) -> bool {
        let state = self.state_deviation();
        let covar = self.covar();
        for i in 0..state.len() {
            let bound = covar[(i, i)].sqrt() * sigma;
            if state[i] > bound || state[i] < -bound {
                return false;
            }
        }
        true
    }
    /// Returns whether this estimate is within 3 sigma, which represent 99.7% for a Normal distribution
    fn within_3sigma(&self) -> bool {
        self.within_sigma(3.0)
    }
    /// Returns the header
    fn header(epoch_fmt: EpochFormat, covar_fmt: CovarFormat) -> Vec<String> {
        let dim = <T as State>::Size::dim();
        let mut hdr_v = Vec::with_capacity(3 * dim + 1);
        hdr_v.push(format!("{epoch_fmt}"));
        for i in 0..dim {
            hdr_v.push(format!("state_{i}"));
        }
        // Serialize the covariance
        for i in 0..dim {
            for j in 0..dim {
                hdr_v.push(format!("{covar_fmt}_{i}_{j}"));
            }
        }
        hdr_v
    }
    /// Returns the default header
    fn default_header() -> Vec<String> {
        Self::header(EpochFormat::GregorianUtc, CovarFormat::Sqrt)
    }

    /// Returns the covariance element at position (i, j) formatted with this estimate's covariance formatter
    fn covar_ij(&self, i: usize, j: usize) -> f64 {
        match self.covar_fmt() {
            CovarFormat::Sqrt => self.covar()[(i, j)].sqrt(),
            CovarFormat::Sigma1 => self.covar()[(i, j)],
            CovarFormat::Sigma3 => self.covar()[(i, j)] * 3.0,
            CovarFormat::MulSigma(x) => self.covar()[(i, j)] * x,
        }
    }
}

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
    /// The Epoch format upon serialization
    pub epoch_fmt: EpochFormat,
    /// The covariance format upon serialization
    pub covar_fmt: CovarFormat,
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
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
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
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
        }
    }
}

impl KfEstimate<Orbit> {
    /// Generates an initial Kalman filter state estimate dispersed from the nominal state using the provided standard deviation parameters.
    ///
    /// The resulting estimate will have a diagonal covariance matrix constructed from the variances of each parameter.
    /// *Limitation:* This method incorrectly assumes all parameters are statistically independent.
    pub fn disperse_from_diag(
        nominal_state: Orbit,
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
        let delta = nominal_state - dispersed_state.state;

        // Build the covariance as three times the absolute value of the error, squared.

        let diag = Vector6::new(
            (3.0 * delta.x_km.abs()).powi(2),
            (3.0 * delta.y_km.abs()).powi(2),
            (3.0 * delta.z_km.abs()).powi(2),
            (3.0 * delta.vx_km_s.abs()).powi(2),
            (3.0 * delta.vy_km_s.abs()).powi(2),
            (3.0 * delta.vz_km_s.abs()).powi(2),
        );

        // Build the covar from the diagonal
        let covar = Matrix::from_diagonal(&diag);

        Self {
            nominal_state: dispersed_state.state,
            state_deviation: OVector::<f64, U6>::zeros(),
            covar,
            covar_bar: covar,
            predicted: true,
            stm: OMatrix::<f64, U6, U6>::identity(),
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
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
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
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
    fn epoch_fmt(&self) -> EpochFormat {
        self.epoch_fmt
    }
    fn covar_fmt(&self) -> CovarFormat {
        self.covar_fmt
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

impl<T: State> Serialize for KfEstimate<T>
where
    DefaultAllocator: Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<usize, <T as State>::Size>
        + Allocator<f64, <T as State>::VecLength>
        + Allocator<usize, <T as State>::Size, <T as State>::Size>,
    <DefaultAllocator as Allocator<f64, <T as State>::Size>>::Buffer: Copy,
    <DefaultAllocator as Allocator<f64, <T as State>::Size, <T as State>::Size>>::Buffer: Copy,
{
    /// Serializes the estimate
    fn serialize<O>(&self, serializer: O) -> Result<O::Ok, O::Error>
    where
        O: Serializer,
    {
        let dim = <T as State>::Size::dim();
        let mut seq = serializer.serialize_seq(Some(dim * 3 + 1))?;
        match self.epoch_fmt {
            EpochFormat::GregorianUtc => seq.serialize_element(&format!("{}", self.epoch()))?,
            EpochFormat::GregorianTai => seq.serialize_element(&format!("{}", self.epoch()))?,
            EpochFormat::MjdTai => seq.serialize_element(&self.epoch().to_mjd_tai_days())?,
            EpochFormat::MjdTt => seq.serialize_element(&self.epoch().to_mjd_tt_days())?,
            EpochFormat::MjdUtc => seq.serialize_element(&self.epoch().to_mjd_utc_days())?,
            EpochFormat::JdeEt => seq.serialize_element(&self.epoch().to_jde_et_days())?,
            EpochFormat::JdeTai => seq.serialize_element(&self.epoch().to_jde_tai_days())?,
            EpochFormat::JdeTt => seq.serialize_element(&self.epoch().to_jde_tt_days())?,
            EpochFormat::JdeUtc => seq.serialize_element(&self.epoch().to_jde_utc_days())?,
            EpochFormat::TaiSecs(e) => {
                seq.serialize_element(&(self.epoch().to_tai_seconds() - e))?
            }
            EpochFormat::TaiDays(e) => seq.serialize_element(&(self.epoch().to_tai_days() - e))?,
        }
        // Serialize the state
        for i in 0..dim {
            seq.serialize_element(&self.state_deviation[i])?;
        }
        // Serialize the covariance
        for i in 0..dim {
            for j in 0..dim {
                let ser_covar = match self.covar_fmt {
                    CovarFormat::Sqrt => self.covar[(i, j)].sqrt(),
                    CovarFormat::Sigma1 => self.covar[(i, j)],
                    CovarFormat::Sigma3 => self.covar[(i, j)] * 3.0,
                    CovarFormat::MulSigma(x) => self.covar[(i, j)] * x,
                };
                seq.serialize_element(&ser_covar)?;
            }
        }
        seq.end()
    }
}

/// A trait to store a navigation solution, can be used in conjunction with KfEstimate
pub trait NavSolution<T>: Estimate<Orbit>
where
    T: State,
    DefaultAllocator: Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<f64, <T as State>::VecLength>,
{
    fn orbital_state(&self) -> Orbit;
    /// Returns the nominal state as computed by the dynamics
    fn expected_state(&self) -> Orbit;
}

impl NavSolution<Orbit> for KfEstimate<Orbit> {
    fn orbital_state(&self) -> Orbit {
        self.state()
    }
    fn expected_state(&self) -> Orbit {
        self.nominal_state()
    }
}

impl EstimateFrom<Spacecraft, StdMeasurement> for Orbit {
    fn extract(from: Spacecraft) -> Self {
        from.orbit
    }

    fn sensitivity(
        msr: &StdMeasurement,
        receiver: Self,
        transmitter: Self,
    ) -> OMatrix<f64, <StdMeasurement as Measurement>::MeasurementSize, Self::Size>
    where
        DefaultAllocator:
            Allocator<f64, <StdMeasurement as Measurement>::MeasurementSize, Self::Size>,
    {
        <Orbit as EstimateFrom<Orbit, StdMeasurement>>::sensitivity(msr, receiver, transmitter)
    }
}

impl EstimateFrom<Orbit, StdMeasurement> for Orbit {
    fn extract(from: Orbit) -> Self {
        from
    }

    fn sensitivity(
        msr: &StdMeasurement,
        receiver: Self,
        transmitter: Self,
    ) -> OMatrix<f64, <StdMeasurement as Measurement>::MeasurementSize, Self::Size>
    where
        DefaultAllocator:
            Allocator<f64, <StdMeasurement as Measurement>::MeasurementSize, Self::Size>,
    {
        let delta_r = receiver.radius() - transmitter.radius();
        let delta_v = receiver.velocity() - transmitter.velocity();
        let ρ = msr.observation()[0];
        let ρ_dot = msr.observation()[1];
        let m11 = delta_r.x / ρ;
        let m12 = delta_r.y / ρ;
        let m13 = delta_r.z / ρ;
        let m21 = delta_v.x / ρ - ρ_dot * delta_r.x / ρ.powi(2);
        let m22 = delta_v.y / ρ - ρ_dot * delta_r.y / ρ.powi(2);
        let m23 = delta_v.z / ρ - ρ_dot * delta_r.z / ρ.powi(2);

        Matrix2x6::new(m11, m12, m13, 0.0, 0.0, 0.0, m21, m22, m23, m11, m12, m13)
    }
}

impl EstimateFrom<Spacecraft, StdMeasurement> for Spacecraft {
    fn extract(from: Spacecraft) -> Self {
        from
    }

    fn sensitivity(
        _msr: &StdMeasurement,
        _receiver: Self,
        _transmitter: Orbit,
    ) -> OMatrix<f64, <StdMeasurement as Measurement>::MeasurementSize, Self::Size>
    where
        DefaultAllocator:
            Allocator<f64, <StdMeasurement as Measurement>::MeasurementSize, Self::Size>,
    {
        todo!()
    }
}

#[test]
fn test_estimate_from_disp() {
    use crate::cosmic::Cosm;
    use crate::utils::rss_orbit_errors;

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

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

    let (init_rss_pos_km, init_rss_vel_km_s) = rss_orbit_errors(&initial_state, &initial_state_dev);

    let delta = initial_state - initial_state_dev;

    println!("Truth initial state:\n{initial_state}\n{initial_state:x}");
    println!("Filter initial state:\n{initial_state_dev}\n{initial_state_dev:x}");
    println!("Initial state dev:\t{init_rss_pos_km:.3} km\t{init_rss_vel_km_s:.3} km/s\n{delta}",);
    println!("covariance: {}", initial_estimate.covar);

    // Check that the error is in the square root of the covariance
    assert!(delta.x_km < initial_estimate.covar[(0, 0)].sqrt());
    assert!(delta.y_km < initial_estimate.covar[(1, 1)].sqrt());
    assert!(delta.z_km < initial_estimate.covar[(2, 2)].sqrt());
    assert!(delta.vx_km_s < initial_estimate.covar[(3, 3)].sqrt());
    assert!(delta.vy_km_s < initial_estimate.covar[(4, 4)].sqrt());
    assert!(delta.vz_km_s < initial_estimate.covar[(5, 5)].sqrt());
}
