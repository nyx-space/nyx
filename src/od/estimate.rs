/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::serde::ser::SerializeSeq;
use super::serde::{Serialize, Serializer};
use super::{CovarFormat, EpochFormat};
use super::{EstimateFrom, State};
use crate::cosmic::Orbit;
use crate::hifitime::Epoch;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector};
use crate::Spacecraft;
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
        hdr_v.push(format!("{}", epoch_fmt));
        for i in 0..dim {
            hdr_v.push(format!("state_{}", i));
        }
        // Serialize the covariance
        for i in 0..dim {
            for j in 0..dim {
                hdr_v.push(format!("{}_{}_{}", covar_fmt, i, j));
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
#[derive(Debug, Clone, PartialEq)]
pub struct KfEstimate<T: State>
where
    DefaultAllocator: Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<usize, <T as State>::Size>
        + Allocator<f64, <T as State>::VecLength>
        + Allocator<usize, <T as State>::Size, <T as State>::Size>,
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
{
    pub fn from_covar(
        nominal_state: T,
        covar: OMatrix<f64, <T as State>::Size, <T as State>::Size>,
    ) -> Self {
        Self {
            nominal_state,
            state_deviation: OVector::<f64, <T as State>::Size>::zeros(),
            covar: covar.clone(),
            covar_bar: covar,
            predicted: true,
            stm: OMatrix::<f64, <T as State>::Size, <T as State>::Size>::identity(),
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
        self.state_deviation.clone()
    }

    fn covar(&self) -> OMatrix<f64, <T as State>::Size, <T as State>::Size> {
        self.covar.clone()
    }

    fn predicted_covar(&self) -> OMatrix<f64, <T as State>::Size, <T as State>::Size> {
        self.covar_bar.clone()
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
            &self.epoch().as_gregorian_utc_str(),
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
{
    /// Serializes the estimate
    fn serialize<O>(&self, serializer: O) -> Result<O::Ok, O::Error>
    where
        O: Serializer,
    {
        let dim = <T as State>::Size::dim();
        let mut seq = serializer.serialize_seq(Some(dim * 3 + 1))?;
        match self.epoch_fmt {
            EpochFormat::GregorianUtc => {
                seq.serialize_element(&self.epoch().as_gregorian_utc_str())?
            }
            EpochFormat::GregorianTai => {
                seq.serialize_element(&self.epoch().as_gregorian_tai_str())?
            }
            EpochFormat::MjdTai => seq.serialize_element(&self.epoch().as_mjd_tai_days())?,
            EpochFormat::MjdTt => seq.serialize_element(&self.epoch().as_mjd_tt_days())?,
            EpochFormat::MjdUtc => seq.serialize_element(&self.epoch().as_mjd_utc_days())?,
            EpochFormat::JdeEt => seq.serialize_element(&self.epoch().as_jde_et_days())?,
            EpochFormat::JdeTai => seq.serialize_element(&self.epoch().as_jde_tai_days())?,
            EpochFormat::JdeTt => seq.serialize_element(&self.epoch().as_jde_tt_days())?,
            EpochFormat::JdeUtc => seq.serialize_element(&self.epoch().as_jde_utc_days())?,
            EpochFormat::TaiSecs(e) => {
                seq.serialize_element(&(self.epoch().as_tai_seconds() - e))?
            }
            EpochFormat::TaiDays(e) => seq.serialize_element(&(self.epoch().as_tai_days() - e))?,
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

impl EstimateFrom<Spacecraft> for Orbit {
    fn extract(from: Spacecraft) -> Self {
        from.orbit
    }
}

impl EstimateFrom<Orbit> for Orbit {
    fn extract(from: Orbit) -> Self {
        from
    }
}

impl EstimateFrom<Spacecraft> for Spacecraft {
    fn extract(from: Spacecraft) -> Self {
        from
    }
}
