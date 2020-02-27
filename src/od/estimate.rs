extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, DimName, MatrixMN, VectorN};
use super::serde::ser::SerializeSeq;
use super::serde::{Serialize, Serializer};
use super::{CovarFormat, EpochFormat};
use crate::hifitime::Epoch;
use std::fmt;

/// Stores an Estimate, as the result of a `time_update` or `measurement_update`.
pub trait Estimate<S>
where
    Self: Clone + PartialEq + Sized,
    S: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    /// An empty estimate. This is useful if wanting to store an estimate outside the scope of a filtering loop.
    fn zeros() -> Self;
    /// Date time of this Estimate
    fn dt(&self) -> Epoch;
    /// The estimated state, or state deviation (check filter docs).
    fn state(&self) -> &VectorN<f64, S>;
    /// The Covariance of this estimate
    fn covar(&self) -> &MatrixMN<f64, S, S>;
    /// The estimated state, or state deviation (check filter docs).
    fn set_state(&mut self, new_state: VectorN<f64, S>);
    /// The Covariance of this estimate
    fn set_covar(&mut self, new_covar: MatrixMN<f64, S, S>);
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    fn predicted(&self) -> bool;
    /// The STM used to compute this Estimate
    fn stm(&self) -> &MatrixMN<f64, S, S>;
    /// The Epoch format upon serialization
    fn epoch_fmt(&self) -> EpochFormat;
    /// The covariance format upon serialization
    fn covar_fmt(&self) -> CovarFormat;
    /// Returns the header
    fn header(epoch_fmt: EpochFormat, covar_fmt: CovarFormat) -> Vec<String> {
        let mut hdr_v = Vec::with_capacity(3 * S::dim() + 1);
        hdr_v.push(format!("{}", epoch_fmt));
        for i in 0..S::dim() {
            hdr_v.push(format!("state_{}", i));
        }
        // Serialize the covariance
        for i in 0..S::dim() {
            for j in 0..S::dim() {
                hdr_v.push(format!("{}_{}_{}", covar_fmt, i, j));
            }
        }
        hdr_v
    }
    /// Returns the default header
    fn default_header() -> Vec<String> {
        Self::header(EpochFormat::GregorianUtc, CovarFormat::Sqrt)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct KfEstimate<S>
where
    S: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    /// Date time of this Estimate
    pub dt: Epoch,
    /// The estimated state, or state deviation (check filter docs).
    pub state: VectorN<f64, S>,
    /// The Covariance of this estimate
    pub covar: MatrixMN<f64, S, S>,
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    pub predicted: bool,
    /// The STM used to compute this Estimate
    pub stm: MatrixMN<f64, S, S>,
    /// The Epoch format upon serialization
    pub epoch_fmt: EpochFormat,
    /// The covariance format upon serialization
    pub covar_fmt: CovarFormat,
}

impl<S> KfEstimate<S>
where
    S: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    pub fn from_covar(dt: Epoch, covar: MatrixMN<f64, S, S>) -> Self {
        Self {
            dt,
            state: VectorN::<f64, S>::zeros(),
            covar,
            predicted: true,
            stm: MatrixMN::<f64, S, S>::zeros(),
            epoch_fmt: EpochFormat::MjdTai,
            covar_fmt: CovarFormat::Sqrt,
        }
    }
}

impl<S> Estimate<S> for KfEstimate<S>
where
    S: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    fn zeros() -> Self {
        Self {
            dt: Epoch::from_tai_seconds(0.0),
            state: VectorN::<f64, S>::zeros(),
            covar: MatrixMN::<f64, S, S>::zeros(),
            predicted: true,
            stm: MatrixMN::<f64, S, S>::zeros(),
            epoch_fmt: EpochFormat::MjdTai,
            covar_fmt: CovarFormat::Sqrt,
        }
    }

    fn dt(&self) -> Epoch {
        self.dt
    }

    fn state(&self) -> &VectorN<f64, S> {
        &self.state
    }

    fn covar(&self) -> &MatrixMN<f64, S, S> {
        &self.covar
    }

    fn predicted(&self) -> bool {
        self.predicted
    }
    fn stm(&self) -> &MatrixMN<f64, S, S> {
        &self.stm
    }
    fn epoch_fmt(&self) -> EpochFormat {
        self.epoch_fmt
    }
    fn covar_fmt(&self) -> CovarFormat {
        self.covar_fmt
    }
    fn set_state(&mut self, new_state: VectorN<f64, S>) {
        self.state = new_state;
    }
    fn set_covar(&mut self, new_covar: MatrixMN<f64, S, S>) {
        self.covar = new_covar;
    }
}

impl<S> fmt::Display for KfEstimate<S>
where
    S: DimName,
    DefaultAllocator:
        Allocator<f64, S> + Allocator<f64, S, S> + Allocator<usize, S> + Allocator<usize, S, S>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "=== PREDICTED: {} ===\nEstState {} Covariance {}\n=====================",
            &self.predicted, &self.state, &self.covar
        )
    }
}

impl<S> fmt::LowerExp for KfEstimate<S>
where
    S: DimName,
    DefaultAllocator:
        Allocator<f64, S> + Allocator<f64, S, S> + Allocator<usize, S> + Allocator<usize, S, S>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "=== PREDICTED: {} ===\nEstState {:e} Covariance {:e}\n=====================",
            &self.predicted, &self.state, &self.covar
        )
    }
}

impl<S> Serialize for KfEstimate<S>
where
    S: DimName,
    DefaultAllocator:
        Allocator<f64, S> + Allocator<f64, S, S> + Allocator<usize, S> + Allocator<usize, S, S>,
{
    /// Serializes the estimate
    fn serialize<O>(&self, serializer: O) -> Result<O::Ok, O::Error>
    where
        O: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(S::dim() * 3 + 1))?;
        match self.epoch_fmt {
            EpochFormat::GregorianUtc => seq.serialize_element(&self.dt.as_gregorian_utc_str())?,
            EpochFormat::GregorianTai => seq.serialize_element(&self.dt.as_gregorian_utc_tai())?,
            EpochFormat::MjdTai => seq.serialize_element(&self.dt.as_mjd_tai_days())?,
            EpochFormat::MjdTt => seq.serialize_element(&self.dt.as_mjd_tt_days())?,
            EpochFormat::MjdUtc => seq.serialize_element(&self.dt.as_mjd_utc_days())?,
            EpochFormat::JdeEt => seq.serialize_element(&self.dt.as_jde_et_days())?,
            EpochFormat::JdeTai => seq.serialize_element(&self.dt.as_jde_tai_days())?,
            EpochFormat::JdeTt => seq.serialize_element(&self.dt.as_jde_tt_days())?,
            EpochFormat::JdeUtc => seq.serialize_element(&self.dt.as_jde_utc_days())?,
            EpochFormat::TaiSecs(e) => seq.serialize_element(&(self.dt.as_tai_seconds() - e))?,
            EpochFormat::TaiDays(e) => seq.serialize_element(&(self.dt.as_tai_days() - e))?,
        }
        // Serialize the state
        for i in 0..S::dim() {
            seq.serialize_element(&self.state[(i, 0)])?;
        }
        // Serialize the covariance
        for i in 0..S::dim() {
            for j in 0..S::dim() {
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
