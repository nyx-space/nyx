extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, DimName, MatrixMN, VectorN};
use super::serde::ser::SerializeSeq;
use super::serde::{Serialize, Serializer};
use crate::hifitime::Epoch;
use std::fmt;

/// Specifies the format of the Epoch during serialization
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EpochFormat {
    /// Default is MJD TAI, as defined in [hifitime](https://docs.rs/hifitime/).
    MjdTai,
    MjdTt,
    MjdUtc,
    JdeEt,
    JdeTai,
    JdeTt,
    JdeUtc,
    /// Seconds past TAI Epoch
    TaiSecs,
    /// Days past TAI Epoch
    TaiDays,
}

impl fmt::Display for EpochFormat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            EpochFormat::MjdTai => write!(f, "MJD TAI"),
            EpochFormat::MjdTt => write!(f, "MJD TT"),
            EpochFormat::MjdUtc => write!(f, "MJD UTC"),
            EpochFormat::JdeEt => write!(f, "JDE ET"),
            EpochFormat::JdeTai => write!(f, "JDE TAI"),
            EpochFormat::JdeTt => write!(f, "JDE TT"),
            EpochFormat::JdeUtc => write!(f, "JDE UTC"),
            EpochFormat::TaiSecs => write!(f, "TAI+ s"),
            EpochFormat::TaiDays => write!(f, "TAI+ days"),
        }
    }
}

/// Specifies the format of the covariance during serialization
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CovarFormat {
    /// Default: allows plotting the variance of the elements instead of the covariance
    Sqrt,
    /// Keeps the covariance as computed, i.e. one sigma (~68%), causes e.g. positional elements in km^2.
    Sigma1,
    /// Three sigma covers about 99.7% of the distribution
    Sigma3,
    /// Allows specifying a custom multiplication factor of each element of the covariance.
    MulSigma(f64),
}

impl fmt::Display for CovarFormat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            CovarFormat::Sqrt => write!(f, "exptd_val_"),
            CovarFormat::Sigma1 => write!(f, "covar_"),
            CovarFormat::Sigma3 => write!(f, "3sig_covar"),
            CovarFormat::MulSigma(x) => write!(f, "{}sig_covar", x),
        }
    }
}

/// Stores an Estimate, as the result of a `time_update` or `measurement_update`.
#[derive(Debug, Clone, PartialEq)]
pub struct Estimate<S>
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

impl<S> Estimate<S>
where
    S: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    /// An empty estimate. This is useful if wanting to store an estimate outside the scope of a filtering loop.
    pub fn empty() -> Estimate<S> {
        Estimate {
            dt: Epoch::from_tai_seconds(0.0),
            state: VectorN::<f64, S>::zeros(),
            covar: MatrixMN::<f64, S, S>::zeros(),
            predicted: true,
            stm: MatrixMN::<f64, S, S>::zeros(),
            epoch_fmt: EpochFormat::MjdTai,
            covar_fmt: CovarFormat::Sqrt,
        }
    }

    pub fn header(epoch_fmt: EpochFormat, covar_fmt: CovarFormat) -> Vec<String> {
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

    pub fn from_covar(dt: Epoch, covar: MatrixMN<f64, S, S>) -> Estimate<S> {
        Estimate {
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

impl<S> fmt::Display for Estimate<S>
where
    S: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S> + Allocator<usize, S> + Allocator<usize, S, S>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "=== PREDICTED: {} ===\nEstState {} Covariance {}\n=====================",
            &self.predicted, &self.state, &self.covar
        )
    }
}

impl<S> Serialize for Estimate<S>
where
    S: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S> + Allocator<usize, S> + Allocator<usize, S, S>,
{
    /// Serializes the estimate
    fn serialize<O>(&self, serializer: O) -> Result<O::Ok, O::Error>
    where
        O: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(S::dim() * 3 + 1))?;
        match self.epoch_fmt {
            EpochFormat::MjdTai => seq.serialize_element(&self.dt.as_mjd_tai_days())?,
            EpochFormat::MjdTt => seq.serialize_element(&self.dt.as_mjd_tt_days())?,
            EpochFormat::MjdUtc => seq.serialize_element(&self.dt.as_mjd_utc_days())?,
            EpochFormat::JdeEt => seq.serialize_element(&self.dt.as_jde_et_days())?,
            EpochFormat::JdeTai => seq.serialize_element(&self.dt.as_jde_tai_days())?,
            EpochFormat::JdeTt => seq.serialize_element(&self.dt.as_jde_tt_days())?,
            EpochFormat::JdeUtc => seq.serialize_element(&self.dt.as_jde_utc_days())?,
            EpochFormat::TaiSecs => seq.serialize_element(&self.dt.as_tai_seconds())?,
            EpochFormat::TaiDays => seq.serialize_element(&self.dt.as_tai_days())?,
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
