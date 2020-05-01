use super::serde::ser::SerializeSeq;
use super::serde::{Serialize, Serializer};
use super::{CovarFormat, EpochFormat};
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, DimName, MatrixMN, VectorN};
use crate::hifitime::Epoch;
use std::f64::INFINITY;
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
    fn state(&self) -> VectorN<f64, S>;
    /// The Covariance of this estimate
    fn covar(&self) -> MatrixMN<f64, S, S>;
    // Sets the epoch
    fn set_dt(&mut self, dt: Epoch);
    /// Sets the estimated state, or state deviation (check filter docs).
    fn set_state(&mut self, new_state: VectorN<f64, S>);
    /// Sets the Covariance of this estimate
    fn set_covar(&mut self, new_covar: MatrixMN<f64, S, S>);
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    fn predicted(&self) -> bool;
    /// The STM used to compute this Estimate
    fn stm(&self) -> MatrixMN<f64, S, S>;
    /// The Epoch format upon serialization
    fn epoch_fmt(&self) -> EpochFormat;
    /// The covariance format upon serialization
    fn covar_fmt(&self) -> CovarFormat;
    /// Returns whether this estimate is within some bound
    /// The 68-95-99.7 rule is a good way to assess whether the filter is operating normally
    fn within_sigma(&self, sigma: f64) -> bool {
        let state = self.state();
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

/// Kalman filter Estimate
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
            epoch_fmt: EpochFormat::GregorianUtc,
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
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
        }
    }

    fn dt(&self) -> Epoch {
        self.dt
    }

    fn state(&self) -> VectorN<f64, S> {
        self.state.clone()
    }

    fn covar(&self) -> MatrixMN<f64, S, S> {
        self.covar.clone()
    }

    fn predicted(&self) -> bool {
        self.predicted
    }
    fn stm(&self) -> MatrixMN<f64, S, S> {
        self.stm.clone()
    }
    fn epoch_fmt(&self) -> EpochFormat {
        self.epoch_fmt
    }
    fn covar_fmt(&self) -> CovarFormat {
        self.covar_fmt
    }
    fn set_dt(&mut self, dt: Epoch) {
        self.dt = dt;
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
        let word = if self.predicted {
            "Prediction"
        } else {
            "Estimate"
        };
        write!(
            f,
            "=== {} @ {} UTC -- within 3 sigma: {} ===\nEstState {} Covariance {}\n=====================",
            word,
            &self.dt.as_gregorian_utc_str(),
            self.within_3sigma(),
            &self.state,
            &self.covar
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
            EpochFormat::GregorianTai => seq.serialize_element(&self.dt.as_gregorian_tai_str())?,
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

/// Information filter Estimate
#[derive(Debug, Clone, PartialEq)]
pub struct IfEstimate<S>
where
    S: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    /// Date time of this Estimate
    pub dt: Epoch,
    /// The information state
    pub info_state: VectorN<f64, S>,
    /// The information matrix, which is the inverse of the covariance
    pub info_mat: MatrixMN<f64, S, S>,
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    pub predicted: bool,
    /// The STM used to compute this Estimate
    pub stm: MatrixMN<f64, S, S>,
    /// The Epoch format upon serialization
    pub epoch_fmt: EpochFormat,
    /// The covariance format upon serialization
    pub covar_fmt: CovarFormat,
}

impl<S> IfEstimate<S>
where
    S: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    pub fn from_covar(dt: Epoch, covar: MatrixMN<f64, S, S>) -> Self {
        let mut info_mat = covar;
        if !info_mat.try_inverse_mut() {
            panic!("provided covariance is singular");
        }

        Self {
            dt,
            info_state: VectorN::<f64, S>::zeros(),
            info_mat,
            predicted: true,
            stm: MatrixMN::<f64, S, S>::zeros(),
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
        }
    }

    /// Returns the covariance, if there is enough information to invert the information matrix
    pub fn try_covar(&self) -> Option<MatrixMN<f64, S, S>> {
        let mut covar = self.info_mat.clone();
        if !covar.try_inverse_mut() {
            None
        } else {
            Some(&covar * &covar.transpose())
        }
    }
}

impl<S> Estimate<S> for IfEstimate<S>
where
    S: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    fn zeros() -> Self {
        let mut info_state = VectorN::<f64, S>::zeros();
        let mut info_mat = MatrixMN::<f64, S, S>::zeros();
        // Initialize everything to infinity
        for i in 0..S::dim() {
            info_state[i] = INFINITY;
            info_mat[(i, i)] = INFINITY;
        }
        Self {
            dt: Epoch::from_tai_seconds(0.0),
            info_state,
            info_mat,
            predicted: true,
            stm: MatrixMN::<f64, S, S>::zeros(),
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
        }
    }

    fn dt(&self) -> Epoch {
        self.dt
    }

    /// Will panic if the information matrix inversion fails
    fn state(&self) -> VectorN<f64, S> {
        &self.covar() * &self.info_state
    }

    /// Will panic if the information matrix inversion fails
    fn covar(&self) -> MatrixMN<f64, S, S> {
        self.try_covar().unwrap()
    }

    fn predicted(&self) -> bool {
        self.predicted
    }
    fn stm(&self) -> MatrixMN<f64, S, S> {
        self.stm.clone()
    }
    fn epoch_fmt(&self) -> EpochFormat {
        self.epoch_fmt
    }
    fn covar_fmt(&self) -> CovarFormat {
        self.covar_fmt
    }
    fn set_dt(&mut self, dt: Epoch) {
        self.dt = dt;
    }
    /// WARNING: This sets the information state, not the filter state
    fn set_state(&mut self, new_info_state: VectorN<f64, S>) {
        self.info_state = new_info_state;
    }
    /// WARNING: This sets the information matrix
    fn set_covar(&mut self, new_info_mat: MatrixMN<f64, S, S>) {
        self.info_mat = new_info_mat;
    }
}

impl<S> fmt::Display for IfEstimate<S>
where
    S: DimName,
    DefaultAllocator:
        Allocator<f64, S> + Allocator<f64, S, S> + Allocator<usize, S> + Allocator<usize, S, S>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.try_covar() {
            Some(covar) => write!(
                f,
                "=== PREDICTED: {} ===\nEstState {} Covariance {}\n=====================",
                &self.predicted,
                self.state(),
                covar
            ),
            None => write!(f, "=== PREDICTED: {} === Not invertible", &self.predicted),
        }
    }
}

impl<S> fmt::LowerExp for IfEstimate<S>
where
    S: DimName,
    DefaultAllocator:
        Allocator<f64, S> + Allocator<f64, S, S> + Allocator<usize, S> + Allocator<usize, S, S>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.try_covar() {
            Some(covar) => write!(
                f,
                "=== PREDICTED: {} ===\nEstState {:e} Covariance {:e}\n=====================",
                &self.predicted,
                self.state(),
                covar
            ),
            None => write!(f, "=== PREDICTED: {} === Not invertible", &self.predicted),
        }
    }
}

impl<S> Serialize for IfEstimate<S>
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
            EpochFormat::GregorianTai => seq.serialize_element(&self.dt.as_gregorian_tai_str())?,
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
        match self.try_covar() {
            Some(covar) => {
                let state = self.state();
                // Serialize the state
                for i in 0..S::dim() {
                    seq.serialize_element(&state[(i, 0)])?;
                }
                // Serialize the covariance
                for i in 0..S::dim() {
                    for j in 0..S::dim() {
                        let ser_covar = match self.covar_fmt {
                            CovarFormat::Sqrt => covar[(i, j)].sqrt(),
                            CovarFormat::Sigma1 => covar[(i, j)],
                            CovarFormat::Sigma3 => covar[(i, j)] * 3.0,
                            CovarFormat::MulSigma(x) => covar[(i, j)] * x,
                        };
                        seq.serialize_element(&ser_covar)?;
                    }
                }
            }
            None => {
                // Set all of the numbers to 1e32
                for _ in 0..S::dim() {
                    seq.serialize_element(&1e32)?;
                }
                // Serialize the covariance
                for _ in 0..S::dim() {
                    for _ in 0..S::dim() {
                        seq.serialize_element(&1e32)?;
                    }
                }
            }
        }
        seq.end()
    }
}
