use super::serde::ser::SerializeSeq;
use super::serde::{Serialize, Serializer};
use super::State;
use super::{CovarFormat, EpochFormat};
use crate::celestia::Orbit;
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, DimName, MatrixMN, VectorN, U42, U6};
use crate::hifitime::Epoch;
use crate::SpacecraftState;
use std::cmp::PartialEq;
use std::f64::INFINITY;
use std::fmt;

/// Stores an Estimate, as the result of a `time_update` or `measurement_update`.
pub trait Estimate<S, P, T: State<S, P>>
where
    Self: Clone + PartialEq + Sized,
    S: DimName,
    P: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, P> + Allocator<f64, S, S>,
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
        self.nominal_state() + self.state_deviation()
        // let mut state = self.nominal_state();
        // state.set(
        //     state.epoch(),
        //     &(state.as_vector().unwrap() + self.state_deviation()),
        // );
        // state
    }
    /// The state deviation as computed by the filter.
    fn state_deviation(&self) -> VectorN<f64, S>;
    /// The nominal state as reported by the filter dynamics
    fn nominal_state(&self) -> T;
    /// The Covariance of this estimate. Will return the predicted covariance if this is a time update/prediction.
    fn covar(&self) -> MatrixMN<f64, S, S>;
    /// The predicted covariance of this estimate from the time update
    fn predicted_covar(&self) -> MatrixMN<f64, S, S>;
    /// Sets the state deviation.
    fn set_state_deviation(&mut self, new_state: VectorN<f64, S>);
    /// Sets the Covariance of this estimate
    fn set_covar(&mut self, new_covar: MatrixMN<f64, S, S>);
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    fn predicted(&self) -> bool;
    /// The STM used to compute this Estimate
    fn stm(&self) -> &MatrixMN<f64, S, S>;
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
pub struct KfEstimate<S, P, T: State<S, P>>
where
    S: DimName,
    P: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    /// The estimated state
    pub nominal_state: T,
    /// The state deviation
    pub state_deviation: VectorN<f64, S>,
    /// The Covariance of this estimate
    pub covar: MatrixMN<f64, S, S>,
    /// The predicted covariance of this estimate
    pub covar_bar: MatrixMN<f64, S, S>,
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    pub predicted: bool,
    /// The STM used to compute this Estimate
    pub stm: MatrixMN<f64, S, S>,
    /// The Epoch format upon serialization
    pub epoch_fmt: EpochFormat,
    /// The covariance format upon serialization
    pub covar_fmt: CovarFormat,
}

impl<S, P, T: State<S, P>> KfEstimate<S, P, T>
where
    S: DimName,
    P: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S> + Allocator<f64, P>,
{
    pub fn from_covar(nominal_state: T, covar: MatrixMN<f64, S, S>) -> Self {
        Self {
            nominal_state,
            state_deviation: VectorN::<f64, S>::zeros(),
            covar: covar.clone(),
            covar_bar: covar,
            predicted: true,
            stm: MatrixMN::<f64, S, S>::zeros(),
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
        }
    }
}

impl<S, P, T: State<S, P>> Estimate<S, P, T> for KfEstimate<S, P, T>
where
    S: DimName,
    P: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    fn zeros(nominal_state: T) -> Self {
        Self {
            nominal_state,
            state_deviation: VectorN::<f64, S>::zeros(),
            covar: MatrixMN::<f64, S, S>::zeros(),
            covar_bar: MatrixMN::<f64, S, S>::zeros(),
            predicted: true,
            stm: MatrixMN::<f64, S, S>::zeros(),
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
        }
    }

    fn nominal_state(&self) -> T {
        self.nominal_state.clone()
    }

    fn state_deviation(&self) -> VectorN<f64, S> {
        self.state_deviation.clone()
    }

    fn covar(&self) -> MatrixMN<f64, S, S> {
        self.covar.clone()
    }

    fn predicted_covar(&self) -> MatrixMN<f64, S, S> {
        self.covar_bar.clone()
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
    fn set_state_deviation(&mut self, new_state: VectorN<f64, S>) {
        self.state_deviation = new_state;
    }
    fn set_covar(&mut self, new_covar: MatrixMN<f64, S, S>) {
        self.covar = new_covar;
    }
}

impl<S, P, T: State<S, P>> fmt::Display for KfEstimate<S, P, T>
where
    S: DimName,
    P: DimName,
    DefaultAllocator: Allocator<f64, S>
        + Allocator<f64, P>
        + Allocator<f64, S, S>
        + Allocator<usize, S>
        + Allocator<usize, S, S>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let word = if self.predicted {
            "Prediction"
        } else {
            "Estimate"
        };
        let mut fmt_cov = Vec::with_capacity(S::dim());
        for i in 0..S::dim() {
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

impl<S, P, T: State<S, P>> fmt::LowerExp for KfEstimate<S, P, T>
where
    S: DimName,
    P: DimName,
    DefaultAllocator: Allocator<f64, S>
        + Allocator<f64, P>
        + Allocator<f64, S, S>
        + Allocator<usize, S>
        + Allocator<usize, S, S>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "=== PREDICTED: {} ===\nEstState {:e} Covariance {:e}\n=====================",
            &self.predicted, &self.state_deviation, &self.covar
        )
    }
}

impl<S, P, T: State<S, P>> Serialize for KfEstimate<S, P, T>
where
    S: DimName,
    P: DimName,
    DefaultAllocator: Allocator<f64, S>
        + Allocator<f64, P>
        + Allocator<f64, S, S>
        + Allocator<usize, S>
        + Allocator<usize, S, S>,
{
    /// Serializes the estimate
    fn serialize<O>(&self, serializer: O) -> Result<O::Ok, O::Error>
    where
        O: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(S::dim() * 3 + 1))?;
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
        for i in 0..S::dim() {
            seq.serialize_element(&self.state_deviation[i])?;
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
pub struct IfEstimate<S, P, T: State<S, P>>
where
    S: DimName,
    P: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, P> + Allocator<f64, S, S>,
{
    /// The nominal state
    pub nominal_state: T,
    /// The information state
    pub info_state: VectorN<f64, S>,
    /// The information matrix, which is the inverse of the covariance
    pub info_mat: MatrixMN<f64, S, S>,
    /// The predicted information matrix, which is the inverse of the covariance
    pub info_mat_bar: MatrixMN<f64, S, S>,
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    pub predicted: bool,
    /// The STM used to compute this Estimate
    pub stm: MatrixMN<f64, S, S>,
    /// The Epoch format upon serialization
    pub epoch_fmt: EpochFormat,
    /// The covariance format upon serialization
    pub covar_fmt: CovarFormat,
}

impl<S, P, T: State<S, P>> IfEstimate<S, P, T>
where
    S: DimName,
    P: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, P> + Allocator<f64, S, S>,
{
    pub fn from_covar(nominal_state: T, covar: MatrixMN<f64, S, S>) -> Self {
        let mut info_mat = covar;
        if !info_mat.try_inverse_mut() {
            panic!("provided covariance is singular");
        }

        Self {
            nominal_state,
            info_state: VectorN::<f64, S>::zeros(),
            info_mat: info_mat.clone(),
            info_mat_bar: info_mat,
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

    /// Returns the covariance, if there is enough information to invert the information matrix
    pub fn try_predicted_covar(&self) -> Option<MatrixMN<f64, S, S>> {
        let mut covar = self.info_mat_bar.clone();
        if !covar.try_inverse_mut() {
            None
        } else {
            Some(&covar * &covar.transpose())
        }
    }
}

impl<S, P, T: State<S, P>> Estimate<S, P, T> for IfEstimate<S, P, T>
where
    S: DimName,
    P: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, P> + Allocator<f64, S, S>,
{
    fn zeros(nominal_state: T) -> Self {
        let mut info_state = VectorN::<f64, S>::zeros();
        let mut info_mat = MatrixMN::<f64, S, S>::zeros();
        // Initialize everything to infinity
        for i in 0..S::dim() {
            info_state[i] = INFINITY;
            info_mat[(i, i)] = INFINITY;
        }
        Self {
            nominal_state,
            info_state,
            info_mat: info_mat.clone(),
            info_mat_bar: info_mat,
            predicted: true,
            stm: MatrixMN::<f64, S, S>::zeros(),
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
        }
    }

    fn nominal_state(&self) -> T {
        self.nominal_state.clone()
    }

    /// Will panic if the information matrix inversion fails
    fn state_deviation(&self) -> VectorN<f64, S> {
        &self.covar() * &self.info_state
    }

    /// Will panic if the information matrix inversion fails
    fn covar(&self) -> MatrixMN<f64, S, S> {
        self.try_covar().unwrap()
    }

    /// Will panic if the information matrix inversion fails
    fn predicted_covar(&self) -> MatrixMN<f64, S, S> {
        self.try_predicted_covar().unwrap()
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
    /// WARNING: This sets the information state, not the filter state
    fn set_state_deviation(&mut self, new_info_state: VectorN<f64, S>) {
        self.info_state = new_info_state;
    }
    /// WARNING: This sets the information matrix
    fn set_covar(&mut self, new_info_mat: MatrixMN<f64, S, S>) {
        self.info_mat = new_info_mat;
    }
}

impl<S, P, T: State<S, P>> fmt::Display for IfEstimate<S, P, T>
where
    S: DimName,
    P: DimName,
    DefaultAllocator: Allocator<f64, S>
        + Allocator<f64, P>
        + Allocator<f64, S, S>
        + Allocator<usize, S>
        + Allocator<usize, S, S>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.try_covar() {
            Some(covar) => {
                let mut fmt_cov = Vec::with_capacity(S::dim());
                for i in 0..S::dim() {
                    fmt_cov.push(format!("{:e}", covar[(i, i)]));
                }
                write!(
                    f,
                    "=== ESTIMATE @ {} -- within 3 sigma: {} ===\nstate {}\nsigmas [{}]\n",
                    &self.epoch().as_gregorian_utc_str(),
                    self.within_3sigma(),
                    &self.state(),
                    fmt_cov.join(",")
                )
            }
            None => write!(
                f,
                "=== PREDICTION @ {} === Not invertible",
                &self.epoch().as_gregorian_utc_str(),
            ),
        }
    }
}

impl<S, P, T: State<S, P>> fmt::LowerExp for IfEstimate<S, P, T>
where
    S: DimName,
    P: DimName,
    DefaultAllocator: Allocator<f64, S>
        + Allocator<f64, P>
        + Allocator<f64, S, S>
        + Allocator<usize, S>
        + Allocator<usize, S, S>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.try_covar() {
            Some(covar) => write!(
                f,
                "=== PREDICTED: {} ===\nEstState {:e} Covariance {:e}\n=====================",
                &self.predicted,
                self.state_deviation(),
                covar
            ),
            None => write!(f, "=== PREDICTED: {} === Not invertible", &self.predicted),
        }
    }
}

impl<S, P, T: State<S, P>> Serialize for IfEstimate<S, P, T>
where
    S: DimName,
    P: DimName,
    DefaultAllocator: Allocator<f64, S>
        + Allocator<f64, P>
        + Allocator<f64, S, S>
        + Allocator<usize, S>
        + Allocator<usize, S, S>,
{
    /// Serializes the estimate
    fn serialize<O>(&self, serializer: O) -> Result<O::Ok, O::Error>
    where
        O: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(S::dim() * 3 + 1))?;
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
        match self.try_covar() {
            Some(covar) => {
                let state = self.state_deviation();
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

/// A trait to store a navigation solution, can be used in conjunction with KfEstimate or IfEstimate
pub trait NavSolution<T>: Estimate<U6, U42, T>
where
    T: State<U6, U42>,
{
    fn orbital_state(&self) -> Orbit;
    /// Returns the nominal state as computed by the dynamics
    fn expected_state(&self) -> Orbit;
}

impl NavSolution<Orbit> for KfEstimate<U6, U42, Orbit> {
    fn orbital_state(&self) -> Orbit {
        self.state()
    }
    fn expected_state(&self) -> Orbit {
        self.nominal_state()
    }
}

impl NavSolution<Orbit> for IfEstimate<U6, U42, Orbit> {
    fn orbital_state(&self) -> Orbit {
        self.state()
    }
    fn expected_state(&self) -> Orbit {
        self.nominal_state()
    }
}

impl NavSolution<SpacecraftState> for KfEstimate<U6, U42, SpacecraftState> {
    fn orbital_state(&self) -> Orbit {
        self.state().orbit
    }
    fn expected_state(&self) -> Orbit {
        self.nominal_state().orbit
    }
}

impl NavSolution<SpacecraftState> for IfEstimate<U6, U42, SpacecraftState> {
    fn orbital_state(&self) -> Orbit {
        self.state().orbit
    }
    fn expected_state(&self) -> Orbit {
        self.nominal_state().orbit
    }
}
