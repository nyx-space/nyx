use super::serde::ser::SerializeSeq;
use super::serde::{Serialize, Serializer};
use super::{CovarFormat, EpochFormat};
use super::{EstimateFrom, State};
use crate::celestia::Orbit;
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, DimName, MatrixMN, VectorN};
use crate::hifitime::Epoch;
use crate::SpacecraftState;
use std::cmp::PartialEq;
use std::f64::INFINITY;
use std::fmt;

/// Stores an Estimate, as the result of a `time_update` or `measurement_update`.
pub trait Estimate<T: State>
where
    Self: Clone + PartialEq + Sized,
    DefaultAllocator:
        Allocator<f64, <T as State>::Size> + Allocator<f64, <T as State>::Size, <T as State>::Size>,
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
    fn state_deviation(&self) -> VectorN<f64, <T as State>::Size>;
    /// The nominal state as reported by the filter dynamics
    fn nominal_state(&self) -> T;
    /// The Covariance of this estimate. Will return the predicted covariance if this is a time update/prediction.
    fn covar(&self) -> MatrixMN<f64, <T as State>::Size, <T as State>::Size>;
    /// The predicted covariance of this estimate from the time update
    fn predicted_covar(&self) -> MatrixMN<f64, <T as State>::Size, <T as State>::Size>;
    /// Sets the state deviation.
    fn set_state_deviation(&mut self, new_state: VectorN<f64, <T as State>::Size>);
    /// Sets the Covariance of this estimate
    fn set_covar(&mut self, new_covar: MatrixMN<f64, <T as State>::Size, <T as State>::Size>);
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    fn predicted(&self) -> bool;
    /// The STM used to compute this Estimate
    fn stm(&self) -> &MatrixMN<f64, <T as State>::Size, <T as State>::Size>;
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
    DefaultAllocator:
        Allocator<f64, <T as State>::Size> + Allocator<f64, <T as State>::Size, <T as State>::Size>,
{
    /// The estimated state
    pub nominal_state: T,
    /// The state deviation
    pub state_deviation: VectorN<f64, <T as State>::Size>,
    /// The Covariance of this estimate
    pub covar: MatrixMN<f64, <T as State>::Size, <T as State>::Size>,
    /// The predicted covariance of this estimate
    pub covar_bar: MatrixMN<f64, <T as State>::Size, <T as State>::Size>,
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    pub predicted: bool,
    /// The STM used to compute this Estimate
    pub stm: MatrixMN<f64, <T as State>::Size, <T as State>::Size>,
    /// The Epoch format upon serialization
    pub epoch_fmt: EpochFormat,
    /// The covariance format upon serialization
    pub covar_fmt: CovarFormat,
}

impl<T: State> KfEstimate<T>
where
    DefaultAllocator:
        Allocator<f64, <T as State>::Size> + Allocator<f64, <T as State>::Size, <T as State>::Size>,
{
    pub fn from_covar(
        nominal_state: T,
        covar: MatrixMN<f64, <T as State>::Size, <T as State>::Size>,
    ) -> Self {
        Self {
            nominal_state,
            state_deviation: VectorN::<f64, <T as State>::Size>::zeros(),
            covar: covar.clone(),
            covar_bar: covar,
            predicted: true,
            stm: MatrixMN::<f64, <T as State>::Size, <T as State>::Size>::zeros(),
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
        }
    }
}

impl<T: State> Estimate<T> for KfEstimate<T>
where
    DefaultAllocator:
        Allocator<f64, <T as State>::Size> + Allocator<f64, <T as State>::Size, <T as State>::Size>,
{
    fn zeros(nominal_state: T) -> Self {
        Self {
            nominal_state,
            state_deviation: VectorN::<f64, <T as State>::Size>::zeros(),
            covar: MatrixMN::<f64, <T as State>::Size, <T as State>::Size>::zeros(),
            covar_bar: MatrixMN::<f64, <T as State>::Size, <T as State>::Size>::zeros(),
            predicted: true,
            stm: MatrixMN::<f64, <T as State>::Size, <T as State>::Size>::zeros(),
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
        }
    }

    fn nominal_state(&self) -> T {
        self.nominal_state.clone()
    }

    fn state_deviation(&self) -> VectorN<f64, <T as State>::Size> {
        self.state_deviation.clone()
    }

    fn covar(&self) -> MatrixMN<f64, <T as State>::Size, <T as State>::Size> {
        self.covar.clone()
    }

    fn predicted_covar(&self) -> MatrixMN<f64, <T as State>::Size, <T as State>::Size> {
        self.covar_bar.clone()
    }

    fn predicted(&self) -> bool {
        self.predicted
    }
    fn stm(&self) -> &MatrixMN<f64, <T as State>::Size, <T as State>::Size> {
        &self.stm
    }
    fn epoch_fmt(&self) -> EpochFormat {
        self.epoch_fmt
    }
    fn covar_fmt(&self) -> CovarFormat {
        self.covar_fmt
    }
    fn set_state_deviation(&mut self, new_state: VectorN<f64, <T as State>::Size>) {
        self.state_deviation = new_state;
    }
    fn set_covar(&mut self, new_covar: MatrixMN<f64, <T as State>::Size, <T as State>::Size>) {
        self.covar = new_covar;
    }
}

impl<T: State> fmt::Display for KfEstimate<T>
where
    DefaultAllocator: Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<usize, <T as State>::Size>
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

/// Information filter Estimate
#[derive(Debug, Clone, PartialEq)]
pub struct IfEstimate<T: State>
where
    DefaultAllocator:
        Allocator<f64, <T as State>::Size> + Allocator<f64, <T as State>::Size, <T as State>::Size>,
{
    /// The nominal state
    pub nominal_state: T,
    /// The information state
    pub info_state: VectorN<f64, <T as State>::Size>,
    /// The information matrix, which is the inverse of the covariance
    pub info_mat: MatrixMN<f64, <T as State>::Size, <T as State>::Size>,
    /// The predicted information matrix, which is the inverse of the covariance
    pub info_mat_bar: MatrixMN<f64, <T as State>::Size, <T as State>::Size>,
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    pub predicted: bool,
    /// The STM used to compute this Estimate
    pub stm: MatrixMN<f64, <T as State>::Size, <T as State>::Size>,
    /// The Epoch format upon serialization
    pub epoch_fmt: EpochFormat,
    /// The covariance format upon serialization
    pub covar_fmt: CovarFormat,
}

impl<T: State> IfEstimate<T>
where
    DefaultAllocator:
        Allocator<f64, <T as State>::Size> + Allocator<f64, <T as State>::Size, <T as State>::Size>,
{
    pub fn from_covar(
        nominal_state: T,
        covar: MatrixMN<f64, <T as State>::Size, <T as State>::Size>,
    ) -> Self {
        let mut info_mat = covar;
        if !info_mat.try_inverse_mut() {
            panic!("provided covariance is singular");
        }

        Self {
            nominal_state,
            info_state: VectorN::<f64, <T as State>::Size>::zeros(),
            info_mat: info_mat.clone(),
            info_mat_bar: info_mat,
            predicted: true,
            stm: MatrixMN::<f64, <T as State>::Size, <T as State>::Size>::zeros(),
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
        }
    }

    /// Returns the covariance, if there is enough information to invert the information matrix
    pub fn try_covar(&self) -> Option<MatrixMN<f64, <T as State>::Size, <T as State>::Size>> {
        let mut covar = self.info_mat.clone();
        if !covar.try_inverse_mut() {
            None
        } else {
            Some(&covar * &covar.transpose())
        }
    }

    /// Returns the covariance, if there is enough information to invert the information matrix
    pub fn try_predicted_covar(
        &self,
    ) -> Option<MatrixMN<f64, <T as State>::Size, <T as State>::Size>> {
        let mut covar = self.info_mat_bar.clone();
        if !covar.try_inverse_mut() {
            None
        } else {
            Some(&covar * &covar.transpose())
        }
    }
}

impl<T: State> Estimate<T> for IfEstimate<T>
where
    DefaultAllocator:
        Allocator<f64, <T as State>::Size> + Allocator<f64, <T as State>::Size, <T as State>::Size>,
{
    fn zeros(nominal_state: T) -> Self {
        let mut info_state = VectorN::<f64, <T as State>::Size>::zeros();
        let mut info_mat = MatrixMN::<f64, <T as State>::Size, <T as State>::Size>::zeros();
        // Initialize everything to infinity
        for i in 0..<T as State>::Size::dim() {
            info_state[i] = INFINITY;
            info_mat[(i, i)] = INFINITY;
        }
        Self {
            nominal_state,
            info_state,
            info_mat: info_mat.clone(),
            info_mat_bar: info_mat,
            predicted: true,
            stm: MatrixMN::<f64, <T as State>::Size, <T as State>::Size>::zeros(),
            epoch_fmt: EpochFormat::GregorianUtc,
            covar_fmt: CovarFormat::Sqrt,
        }
    }

    fn nominal_state(&self) -> T {
        self.nominal_state.clone()
    }

    /// Will panic if the information matrix inversion fails
    fn state_deviation(&self) -> VectorN<f64, <T as State>::Size> {
        &self.covar() * &self.info_state
    }

    /// Will panic if the information matrix inversion fails
    fn covar(&self) -> MatrixMN<f64, <T as State>::Size, <T as State>::Size> {
        self.try_covar().unwrap()
    }

    /// Will panic if the information matrix inversion fails
    fn predicted_covar(&self) -> MatrixMN<f64, <T as State>::Size, <T as State>::Size> {
        self.try_predicted_covar().unwrap()
    }

    fn predicted(&self) -> bool {
        self.predicted
    }
    fn stm(&self) -> &MatrixMN<f64, <T as State>::Size, <T as State>::Size> {
        &self.stm
    }
    fn epoch_fmt(&self) -> EpochFormat {
        self.epoch_fmt
    }
    fn covar_fmt(&self) -> CovarFormat {
        self.covar_fmt
    }
    /// WARNING: This sets the information state, not the filter state
    fn set_state_deviation(&mut self, new_info_state: VectorN<f64, <T as State>::Size>) {
        self.info_state = new_info_state;
    }
    /// WARNING: This sets the information matrix
    fn set_covar(&mut self, new_info_mat: MatrixMN<f64, <T as State>::Size, <T as State>::Size>) {
        self.info_mat = new_info_mat;
    }
}

impl<T: State> fmt::Display for IfEstimate<T>
where
    DefaultAllocator: Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<usize, <T as State>::Size>
        + Allocator<usize, <T as State>::Size, <T as State>::Size>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self.try_covar() {
            Some(covar) => {
                let dim = <T as State>::Size::dim();
                let mut fmt_cov = Vec::with_capacity(dim);
                for i in 0..dim {
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

impl<T: State> fmt::LowerExp for IfEstimate<T>
where
    DefaultAllocator: Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<usize, <T as State>::Size>
        + Allocator<usize, <T as State>::Size, <T as State>::Size>,
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

impl<T: State> Serialize for IfEstimate<T>
where
    DefaultAllocator: Allocator<f64, <T as State>::Size>
        + Allocator<f64, <T as State>::Size, <T as State>::Size>
        + Allocator<usize, <T as State>::Size>
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
        match self.try_covar() {
            Some(covar) => {
                let state = self.state_deviation();
                // Serialize the state
                for i in 0..dim {
                    seq.serialize_element(&state[(i, 0)])?;
                }
                // Serialize the covariance
                for i in 0..dim {
                    for j in 0..dim {
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
                for _ in 0..dim {
                    seq.serialize_element(&1e32)?;
                }
                // Serialize the covariance
                for _ in 0..dim {
                    for _ in 0..dim {
                        seq.serialize_element(&1e32)?;
                    }
                }
            }
        }
        seq.end()
    }
}

/// A trait to store a navigation solution, can be used in conjunction with KfEstimate or IfEstimate
pub trait NavSolution<T>: Estimate<Orbit>
where
    T: State,
    DefaultAllocator:
        Allocator<f64, <T as State>::Size> + Allocator<f64, <T as State>::Size, <T as State>::Size>,
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

impl NavSolution<Orbit> for IfEstimate<Orbit> {
    fn orbital_state(&self) -> Orbit {
        self.state()
    }
    fn expected_state(&self) -> Orbit {
        self.nominal_state()
    }
}

// impl NavSolution<SpacecraftState> for KfEstimate<SpacecraftState> {
//     fn orbital_state(&self) -> Orbit {
//         self.state().orbit
//     }
//     fn expected_state(&self) -> Orbit {
//         self.nominal_state().orbit
//     }
// }

// impl NavSolution<SpacecraftState> for IfEstimate<SpacecraftState> {
//     fn orbital_state(&self) -> Orbit {
//         self.state().orbit
//     }
//     fn expected_state(&self) -> Orbit {
//         self.nominal_state().orbit
//     }
// }

impl EstimateFrom<SpacecraftState> for Orbit {
    fn extract(from: &SpacecraftState) -> Self {
        from.orbit
    }

    fn add_dev(to: &SpacecraftState, dev: VectorN<f64, Self::Size>) -> SpacecraftState {
        *to + dev
    }
}

impl EstimateFrom<Orbit> for Orbit {
    fn extract(from: &Orbit) -> Self {
        *from
    }

    fn add_dev(to: &Orbit, dev: VectorN<f64, Self::Size>) -> Orbit {
        *to + dev
    }
}
