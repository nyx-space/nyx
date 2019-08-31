extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, DimName, MatrixMN, VectorN};
use super::serde::ser::SerializeSeq;
use super::serde::{Serialize, Serializer};
use std::fmt;

/// Defines both a Classical and an Extended Kalman filter (CKF and EKF)
#[derive(Debug, Clone)]
pub struct KF<S, M>
where
    S: DimName,
    M: DimName,
    DefaultAllocator:
        Allocator<f64, M> + Allocator<f64, S> + Allocator<f64, M, M> + Allocator<f64, M, S> + Allocator<f64, S, S>,
{
    /// The previous estimate used in the KF computations.
    pub prev_estimate: Estimate<S>,
    /// Sets the Measurement noise (usually noted R)
    pub measurement_noise: MatrixMN<f64, M, M>,
    /// Determines whether this KF should operate as a Conventional/Classical Kalman filter or an Extended Kalman Filter.
    /// Recall that one should switch to an Extended KF only once the estimate is good (i.e. after a few good measurement updates on a CKF).
    pub ekf: bool,
    h_tilde: MatrixMN<f64, M, S>,
    stm: MatrixMN<f64, S, S>,
    stm_updated: bool,
    h_tilde_updated: bool,
}

impl<S, M> KF<S, M>
where
    S: DimName,
    M: DimName,
    DefaultAllocator: Allocator<f64, M>
        + Allocator<f64, S>
        + Allocator<f64, M, M>
        + Allocator<f64, M, S>
        + Allocator<f64, S, M>
        + Allocator<f64, S, S>,
{
    /// Initializes this KF with an initial estimate and measurement noise.
    pub fn initialize(initial_estimate: Estimate<S>, measurement_noise: MatrixMN<f64, M, M>) -> KF<S, M> {
        KF {
            prev_estimate: initial_estimate,
            measurement_noise,
            ekf: false,
            h_tilde: MatrixMN::<f64, M, S>::zeros(),
            stm: MatrixMN::<f64, S, S>::identity(),
            stm_updated: false,
            h_tilde_updated: false,
        }
    }
    /// Update the State Transition Matrix (STM). This function **must** be called in between each
    /// call to `time_update` or `measurement_update`.
    pub fn update_stm(&mut self, new_stm: MatrixMN<f64, S, S>) {
        self.stm = new_stm;
        self.stm_updated = true;
    }

    /// Update the sensitivity matrix (or "H tilde"). This function **must** be called prior to each
    /// call to `measurement_update`.
    pub fn update_h_tilde(&mut self, h_tilde: MatrixMN<f64, M, S>) {
        self.h_tilde = h_tilde;
        self.h_tilde_updated = true;
    }

    /// Computes a time update (i.e. advances the filter estimate with the updated STM).
    ///
    /// May return a FilterError if the STM was not updated.
    pub fn time_update(&mut self) -> Result<Estimate<S>, FilterError> {
        if !self.stm_updated {
            return Err(FilterError::StateTransitionMatrixNotUpdated);
        }
        let covar_bar = self.stm.clone() * self.prev_estimate.covar.clone() * self.stm.transpose();
        let state_bar = if self.ekf {
            VectorN::<f64, S>::zeros()
        } else {
            self.stm.clone() * self.prev_estimate.state.clone()
        };
        let estimate = Estimate {
            state: state_bar,
            covar: covar_bar,
            stm: self.stm.clone(),
            predicted: true,
        };
        self.stm_updated = false;
        self.prev_estimate = estimate.clone();
        Ok(estimate)
    }

    /// Computes the measurement update with a provided real observation and computed observation.
    ///
    /// May return a FilterError if the STM or sensitivity matrices were not updated.
    pub fn measurement_update(
        &mut self,
        real_obs: VectorN<f64, M>,
        computed_obs: VectorN<f64, M>,
    ) -> Result<Estimate<S>, FilterError> {
        if !self.stm_updated {
            return Err(FilterError::StateTransitionMatrixNotUpdated);
        }
        if !self.h_tilde_updated {
            return Err(FilterError::SensitivityNotUpdated);
        }
        // Compute Kalman gain
        let covar_bar = self.stm.clone() * self.prev_estimate.covar.clone() * self.stm.transpose();
        let mut h_tilde_t = MatrixMN::<f64, S, M>::zeros();
        self.h_tilde.transpose_to(&mut h_tilde_t);
        let mut invertible_part = self.h_tilde.clone() * covar_bar.clone() * h_tilde_t.clone() + self.measurement_noise.clone();
        if !invertible_part.try_inverse_mut() {
            return Err(FilterError::GainIsSingular);
        }
        let gain = covar_bar.clone() * h_tilde_t * invertible_part;

        // Compute observation deviation (usually marked as y_i)
        let delta_obs = real_obs - computed_obs;

        // Compute the state estimate
        let state_hat = if self.ekf {
            gain.clone() * delta_obs
        } else {
            // Must do a time update first
            let state_bar = self.stm.clone() * self.prev_estimate.state.clone();
            let innovation = delta_obs - (self.h_tilde.clone() * state_bar.clone());
            state_bar + gain.clone() * innovation
        };

        // Compute the covariance (Jacobi formulation)
        let covar = (MatrixMN::<f64, S, S>::identity() - gain * self.h_tilde.clone()) * covar_bar;

        // And wrap up
        let estimate = Estimate {
            state: state_hat,
            covar,
            stm: self.stm.clone(),
            predicted: false,
        };
        self.stm_updated = false;
        self.h_tilde_updated = false;
        self.prev_estimate = estimate.clone();
        Ok(estimate)
    }
}

/// Stores an Estimate, as the result of a `time_update` or `measurement_update`.
#[derive(Debug, Clone, PartialEq)]
pub struct Estimate<S>
where
    S: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    /// The estimated state
    pub state: VectorN<f64, S>,
    /// The Covariance of this estimate
    pub covar: MatrixMN<f64, S, S>,
    /// Whether or not this is a predicted estimate from a time update, or an estimate from a measurement
    pub predicted: bool,
    /// The STM used to compute this Estimate
    pub stm: MatrixMN<f64, S, S>,
}

impl<S> Estimate<S>
where
    S: DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    /// An empty estimate. This is useful if wanting to store an estimate outside the scope of a filtering loop.
    pub fn empty() -> Estimate<S> {
        Estimate {
            state: VectorN::<f64, S>::zeros(),
            covar: MatrixMN::<f64, S, S>::zeros(),
            predicted: true,
            stm: MatrixMN::<f64, S, S>::zeros(),
        }
    }

    pub fn header() -> Vec<String> {
        let mut hdr_v = Vec::with_capacity(3 * S::dim());
        for i in 0..S::dim() {
            hdr_v.push(format!("state_{}", i));
        }
        // Serialize the covariance
        for i in 0..S::dim() {
            for j in 0..S::dim() {
                hdr_v.push(format!("covar_{}_{}", i, j));
            }
        }
        hdr_v
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
        let mut seq = serializer.serialize_seq(Some(S::dim() * 3))?;
        // Serialize the state
        for i in 0..S::dim() {
            seq.serialize_element(&self.state[(i, 0)])?;
        }
        // Serialize the covariance
        for i in 0..S::dim() {
            for j in 0..S::dim() {
                seq.serialize_element(&self.covar[(i, j)])?;
            }
        }
        seq.end()
    }
}

/// Stores the different kinds of filter errors.
#[derive(Debug, PartialEq)]
pub enum FilterError {
    StateTransitionMatrixNotUpdated,
    SensitivityNotUpdated,
    GainIsSingular,
}

impl fmt::Display for FilterError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            FilterError::StateTransitionMatrixNotUpdated => {
                write!(f, "STM was not updated prior to time or measurement update")
            }
            FilterError::SensitivityNotUpdated => write!(
                f,
                "The measurement matrix H_tilde was not updated prior to measurement update"
            ),
            FilterError::GainIsSingular => write!(f, "Gain could not be computed because H*P_bar*H + R is singular"),
        }
    }
}
