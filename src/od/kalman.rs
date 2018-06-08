extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, Dim, DimName, MatrixMN, VectorN};
use std::fmt;

#[derive(Debug, Clone)]
pub struct KF<S, M>
where
    S: Dim + DimName,
    M: Dim + DimName,
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
    S: Dim + DimName,
    M: Dim + DimName,
    DefaultAllocator: Allocator<f64, M>
        + Allocator<f64, S>
        + Allocator<f64, M, M>
        + Allocator<f64, M, S>
        + Allocator<f64, S, M>
        + Allocator<f64, S, S>,
{
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
    pub fn update_stm(&mut self, new_stm: MatrixMN<f64, S, S>) {
        self.stm = new_stm;
        self.stm_updated = true;
    }

    pub fn update_h_tilde(&mut self, h_tilde: MatrixMN<f64, M, S>) {
        self.h_tilde = h_tilde;
        self.h_tilde_updated = true;
    }

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
            covar: covar,
            stm: self.stm.clone(),
            predicted: false,
        };
        self.stm_updated = false;
        self.h_tilde_updated = false;
        self.prev_estimate = estimate.clone();
        Ok(estimate)
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Estimate<S>
where
    S: Dim + DimName,
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
    S: Dim + DimName,
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    pub fn empty() -> Estimate<S> {
        Estimate {
            state: VectorN::<f64, S>::zeros(),
            covar: MatrixMN::<f64, S, S>::zeros(),
            predicted: true,
            stm: MatrixMN::<f64, S, S>::zeros(),
        }
    }
}

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
