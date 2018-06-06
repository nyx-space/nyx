extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, Dim, DimName, MatrixMN, VectorN};
use std::fmt;
use std::ops::Mul;
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
    DefaultAllocator:
        Allocator<f64, M> + Allocator<f64, S> + Allocator<f64, M, M> + Allocator<f64, M, S> + Allocator<f64, S, S>,
{
    pub fn update_stm(&mut self, new_stm: MatrixMN<f64, S, S>) {
        self.stm = new_stm;
        self.stm_updated = true;
    }

    pub fn update_h_tilde(&mut self, h_tilde: MatrixMN<f64, S, S>) {
        self.h_tilde = h_tilde;
        self.h_tilde_updated = true;
    }

    pub fn time_update(&mut self) -> Result<Estimate<S>, FilterError> {
        if !self.stm_updated {
            return Err(FilterError::STM_NOT_UPDATED);
        }
        let covar_bar = self.stm * self.prev_estimate.covar * self.stm.transpose();
        let state_bar = if self.ekf {
            VectorN::zeros()
        } else {
            self.stm * self.prev_estimate.state()
        };
        let estimate = Estimate {
            state: state_bar,
            covar: covar_bar,
            stm: self.stm,
            predicted: true,
        };
        self.stm_updated = false;
        self.prev_estimate = estimate;
        Ok(estimate)
    }

    pub fn measurement_update(
        &mut self,
        real_obs: VectorN<f64, M>,
        computed_obs: VectorN<f64, M>,
    ) -> Result<Estimate<S>, FilterError> {
        if !self.stm_updated {
            return Err(FilterError::STM_NOT_UPDATED);
        }
        if !self.h_tilde_updated {
            return Err(FilterError::H_TILDE_NOT_UPDATED);
        }
        // Compute Kalman gain
        let covar_bar = self.stm * self.prev_estimate.covar * self.stm.transpose();
        let mut invertible_part = self.h_tilde * covar_bar * self.h_tilde.transpose() + self.measurement_noise;
        if !invertible_part.try_inverse_mut() {
            return Err(FilterError::GAIN_COMP_SINGULAR);
        }
        let gain = covar_bar * self.h_tilde.transpose() * invertible_part;

        // Compute observation deviation (usually marked as y_i)
        let delta_obs = real_obs - computed_obs;

        // Compute the state estimate
        let state_hat = if self.ekf {
            gain * delta_obs
        } else {
            // Must do a time update first
            let state_bar = self.stm * self.prev_estimate.state();
            let innovation = delta_obs - (self.h_tilde * state_bar);
            state_bar + gain * innovation
        };

        // Compute the covariance (Jacobi formulation)
        let covar = (MatrixMN::identity() - gain * self.h_tilde) * covar_bar;

        // And wrap up
        let estimate = Estimate {
            state: state_hat,
            covar: covar,
            stm: self.stm,
            predicted: false,
        };
        self.stm_updated = false;
        self.h_tilde_updated = false;
        self.prev_estimate = estimate;
        Ok(estimate)
    }
}

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

pub enum FilterError {
    STM_NOT_UPDATED,
    H_TILDE_NOT_UPDATED,
    GAIN_COMP_SINGULAR,
}

impl fmt::Display for FilterError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            FilterError::STM_NOT_UPDATED => write!(f, "STM was not updated prior to time or measurement update"),
            FilterError::H_TILDE_NOT_UPDATED => write!(
                f,
                "The measurement matrix H_tilde was not updated prior to measurement update"
            ),
            FilterError::GAIN_COMP_SINGULAR => write!(f, "Gain could not be computed because H*P_bar*H + R is singular"),
        }
    }
}
