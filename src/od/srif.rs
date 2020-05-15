extern crate nalgebra as na;

use self::na::linalg::QR;
use crate::dimensions::allocator::Allocator;
use crate::dimensions::dimension::{DimMin, DimMinimum, DimNameAdd, DimNameSum};
use crate::dimensions::{DefaultAllocator, DimName, MatrixMN, VectorN, U1};

pub use super::estimate::{Estimate, IfEstimate};
pub use super::residual::Residual;
use super::{CovarFormat, EpochFormat, EstimableState, Filter, FilterError};

/// Defines both a Classical and an Extended Kalman filter (CKF and EKF)
#[derive(Debug, Clone)]
pub struct SRIF<S, A, M, T>
where
    S: DimName,
    A: DimName,
    M: DimName,
    T: EstimableState<S>,
    DefaultAllocator: Allocator<f64, S>
        + Allocator<f64, M, M>
        + Allocator<f64, M, S>
        + Allocator<f64, S, S>
        + Allocator<f64, A, A>,
{
    /// The previous estimate used in the KF computations.
    pub prev_estimate: IfEstimate<S, T>,
    /// Sets the Measurement noise (usually noted R)
    pub inv_measurement_noise: MatrixMN<f64, M, M>,
    /// Sets the process noise (usually noted Q) in the frame of the estimated state
    pub process_noise: Option<MatrixMN<f64, A, A>>,
    /// Enables state noise compensation (process noise) only be applied if the time between measurements is less than the process_noise_dt amount in seconds
    pub process_noise_dt: Option<f64>,
    /// Determines whether this KF should operate as a Conventional/Classical Kalman filter or an Extended Kalman Filter.
    /// Recall that one should switch to an Extended KF only once the estimate is good (i.e. after a few good measurement updates on a CKF).
    pub ekf: bool,
    h_tilde: MatrixMN<f64, M, S>,
    stm: MatrixMN<f64, S, S>,
    stm_updated: bool,
    h_tilde_updated: bool,
    epoch_fmt: EpochFormat, // Stored here only for simplification, kinda ugly
    covar_fmt: CovarFormat, // Idem
}

impl<S, A, M, T> SRIF<S, A, M, T>
where
    S: DimName + DimNameAdd<M> + DimMin<M>,
    A: DimName,
    M: DimName + DimNameAdd<S>,
    T: EstimableState<S>,
    DefaultAllocator: Allocator<f64, S>
        + Allocator<f64, M, M>
        + Allocator<f64, M, S>
        + Allocator<f64, DimNameSum<S, M>, S>
        + Allocator<f64, S, S>
        + Allocator<f64, A, A>,
{
    /// Initializes this KF with an initial estimate and measurement noise.
    pub fn initialize(
        initial_estimate: IfEstimate<S, T>,
        process_noise: MatrixMN<f64, A, A>,
        measurement_noise: MatrixMN<f64, M, M>,
        process_noise_dt: Option<f64>,
    ) -> Self {
        let inv_measurement_noise = measurement_noise
            .try_inverse()
            .expect("measurement noise singular");

        let epoch_fmt = initial_estimate.epoch_fmt;
        let covar_fmt = initial_estimate.covar_fmt;

        Self {
            prev_estimate: initial_estimate,
            inv_measurement_noise,
            process_noise: Some(process_noise),
            process_noise_dt,
            ekf: false,
            h_tilde: MatrixMN::<f64, M, S>::zeros(),
            stm: MatrixMN::<f64, S, S>::identity(),
            stm_updated: false,
            h_tilde_updated: false,
            epoch_fmt,
            covar_fmt,
        }
    }
}

impl<S, A, M, T> Filter<S, A, M, T> for SRIF<S, A, M, T>
where
    S: DimName + DimNameAdd<M> + DimNameAdd<S> + DimNameAdd<U1> + DimMin<U1>,
    A: DimName,
    M: DimName + DimNameAdd<S> + DimNameAdd<M> + DimNameAdd<U1>,
    DimNameSum<S, M>: DimMin<DimNameSum<S, U1>>,
    T: EstimableState<S>,
    DefaultAllocator: Allocator<f64, M>
        + Allocator<f64, S>
        + Allocator<f64, M, M>
        + Allocator<f64, M, S>
        + Allocator<f64, DimNameSum<S, M>, DimNameSum<S, U1>>
        + Allocator<f64, DimNameSum<S, U1>, DimNameSum<S, M>>
        + Allocator<f64, DimNameSum<S, M>>
        + Allocator<f64, DimNameSum<S, U1>>
        + Allocator<f64, DimMinimum<DimNameSum<S, M>, DimNameSum<S, U1>>>
        + Allocator<f64, DimMinimum<DimNameSum<S, M>, DimNameSum<S, U1>>, DimNameSum<S, U1>>
        + Allocator<f64, S, M>
        + Allocator<f64, S, S>
        + Allocator<f64, A, A>
        + Allocator<f64, S, A>
        + Allocator<f64, A, S>
        + Allocator<f64, S, U1>,
{
    type Estimate = IfEstimate<S, T>;

    /// Returns the previous estimate
    fn previous_estimate(&self) -> &Self::Estimate {
        &self.prev_estimate
    }

    /// Update the State Transition Matrix (STM). This function **must** be called in between each
    /// call to `time_update` or `measurement_update`.
    fn update_stm(&mut self, new_stm: MatrixMN<f64, S, S>) {
        self.stm = new_stm;
        self.stm_updated = true;
    }

    /// Update the sensitivity matrix (or "H tilde"). This function **must** be called prior to each
    /// call to `measurement_update`.
    fn update_h_tilde(&mut self, h_tilde: MatrixMN<f64, M, S>) {
        self.h_tilde = h_tilde;
        self.h_tilde_updated = true;
    }

    /// Computes a time update/prediction (i.e. advances the filter estimate with the updated STM).
    ///
    /// May return a FilterError if the STM was not updated.
    fn time_update(&mut self, nominal_state: T) -> Result<Self::Estimate, FilterError> {
        if !self.stm_updated {
            return Err(FilterError::StateTransitionMatrixNotUpdated);
        }

        let stm_inv = self
            .stm
            .clone()
            .try_inverse()
            .expect("state transition matrix singular");

        let r_bar = (&self.prev_estimate.info_mat * &stm_inv).abs();

        let state_bar = if self.ekf {
            VectorN::<f64, S>::zeros()
        } else {
            &self.stm * &self.prev_estimate.state_deviation()
        };

        let b_bar = &r_bar * &state_bar;

        let estimate = IfEstimate {
            nominal_state,
            info_state: b_bar,
            info_mat: r_bar,
            stm: self.stm.clone(),
            predicted: true,
            epoch_fmt: self.epoch_fmt,
            covar_fmt: self.covar_fmt,
        };

        self.stm_updated = false;
        self.prev_estimate = estimate.clone();
        Ok(estimate)
    }

    /// Computes the measurement update with a provided real observation and computed observation.
    ///
    /// May return a FilterError if the STM or sensitivity matrices were not updated.
    fn measurement_update(
        &mut self,
        nominal_state: T,
        real_obs: VectorN<f64, M>,
        computed_obs: VectorN<f64, M>,
    ) -> Result<(Self::Estimate, Residual<M>), FilterError> {
        if !self.h_tilde_updated {
            return Err(FilterError::SensitivityNotUpdated);
        }
        // Compute the time update (copy/paste to keep the `_bar` markings)
        if !self.stm_updated {
            return Err(FilterError::StateTransitionMatrixNotUpdated);
        }
        let stm_inv = self
            .stm
            .clone()
            .try_inverse()
            .expect("state transition matrix singular");

        let mut r_bar = (&self.prev_estimate.info_mat * &stm_inv).abs();
        if let Some(pcr_dt) = self.process_noise_dt {
            let delta_t = nominal_state.epoch() - self.prev_estimate.epoch();

            if delta_t <= pcr_dt {
                // Let's compute the Gamma matrix, an approximation of the time integral
                // which assumes that the acceleration is constant between these two measurements.
                let mut gamma = MatrixMN::<f64, S, A>::zeros();
                for i in 0..A::dim() {
                    gamma[(i, i)] = delta_t.powi(2);
                    gamma[(i + A::dim(), i)] = delta_t;
                }
                // Let's add the process noise
                let mut covar_prc = delta_t.powi(2)
                    * (&gamma * self.process_noise.as_ref().unwrap() * &gamma.transpose());
                if !covar_prc.try_inverse_mut() {
                    panic!("process noise singular");
                }
                r_bar += covar_prc;
            }
        }

        let state_bar = if self.ekf {
            VectorN::<f64, S>::zeros()
        } else {
            &self.stm * &self.prev_estimate.state_deviation()
        };

        let b_bar = &r_bar * &state_bar;

        // Done with time update

        // Compute observation deviation (usually marked as y_i)
        let prefit = &self.inv_measurement_noise * &(real_obs - computed_obs);
        // Whiten the h tilde and y
        let h_tilde = &self.inv_measurement_noise * &self.h_tilde;

        // Householder transformation
        let mut stacked = MatrixMN::<f64, DimNameSum<S, M>, DimNameSum<S, U1>>::zeros();
        stacked.fixed_slice_mut::<S, S>(0, 0).copy_from(&r_bar);
        stacked
            .fixed_slice_mut::<S, U1>(0, S::dim())
            .copy_from(&b_bar);
        stacked
            .fixed_slice_mut::<M, S>(S::dim(), 0)
            .copy_from(&h_tilde);
        stacked
            .fixed_slice_mut::<M, U1>(S::dim(), S::dim())
            .copy_from(&prefit);

        // nalgebra uses a Householder transformation to compute the Q, R
        // The householder transform used in the SRIF implementation seems to be opposite of the R computed here (as per gokalman tests).
        let qr = QR::new(stacked);
        let hh_rtn = -qr.r();

        // Extract the data
        let rk = hh_rtn.fixed_slice::<S, S>(0, 0).into_owned();
        let bk = hh_rtn.fixed_slice::<S, U1>(0, S::dim()).to_owned();
        // Get the error estimate -- Don't know where I'd store this for now
        // let ek = hh_rtn[(S::dim() + 1, S::dim() + 1)];

        let mut info_state = VectorN::<f64, S>::zeros();
        for i in 0..S::dim() {
            info_state[i] = bk[(i, 0)];
        }

        let postfit = &prefit - (&self.h_tilde * &state_bar);

        // And wrap up
        let res = Residual::new(nominal_state.epoch(), prefit, postfit);
        let estimate = IfEstimate {
            nominal_state,
            info_state,
            info_mat: rk,
            stm: self.stm.clone(),
            predicted: false,
            epoch_fmt: self.epoch_fmt,
            covar_fmt: self.covar_fmt,
        };
        self.stm_updated = false;
        self.h_tilde_updated = false;
        self.prev_estimate = estimate.clone();
        Ok((estimate, res))
    }

    fn is_extended(&self) -> bool {
        false
    }

    fn set_extended(&mut self, _: bool) {
        unimplemented!();
    }

    fn set_process_noise(&mut self, prc: MatrixMN<f64, A, A>) {
        self.process_noise = Some(prc);
    }
}
