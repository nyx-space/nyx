extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::dimension::{DimMin, DimMinimum, DimNameAdd, DimNameSum};
use self::na::linalg::QR;
use self::na::{DefaultAllocator, DimName, MatrixMN, VectorN, U1};
use crate::hifitime::Epoch;

pub use super::estimate::{Estimate, IfEstimate};
pub use super::residual::Residual;
use super::{CovarFormat, EpochFormat, Filter, FilterError};

/// Defines both a Classical and an Extended Kalman filter (CKF and EKF)
#[derive(Debug, Clone)]
pub struct SRIF<S, A, M>
where
    S: DimName,
    A: DimName,
    M: DimName,
    DefaultAllocator: Allocator<f64, S>
        + Allocator<f64, M, M>
        + Allocator<f64, M, S>
        + Allocator<f64, S, S>
        + Allocator<f64, A, A>,
{
    /// The previous estimate used in the KF computations.
    pub prev_estimate: IfEstimate<S>,
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

impl<S, A, M> SRIF<S, A, M>
where
    S: DimName + DimNameAdd<M> + DimMin<M>,
    A: DimName,
    M: DimName + DimNameAdd<S>,
    DefaultAllocator: Allocator<f64, S>
        + Allocator<f64, M, M>
        + Allocator<f64, M, S>
        + Allocator<f64, DimNameSum<S, M>, S>
        + Allocator<f64, S, S>
        + Allocator<f64, A, A>,
{
    /// Initializes this KF with an initial estimate and measurement noise.
    pub fn initialize(
        initial_estimate: IfEstimate<S>,
        process_noise: MatrixMN<f64, A, A>,
        measurement_noise: MatrixMN<f64, M, M>,
        process_noise_dt: Option<f64>,
    ) -> Self {
        let inv_measurement_noise = measurement_noise
            .try_inverse()
            .expect("measurement noise singular");

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
            epoch_fmt: EpochFormat::MjdTai,
            covar_fmt: CovarFormat::Sqrt,
        }
    }
}

impl<S, A, M> Filter<S, A, M> for SRIF<S, A, M>
where
    S: DimName + DimNameAdd<M> + DimNameAdd<S> + DimNameAdd<U1>,
    A: DimName,
    M: DimName + DimNameAdd<S> + DimNameAdd<M> + DimNameAdd<U1>,
    DefaultAllocator: Allocator<f64, M>
        + Allocator<f64, S>
        + Allocator<f64, M, M>
        + Allocator<f64, M, S>
        + Allocator<f64, DimNameSum<S, M>, DimNameSum<S, U1>>
        + Allocator<f64, S, M>
        + Allocator<f64, S, S>
        + Allocator<f64, A, A>
        + Allocator<f64, S, A>
        + Allocator<f64, A, S>
        + Allocator<f64, S, U1>,
{
    type Estimate = IfEstimate<S>;

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
    fn time_update(&mut self, dt: Epoch) -> Result<Self::Estimate, FilterError> {
        if !self.stm_updated {
            return Err(FilterError::StateTransitionMatrixNotUpdated);
        }

        let stm_inv = self
            .stm
            .clone()
            .try_inverse()
            .expect("state transition matrix singular");

        let r_bar = &self.prev_estimate.info_mat * &stm_inv;

        let state_bar = if self.ekf {
            VectorN::<f64, S>::zeros()
        } else {
            &self.stm * &self.prev_estimate.state()
        };

        let b_bar = &r_bar * &state_bar;

        let estimate = IfEstimate {
            dt,
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
        dt: Epoch,
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

        let r_bar = &self.prev_estimate.info_mat * &stm_inv;

        let state_bar = if self.ekf {
            VectorN::<f64, S>::zeros()
        } else {
            &self.stm * &self.prev_estimate.state()
        };

        let b_bar = &r_bar * &state_bar;

        // Done with time update

        // Compute observation deviation (usually marked as y_i)
        let prefit = &self.inv_measurement_noise * &(real_obs - computed_obs);
        // Whiten the h tilde and y
        let h_tilde = &self.inv_measurement_noise * &self.h_tilde;
        /*
        dynamics::na::Matrix<f64,
            <S as dynamics::na::DimNameAdd<M>>::Output,
            <S as dynamics::na::DimNameAdd<dynamics::na::U1>>::Output,
            <dynamics::na::DefaultAllocator as
                od::hyperdual::Allocator<f64,
                    <S as dynamics::na::DimNameAdd<M>>::Output,
                    <S as dynamics::na::DimNameAdd<dynamics::na::U1>>::Output
                >
            >::Buffer>
        */

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
        let qr = QR::new(stacked);

        // The householder transform used in the SRIF implementation seems to be opposite of the R computed here (as per gokalman tests).
        // let hh_rtn = -qr.r();
        let hh_rtn = stacked.clone();

        // Extract the data
        let rk = hh_rtn.fixed_slice::<S, S>(0, 0).into_owned();
        let bk = hh_rtn.fixed_slice::<S, U1>(0, S::dim()).to_owned();
        // 	bkMat := A.Slice(0, n, n, n+1)
        //	ekMat := A.Slice(n, n+m, n, n+1)
        let ek = hh_rtn.fixed_slice::<M, U1>(S::dim(), S::dim()).to_owned();

        let mut info_state = VectorN::<f64, S>::zeros();
        let mut postfit = VectorN::<f64, M>::zeros();
        // Can probably be simplified with a from_slice call
        for i in 0..S::dim() {
            info_state[i] = bk[(i, 0)];
        }

        for i in 0..M::dim() {
            postfit[i] = ek[(i, 0)];
        }

        /*
            tmpEst := NewSRIFEstimate(kf.Φ, bk, realObservation, &y, Rk, &RBar)
        est = &tmpEst
        kf.prevEst = est.(*SRIFEstimate)
        kf.step++
        kf.locked = true
        */

        // And wrap up
        let estimate = IfEstimate {
            dt,
            info_state,
            info_mat: rk,
            stm: self.stm.clone(),
            predicted: false,
            epoch_fmt: self.epoch_fmt,
            covar_fmt: self.covar_fmt,
        };
        let res = Residual::new(dt, prefit, postfit);
        self.stm_updated = false;
        self.h_tilde_updated = false;
        self.prev_estimate = estimate.clone();
        Ok((estimate, res))
    }

    fn is_extended(&self) -> bool {
        self.ekf
    }

    fn set_extended(&mut self, status: bool) {
        self.ekf = status;
    }

    fn set_process_noise(&mut self, prc: MatrixMN<f64, A, A>) {
        self.process_noise = Some(prc);
    }
}

#[test]
fn test_hh2() {
    use self::na::{Matrix2, Vector1, VectorN, U1, U2, U3};

    let r_bar = Matrix2::<f64>::new(1.0, -1.0, 1.0, -1.0);
    let b_bar = VectorN::<f64, U2>::new(2.0, 2.0);
    let h_tilde = MatrixMN::<f64, U1, U2>::new(3.0, -3.0);
    let prefit = Vector1::new(4.0);

    let mut stacked = MatrixMN::<f64, U3, U3>::zeros();
    stacked.fixed_slice_mut::<U2, U2>(0, 0).copy_from(&r_bar);
    stacked
        .fixed_slice_mut::<U2, U1>(0, U2::dim())
        .copy_from(&b_bar);
    stacked
        .fixed_slice_mut::<U1, U2>(U2::dim(), 0)
        .copy_from(&h_tilde);
    stacked
        .fixed_slice_mut::<U1, U1>(U2::dim(), U2::dim())
        .copy_from(&prefit);

    println!("{}", &stacked);

    let qr = QR::new(stacked);

    println!("Q {} R {}", qr.q(), qr.r());
}
