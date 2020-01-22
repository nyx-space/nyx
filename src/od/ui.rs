extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, DimName};

pub use super::estimate::*;
pub use super::kalman::*;
pub use super::ranging::*;
pub use super::residual::*;
pub use super::*;

use crate::propagators::error_ctrl::ErrorCtrl;
use crate::propagators::Propagator;

pub struct ODProcess<
    'a,
    D: Estimable<N::MeasurementInput, LinStateSize = M::StateSize>,
    E: ErrorCtrl,
    M: Measurement,
    N: MeasurementDevice<M>,
> where
    DefaultAllocator: Allocator<f64, D::StateSize>
        + Allocator<f64, M::MeasurementSize>
        + Allocator<f64, M::MeasurementSize, M::StateSize>
        + Allocator<f64, D::LinStateSize>
        + Allocator<f64, M::MeasurementSize, M::MeasurementSize>
        + Allocator<f64, M::MeasurementSize, D::LinStateSize>
        + Allocator<f64, D::LinStateSize, M::MeasurementSize>
        + Allocator<f64, D::LinStateSize, D::LinStateSize>,
{
    /// Propagator used for the estimation
    pub prop: &'a mut Propagator<'a, D, E>,
    /// Kalman filter itself
    pub kf: &'a mut KF<D::LinStateSize, M::MeasurementSize>,
    /// List of measurement devices used
    pub devices: &'a [N],
    /// Whether or not these devices can make simultaneous measurements of the spacecraft
    pub simultaneous_msr: bool,
    /// Vector of estimates available after a pass
    pub estimates: Vec<Estimate<D::LinStateSize>>,
    /// Vector of residuals available after a pass
    pub residuals: Vec<Residual<M::MeasurementSize>>,
}

impl<
        'a,
        D: Estimable<N::MeasurementInput, LinStateSize = M::StateSize>,
        E: ErrorCtrl,
        M: Measurement,
        N: MeasurementDevice<M>,
    > ODProcess<'a, D, E, M, N>
where
    DefaultAllocator: Allocator<f64, D::StateSize>
        + Allocator<f64, M::MeasurementSize>
        + Allocator<f64, M::MeasurementSize, M::StateSize>
        + Allocator<f64, D::LinStateSize>
        + Allocator<f64, M::MeasurementSize, M::MeasurementSize>
        + Allocator<f64, M::MeasurementSize, D::LinStateSize>
        + Allocator<f64, D::LinStateSize, M::MeasurementSize>
        + Allocator<f64, D::LinStateSize, D::LinStateSize>,
{
    /// Allows to smooth the provided estimates. Returns an array of smoothed estimates.
    ///
    /// Estimates must be ordered in chronological order. This function will smooth the
    /// estimates from the last in the list to the first one.
    pub fn smooth(&mut self) -> Result<Vec<Estimate<D::LinStateSize>>, FilterError> {
        debug!("Smoothing {} estimates", self.estimates.len());
        let mut smoothed = Vec::with_capacity(self.estimates.len());

        for estimate in self.estimates.iter().rev() {
            let mut sm_est = estimate.clone();
            // TODO: Ensure that SNC was _not_ enabled
            let mut stm_inv = estimate.stm.clone();
            if !stm_inv.try_inverse_mut() {
                return Err(FilterError::StateTransitionMatrixSingular);
            }
            sm_est.covar = &stm_inv * &estimate.covar * &stm_inv.transpose();
            sm_est.state = &stm_inv * &estimate.state;
            smoothed.push(sm_est);
        }

        // And reverse to maintain order
        smoothed.reverse();
        Ok(smoothed)
    }

    /// Allows processing all measurements without covariance mapping. Only works for CKFs.
    pub fn process_measurements(&mut self, measurements: &[(Epoch, M)]) -> Option<FilterError> {
        info!("Processing {} measurements", measurements.len());

        let mut prev_dt = self.kf.prev_estimate.dt;

        for (next_epoch, real_meas) in measurements.iter() {
            // Propagate the dynamics to the measurement, and then start the filter.
            let delta_time = *next_epoch - prev_dt;
            prev_dt = *next_epoch; // Update the epoch for the next computation
            self.prop.until_time_elapsed(delta_time);
            // Update the STM of the KF
            self.kf.update_stm(self.prop.dynamics.stm());
            let (dt, meas_input) = self
                .prop
                .dynamics
                .to_measurement(&self.prop.dynamics.state());
            // Get the computed observations
            for device in self.devices.iter() {
                let computed_meas: M = device.measure(&meas_input);
                if computed_meas.visible() {
                    self.kf.update_h_tilde(computed_meas.sensitivity());
                    match self.kf.measurement_update(
                        dt,
                        real_meas.observation(),
                        computed_meas.observation(),
                    ) {
                        Ok((est, res)) => {
                            self.estimates.push(est);
                            self.residuals.push(res);
                        }
                        Err(e) => return Some(e),
                    }
                    if !self.simultaneous_msr {
                        break;
                    }
                }
            }
        }

        None
    }
}
