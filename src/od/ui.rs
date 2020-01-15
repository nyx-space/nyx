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

/// Allows to smooth the provided estimates. Returns an array of smoothed estimates.
///
/// Estimates must be ordered in chronological order. This function will smooth the
/// estimates from the last in the list to the first one.
pub fn smooth<S: DimName>(estimates: &[Estimate<S>]) -> Result<Vec<Estimate<S>>, FilterError>
where
    DefaultAllocator: Allocator<f64, S> + Allocator<f64, S, S>,
{
    debug!("Smoothing {} estimates", estimates.len());
    let mut smoothed = Vec::with_capacity(estimates.len());

    for estimate in estimates.iter().rev() {
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
pub fn process_measurements<
    D: Estimable<N::MeasurementInput, LinStateSize = M::StateSize>,
    E: ErrorCtrl,
    M: Measurement,
    N: MeasurementDevice<M>,
>(
    kf: &mut KF<D::LinStateSize, M::MeasurementSize>,
    prop: &mut Propagator<D, E>,
    measurements: Vec<(Epoch, M)>,
    devices: Vec<N>,
) -> Result<
    (
        Vec<Estimate<D::LinStateSize>>,
        Vec<Residual<M::MeasurementSize>>,
    ),
    FilterError,
>
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
    info!("Processing {} measurements", measurements.len());

    let mut estimates = Vec::with_capacity(measurements.len());
    let mut residuals = Vec::with_capacity(measurements.len());

    let mut prev_dt = kf.prev_estimate.dt;

    for (next_epoch, real_meas) in measurements.iter() {
        // Propagate the dynamics to the measurement, and then start the filter.
        let delta_time = *next_epoch - prev_dt;
        prev_dt = *next_epoch; // Update the epoch for the next computation
        prop.until_time_elapsed(delta_time);
        // Update the STM of the KF
        kf.update_stm(prop.dynamics.stm());
        let (dt, meas_input) = prop.dynamics.to_measurement(&prop.dynamics.state());
        // Get the computed observations
        for device in devices.iter() {
            let computed_meas: M = device.measure(&meas_input);
            if computed_meas.visible() {
                kf.update_h_tilde(computed_meas.sensitivity());
                let (est, res) = kf.measurement_update(
                    dt,
                    real_meas.observation(),
                    computed_meas.observation(),
                )?;
                estimates.push(est);
                residuals.push(res);
                break; // We know that only one station is in visibility at each time.
            }
        }
    }

    Ok((estimates, residuals))
}
