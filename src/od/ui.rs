extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, DimName};

pub use super::estimate::*;
pub use super::kalman::*;
pub use super::ranging::*;
pub use super::residual::*;
pub use super::*;

/// Allows to smooth the provided estimates. Returns an array of smoothed estimates.
///
/// Estimates must be ordered in chronological order. This function will smooth the
/// estimates from the last in the list to the first one.
pub fn smooth<S: DimName>(estimates: &Vec<Estimate<S>>) -> Result<Vec<Estimate<S>>, FilterError>
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
