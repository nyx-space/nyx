extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, DimName, Vector3, Vector6, VectorN, U3};

// This determines when to take into consideration the magnitude of the state_delta and
// prevents dividing by too small of a number.
const REL_ERR_THRESH: f64 = 0.1;

/// The Error Control trait manages how a propagator computes the error in the current step.
pub trait ErrorCtrl {
    /// Computes the actual error of the current step.
    ///
    /// The `error_est` is the estimated error computed from the difference in the two stages of
    /// of the RK propagator. The `candidate` variable is the candidate state, and `cur_state` is
    /// the current state. This function must return the error.
    fn estimate<N: DimName>(error_est: &VectorN<f64, N>, candidate: &VectorN<f64, N>, cur_state: &VectorN<f64, N>) -> f64
    where
        DefaultAllocator: Allocator<f64, N>;
}

/// A largest error control which effectively computes the largest error at each component
///
/// This is a standard error computation algorithm, but it's argubly bad if the state's components have different units.
/// It calculates the largest local estimate of the error from the integration (`error_est`)
/// given the difference in the candidate state and the previous state (`state_delta`).
/// This error estimator is from the physical model estimator of GMAT
/// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/PhysicalModel.cpp#L987]
pub struct LargestError;
impl ErrorCtrl for LargestError {
    fn estimate<N: DimName>(error_est: &VectorN<f64, N>, candidate: &VectorN<f64, N>, cur_state: &VectorN<f64, N>) -> f64
    where
        DefaultAllocator: Allocator<f64, N>,
    {
        let state_delta = candidate - cur_state;
        let mut max_err = 0.0;
        for (i, prop_err_i) in error_est.iter().enumerate() {
            let err = if state_delta[i] > REL_ERR_THRESH {
                (prop_err_i / state_delta[i]).abs()
            } else {
                prop_err_i.abs()
            };
            if err > max_err {
                max_err = err;
            }
        }
        max_err
    }
}

/// A largest step error control which effectively computes the L1 norm of the provided Vector of size 3
///
/// Note that this error controller should be preferrably be used only with slices of a state with the same units.
/// For example, one should probably use this for position independently of using it for the velocity.
/// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/ODEModel.cpp#L3033]
pub struct LargestStep;
impl ErrorCtrl for LargestStep {
    fn estimate<N: DimName>(error_est: &VectorN<f64, N>, candidate: &VectorN<f64, N>, cur_state: &VectorN<f64, N>) -> f64
    where
        DefaultAllocator: Allocator<f64, N>,
    {
        let state_delta = candidate - cur_state;
        let mut mag = 0.0f64;
        let mut err = 0.0f64;
        for i in 0..N::dim() {
            mag += state_delta[i].abs();
            err += error_est[i].abs();
        }

        if mag > REL_ERR_THRESH {
            err / mag
        } else {
            err
        }
    }
}

/// A largest state error control
///
/// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/ODEModel.cpp#L3018]
pub fn largest_state(error_est: &Vector3<f64>, candidate: &Vector3<f64>, cur_state: &Vector3<f64>) -> f64 {
    let sum_state = candidate + cur_state;
    let mag = (sum_state[(0, 0)].abs() + sum_state[(1, 0)].abs() + sum_state[(2, 0)].abs()) * 0.5;
    let err = error_est[(0, 0)].abs() + error_est[(1, 0)].abs() + error_est[(2, 0)].abs();
    if mag > REL_ERR_THRESH {
        err / mag
    } else {
        err
    }
}

/// An RSS step error control which effectively computes the L2 norm of the provided Vector of size 3
///
/// Note that this error controller should be preferrably be used only with slices of a state with the same units.
/// For example, one should probably use this for position independently of using it for the velocity.
/// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/ODEModel.cpp#L3045]
pub fn rss_step(error_est: &Vector3<f64>, candidate: &Vector3<f64>, cur_state: &Vector3<f64>) -> f64 {
    let mag = (candidate - cur_state).norm();
    let err = error_est.norm();
    if mag > REL_ERR_THRESH {
        err / mag
    } else {
        err
    }
}

/// An RSS state error control: when in doubt, use this error controller, especially for high accurracy.
///
/// Here is the warning from GMAT R2016a on this error controller:
/// > This is a more stringent error control method than [`rss_step`] that is often used as the default in other software such as STK.
/// > If you set [the] accuracy to a very small number, 1e-13 for example, and set  the error control to [`rss_step`], integrator
/// > performance will be poor, for little if any improvement in the accuracy of the orbit integration.
/// For more best practices of these integrators (which clone those in GMAT), please refer to the
/// [GMAT reference](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/doc/help/src/Resource_NumericalIntegrators.xml#L1292).
/// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/ODEModel.cpp#L3004]
pub fn rss_state(error_est: &Vector3<f64>, candidate: &Vector3<f64>, cur_state: &Vector3<f64>) -> f64 {
    let mag = 0.5 * (candidate + cur_state).norm();
    let err = error_est.norm();
    if mag > REL_ERR_THRESH {
        err / mag
    } else {
        err
    }
}

/// A largest step error control which effectively computes the L1 norm of the provided vector
/// composed of two vectors of the same unit, both of size 3 (e.g. position + velocity).
pub fn largest_step_pos_vel(error_est: &Vector6<f64>, candidate: &Vector6<f64>, cur_state: &Vector6<f64>) -> f64 {
    let err_radius = LargestStep::estimate(
        &error_est.fixed_rows::<U3>(0).into_owned(),
        &candidate.fixed_rows::<U3>(3).into_owned(),
        &cur_state.fixed_rows::<U3>(0).into_owned(),
    );
    let err_velocity = LargestStep::estimate(
        &error_est.fixed_rows::<U3>(3).into_owned(),
        &candidate.fixed_rows::<U3>(3).into_owned(),
        &cur_state.fixed_rows::<U3>(3).into_owned(),
    );

    if err_radius > err_velocity {
        err_radius
    } else {
        err_velocity
    }
}

/// An RSS state error control which effectively for the provided vector
/// composed of two vectors of the same unit, both of size 3 (e.g. position + velocity).
pub fn rss_state_pos_vel(error_est: &Vector6<f64>, candidate: &Vector6<f64>, cur_state: &Vector6<f64>) -> f64 {
    let err_radius = rss_state(
        &error_est.fixed_rows::<U3>(0).into_owned(),
        &candidate.fixed_rows::<U3>(0).into_owned(),
        &cur_state.fixed_rows::<U3>(0).into_owned(),
    );
    let err_velocity = rss_state(
        &error_est.fixed_rows::<U3>(3).into_owned(),
        &candidate.fixed_rows::<U3>(3).into_owned(),
        &cur_state.fixed_rows::<U3>(3).into_owned(),
    );

    if err_radius > err_velocity {
        err_radius
    } else {
        err_velocity
    }
}

/// A largest step error control which effectively computes the L1 norm of the provided vector
/// composed of two vectors of the same unit, both of size 3 (e.g. position + velocity).
pub fn rss_step_pos_vel(error_est: &Vector6<f64>, candidate: &Vector6<f64>, cur_state: &Vector6<f64>) -> f64 {
    let err_radius = rss_step(
        &error_est.fixed_rows::<U3>(0).into_owned(),
        &candidate.fixed_rows::<U3>(0).into_owned(),
        &cur_state.fixed_rows::<U3>(0).into_owned(),
    );
    let err_velocity = rss_step(
        &error_est.fixed_rows::<U3>(3).into_owned(),
        &candidate.fixed_rows::<U3>(3).into_owned(),
        &cur_state.fixed_rows::<U3>(3).into_owned(),
    );

    if err_radius > err_velocity {
        err_radius
    } else {
        err_velocity
    }
}
