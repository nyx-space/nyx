extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, Dim, DimName, U3, Vector3, Vector6, VectorN};

// This determines when to take into consideration the magnitude of the state_delta -- prevents dividing by too small of a number.
const REL_ERR_THRESH: f64 = 0.1;

/// A largest error control which effectively computes the largest error at each component
///
/// This is a standard error computation algorithm, but it's argubly bad if the state's components have different units.
/// It calculates the largest local estimate of the error from the integration (`prop_err`)
/// given the difference in the candidate state and the previous state (`state_delta`).
/// This error estimator is from the physical model estimator of GMAT
/// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/PhysicalModel.cpp#L987]
pub fn largest_error<N: Dim + DimName>(
    prop_err: &VectorN<f64, N>,
    candidate: &VectorN<f64, N>,
    cur_state: &VectorN<f64, N>,
) -> f64
where
    DefaultAllocator: Allocator<f64, N>,
{
    let state_delta = candidate - cur_state;
    let mut max_err = 0.0;
    for (i, prop_err_i) in prop_err.iter().enumerate() {
        let err = if state_delta[(i, 0)] > REL_ERR_THRESH {
            (prop_err_i / state_delta[(i, 0)]).abs()
        } else {
            prop_err_i.abs()
        };
        if err > max_err {
            max_err = err;
        }
    }
    max_err
}

/// A largest step error control which effectively computes the L1 norm of the provided Vector of size 3
///
/// Note that this error controller should be preferrably be used only with slices of a state with the same units.
/// For example, one should probably use this for position independently of using it for the velocity.
/// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/ODEModel.cpp#L3033]
pub fn largest_step(prop_err: &Vector3<f64>, candidate: &Vector3<f64>, cur_state: &Vector3<f64>) -> f64 {
    let state_delta = candidate - cur_state;
    let mag = state_delta[(0, 0)].abs() + state_delta[(1, 0)].abs() + state_delta[(2, 0)].abs();
    let err = prop_err[(0, 0)].abs() + prop_err[(1, 0)].abs() + prop_err[(2, 0)].abs();
    if mag > REL_ERR_THRESH {
        err / mag
    } else {
        err
    }
}

/// A largest state error control
///
/// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/ODEModel.cpp#L3018]
pub fn largest_state(prop_err: &Vector3<f64>, candidate: &Vector3<f64>, cur_state: &Vector3<f64>) -> f64 {
    let sum_state = candidate + cur_state;
    let mag = (sum_state[(0, 0)].abs() + sum_state[(1, 0)].abs() + sum_state[(2, 0)].abs()) * 0.5;
    let err = prop_err[(0, 0)].abs() + prop_err[(1, 0)].abs() + prop_err[(2, 0)].abs();
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
pub fn rss_step(prop_err: &Vector3<f64>, candidate: &Vector3<f64>, cur_state: &Vector3<f64>) -> f64 {
    let mag = (candidate - cur_state).norm();
    let err = prop_err.norm();
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
pub fn rss_state(prop_err: &Vector3<f64>, candidate: &Vector3<f64>, cur_state: &Vector3<f64>) -> f64 {
    let mag = 0.5 * (candidate + cur_state).norm();
    let err = prop_err.norm();
    if mag > REL_ERR_THRESH {
        err / mag
    } else {
        err
    }
}

/// A largest step error control which effectively computes the L1 norm of the provided vector
/// composed of two vectors of the same unit, both of size 3 (e.g. position + velocity).
pub fn largest_step_pos_vel(prop_err: &Vector6<f64>, candidate: &Vector6<f64>, cur_state: &Vector6<f64>) -> f64 {
    let err_radius = largest_step(
        &prop_err.fixed_rows::<U3>(0).into_owned(),
        &candidate.fixed_rows::<U3>(3).into_owned(),
        &cur_state.fixed_rows::<U3>(0).into_owned(),
    );
    let err_velocity = largest_step(
        &prop_err.fixed_rows::<U3>(3).into_owned(),
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
pub fn rss_state_pos_vel(prop_err: &Vector6<f64>, candidate: &Vector6<f64>, cur_state: &Vector6<f64>) -> f64 {
    let err_radius = rss_state(
        &prop_err.fixed_rows::<U3>(0).into_owned(),
        &candidate.fixed_rows::<U3>(3).into_owned(),
        &cur_state.fixed_rows::<U3>(0).into_owned(),
    );
    let err_velocity = rss_state(
        &prop_err.fixed_rows::<U3>(3).into_owned(),
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
pub fn rss_step_pos_vel(prop_err: &Vector6<f64>, candidate: &Vector6<f64>, cur_state: &Vector6<f64>) -> f64 {
    let err_radius = rss_step(
        &prop_err.fixed_rows::<U3>(0).into_owned(),
        &candidate.fixed_rows::<U3>(3).into_owned(),
        &cur_state.fixed_rows::<U3>(0).into_owned(),
    );
    let err_velocity = rss_step(
        &prop_err.fixed_rows::<U3>(3).into_owned(),
        &candidate.fixed_rows::<U3>(3).into_owned(),
        &cur_state.fixed_rows::<U3>(3).into_owned(),
    );

    if err_radius > err_velocity {
        err_radius
    } else {
        err_velocity
    }
}
