extern crate nalgebra as na;

use self::na::{DefaultAllocator, Dim, DimName, U1, U3, U6, Vector3, VectorN};
use self::na::allocator::Allocator;

// This determines when to take into consideration the magnitude of the state_delta -- prevents dividing by too small of a number.
const REL_ERR_THRESH: f64 = 0.1;

/// A largest error control which effectively computes the largest error at each component
///
/// This is a standard error computation algorithm, but it's argubly bad if the state's components have different units.
/// It calculates the largest local estimate of the error from the integration (`prop_err`)
/// given the difference in the candidate state and the previous state (`state_delta`).
/// This error estimator is from the physical model estimator of GMAT
/// (Source)[https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/forcemodel/PhysicalModel.cpp#L987]
pub fn largest_error<N: Dim + DimName>(prop_err: &VectorN<f64, N>, state_delta: &VectorN<f64, N>) -> f64
where
    DefaultAllocator: Allocator<f64, N>,
{
    let mut max_err = 0.0;
    let rel_threshold = 0.1;
    for (i, prop_err_i) in prop_err.iter().enumerate() {
        let err = if state_delta[(i, 0)] > rel_threshold {
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
pub fn largest_step(prop_err: &VectorN<f64, U3>, state_delta: &VectorN<f64, U3>) -> f64
where
    DefaultAllocator: Allocator<f64, U3>,
{
    let mag = state_delta[(0, 0)].abs() + state_delta[(1, 0)].abs() + state_delta[(2, 0)].abs();
    let err = prop_err[(0, 0)].abs() + prop_err[(1, 0)].abs() + prop_err[(2, 0)].abs();
    if mag > REL_ERR_THRESH {
        err / mag
    } else {
        err
    }
}

/// A largest step error control which effectively computes the L1 norm of the provided vector
/// composed of two vectors of the same unit, both of size 3 (e.g. position + velocity).
pub fn largest_step_pos_vel(prop_err: &VectorN<f64, U6>, state_delta: &VectorN<f64, U6>) -> f64
where
    DefaultAllocator: Allocator<f64, U6>,
{
    // HACK: I'm sure there's an easier way than this, cf. [nalgebra#issue-341](https://github.com/sebcrozet/nalgebra/issues/341).
    let mut prop_err_radius = [0.0; 3];
    let mut prop_err_velocity = [0.0; 3];
    for (i, val) in prop_err.iter().enumerate() {
        if i < 3 {
            prop_err_radius[i] = *val;
        } else {
            prop_err_velocity[i - 3] = *val;
        }
    }

    let mut state_delta_radius = [0.0; 3];
    let mut state_delta_velocity = [0.0; 3];
    for (i, val) in state_delta.iter().enumerate() {
        if i < 3 {
            state_delta_radius[i] = *val;
        } else {
            state_delta_velocity[i - 3] = *val;
        }
    }

    let err_radius = largest_step(
        &Vector3::from_row_slice(&prop_err_radius),
        &Vector3::from_row_slice(&state_delta_radius),
    );
    let err_velocity = largest_step(
        &Vector3::from_row_slice(&prop_err_velocity),
        &Vector3::from_row_slice(&state_delta_velocity),
    );
    if err_radius > err_velocity {
        err_radius
    } else {
        err_velocity
    }
}
