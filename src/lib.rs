/// The `propagators` includes all the propagators / integrators available in `nyx`.
///
/// # Full example
/// ```
/// extern crate nalgebra;
/// extern crate nyx;
/// use nalgebra::{U1, U3, Vector6};
///
/// fn two_body_dynamics(_t: f64, state: &Vector6<f64>) -> Vector6<f64> {
///     let radius = state.fixed_slice::<U3, U1>(0, 0);
///     let velocity = state.fixed_slice::<U3, U1>(3, 0);
///     let body_acceleration = (-398_600.441500000015366822 / radius.norm().powi(3)) * radius;
///     Vector6::from_iterator(velocity.iter().chain(body_acceleration.iter()).cloned())
/// }
///
/// fn main() {
///     use std::f64;
///     use nyx::propagators::{Dormand45, Options, Propagator};
///     // Initial spacecraft state
///     let mut state =
///         Vector6::from_row_slice(&[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);
///     // Final expected spaceraft state
///     let rslt = Vector6::from_row_slice(&[
///         -5971.194191668567,
///         3945.5066531626767,
///         2864.636618498951,
///         0.04909695770740798,
///         -4.185093318527218,
///         5.848940867713008,
///     ]);
///
///     let mut cur_t = 0.0;
///     let mut iterations = 0;
///     let mut prop = Propagator::new::<Dormand45>(&Options::with_adaptive_step(0.1, 30.0, 1e-12));
///     loop {
///         let (t, state_t) = prop.derive(cur_t, &state, two_body_dynamics);
///         iterations += 1;
///         cur_t = t;
///         state = state_t;
///         if cur_t >= 3600.0 * 24.0 {
///             let details = prop.clone().latest_details();
///             if details.error > 1e-2 {
///                 assert!(
///                              details.step - 1e-1 < f64::EPSILON,
///                              "step size should be at its minimum because error is higher than tolerance: {:?}",
///                              details
///                          );
///             }
///             println!("{:?}", prop.latest_details());
///             assert_eq!(state, rslt, "geo prop failed");
///             assert_eq!(iterations, 864_000, "wrong number of iterations");
///             break;
///         }
///     }
/// }
/// ```
pub mod propagators;
