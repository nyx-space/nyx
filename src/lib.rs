//! # nyx-space
//!
//! [Nyx](https://en.wikipedia.org/wiki/Nyx) is a high fidelity, fast, reliable and validated astrodynamical toolkit library written in Rust.
//! It will _eventually_ provide most functionality in Python for rapid prototyping.
//!
//! The target audience is researchers and astrodynamics engineers. The rationale for using Rust is to allow for very fast computations, guaranteed thread safety,
//! and portability to all platforms supported by [Rust](https://forge.rust-lang.org/platform-support.html).
//!
//! To some extend, the ultimate goal of this library is to retire [SPICE Toolkit](https://naif.jpl.nasa.gov/naif/toolkit.html).
//!
//! NOTE: It is recommended to compile all code in `nyx` with the `--release` flag. A lot of heavy
//! computation is done in this library, and no one likes waiting for production code to run.
//! ## Features
//!
//!  * Propagators / Integrators of equations of motions (cf. the `propagators` module)
//!  * Two Body dynamics with planets defined as in GMAT / STK.
//!  * Angular momentum dynamics for a rigid body
//!  * Convenient and explicit definition of the dynamics for a simulation (cf. the [dynamics documentation](./dynamics/index.html))
//!
//! ## Usage
//!
//! Put this in your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! nyx-space = "0.0.2"
//! ```
//!
//! And add the following to your crate root:
//!
//! ```rust
//! extern crate nyx_space as nyx;
//! ```

/// Provides all the propagators / integrators available in `nyx`.
///
/// # Full example
/// ```
/// extern crate nalgebra;
/// extern crate nyx_space as nyx;
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
///             let details = prop.latest_details();
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

/// Provides several dynamics used for orbital mechanics and attitude dynamics, which can be elegantly combined.
///
/// # Combining dynamics in a full spacecraft model.
/// ```
/// extern crate nalgebra as na;
/// extern crate nyx_space as nyx;
/// // Warning: this is arguably a bad example: attitude dynamics very significantly
/// // faster than orbital mechanics. Hence you really should use different propagators
/// // for the attitude and orbital position and velocity.
/// use self::nyx::dynamics::Dynamics;
/// use self::nyx::dynamics::celestial::TwoBody;
/// use self::nyx::dynamics::momentum::AngularMom;
/// use self::nyx::celestia::EARTH;
/// use self::nyx::propagators::{CashKarp45, Options, Propagator};
/// use self::na::{Matrix3, U9, Vector3, Vector6, VectorN};
///
/// // In the following struct, we only store the dynamics because this is only a proof
/// // of concept. An engineer could add more useful information to this struct, such
/// // as a short cut to the position or an attitude.
/// #[derive(Copy, Clone)]
/// pub struct PosVelAttMom {
///     pub twobody: TwoBody,
///     pub momentum: AngularMom,
/// }
///
/// impl Dynamics for PosVelAttMom {
///     type StateSize = U9;
///     fn time(&self) -> f64 {
///         // Both dynamical models have the same time because they share the propagator.
///         self.twobody.time()
///     }
///
///     fn state(&self) -> VectorN<f64, Self::StateSize> {
///         let twobody_state = self.twobody.state();
///         let momentum_state = self.momentum.state();
///         // We're channing both states to create a combined state.
///         // The most important part here is make sure that the `state` and `set_state` handle the state in the same order.
///         <VectorN<f64, U9>>::from_iterator(
///             twobody_state.iter().chain(momentum_state.iter()).cloned(),
///         )
///     }
///
///     fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
///         // HACK: Reconstructing the Vector6 from scratch because for some reason it isn't the correct type when using `fixed_slice`.
///         // No doubt there's a more clever way to handle this, I just haven't figured it out yet.
///         let mut pos_vel_vals = [0.0; 6];
///         let mut mom_vals = [0.0; 3];
///         for (i, val) in new_state.iter().enumerate() {
///             if i < 6 {
///                 pos_vel_vals[i] = *val;
///             } else {
///                 mom_vals[i - 6] = *val;
///             }
///         }
///         self.twobody
///             .set_state(new_t, &Vector6::from_row_slice(&pos_vel_vals));
///         self.momentum
///             .set_state(new_t, &Vector3::from_row_slice(&mom_vals));
///     }
///
///     fn eom(
///         &self,
///         _t: f64,
///         state: &VectorN<f64, Self::StateSize>,
///     ) -> VectorN<f64, Self::StateSize> {
///         // Same issue as in `set_state`.
///         let mut pos_vel_vals = [0.0; 6];
///         let mut mom_vals = [0.0; 3];
///         for (i, val) in state.iter().enumerate() {
///             if i < 6 {
///                 pos_vel_vals[i] = *val;
///             } else {
///                 mom_vals[i - 6] = *val;
///             }
///         }
///         let dpos_vel_dt = self.twobody
///             .eom(_t, &Vector6::from_row_slice(&pos_vel_vals));
///         let domega_dt = self.momentum.eom(_t, &Vector3::from_row_slice(&mom_vals));
///         <VectorN<f64, U9>>::from_iterator(dpos_vel_dt.iter().chain(domega_dt.iter()).cloned())
///     }
/// }
///
/// // Let's initialize our combined dynamics.
///
/// fn main(){
///     let dyn_twobody = TwoBody::around(
///         &Vector6::from_row_slice(&[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]),
///         &EARTH,
///     );
///
///     let omega = Vector3::new(0.1, 0.4, -0.2);
///     let tensor = Matrix3::new(10.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 2.0);
///     let dyn_mom = AngularMom::from_tensor_matrix(&tensor, &omega);
///
///     let mut full_model = PosVelAttMom {
///         twobody: dyn_twobody,
///         momentum: dyn_mom,
///     };
///
///     let init_momentum = full_model.momentum.momentum().norm();
///     let mom_tolerance = 1e-8;
///
///     // And now let's define the propagator and propagate for a short amount of time.
///     let mut prop = Propagator::new::<CashKarp45>(&Options::with_adaptive_step(0.01, 30.0, 1e-12));
///
///     // And propagate
///     loop {
///         let (t, state) = prop.derive(
///             full_model.time(),
///             &full_model.state(),
///             |t_: f64, state_: &VectorN<f64, U9>| full_model.eom(t_, state_),
///         );
///         full_model.set_state(t, &state);
///         if full_model.time() >= 3600.0 {
///             println!("{:?}", prop.latest_details());
///             println!("{}", full_model.state());
///             let delta_mom =
///                 ((full_model.momentum.momentum().norm() - init_momentum) / init_momentum).abs();
///             if delta_mom > mom_tolerance {
///                 panic!(
///                     "angular momentum prop failed: momentum changed by {:e} (> {:e})",
///                     delta_mom, mom_tolerance
///                 );
///             }
///             break;
///         }
///     }
/// }
/// ```
pub mod dynamics;

/// Provides the solar system planets, and (eventually) ephemeride management.
pub mod celestia;

/// Include utility functions shared by different modules, and which may be useful to engineers.
pub mod utils;

#[macro_use]
extern crate lazy_static;
