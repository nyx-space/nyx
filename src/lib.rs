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
//! nyx-space = "0.0.3"
//! ```
//!
//! And add the following to your crate root:
//!
//! ```rust
//! extern crate nyx_space as nyx;
//! ```

/// Provides all the propagators / integrators available in `nyx`.
///
/// # Custom derivative function example
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
///     use nyx::propagators::{Dormand45, error_ctrl, Options, Propagator};
///     let prop_time = 24.0 * 3_600.0;
///     let accuracy = 1e-12;
///     let min_step = 0.1;
///     let max_step = 30.0;
///     // Initial spacecraft state
///     let mut init_state =
///         Vector6::from_row_slice(&[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);
///     // Final expected spaceraft state
///     let rslt = Vector6::from_row_slice(&[
///         -5971.194191821826,
///         3945.506657649147,
///         2864.636612371127,
///         0.049096952217479194,
///         -4.1850933148636145,
///         5.848940870294863,
///     ]);
///
///     let mut cur_t = 0.0;
///     let mut prop = Propagator::new::<Dormand45>(&Options::with_adaptive_step(0.1, 30.0, 1e-12));
///     loop {
///         let (t, state) = prop.derive(
///             cur_t,
///             &init_state,
///             two_body_dynamics,
///             error_ctrl::rss_state_pos_vel,
///         );
///         if t < prop_time {
///             // We haven't passed the time based stopping condition.
///             cur_t = t;
///             init_state = state;
///         } else {
///             // At this point, we've passed the condition, so let's switch to a fixed step of _exactly_ the
///             // previous time step minus the amount by which we overshot. This allows us to propagate in time for
///             // _exactly_ the time we want to propagate for.
///             let prev_details = prop.latest_details().clone();
///             let overshot = t - prop_time;
///             prop.set_fixed_step(prev_details.step - overshot);
///             // Take one final step
///             let (t, state) = prop.derive(
///                 cur_t,
///                 &init_state,
///                 two_body_dynamics,
///                 error_ctrl::rss_state_pos_vel,
///             );
///
///             assert!(
///                 (t - prop_time).abs() < 1e-12,
///                 "propagated for {} instead of {}",
///                 t,
///                 prop_time
///             );
///
///             // Let's check that, prior to the refined step, we either hit the accuracy wanted,
///             // or we are using the minimum step size.
///             if prev_details.error > accuracy {
///                 assert!(prev_details.step - min_step < f64::EPSILON);
///             }
///
///             assert_eq!(state, rslt, "leo prop failed");
///             break;
///         }
///     }
/// }
/// ```
pub mod propagators;

/// Provides several dynamics used for orbital mechanics and attitude dynamics, which can be elegantly combined.
///
/// # Simple two body propagation
/// ```
/// extern crate nalgebra as na;
/// extern crate hifitime;
/// extern crate nyx_space as nyx;
///
/// fn main() {
///     use nyx::propagators::*;
///     use nyx::celestia::{State, EARTH};
///     use nyx::dynamics::Dynamics;
///     use nyx::dynamics::celestial::TwoBody;
///     use self::na::Vector6;
///     use std::f64;
///     use hifitime::SECONDS_PER_DAY;
///     use hifitime::julian::ModifiedJulian;
///
///     let prop_time = 24.0 * 3_600.0;
///     let accuracy = 1e-12;
///     let min_step = 0.1;
///     let max_step = 60.0;
///
///     let dt = ModifiedJulian { days: 21545.0 };
///     let initial_state = State::from_cartesian::<EARTH>(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0, dt);
///
///     println!("Initial state:\n{0}\n{0:o}\n", initial_state);
///
///     let rslt = State::from_cartesian::<EARTH>(
///         -5971.1941916712285,
///         3945.5066532419537,
///         2864.636618390466,
///         0.04909695760948815,
///         -4.1850933184621315,
///         5.848940867758592,
///         ModifiedJulian { days: 21546.0 }
///     );
///
///     let mut prop = Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy));
///     let mut dyn = TwoBody::from_state_vec::<EARTH>(&initial_state.to_cartesian_vec());
///     let (final_t, final_state_vec) = prop.until_time_elapsed(prop_time, dyn, error_ctrl::rss_step_pos_vel);
///     dyn.set_state(final_t, &final_state_vec);
///
///     let final_dt = ModifiedJulian {
///         days: dt.days + final_t / SECONDS_PER_DAY,
///     };
///     let final_state = State::from_cartesian_vec::<EARTH>(&dyn.state(), final_dt);
///     assert_eq!(final_state, rslt, "two body prop failed",);
///
///     println!("Final state:\n{0}\n{0:o}", final_state);
/// }
/// ```
///
/// # Combining dynamics in a full spacecraft model.
/// ```
/// extern crate nalgebra as na;
/// extern crate hifitime;
/// extern crate nyx_space as nyx;
///
/// // Warning: this is arguably a bad example: attitude dynamics very significantly
/// // faster than orbital mechanics. Hence you really should use different propagators
/// // for the attitude and orbital position and velocity.
/// // If you run this example, you'll notice that the step size used ends up being absolutely tiny
/// // which is needed for the attitude, but not for the astrodynamics (on my machine I'm at 0.04376 seconds).
/// use self::nyx::dynamics::Dynamics;
/// use self::nyx::dynamics::celestial::TwoBody;
/// use self::nyx::dynamics::momentum::AngularMom;
/// use self::nyx::celestia::{State, EARTH};
/// use self::nyx::propagators::{error_ctrl, CashKarp45, Options, Propagator};
/// use self::na::{Matrix3, U3, U6, U9, Vector3, Vector6, VectorN};
/// use std::f64;
/// use hifitime::julian::ModifiedJulian;
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
///         self.twobody
///             .set_state(new_t, &new_state.fixed_rows::<U6>(0).into_owned());
///         self.momentum
///             .set_state(new_t, &new_state.fixed_rows::<U3>(6).into_owned());
///     }
///
///     fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
///         let dpos_vel_dt = self.twobody
///             .eom(_t, &state.fixed_rows::<U6>(0).into_owned());
///         let domega_dt = self.momentum
///             .eom(_t, &state.fixed_rows::<U3>(6).into_owned());
///         <VectorN<f64, U9>>::from_iterator(dpos_vel_dt.iter().chain(domega_dt.iter()).cloned())
///     }
/// }
///
/// // Let's initialize our combined dynamics.
/// fn main() {
///     let prop_time = 3_600.0;
///     let accuracy = 1e-13;
///     let min_step = 0.01;
///     let max_step = 60.0;
///
///     let dyn_twobody = TwoBody::from_state_vec::<EARTH>(&Vector6::new(
///         -2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0,
///     ));
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
///     let mut prop =
///         Propagator::new::<CashKarp45>(&Options::with_adaptive_step(min_step, max_step, accuracy));
///
///     let (final_t, final_state) =
///         prop.until_time_elapsed(prop_time, full_model, error_ctrl::largest_error::<U9>);
///     full_model.set_state(final_t, &final_state);
///
///     let prev_details = prop.latest_details().clone();
///     println!("{:?}", prev_details);
///     if prev_details.error > accuracy {
///         assert!(
///             prev_details.step - min_step < f64::EPSILON,
///             "step size should be at its minimum because error is higher than tolerance: {:?}",
///             prev_details
///         );
///     }
///
///     let delta_mom = ((full_model.momentum.momentum().norm() - init_momentum) / init_momentum).abs();
///     if delta_mom > mom_tolerance {
///         panic!(
///             "angular momentum prop failed: momentum changed by {:e} (> {:e})",
///             delta_mom, mom_tolerance
///         );
///     }
///
///     println!("Final momentum: {:?}", full_model.momentum.momentum());
///
///     println!(
///         "Final orbital state:\n{0}\n{0:o}",
///         State::from_cartesian_vec::<EARTH>(
///             &full_model.twobody.state(),
///             ModifiedJulian { days: 21545.0 }
///         )
///     );
/// }
/// ```
pub mod dynamics;

/// Provides the solar system planets, and state and (later) ephemeride management.
///
/// # State creation and management
/// ```
/// extern crate hifitime;
/// extern crate nyx_space as nyx;
///
/// fn main(){
///     use hifitime::julian::ModifiedJulian;
///     use nyx::celestia::{State, EARTH};
///     let dt = ModifiedJulian { days: 21545.0 };
///     // The parameter is anything which implements `CelestialBody`.
///     // In this case, we're creating these states around Earth.
///     let cart = State::from_cartesian::<EARTH>(
///         5946.673548288958,
///         1656.154606023661,
///         2259.012129598249,
///         -3.098683050943824,
///         4.579534132135011,
///         6.246541551539432,
///         dt,
///     );
///     let kep = State::from_keplerian::<EARTH>(
///         7712.186117895041,
///         0.15899999999999995,
///         53.75369,
///         1.99863286421117e-05,
///         359.787880000004,
///         25.434003407751188,
///         dt,
///     );
///     // We can check whether two states are equal.
///     if cart != kep {
///         panic!("This won't happen");
///     }
///     // Of more interest, we can fetch specific orbital elements.
///     println!("sma = {} km   inc = {} degrees", cart.sma(), cart.inc());
///     // Note that the state data is stored as X, Y, Z, VX, VY, VZ.
///     // Hence, the following print statement may display some rounded values despite
///     // being created with fixed values. GMAT has the same "issue"
///     // (but `nyx` won't change your script).
///     println!("ecc = {} km   RAAN = {} degrees", kep.ecc(), cart.raan());
/// }
/// ```
pub mod celestia;

/// Include utility functions shared by different modules, and which may be useful to engineers.
pub mod utils;

/// Provides all the input/output needs for this library, including loading of SPICE kernels, and gravity potential files.
pub mod io;

/// Provides all the orbital determination tools.
pub mod od;

#[macro_use]
extern crate log;
