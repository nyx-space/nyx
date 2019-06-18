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
//!  * Orbital state definition with transformations to other frames
//!
//! ## Usage
//!
//! Put this in your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! nyx-space = "0.0.7"
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
/// extern crate nalgebra as na;
/// extern crate nyx_space as nyx;
/// use self::na::Vector6;
/// use nyx::celestia::EARTH;
/// use nyx::dynamics::celestial::TwoBody;
/// use nyx::dynamics::Dynamics;
/// use nyx::propagators::error_ctrl::RSSStepPV;
/// use nyx::propagators::*;
///
/// fn main() {
///     let prop_time = 24.0 * 3_600.0;
///     let accuracy = 1e-12;
///     let min_step = 0.1;
///     let max_step = 60.0;
///
///     let init = Vector6::new(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0);
///
///     let rslt = Vector6::from_row_slice(&[
///         -5971.1941916712285,
///         3945.5066532419537,
///         2864.636618390466,
///         0.04909695760948815,
///         -4.1850933184621315,
///         5.848940867758592,
///     ]);
///
///     let mut dyn = TwoBody::from_state_vec::<EARTH>(init);
///     let mut prop = Propagator::new::<RK89>(
///         &mut dyn,
///         &PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSStepPV {}),
///     );
///     prop.until_time_elapsed(prop_time);
///     assert_eq!(prop.state(), rslt, "two body prop failed");
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
/// use hifitime::julian::ModifiedJulian;
/// use hifitime::SECONDS_PER_DAY;
/// use nyx::celestia::{State, EARTH, ECI};
/// use nyx::dynamics::celestial::TwoBody;
/// use nyx::dynamics::Dynamics;
/// use nyx::propagators::error_ctrl::RSSStepPV;
/// use nyx::propagators::{PropOpts, Propagator, RK89};
///
/// fn main() {
///     let dt = ModifiedJulian { days: 21545.0 };
///     let initial_state = State::from_cartesian_eci(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0, dt);
///
///     println!("Initial state:\n{0}\n{0:o}\n", initial_state);
///
///     let prop_time = 24.0 * 3_600.0;
///     let accuracy = 1e-12;
///     let min_step = 0.1;
///     let max_step = 60.0;
///
///     let rslt = State::from_cartesian_eci(
///         -5971.1941916712285,
///         3945.5066532419537,
///         2864.636618390466,
///         0.04909695760948815,
///         -4.1850933184621315,
///         5.848940867758592,
///         ModifiedJulian { days: 21546.0 },
///     );
///
///     let mut dyn = TwoBody::from_state_vec::<EARTH>(initial_state.to_cartesian_vec());
///     let mut prop = Propagator::new::<RK89>(
///         &mut dyn,
///         &PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSStepPV {}),
///     );
///     let (final_t, final_state0) = prop.until_time_elapsed(prop_time);
///
///     let final_dt = ModifiedJulian {
///         days: dt.days + final_t / SECONDS_PER_DAY,
///     };
///     let final_state = State::from_cartesian_vec::<EARTH, ModifiedJulian>(&prop.state(), final_dt, ECI {});
///     assert_eq!(final_state, rslt, "two body prop failed");
///     assert_eq!(prop.state(), final_state0, "until_time_elapsed returns the wrong value");
///
///     println!("Final state:\n{0}\n{0:o}", final_state);
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
///     use nyx::celestia::{State, EARTH, ECI};
///     let dt = ModifiedJulian { days: 21545.0 };
///     // The parameter is anything which implements `CelestialBody`.
///     // In this case, we're creating these states around Earth.
///     let cart = State::from_cartesian::<EARTH, ModifiedJulian>(
///         5946.673548288958,
///         1656.154606023661,
///         2259.012129598249,
///         -3.098683050943824,
///         4.579534132135011,
///         6.246541551539432,
///         dt,
///         ECI {},
///     );
///     let cart_simple = State::from_cartesian_eci(
///         5946.673548288958,
///         1656.154606023661,
///         2259.012129598249,
///         -3.098683050943824,
///         4.579534132135011,
///         6.246541551539432,
///         dt,
///     );
///     let kep = State::from_keplerian::<EARTH, ModifiedJulian>(
///         7712.186117895041,
///         0.15899999999999995,
///         53.75369,
///         1.99863286421117e-05,
///         359.787880000004,
///         25.434003407751188,
///         dt,
///         ECI {},
///     );
///     // We can check whether two states are equal.
///     if cart != kep {
///         panic!("This won't happen");
///     }
///     if cart != cart_simple {
///         panic!("This won't happen either");
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
#[macro_use]
extern crate prost_derive;
