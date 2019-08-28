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
//!  * Multi body dynamics (known bug for heliocentric propagation: https://gitlab.com/chrisrabotin/nyx/issues/61)
//!  * Multi body dynamics estimation (i.e. state transition matrix computation, using hyperdual numbers, cf. conf. paper AAS 19-716)
//!
//! ## Usage
//!
//! Put this in your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! nyx-space = "0.0.11"
//! ```
//!
//! And add the following to your crate root:
//!
//! ```rust
//! extern crate nyx_space as nyx;
//! ```

/// Provides all the propagators / integrators available in `nyx`.
pub mod propagators;

/// Provides several dynamics used for orbital mechanics and attitude dynamics, which can be elegantly combined.
///
/// # Simple two body propagation
/// ```
/// extern crate nalgebra as na;
/// extern crate hifitime;
/// extern crate nyx_space as nyx;
/// use hifitime::{Epoch, SECONDS_PER_DAY};
/// use nyx::celestia::{Cosm, Geoid, State};
/// use nyx::dynamics::celestial::CelestialDynamics;
/// use nyx::dynamics::Dynamics;
/// use nyx::propagators::error_ctrl::RSSStepPV;
/// use nyx::propagators::{PropOpts, Propagator, RK89};
///
/// fn main() {
///     let cosm = Cosm::from_xb("./de438s");
///     let earth_geoid = cosm.geoid_from_id(3);
///
///     let dt = Epoch::from_mjd_tai(21_545.0);
///     let initial_state = State::<Geoid>::from_cartesian(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0, dt, earth_geoid);
///
///     println!("Initial state:\n{0}\n{0:o}\n", initial_state);
///
///     let prop_time = 24.0 * 3_600.0;
///     let accuracy = 1e-12;
///     let min_step = 0.1;
///     let max_step = 60.0;
///
///     let rslt = State::<Geoid>::from_cartesian(
///             -5_971.194_376_797_643,
///             3_945.517_912_574_178_4,
///             2_864.620_957_744_429_2,
///             0.049_083_101_605_507_95,
///             -4.185_084_125_817_658,
///             5.848_947_462_472_877,
///             Epoch::from_mjd_tai(21_546.0),
///             earth_geoid,
///     );
///
///     let mut dynamics = CelestialDynamics::two_body(initial_state);
///     let mut prop = Propagator::new::<RK89>(
///         &mut dynamics,
///         &PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSStepPV {}),
///     );
///     prop.until_time_elapsed(prop_time);
///
///     assert_eq!(prop.dynamics.state, rslt, "two body prop failed");
///
///     println!("Final state:\n{0}\n{0:o}", prop.dynamics.state);
/// }
/// ```
///
/// # Multibody propagation of a Halo orbit
/// Multibody propagation is **an order of magnitude faster** in nyx than in GMAT.
/// In nyx, the following function is executed in 0.14 seconds in release mode.
/// ```
/// extern crate nalgebra as na;
/// extern crate hifitime;
/// extern crate nyx_space as nyx;
/// fn main() {
///     use hifitime::Epoch;
///     use na::Vector6;
///     use nyx::celestia::{bodies, Cosm, Geoid, State};
///     use nyx::dynamics::celestial::CelestialDynamics;
///     use nyx::propagators::*;
///     use nyx::utils::rss_state_errors;
///
///     let prop_time = 24.0 * 3_600.0;
///
///     let cosm = Cosm::from_xb("./de438s");
///     let earth_geoid = cosm.geoid_from_id(bodies::EARTH);
///
///     let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
///
///     let halo_rcvr = State::<Geoid>::from_cartesian(
///         333_321.004_516,
///         -76_134.198_887,
///         -20_873.831_939,
///         0.257_153_712,
///         0.930_284_066,
///         0.346_177,
///         start_time,
///         earth_geoid,
///     );
///
///     // GMAT data
///     let rslt = Vector6::new(
///         345_350.664_030_479,
///         5_930.672_047_088,
///         7_333.283_779_286,
///         2.129_819_943e-2,
///         9.566_789_568e-1,
///         3.028_175_811e-1,
///     );
///
///     let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
///     let mut dynamics = CelestialDynamics::new(halo_rcvr, bodies, &cosm);
///
///     let mut prop = Propagator::new::<RK89>(&mut dynamics, &PropOpts::default());
///     prop.until_time_elapsed(prop_time);
///     let (err_r, err_v) = rss_state_errors(&prop.state(), &rslt);
///
///     println!(
///         "RSS errors:\tpos = {:.5e} km\tvel = {:.5e} km/s\ninit\t{}\nfinal\t{}",
///         err_r, err_v, halo_rcvr, prop.dynamics.state
///     );
///     assert!(err_r < 1e-3, format!("multi body failed in position: {:.5e}", err_r));
///     assert!(err_v < 1e-6, format!("multi body failed in velocity: {:.5e}", err_v));
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
///     use hifitime::Epoch;
///     use nyx::celestia::{Cosm, Geoid, State};
///     let cosm = Cosm::from_xb("./de438s");
///     // In this case, we're creating these states around a Geoid which is Earth.
///     // But for simplicity, we're actually going to use the GMAT value for Earth GM (de438s has a slightly different value).
///     let mut earth_geoid = cosm.geoid_from_id(399);
///     earth_geoid.gm = 398_600.441_5;
///     let dt = Epoch::from_mjd_tai(21545.0);
///     let cart = State::<Geoid>::from_cartesian(
///             5_946.673_548_288_958,
///             1_656.154_606_023_661,
///             2_259.012_129_598_249,
///             -3.098_683_050_943_824,
///             4.579_534_132_135_011,
///             6.246_541_551_539_432,
///             dt,
///             earth_geoid,
///     );
///
///     let kep = State::<Geoid>::from_keplerian(
///            7_712.186_117_895_041,
///            0.158_999_999_999_999_95,
///            53.75369,
///            1.998_632_864_211_17e-5,
///            359.787_880_000_004,
///            25.434_003_407_751_188,
///            dt,
///            earth_geoid
///     );
///     // We can check whether two states are equal.
///     if cart != kep {
///         dbg!("{:?}", cart-kep);
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
#[macro_use]
extern crate prost_derive;
extern crate hifitime;
