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
//!  * Multi body dynamics
//!
//! ## Usage
//!
//! Put this in your `Cargo.toml`:
//!
//! ```toml
//! [dependencies]
//! nyx-space = "0.0.10"
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
/// use hifitime::julian::ModifiedJulian;
/// use hifitime::SECONDS_PER_DAY;
/// use nyx::celestia::{Cosm, Geoid, State};
/// use nyx::dynamics::celestial::CelestialDynamics;
/// use nyx::dynamics::Dynamics;
/// use nyx::propagators::error_ctrl::RSSStepPV;
/// use nyx::propagators::{PropOpts, Propagator, RK89};
///
/// fn main() {
///     let cosm = Cosm::from_xb("./de438s");
///     let earth_geoid = cosm.geoid_from_id(3).unwrap();
///
///     let dt = ModifiedJulian { days: 21545.0 };
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
///             -5_971.194_376_784_884,
///             3_945.517_912_191_541,
///             2_864.620_958_267_658_4,
///             0.049_083_102_073_914_83,
///             -4.185_084_126_130_087_5,
///             5.848_947_462_252_259_5,
///             ModifiedJulian { days: 21546.0 },
///             earth_geoid,
///     );
///
///     let mut dyn = CelestialDynamics::two_body(initial_state);
///     let mut prop = Propagator::new::<RK89>(
///         &mut dyn,
///         &PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSStepPV {}),
///     );
///     let (final_t, final_state0) = prop.until_time_elapsed(prop_time);
///
///     let final_dt = ModifiedJulian {
///         days: dt.days + final_t / SECONDS_PER_DAY,
///     };
///     let final_state = State::from_cartesian_vec(&prop.state(), final_dt, earth_geoid);
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
///     use nyx::celestia::{Cosm, Geoid, State};
///     let cosm = Cosm::from_xb("./de438s");
///     // In this case, we're creating these states around a Geoid which is Earth.
///     // But for simplicity, we're actually going to use the GMAT value for Earth GM (de438s has a slightly different value).
///     let mut earth_geoid = cosm.geoid_from_id(399).unwrap();
///     earth_geoid.gm = 398_600.441_5;
///     let dt = ModifiedJulian { days: 21545.0 };
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
