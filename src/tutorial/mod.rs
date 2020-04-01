/*!

This tutorial covers some of the basics of using nyx and its builtin features.

The main advantage of nyx is it's incredible execution speed. It's therefore very well suited to run
Monte Carlos analyzes, i.e. dozens of variations of a given scenario.

We'll start by setting up a test project and play around with creating some states.
Then we'll look into setting up a two body propagation with the different propagator/integrators available.
At that time we'll see how to save the states of the propagator into a file, or how to setup
other processing of these states, like searching for eclipses. We'll subsequently dive into the few
dynamical models provided by nyx (multibody dynamics, solar radiation pressure, etc.), including how to setup
thrusting profiles. Finally, we'll put all of that together and setup an orbit determination scenario with
simulated range and Doppler measurements, and processing those in different filters (Kalman and SRIF).

Nyx is highly configurable, and used at Advanced Space for much more than what's being shown here.
For example, it is relatively easy to create new dynamics which leverage those available in nyx with new equations of motion;
or the definition of new `MeasurementDevice`s which account for specific clock errors when running
specific high fidelity orbit determination scenarios.

**How to use this tutorial:**

+ It's important to read the comments in the listings: that's where most of the code is explained.
+ Copy and paste each code listing into the `fn main(){ ... }` function and run the code. In general, a program in nyx a few times faster
than GMAT (and a few orders of magnitude faster for some stuff like eclipse locators).
+ It's also a good idea to have a set of small toy problems you've solved with other tools before, and want to try solving with nyx.
That way, you can compare how you would solve them in your usual tool, and see how easy or complicated they are to solve with nyx.
Note that nyx is still going through rapid iteration, and every feedback on how to use it is valuable. Let us know on the
[issues](https://gitlab.com/chrisrabotin/nyx/issues) page, or by email at christopher.rabotin [@] gmail.com or rabotin [@] advanced-space.com.

**Known problems:**

+ The TAI time printed is incorrect, cf. [hifitime#67](https://github.com/ChristopherRabotin/hifitime/issues/67)
+ Frame transformations do not account for orientation of the frame, cf. [nyx#93](https://gitlab.com/chrisrabotin/nyx/-/issues/93)

# Table of contents
0. [Setup](#setup)
1. [Celestial computations](#celestial-computations)
    - [Orbital state manipulation](#orbital-state-manipulation)
    - [Planetary state computation](#planetary-state-computation)
    - [Reference frame changes](#reference-frame-changes)
    - [Visibility computation](#visibility-computation)
2. [Propagation](#propagation)
    - [Propagation with different Runge Kutta methods](#propagation-with-different-runge-kutta-methods)
    - [Propagation to different stopping conditions](#propagation-to-different-stopping-conditions)
    - [Building new events](#building-new-events)
    - [Saving all the states](#saving-all-the-states)
    - [Eclipse locator](#eclipse-locators)
3. [Dynamical models](#dynamical-models)
    - [Multibody dynamics](#multibody-dynamics)
    - [Basic attitude dynamics](#basic-attitude-dynamics)
    - [Finite burns with fuel depletion](#finite-burns-with-fuel-depletion)
    - [Sub-Optimal Control of continuous thrust](#sub-optimal-control-of-continuous-thrust)
    - Solar radiation pressure modeling
    - Basic drag models (cannonball)
    - Defining new dynamical models
4. [Orbit determination](#orbit-determination)
    - Generating orbit determination measurements
    - Classical and Extended Kalman Filter
    - Saving the estimates
    - State noise compensation (SNC)
    - Smoothing and iterations of CKFs
    - Square Root Information Filer (SRIF)

# Setup

<!-- IMPORTANT -->
<!-- This section is copied from https://docs.rs/crate/csv/1.1.3/source/src/tutorial.rs -->

In this section, we'll get you setup with a simple program that reads CSV data
and prints a "debug" version of each record. This assumes that you have the
[Rust toolchain installed](https://www.rust-lang.org/install.html),
which includes both Rust and Cargo.

We'll start by creating a new Cargo project:

```text
$ cargo new --bin nyxtutorial
$ cd nyxtutorial
```

Once inside `nyxtutorial`, open `Cargo.toml` in your favorite text editor and add
`nyx_space = "0.0.20-alpha.1"` to your `[dependencies]` section.
All high fidelity time management is handled by the hifitime library, and is accessible via `nyx::time`.
The documentation for hifitime is available on [docs.rs](https://docs.rs/hifitime).
Mainly, you'll be using the `Epoch` structure of this library, whose documentation is [here](https://docs.rs/hifitime/1.0.10/hifitime/struct.Epoch.html).

At this point, your
`Cargo.toml` should look something like this:

```text
[package]
name = "nyxtutorial"
version = "0.1.0"
authors = ["Your Name"]

[dependencies]
nyx-space = "0.0.20-alpha.1"
```

Next, let's build your project. Since you added the `nyx_space` crate as a
dependency, Cargo will automatically download it and compile it for you. To
build your project, use Cargo:

```text
$ cargo build
```

This will produce a new binary, `nyxtutorial`, in your `target/debug` directory.
It won't do much at this point, but you can run it:

```text
$ ./target/debug/nyxtutorial
Hello, world!
```

You can also directly run the program with an implicit building step:
```text
$ cargo run
(...)
Hello, world!
```

We'll make our program more useful in the rest of this tutorial. The function which is executed when the program is run is `main()`, located in `src/main.rs`.
Unless specified otherwise, we can replace the contents of `main` function with any of the examples below, execte `cargo run`, and it'll build and run the example.

# Celestial computations
## Orbital state manipulation
As in any astrodynamics toolkit, it's essential to be able to convert a state from one definition to another (e.g. Cartesian state to Keplerian state).
All states are defined with respect to a reference frame. In nyx, that reference frame also contains useful information, such as the center body's gravitational parameter.

Let's start by loading an XB file (which stands for eXchange Binary). Think of these files as a clone of JPL's BSP files (only easier to read and used onboard spacecraft).
Although these files are currently proprietary to Advanced Space, there is an on-going effort to release the specifications and related tools as we think these are much more useful than SPICE's DAF BSP files.

Start by downloading [de438s.exb](https://gitlab.com/chrisrabotin/nyx/-/blob/master/de438.exb) and [de438s.fxb](https://gitlab.com/chrisrabotin/nyx/-/blob/master/de438.fxb)
and placing them in your `nyxtutorial` directory.

Orbital states are stored in an [`State`](../celestia/type.State.html) structure. It has a lot of useful methods, like computing the angular momentum of that orbit, or its period.
It can be initialized from Cartesian coordinates or Keplerian elements, but the data is stored in Cartesian. Therefore, a state initialized from Keplerian elements might not return exactly the same input values you specified.

_Note_: the following code goes inside the `fn main()` function. The same applies for all of the examples in this tutorial. However, it is common in Rust to move the `extern crate` outside of functions.

```
// Import nyx
extern crate nyx_space as nyx;
// Tell Rust that we're using the structures called Epoch, Cosm and State
// from the respective modules. We're also importing _bodies_ which has a mapping from
// the name of the object to its ID in the EXB file.
use nyx::celestia::{bodies, Cosm, State};
use nyx::time::Epoch;

// Load both the de438s.exb and de438s.fxb (resp. the ephemeris and frame information).
let cosm = Cosm::from_xb("./de438s");
// Get a copy of the Earth body and frame information from the cosm.
let earth = cosm.frame_by_id(bodies::EARTH);
// Initialize a Epoch for 31 January 2020 at midnight TAI.
let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 31);

// And initialize a Cartesian state with position, velocity, epoch and center object.
// The position units are in kilometers and the velocity units in kilometers per second.
let cart = State::cartesian(
    -2436.45, -2436.45, 6891.037, // X, Y, Z (km)
    5.088_611, -5.088_611, 0.0, // VX, VY, VZ (km/s)
    dt, earth,
);

// And print this state in Cartesian
println!("{}", cart);
// And print this state in Keplerian
println!("{:o}", cart);

// Now let's recreate the state but from the Keplerian elements

let kep = State::keplerian(
    cart.sma(),
    cart.ecc(),
    cart.inc(),
    cart.raan(),
    cart.aop(),
    cart.ta(),
    dt,
    earth,
);

// This should be very close to zero
println!("{:e}", cart - kep);
```

## Planetary state computation
**IMPORTANT: Update this section after attitude frame transformations have been implemented**

All of the ephemeris computation happens through the `Cosm` structure, whose documentation is [here](../celestia/struct.Cosm.html).

Let's start by getting the position and velocity of the Earth with respect to the Sun at a given time. As of 0.0.20-alpha.1, these are all in J2000 equatorial.

```
extern crate nyx_space as nyx;
use nyx::celestia::{bodies, Cosm, LTCorr};
use nyx::time::Epoch;

// Load the ephemeris
let cosm = Cosm::from_xb("./de438s");
// Define a time
let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

// And finally get the Earth state
// Note the LTCorr:None means we are _not_ correction for light time.
let earth_as_seen_from_sun = cosm.celestial_state(bodies::EARTH, dt, bodies::SUN, LTCorr::None);

// Get the same state with light time correction, but no abberation.
let with_lt = cosm.celestial_state(bodies::EARTH, dt, bodies::SUN, LTCorr::LightTime);

// Get the same state with light time correction, with light time corrections and abberation.
let with_abbr = cosm.celestial_state(bodies::EARTH, dt, bodies::SUN, LTCorr::Abberation);

println!("{}\n{}\n{}", earth_as_seen_from_sun, with_lt, with_abbr);
```

By executing `cargo run` in the `nyxtutorial` folder, the following should be printed to console:

```text
[Geoid 10] 2019-12-31T23:59:23 TAI      position = [-24885932.304049, 133017335.386064, 57663346.118202] km     velocity = [-29.848886, -4.736858, -2.052876] km/s
[Geoid 10] 2019-12-31T23:59:23 TAI      position = [-24871279.186939, 133019660.572113, 57664353.602189] km     velocity = [-29.849413, -4.734134, -2.051694] km/s
[Geoid 10] 2019-12-31T23:59:23 TAI      position = [-24871286.361793, 133019659.365008, 57664353.292134] km     velocity = [-29.849413, -4.734134, -2.051694] km/s
```

## Reference frame changes
We can also use `Cosm` in order to convert `State`s into other frames. In the following example, we take a orbital state around the Earth, compute it a Moon centered frame.

```
extern crate nyx_space as nyx;
use nyx::celestia::{bodies, Cosm, State};
use nyx::time::Epoch;

let cosm = Cosm::from_xb("./de438s");
let earth = cosm.frame_by_id(bodies::EARTH);
let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 31);

// And initialize a Cartesian state with position, velocity, epoch and center object.
let state = State::cartesian(
    -2436.45, -2436.45, 6891.037, // X, Y, Z (km)
    5.088_611, -5.088_611, 0.0, // VX, VY, VZ (km/s)
    dt, earth,
);

// Get a copy of the Moon
let moon = cosm.frame_by_id(bodies::EARTH_MOON);
// And convert that to a Moon centered frame using the Cosm
// Note that the amperstand in front of state is to specify that we're
// taking a reference of the state, instead of copying it.
// This makes things faster and uses less memory.
let seen_from_moon = cosm.frame_chg(&state, moon);

println!("{}\n{}", state, seen_from_moon);
```

Executing this, the output should be:
```text
[Geoid 399] 2020-01-30T23:59:23 TAI     position = [-2436.450000, -2436.450000, 6891.037000] km velocity = [5.088611, -5.088611, 0.000000] km/s
[Geoid 301] 2020-01-30T23:59:23 TAI     position = [-386674.710813, -128504.451448, -8399.855914] km    velocity = [5.393609, -5.923479, -0.379899] km/s
```

## Visibility computation
Now that we know how to setup an State, we can set up several orbit states, and convert them into other frames.

Let's build three states in a circular orbit and check whether they are in line of sight given the position of some celestial object.
```
extern crate nyx_space as nyx;

use nyx::celestia::eclipse::{line_of_sight, EclipseState};
use nyx::celestia::{bodies, Cosm, State};
use nyx::time::Epoch;

// Load the ephemeris
let cosm = Cosm::from_xb("./de438s");
// Get a copy of the Earth
let earth = cosm.frame_by_id(bodies::EARTH);

let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

// Define the semi major axis of these orbits
let sma = earth.equatorial_radius + 300.0;

let sc1 = State::keplerian(sma, 0.001, 0.1, 90.0, 75.0, 0.0, dt, earth);
let sc2 = State::keplerian(sma + 1.0, 0.001, 0.1, 90.0, 75.0, 0.0, dt, earth);
let sc3 = State::keplerian(sma, 0.001, 0.1, 90.0, 75.0, 180.0, dt, earth);

// Both states are out of phase by pi, so the Earth actually prevents both spacecraft
// from being in line of sight of each other.
let sc1_sc3_visibility = line_of_sight(&sc1, &sc3, earth, &cosm);
println!("SC1 <-> SC3: {}", sc1_sc3_visibility);

// This `assert_eq!` is a Rust macro which will cause the program to fail if
// the information on the left is different from the one on the right.
// Therefore, we're expecting the visibility to be "umbra" (shadow).
assert_eq!(sc1_sc3_visibility, EclipseState::Umbra);

// Nearly identical orbits in the same phasing
let sc1_sc2_visibility = line_of_sight(&sc1, &sc2, earth, &cosm);
println!("SC1 <-> SC2: {}", sc1_sc2_visibility);
assert_eq!(sc1_sc2_visibility, EclipseState::Visibilis);
```

# Propagation
Awesome! So far, we've seen how to create a state, so let's take this state and propagate it.
Let's first just use a default propagator with the default options.

```
extern crate nyx_space as nyx;
use nyx::celestia::{bodies, Cosm, State};
use nyx::propagators::{PropOpts, Propagator};
use nyx::time::{Epoch, SECONDS_PER_DAY};

let cosm = Cosm::from_xb("./de438s");
let earth = cosm.frame_by_id(bodies::EARTH);

let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 31);

let state = State::cartesian(
    -2436.45, -2436.45, 6891.037, // X, Y, Z (km)
    5.088_611, -5.088_611, 0.0, // VX, VY, VZ (km/s)
    dt, earth,
);

// Let's initialize two body celestial dynamics.
// (the use statement is here for clarity, but should be moved to the top of the file).
use nyx::dynamics::celestial::CelestialDynamics;
// Note that we're defining this variable as `mut`.
// This means the variable is mutable, i.e. can be changed or modified.
let mut dynamics = CelestialDynamics::two_body(state);

// Let's setup the propagator options.
// The default is to use the same configuration as NASA GMAT:
// An RSS Step step error control, with an adaptive step size.
// The step size
let opts = PropOpts::default();
println!("default options: {}", opts.info());

// Now let's setup the propagator. The default propagator is an RK89.
// We pass it a pointer to the dynamics we've defined above. Specifically,
// we're giving is a mutuable reference to the dynamics: `&` is the reference,
// and `mut` makes it mututable.
let mut prop = Propagator::default(&mut dynamics, &opts);

// Determine the propagation time, in seconds!
let prop_time = SECONDS_PER_DAY;
// Tell the propagator to integrate the trajectory for that specific amount of time.
// You can get the last state from the output of the until_time_elapsed call
let last_state_0 = prop.until_time_elapsed(prop_time);

// Or by calling state() on the dynamics of the propagator.
// Note that in order to call state(), we need to tell Rust that we're using
// the `dynamics` variable as a Dynamics, so we need to `use` dynamics::Dynamics.
use nyx::dynamics::Dynamics;
let last_state_1 = prop.dynamics.state();
println!("{}\n{}", last_state_0, last_state_1);

// We can check that the propagator works well by doing a back propagation
// of that same amount of time, and checking that the value matches
let backprop_state = prop.until_time_elapsed(-prop_time);
println!("{}\n{}", state, backprop_state);

// Finally, you can check the integration step used for the last step
// as follows: (the {:?} specifies that we want to print in Debug mode)
println!("{:?}", prop.latest_details());
```

Execute `cargo run`, and you should get something like this:
```text
default options: [min_step: 1e-3, max_step: 2.7e3, tol: 1e-12, attempts: 50]
[Geoid 399] 2020-01-31T23:59:23 TAI     position = [-5971.194377, 3945.517913, 2864.620958] km  velocity = [0.049083, -4.185084, 5.848947] km/s
[Geoid 399] 2020-01-31T23:59:23 TAI     position = [-5971.194377, 3945.517913, 2864.620958] km  velocity = [0.049083, -4.185084, 5.848947] km/s
[Geoid 399] 2020-01-30T23:59:23 TAI     position = [-2436.450000, -2436.450000, 6891.037000] km velocity = [5.088611, -5.088611, 0.000000] km/s
[Geoid 399] 2020-01-30T23:59:23 TAI     position = [-2436.449999, -2436.450001, 6891.037000] km velocity = [5.088611, -5.088611, -0.000000] km/s
IntegrationDetails { step: -41.788137644246206, error: 0.0000000000003875125525111948, attempts: 1 }
```

Note the difference between the initial state (`state`) and the backward propagated state. That's normal!
It's because of the precision limits of the machine's 64 bit floating point value.

## Propagation with different Runge Kutta methods
There are plenty of [built-in propagators](../propagators/index.html#structs) in nyx. You may want to use CashKarp45 which is 4th order with 5th order correction integrator.

The default propagator works fine for most cases, but sometimes you will want something specific (for validation against another tool for example).

```
extern crate nyx_space as nyx;
use nyx::celestia::{bodies, Cosm, State};
use nyx::dynamics::celestial::CelestialDynamics;
use nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use nyx::time::{Epoch, SECONDS_PER_DAY};

let cosm = Cosm::from_xb("./de438s");
let earth = cosm.frame_by_id(bodies::EARTH);

let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 31);

let state = State::cartesian(
    -2436.45, -2436.45, 6891.037, // X, Y, Z (km)
    5.088_611, -5.088_611, 0.0, // VX, VY, VZ (km/s)
    dt, earth,
);

// Let's initialize two body celestial dynamics.
let mut dynamics = CelestialDynamics::two_body(state);

// Let's setup the propagator options.
// We's setting it to be a with_fixed_step.
let step_size = 10.0; // in seconds
let opts = PropOpts::with_fixed_step(step_size);
println!("propagator options: {}", opts.info());

// Now let's setup the propagator.
// We're using the "turbofish" operator to specify to Rust
// that we want the propagator to be an RK4Fixed type.
let mut prop = Propagator::new::<RK4Fixed>(&mut dynamics, &opts);

let last_state = prop.until_time_elapsed(SECONDS_PER_DAY);

println!("{}", last_state);
```

The output execution should be something like this:
```text
propagator options: [min_step: 1e1, max_step: 1e1, tol: 0e0, attempts: 0]
[Geoid 399] 2020-01-31T23:59:23 TAI     position = [-5971.194374, 3945.517805, 2864.621106] km  velocity = [0.049083, -4.185084, 5.848947] km/s
```

Read more about the [turbofish operator here](https://matematikaadit.github.io/posts/rust-turbofish.html).


## Propagation to different stopping conditions
In many cases, you'll want to propagate until a given condition is met instead of just propagating for an amount of time.

Nyx allows you to propagate until a condition is reached once or more times. We'll look at two such examples.

The condition stopper uses a Brent root solver. Documentation on StopCondition is available [here](../propagators/events/struct.StopCondition.html).

```
extern crate nyx_space as nyx;
use nyx::celestia::{bodies, Cosm, State};
use nyx::dynamics::celestial::CelestialDynamics;
use nyx::propagators::error_ctrl::RSSStepPV;
use nyx::propagators::events::{EventKind, OrbitalEvent, StopCondition};
use nyx::propagators::{PropOpts, Propagator};
use nyx::time::Epoch;

let cosm = Cosm::from_xb("./de438s");
let earth = cosm.frame_by_id(bodies::EARTH);

let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 31);

let state = State::cartesian(
    -2436.45, -2436.45, 6891.037, // X, Y, Z (km)
    5.088_611, -5.088_611, 0.0, // VX, VY, VZ (km/s)
    dt, earth,
);

// Save the period of this orbit.
let period = state.period();
// Define the maximum propagation time to be four times the orbit
let max_prop_time = 4.0 * period;
// Define the tolerance of the stopping condition.

// Define the apoapsis event
let apo_event = OrbitalEvent::new(EventKind::Apoapse);
// And set it up as a condition.
let condition = StopCondition::new(apo_event, max_prop_time, 1e-6);

let mut dynamics = CelestialDynamics::two_body(state);

// Note that the precision of the condition is matched only as good as the
// minimum step of the propagator options.
let opts = PropOpts::with_adaptive_step(1.0, 60.0, 1e-9, RSSStepPV {});

let mut prop = Propagator::default(&mut dynamics, &opts);

match prop.until_event(condition) {
    Err(convergence_error) => println!("Did not converge!\n{:?}", convergence_error),
    Ok(state) => {
        // Print all of the condition crossings / iterations the Brent solver
        println!("{}", prop.event_trackers);
        println!("Final time: {}", state.dt.as_gregorian_utc_str());
        // Or print all of the crossings with this
        println!("All crossings:\n{:?}", prop.event_trackers.found_bounds);
    }
}
```

The output should be something like:
```text
Converged on (3370.134633983629, 3370.134662763675) for event OrbitalEvent { kind: Apoapse, tgt: None, cosm: None }
Final time: 2020-01-31T00:55:33
All crossings:
[[(3360.0, 3420.0), (3420.0, 3390.1780649342677), (3360.0, 3375.0890324671336), (3367.544516233567, 3375.0890324671336), (3375.0890324671336, 3371.31677435035), (3369.430645291958, 3371.31677435035), (3371.31677435035, 3370.373709821154), (3369.902177556556, 3370.373709821154), (3370.373709821154, 3370.137943688855), (3370.0200606227054, 3370.137943688855), (3370.0790021557805, 3370.137943688855), (3370.1084729223176, 3370.137943688855), (3370.1232083055866, 3370.137943688855), (3370.130575997221, 3370.137943688855), (3370.1342598430383, 3370.137943688855), (3370.137943688855, 3370.1361017659465), (3370.1342598430383, 3370.1351808044924), (3370.1342598430383, 3370.1347203237656), (3370.134490083402, 3370.1347203237656), (3370.1346052035838, 3370.1347203237656), (3370.1347203237656, 3370.134662763675), (3370.1346052035838, 3370.134633983629), (3370.134633983629, 3370.134662763675)]]
```

Let's now allow for the apoapse event to happen several times and stop on the third passing.

```
extern crate nyx_space as nyx;
use nyx::celestia::{bodies, Cosm, State};
use nyx::dynamics::celestial::CelestialDynamics;
use nyx::propagators::error_ctrl::RSSStepPV;
use nyx::propagators::events::{EventKind, OrbitalEvent, StopCondition};
use nyx::propagators::{PropOpts, Propagator};
use nyx::time::Epoch;

let cosm = Cosm::from_xb("./de438s");
let earth = cosm.frame_by_id(bodies::EARTH);

let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 31);

let state = State::cartesian(
    -2436.45, -2436.45, 6891.037, // X, Y, Z (km)
    5.088_611, -5.088_611, 0.0, // VX, VY, VZ (km/s)
    dt, earth,
);

let period = state.period();
let max_prop_time = 4.0 * period;

let apo_event = OrbitalEvent::new(EventKind::Apoapse);

// And now let's stop the propagation after three apoapsis crossings.
let condition = StopCondition::after_hits(apo_event, 3, max_prop_time, 1e-6);

let mut dynamics = CelestialDynamics::two_body(state);

let opts = PropOpts::with_adaptive_step(1.0, 60.0, 1e-9, RSSStepPV {});

let mut prop = Propagator::default(&mut dynamics, &opts);

match prop.until_event(condition) {
    Err(convergence_error) => println!("Did not converge!\n{:?}", convergence_error),
    Ok(state) => {
        // Print all of the condition crossings / iterations the Brent solver
        println!("{}", prop.event_trackers);
        println!("Final time: {}", state.dt.as_gregorian_utc_str());
        // Confirm that this is the third apoapse event which is found
        assert!(
            state.dt - dt < 3.0 * period && state.dt - dt >= 2.0 * period,
            "converged on the wrong apoapse"
        );
        assert!(
            (180.0 - state.ta()) < 1e-6,
            "converged, yet convergence critera not met"
        );
    }
}
```

Whose output should be:
```text
Converged on (16850.673179626465, 16850.673294067383) for event OrbitalEvent { kind: Apoapse, tgt: None, cosm: None }
Final time: 2020-01-31T04:40:13
```

## Building new events

Nyx provides a few kinds of events the list is available [here](../propagators/events/trait.Event.html#implementors). Anyone can create new kinds of events by implementing the `Event` trait.
As previous stated, nyx is designed to be very modular and flexible. If you want to track an event which isn't currently supported, you can create your own event.

The best way to build a new event is to copy the code from from [`ScEvent`](../propagators/events.rs.html#178), and modify it to suit your needs in your own project.

## Saving all the states
Oftentimes we want to save the states of a propagation to a file for some kind of post-processing.

In this example, we'll store all of the states of a propagation into a CSV file.

To do this, we'll need to attach "channels" to the propagator. This allows the propagator to send states to another thread.
Open the `Cargo.toml` file and add `csv = "1"` to specify that we're using the [CSV crate](https://docs.rs/csv/), version 1 or above.

```
extern crate nyx_space as nyx;
extern crate csv;
// ^^^ Allows us to the CSV crate
use nyx::celestia::{bodies, Cosm, State};
use nyx::dynamics::celestial::CelestialDynamics;
use nyx::propagators::{PropOpts, Propagator};
use nyx::time::{Epoch, SECONDS_PER_DAY};

let cosm = Cosm::from_xb("./de438s");
let earth = cosm.frame_by_id(bodies::EARTH);

let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 31);

let state = State::cartesian(
    -2436.45, -2436.45, 6891.037, // X, Y, Z (km)
    5.088_611, -5.088_611, 0.0, // VX, VY, VZ (km/s)
    dt, earth,
);

// Let's initialize celestial dynamics with the point masses of the Moon and the Sun
// in addition to the those of the Earth, which are included by default because the
// initial state is defined around the Earth.
let mut dynamics = CelestialDynamics::new(state, vec![bodies::EARTH_MOON, bodies::SUN], &cosm);

let opts = PropOpts::default();
println!("propagator options: {}", opts.info());

// Now let's setup the propagator.
let mut prop = Propagator::default(&mut dynamics, &opts);
// Now, let's define the channels so that the propagator can send and receive the
// output states from the CelestialDynamics.
use std::sync::mpsc;
// ^^^ We're using the standard library (std), so we can just `use` something without
// having to `extern crate` it.
let (tx, rx) = mpsc::channel();
// ^^^ This defines two channels, a sender/transmitter (tx) and receiver (rx).
// They know of each other's existance, we don't need to do anything else to link them.

// Now, let's attach the sending channel to the propagator
prop.tx_chan = Some(&tx);

// And let's prepare the receive the state on another thread.
use std::thread;
thread::spawn(move || {
    // Initialize the CSV file as states.csv
    let mut wtr = csv::Writer::from_path("./states.csv").expect("could not create file");
    // And while more states are published on the channel,
    // we'll take them and serialize them to that CSV file.
    while let Ok(rx_state) = rx.recv() {
        wtr.serialize(rx_state).expect("could not serialize state");
    }
    // No need to worry about closing the file, Rust will take care of that by itself
    // after this thread dies.
});

// And finally, let's propagate for a day.
let last_state = prop.until_time_elapsed(SECONDS_PER_DAY);
println!("{}", last_state);
```

Finally, after running this, we can confirm that the final state reported by the propagator is in the file, and so are all of the other states.
```text
chris@localhost [~/Workspace/nyxtutorial]$ cargo run
(...)
propagator options: [min_step: 1e-3, max_step: 2.7e3, tol: 1e-12, attempts: 50]
[Geoid 399] 2020-01-31T23:59:23 TAI     position = [-5971.186531, 3945.538664, 2864.607799] km  velocity = [0.049052, -4.185089, 5.848946] km/s
 chris@localhost [~/Workspace/nyxtutorial] $ tail -n 1 states.csv
2458880.500372509,-5971.186530720716,3945.538663548477,2864.607798628203,0.04905161464703625,-4.185088589338641,5.8489455915647675
 chris@localhost [~/Workspace/nyxtutorial] $ wc states.csv
  1043   1043 136732 states.csv
```

## Eclipse locators
Here, we're using the propagator channels to feed each state to an [eclipse locator](../eclipse/struct.EclipseLocator.html), and compute the [eclipse state](../celestia/eclipse/enum.EclipseState.html#variants).
In Penumbra, the closer the reported value is, the more light is received by the object.


```
extern crate nyx_space as nyx;
use nyx::celestia::{bodies, Cosm, LTCorr, State};
use nyx::dynamics::celestial::CelestialDynamics;
use nyx::propagators::{PropOpts, Propagator};
use nyx::time::{Epoch, SECONDS_PER_DAY};
use std::sync::mpsc;
use std::thread;

let cosm = Cosm::from_xb("./de438s");
let earth = cosm.frame_by_id(bodies::EARTH);

// GEO are in shadow or near shadow during the equinoxes.
let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 3, 19);

let geo_bird = State::keplerian(42000.0, 0.1, 0.1, 0.0, 0.0, 0.0, start_time, earth);

let (truth_tx, truth_rx) = mpsc::channel();

let bodies = vec![bodies::SUN, bodies::JUPITER_BARYCENTER];

// Let's move the propagation to another thread.
thread::spawn(move || {
    let cosm = Cosm::from_xb("./de438s");
    let mut dynamics = CelestialDynamics::new(geo_bird, bodies, &cosm);
    let mut prop = Propagator::default(&mut dynamics, &PropOpts::with_fixed_step(60.0));
    prop.tx_chan = Some(&truth_tx);
    prop.until_time_elapsed(2.0 * SECONDS_PER_DAY);
});

// Get a copy of the Sun to pass it to the locator.
let sun = cosm.frame_by_id(bodies::SUN);

// Import everything needed for the eclipse locator
use nyx::celestia::eclipse::{EclipseLocator, EclipseState};
// Initialize the EclipseLocator with light time correction
let e_loc = EclipseLocator {
    light_source: sun,
    shadow_bodies: vec![earth],
    cosm: &cosm,
    correction: LTCorr::LightTime,
};

// Receive the states on the main thread.
let mut prev_eclipse_state = EclipseState::Umbra;
let mut cnt_changes = 0;
while let Ok(rx_state) = truth_rx.recv() {
    // For each state, let' use the locator to compute the visibility
    let new_eclipse_state = e_loc.compute(&rx_state);
    if new_eclipse_state != prev_eclipse_state {
        // Notify the user if we've changed states
        println!(
            "{:.6} now in {:?}",
            rx_state.dt.as_jde_tai_days(),
            new_eclipse_state
        );
        prev_eclipse_state = new_eclipse_state;
        cnt_changes += 1;
    }
}

// We should get 62 eclipse state changes with 60 second time steps
// during the Spring solstice for a GEO bird.
assert_eq!(cnt_changes, 62, "wrong number of eclipse state changes");
```

Running this example, we should get the following output, where the time is in JDE days TAI:

```text
2458927.500694 now in Penumbra(0.05597248362684378)
2458927.501389 now in Umbra
2458927.514583 now in Penumbra(0.02019090399088673)
2458927.515278 now in Penumbra(0.069027462776996)
2458927.515972 now in Penumbra(0.13767706147271006)
2458927.516667 now in Penumbra(0.22336196042781828)
2458927.517361 now in Penumbra(0.3237763981611027)
2458927.518056 now in Penumbra(0.4362949211175176)
2458927.518750 now in Penumbra(0.5575536900654162)
2458927.519444 now in Penumbra(0.6830021886277537)
2458927.520139 now in Penumbra(0.8061601587878173)
2458927.520833 now in Penumbra(0.9169007484070181)
2458927.521528 now in Penumbra(0.9945152235541193)
2458927.522222 now in Visibilis
2458928.465972 now in Penumbra(0.9986185347114862)
2458928.466667 now in Penumbra(0.9277404538888158)
2458928.467361 now in Penumbra(0.8192233079898528)
2458928.468056 now in Penumbra(0.6967477897141275)
2458928.468750 now in Penumbra(0.571079475481808)
2458928.469444 now in Penumbra(0.4489997668481214)
2458928.470139 now in Penumbra(0.33523887854028916)
2458928.470833 now in Penumbra(0.23327867943810346)
2458928.471528 now in Penumbra(0.14581960299025368)
2458928.472222 now in Penumbra(0.07519295601630901)
2458928.472917 now in Penumbra(0.024077775064513206)
2458928.473611 now in Umbra
2458928.486111 now in Penumbra(0.0010844822640632245)
2458928.486806 now in Penumbra(0.06920289131705062)
2458928.487500 now in Umbra
2458928.493750 now in Penumbra(0.19451156446859705)
2458928.494444 now in Penumbra(0.04635325320369781)
2458928.495139 now in Umbra
2458928.507639 now in Penumbra(0.0014931188778925764)
2458928.508333 now in Penumbra(0.033537098417180215)
2458928.509028 now in Penumbra(0.08931751245417924)
2458928.509722 now in Penumbra(0.16392736253295587)
2458928.510417 now in Penumbra(0.2548708309066666)
2458928.511111 now in Penumbra(0.3597729360862881)
2458928.511806 now in Penumbra(0.4757890276390606)
2458928.512500 now in Penumbra(0.5991997242634797)
2458928.513194 now in Penumbra(0.7248977916058191)
2458928.513889 now in Penumbra(0.8454378076480954)
2458928.514583 now in Penumbra(0.9485260410749184)
2458928.515278 now in Visibilis
2458929.459722 now in Penumbra(0.9795249988321703)
2458929.460417 now in Penumbra(0.8894095384534675)
2458929.461111 now in Penumbra(0.7738544870636861)
2458929.461806 now in Penumbra(0.6491022932703944)
2458929.462500 now in Penumbra(0.524046157343583)
2458929.463194 now in Penumbra(0.4045660975315516)
2458929.463889 now in Penumbra(0.2948585574837825)
2458929.464583 now in Penumbra(0.1980674290020449)
2458929.465278 now in Penumbra(0.11670472085198429)
2458929.465972 now in Penumbra(0.053115299257025496)
2458929.466667 now in Penumbra(0.010614514304376692)
2458929.467361 now in Umbra
2458929.479861 now in Penumbra(0.015950304712606278)
2458929.480556 now in Penumbra(0.1170710809027578)
2458929.481250 now in Umbra
2458929.487500 now in Penumbra(0.11886667826332224)
2458929.488194 now in Penumbra(0.016638753766634026)
2458929.488889 now in Umbra
```


# Dynamical models
## Multibody dynamics
In almost all cases, we don't want to simply propagate and analyze a two body dynamics scenario (that just isn't how spacecraft fly).

Setting up a celestial dynamics with multibody point masses is quite straightforward.

```
extern crate nyx_space as nyx;
use nyx::celestia::{bodies, Cosm, State};
use nyx::dynamics::celestial::CelestialDynamics;
use nyx::propagators::{PropOpts, Propagator};
use nyx::time::{Epoch, SECONDS_PER_DAY};

let cosm = Cosm::from_xb("./de438s");
let earth = cosm.frame_by_id(bodies::EARTH);

let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 31);

let state = State::cartesian(
    -2436.45, -2436.45, 6891.037, // X, Y, Z (km)
    5.088_611, -5.088_611, 0.0, // VX, VY, VZ (km/s)
    dt, earth,
);

// Define which other masses we want.
let pt_masses = vec![bodies::EARTH_MOON, bodies::SUN];

// Let's initialize celestial dynamics with the extra point masses
// in addition to the those of the Earth, which are included by default because the
// initial state is defined around the Earth.
let mut dynamics = CelestialDynamics::new(state, pt_masses, &cosm);

let opts = PropOpts::default();
println!("propagator options: {}", opts.info());

// Now let's setup the propagator.
let mut prop = Propagator::default(&mut dynamics, &opts);
// And finally, let's propagate for a day.
let last_state = prop.until_time_elapsed(SECONDS_PER_DAY);
println!("{}", last_state);
```

## Basic attitude dynamics
Currently the attitude model is extremely limited. If you need high fidelity attitude
dynamics, we recommend using [Basilisk](http://hanspeterschaub.info/research-basilisk.html),
developed at the University of Colorado at Boulder and used for simulation and flight of
the Emirati Mars Mission.

Regardless, because nyx supports any kind of equation of motion via the Dynamics trait, there
are basic dynamics for computing the momentum over time from the inertia tensor and the angular velocity.

```
extern crate nyx_space as nyx;
use nyx::dimensions::{Matrix3, Vector3};
use nyx::dynamics::momentum::AngularMom;
use nyx::propagators::error_ctrl::LargestStep;
use nyx::propagators::{CashKarp45, PropOpts, Propagator};

// Define the angular velocity of this rigid body
let omega = Vector3::new(0.1, 0.4, -0.2);
// Define the inertia tensor of the spacecraft
let tensor = Matrix3::new(10.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 2.0);

let mut dynamics = AngularMom::from_tensor_matrix(&tensor, &omega);

// Compute the norm of the initial momentum.
// Because the dynamics are defined as `mut`, they will be overwritten
// once we move the variable to the Propagator.
let init_momentum = dynamics.momentum().norm();

// Define a CashKarp 4-5 integrator, and move the dynamics there.
let mut prop = Propagator::new::<CashKarp45>(
    &mut dynamics,
    &PropOpts::with_adaptive_step(0.1, 5.0, 1e-8, LargestStep {}),
);

// Propagate for five seconds.
prop.until_time_elapsed(5.0);

println!("{:?}", prop.latest_details());

// Compute the different in momentum from the start of the simulation to the end.
let delta_mom = ((prop.dynamics.momentum().norm() - init_momentum) / init_momentum).abs();

// Without any external torque, the momentum of the rigid body is conserved.
// So if the relative difference is high, then our propagator is broken (it isn't).
if delta_mom > 1e-8 {
    panic!(
        "angular momentum prop failed: momentum changed by {:e}",
        delta_mom
    );
}
```

## Finite burns with fuel depletion
In this example, we'll see how to setup a thrusting arc. As you will see in the last subsection of this section,
it's relatively easy to code up new dynamical models. Hence, you may want to develop your own thrusting models
better suited to your scenario.

We will need to define a spacercraft, which implement the `Dynamics` trait. However, that spacecraft requires us to define
all of the subsystems of the spacecraft prior to be able to run the propagator.

**Important documentation to read to understand this example:**

+ [Mnvr](../dynamics/thrustctrl/struct.Mnvr.html)
+ [FiniteBurns](../dynamics/thrustctrl/struct.FiniteBurns.html)
+ [ThrustControl](../dynamics/thrustctrl/trait.ThrustControl.html), the trait which defines thrusting control strategy
+ [Spacecraft](../dynamics/spacecraft/struct.Spacecraft.html)
+ [SpacecraftState](../dynamics/spacecraft/struct.SpacecraftState.html)

```
extern crate nyx_space as nyx;
use nyx::celestia::{bodies, Cosm, State};
use nyx::dimensions::Vector3;
use nyx::dynamics::celestial::CelestialDynamics;
use nyx::dynamics::propulsion::{Propulsion, Thruster};
use nyx::dynamics::spacecraft::Spacecraft;
use nyx::dynamics::thrustctrl::{FiniteBurns, Mnvr};
use nyx::dynamics::Dynamics;
use nyx::propagators::{PropOpts, Propagator};
use nyx::time::Epoch;

let cosm = Cosm::from_xb("./de438s");
let earth = cosm.frame_by_id(bodies::EARTH);

let start_time = Epoch::from_gregorian_tai_at_midnight(2002, 1, 1);

let orbit = State::cartesian(
    -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_time, earth,
);

// Define the list of thrusters (hence the `vec![]` call).
// Note that we actually only have one thruster whose thrust level is set to 10.0 Newton
// and specific impulse is set to 300 seconds.
let biprop = vec![Thruster {
    thrust: 10.0,
    isp: 300.0,
}];

// Define the propagation time
let prop_time = 3000.0;

// Define a single maneuver and its schedule (start and end epoch)
// The `thrust_lvl` determines the percentage of thrust we should provide during this maneuver
// The vector defines the orientation of the maneuver in the VNC frame.
let mnvr = Mnvr {
    start: start_time,
    end: start_time + prop_time,
    thrust_lvl: 1.0, // Full thrust
    vector: Vector3::new(1.0, 0.0, 0.0),
};

// Now, let's define a schedule, which expects a vector of maneuvers.
// In this case, we only have one maneuver, but we still need to build a vector.
let schedule = FiniteBurns::from_mnvrs(vec![mnvr], earth);

// Define the spacecraft mass
let dry_mass = 1e3;
let fuel_mass = 756.0;
let fuel_depl = true;

// Now, we can finally define the propulsion subsystem.
// It requires something which implement `ThrustControl` (the schedule),
// the list of thruster (`biprop`), and whether or not we want to
// compute the fuel depletion.
let prop_subsys = Propulsion::new(schedule, biprop, fuel_depl);

// Define the celestial dynamics
let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
let dynamics = CelestialDynamics::new(orbit, bodies, &cosm);

// And finally, we can define the spacecraft, with the celestial dynamics, the
// propulsion subsystem, and the masses.
let mut sc = Spacecraft::with_prop(dynamics, prop_subsys, dry_mass, fuel_mass);

// And setup the propagator as usual.
let mut prop = Propagator::default(&mut sc, &PropOpts::with_fixed_step(10.0));
prop.until_time_elapsed(prop_time);

println!("{}", prop.dynamics.state());
```

The output should be:
```text
[Geoid 399] 2002-01-01T00:49:28 TAI     sma = 7749.128452 km    ecc = 0.003673  inc = 63.434018 deg     raan = 135.000003 deg   aop = 154.631735 deg    ta = 95.445722 deg      1745.8028378702975 kg
```

## Sub-Optimal Control of continuous thrust
Nyx currently only provides the [Ruggiero control law](../dynamics/thrustctrl/struct.Ruggiero.html) as a suboptimal control. This example is very similar to the previous one
with the difference that we add event tracker to see whether the propagator found that we hit those targets (defined as [`Achieve`](../dynamics/thrustctrl/enum.Achieve.html) items).

Also note that we define a set of objectives for the Ruggiero control to know what it should target. The controller will automatically change the thrust direction
based on the osculating orbital state, and therefore do a local optimization of the thrust. A succession of local optimizations is a sub-optimal control.

```
extern crate nyx_space as nyx;
use nyx::celestia::{bodies, Cosm, State};
use nyx::dynamics::celestial::CelestialDynamics;
use nyx::dynamics::propulsion::{Propulsion, Thruster};
use nyx::dynamics::spacecraft::Spacecraft;
use nyx::dynamics::thrustctrl::{Achieve, Ruggiero};
use nyx::dynamics::Dynamics;
use nyx::propagators::events::{EventKind, EventTrackers, OrbitalEvent, SCEvent};
use nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use nyx::time::{Epoch, SECONDS_PER_DAY};
// Source: AAS-2004-5089

let cosm = Cosm::from_xb("./de438s");
let earth = cosm.frame_by_id(bodies::EARTH);

let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

let orbit = State::keplerian(7000.0, 0.01, 0.05, 0.0, 0.0, 1.0, start_time, earth);

let prop_time = 39.91 * SECONDS_PER_DAY;

// Define the dynamics
let dynamics = CelestialDynamics::two_body(orbit);

// Define the thruster
let lowt = vec![Thruster {
    thrust: 1.0,
    isp: 3100.0,
}];

// Define the objectives such that the control law knows what to target.
let objectives = vec![
    Achieve::Sma {
        target: 42000.0,
        tol: 1.0,
    },
    Achieve::Ecc {
        target: 0.01,
        tol: 5e-5,
    },
];

// Track the completion events
let tracker = EventTrackers::from_events(vec![
    SCEvent::orbital(OrbitalEvent::new(EventKind::Sma(42000.0))),
    SCEvent::orbital(OrbitalEvent::new(EventKind::Ecc(0.01))),
]);

// Setup the Ruggiero control law with the objectives and the initial orbit state
let ruggiero = Ruggiero::new(objectives, orbit);

let dry_mass = 1.0;
let fuel_mass = 299.0;

// Define the propulsion subsystem where the control is Ruggiero and the propulsion
// is the 1 Newton EP thruster defined above.
let prop_subsys = Propulsion::new(ruggiero, lowt, true);

let mut sc = Spacecraft::with_prop(dynamics, prop_subsys, dry_mass, fuel_mass);
println!("{:o}", orbit);

let mut prop = Propagator::new::<RK4Fixed>(&mut sc, &PropOpts::with_fixed_step(10.0));
prop.event_trackers = tracker;
prop.until_time_elapsed(prop_time);

println!("{}", prop.event_trackers);
let final_state = prop.dynamics.celestial.state();
let fuel_usage = fuel_mass - sc.fuel_mass;
println!("{:o}", final_state);
println!("fuel usage: {:.3} kg", fuel_usage);

assert!(
    sc.prop.as_ref().unwrap().ctrl.achieved(&final_state),
    "objective not achieved"
);

assert!((fuel_usage - 93.449).abs() < 1.0);
```

Note that the Ruggiero control law is a bit complex to compute, so we recommend that this example
be executed in `release` mode, which will optimize the binary.

```test
$ cargo run --release
   Compiling nyx-space v0.0.20-alpha.1 (/home/chris/Workspace/rust/nyx)
(...)
   Compiling nyxtutorial v0.1.0 (/home/chris/Workspace/rust/nyxtutorial)
    Finished release [optimized] target(s) in 8.89s
     Running `target/release/nyxtutorial`
[Geoid 399] 2019-12-31T23:59:23 TAI     sma = 7000.000000 km    ecc = 0.010000  inc = 0.050000 deg      raan = 0.000000 deg     aop = 360.000000 deg    ta = 1.000000 deg
[ERROR] Event SCEvent { kind: Sma(42000.0), orbital: Some(OrbitalEvent { kind: Sma(42000.0), tgt: None, cosm: None }) } did NOT converge
[ OK  ] Event SCEvent { kind: Ecc(0.01), orbital: Some(OrbitalEvent { kind: Ecc(0.01), tgt: None, cosm: None }) } converged on (2226230, 2226240)
[Geoid 399] 2020-02-09T21:49:47 TAI     sma = 41999.469324 km   ecc = 0.010045  inc = 0.050000 deg      raan = 0.000000 deg     aop = 218.070483 deg    ta = 89.325437 deg
fuel usage: 93.448 kg
```

## Solar radiation pressure modeling
## Basic drag models (cannonball)
## Defining new dynamical models
## New thrusting models
# Orbit determination
## Generating orbit determination measurements
## Classical and Extended Kalman Filter
## Saving the estimates
## State noise compensation (SNC)
## Smoothing and iterations of CKFs
## Square Root Information Filer (SRIF)
 */
