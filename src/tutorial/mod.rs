/*!

This tutorial covers some of the basics of using nyx and its builtin features.
Nyx is highly configurable, and used at Advanced Space for much more than what's being shown here.
For example, it is relatively easy to create new dynamics which leverage those available in nyx with new equations of motion;
or the definition of new `MeasurementDevice`s which account for specific clock errors.

# Table of contents
0. [Setup](#setup)
1. [Celestial computations](#celestial-computations)
    - Orbital state manipulation
    - Planetary and Solar eclipse and visibility computation
    - Light-time corrections and abberations
2. [Propagation](#propagation)
    - Propagation with different Runge Kutta methods
    - Propagation to different stopping conditions
3. [Dynamical models](#dynamical-models)
    - Multibody dynamics using XB files
    - Basic attitude dynamics
    - Finite burns with fuel depletion
    - Sub-Optimal Control of continuous thrust
    - Solar radiation pressure modeling
    - Basic drag models (cannonball)
4. [Orbit determination](#orbit-determination)
    - Classical and Extended Kalman Filter
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
`nyx_space = "0.0.19"` to your `[dependencies]` section. Also add `hifitime = "1"`,
which is the high fidelity time computations library nyx relies on.
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
nyx-space = "0.0.19"
hifitime = "1"
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

Orbital states are stored in an [`OrbitState`](../celestia/type.OrbitState.html) structure. It has a lot of useful methods, like computing the angular momentum of that orbit, or its period.
It can be initialized from Cartesian coordinates or Keplerian elements, but the data is stored in Cartesian. Therefore, a state initialized from Keplerian elements might not return exactly the same input values you specified.


```
// Import hifitime and nyx
extern crate hifitime;
extern crate nyx_space as nyx;
// Tell Rust that we're using the structures called Epoch, Cosm and OrbitState
// from the respective modules. We're also importing _bodies_ which has a mapping from
// the name of the object to its ID in the EXB file.
use hifitime::Epoch;
use nyx::celestia::{bodies, Cosm, OrbitState};

// Load both the de438s.exb and de438s.fxb (resp. the ephemeris and frame information).
let cosm = Cosm::from_xb("./de438s");
// Get a copy of the Earth body and frame information from the cosm.
let earth = cosm.geoid_from_id(bodies::EARTH);
// Initialize a Epoch for 31 January 2020 at midnight TAI.
let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 31);


// And initialize a Cartesian state with position, velocity, epoch and center object.
// The position units are in kilometers and the velocity units in kilometers per second.
let cart = OrbitState::cartesian(
    -2436.45, -2436.45, 6891.037,  // X, Y, Z (km)
    5.088_611, -5.088_611, 0.0,    // VX, VY, VZ (km/s)
    dt, earth,
);

// And print this state in Cartesian
println!("{}", cart);
// And print this state in Keplerian
println!("{:o}", cart);

// Now let's recreate the state but from the Keplerian elements

let kep = OrbitState::keplerian(
    cart.sma(),
    cart.ecc(),
    cart.inc(),
    cart.raan(),
    cart.aop(),
    cart.ta(),
    dt,
    earth
);

// This should be zero
println!("{}", cart - kep);
```

## Planetary and Solar eclipse and visibility computation
## Light-time corrections and abberations
# Propagation
## Propagation with different Runge Kutta methods
## Propagation to different stopping conditions
# Dynamical models
## Multibody dynamics using XB files
## Basic attitude dynamics
## Finite burns with fuel depletion
## Sub-Optimal Control of continuous thrust
## Solar radiation pressure modeling
## Basic drag models (cannonball)
# Orbit determination
## Classical and Extended Kalman Filter
## State noise compensation (SNC)
## Smoothing and iterations of CKFs
## Square Root Information Filer (SRIF)
 */
