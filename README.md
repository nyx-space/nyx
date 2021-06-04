# nyx
[Nyx](https://en.wikipedia.org/wiki/Nyx) is a high fidelity, fast, reliable and **[validated](https://nyxspace.com/MathSpec/)** astrodynamical toolkit library written in Rust.

The target audience is researchers and astrodynamics engineers. The rationale for using Rust is to allow for very fast computations, guaranteed thread safety,
and portability to all platforms supported by [Rust](https://forge.rust-lang.org/platform-support.html).

[![nyx-space on crates.io][cratesio-image]][cratesio]
[![nyx-space on docs.rs][docsrs-image]][docsrs]

[cratesio-image]: https://img.shields.io/crates/v/nyx-space.svg
[cratesio]: https://crates.io/crates/nyx-space
[docsrs-image]: https://docs.rs/nyx-space/badge.svg
[docsrs]: https://docs.rs/nyx-space/

# License
The [AGPLv3 LICENSE](https://nyxspace.com/license/) is strictly enforced.

# Features
Unless specified otherwise in the documentation of specific functions, all vectors and matrices are [statically allocated](https://discourse.nphysics.org/t/statically-typed-matrices-whose-size-is-a-multiple-or-another-one/460/4).

Lots of features are still being worked on, and there currently isn't any guarantee that the API won't change _between_ versions. However, you can be assured that the API will not change for previous versions.
Outstanding mission design features available [here](https://gitlab.com/chrisrabotin/nyx/-/issues?label_name=subsys%3A%3AMD), and orbit determination features [here](https://gitlab.com/chrisrabotin/nyx/-/issues?scope=all&utf8=%E2%9C%93&state=opened&label_name[]=subsys%3A%3AOD).

## Propagation
- [x] Propagation with different Runge Kutta methods (validated in GMAT)
- [x] Convenient and explicit definition of the dynamics for a simulation (cf. [tests/orbitaldyn.rs](tests/orbitaldyn.rs))
- [x] Propagation to different stopping conditions
- [x] Detect orbital events in other frames
## Dynamical models
- [x] Multibody dynamics using XB files
- [x] Finite burns with fuel depletion (including low thrust / ion propulsion)
- [x] Sub-Optimal Control of continuous thrust (e.g. Ruggerio, Petropoulos/Q-law)
- [x] Solar radiation pressure modeling
- [x] Basic drag models (cannonball)
- [x] Spherical harmonics
- [ ] Spacecraft attitude control and some useful optimal control algorithms
## Orbit determination
- [x] Statistical Orbit Determination: Classical and Extended Kalman Filter
- [x] Orbit Determination with multibody dynamics
- [x] Smoothing and iterations of CKFs
- [ ] Square Root Information Filer (SRIF) (removed in version 1.0.0-alpha.1)
- [x] An easy-to-use OD user interface
- [x] Estimation with spherical harmonics enabled
- [ ] Solar radiation pressure (SRP) parameter estimation ([#98](https://gitlab.com/chrisrabotin/nyx/issues/98))
- [x] Covariance mapping and estimate frame transformations
- [x] State noise compensation (SNC)
- [ ] Dynamic model compensation (DMC) ([#86](https://gitlab.com/chrisrabotin/nyx/issues/86))
- [x] High fidelity ground station placement
## Celestial computations
- [x] Orbital state manipulation
- [x] Planetary and Solar eclipse and visibility computation
- [x] Light-time corrections and abberations
- [x] Frame rotations

# Who am I?
An astrodynamics engineer with a heavy background in software. Nyx relies on the drawbacks of
[smd](https://github.com/ChristopherRabotin/smd), a library I wrote in Go while researching at the University
of Colorado at Boulder. I work for Masten Space Systems on the XL-1 Moon Lander ([we do really cool stuff](https://masten.aero/)).

# Examples
Refer to the [showcase](https://nyxspace.com/showcase/).