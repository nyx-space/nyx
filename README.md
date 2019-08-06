# nyx
[Nyx](https://en.wikipedia.org/wiki/Nyx) is a high fidelity, fast, reliable and validated astrodynamical toolkit library written in Rust.
It will _eventually_ provide most functionality in Python for rapid prototyping.

The target audience is researchers and astrodynamics engineers. The rationale for using Rust is to allow for very fast computations, guaranteed thread safety,
and portability to all platforms supported by [Rust](https://forge.rust-lang.org/platform-support.html).

[![nyx-space on crates.io][cratesio-image]][cratesio]
[![nyx-space on docs.rs][docsrs-image]][docsrs]

[cratesio-image]: https://img.shields.io/crates/v/nyx-space.svg
[cratesio]: https://crates.io/crates/nyx-space
[docsrs-image]: https://docs.rs/nyx-space/badge.svg
[docsrs]: https://docs.rs/nyx-space/

# License
The [LICENSE](https://github.com/ChristopherRabotin/nyx/blob/master/LICENSE) will be strictly enforced once/if this toolkit
reaches production-level quality.

# Features
- [x] Propagation with different Runge Kutta methods (validated in GMAT)
- [x] Convenient and explicit definition of the dynamics for a simulation (cf. [tests/lib.rs](tests/lib.rs))
- [x] Orbital state manipulation (from GMAT source code and validated in GMAT)
- [x] Statistical Orbit Determination: Classical and Extended Kalman Filter
- [x] Multibody dynamics using XB files (caveat:https://gitlab.com/chrisrabotin/nyx/issues/61)
- [ ] Orbit Determination with multibody dynamics
- [ ] Planetary and Solar eclipse and visibility computation
- [ ] Finite burns with fuel depletion (including low thrust / ion propulsion)
- [ ] Monte Carlo simulations on different parameters
- [ ] Sub-Optimal Control of continuous thrust (Ruggerio, Naasz, Petropoulos)
- [ ] Link budget computations
- [ ] Spacecraft attitude control and some useful optimal control algorithms

_Note:_ Some of these features may only be made available only through a commercial license in the future.

# Who am I?
A astrodynamics engineer with a heavy background in software. Nyx relies on the fallbacks of
[smd](https://github.com/ChristopherRabotin/smd), a library I wrote in Go while researching at the University
of Colorado at Boulder. I work for Advanced Space ([we do cool stuff](http://advanced-space.com/)), but this code is mostly developed on my leisure time.
