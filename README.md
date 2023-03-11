# nyx
[Nyx](https://nyxspace.com) is a high fidelity, fast, reliable and **[validated]([https://nyxspace.com/MathSpec/](https://nyxspace.com/nyxspace/MathSpec/))** astrodynamics toolkit library written in Rust. It provides convenient interfaces to both Rust and Python.

The target audience is mission designers, astrodynamics engineers/hobbyists, and GNC engineers. The rationale for using Rust is to allow for very fast computations, guaranteed thread safety,
and portability to all platforms supported by [Rust](https://forge.rust-lang.org/platform-support.html).

[![nyx-space on crates.io][cratesio-image]][cratesio]
[![nyx-space on docs.rs][docsrs-image]][docsrs]
[![LoC](https://tokei.rs/b1/gitlab/nyx-space/nyx?category=lines)](https://github.com/nyx-space/nyx).

[cratesio-image]: https://img.shields.io/crates/v/nyx-space.svg
[cratesio]: https://crates.io/crates/nyx-space
[docsrs-image]: https://docs.rs/nyx-space/badge.svg
[docsrs]: https://docs.rs/nyx-space/

# License
The [AGPLv3 LICENSE](https://nyxspace.com/license/) is enforced.

# Features

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
- [x] Light-time corrections and aberrations
- [x] Frame rotations

# Who am I?
An astrodynamics engineer with a heavy background in software. I currently work for Rocket Lab USA on the GNC of the Blue Ghost lunar lander.

# Examples
Refer to the [showcase](https://nyxspace.com/showcase/).

# Python todo

The "[maturin](https://crates.io/crates/maturin)" python package is used to build the python bindings.

```
pip install maturin
```

Build the python bindings using the following command.
```
maturin build --cargo-extra-args="--features python"
```

This creates a wheel file in `./target/wheels/` which can be installed using `pip install <filename.whl>`.

For development mode, the following command may be used that automatically installs the python module

```
maturin develop --cargo-extra-args="--features python"
```

## Python code example

This minimal example runs the scenario defined in `data/simple-scenario.toml` using the Python bindings.

```
from nyx_space import md
from nyx_space import io
from nyx_space import cosmic

# Initialize the cosm which stores the ephemeris
cosm = cosmic.Cosm.de438()

with open('data/simple-scenario.toml', 'r') as f:
    scen_data = f.read()

scenario = io.ScenarioSerde.from_toml_str(scen_data)

md.MDProcess.execute_all_in_scenario(scenario, cosm)

```
