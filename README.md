# nyx
[Nyx](https://en.wikipedia.org/wiki/Nyx) is a high fidelity astrodynamical toolkit library written in Rust.
It will provide most functionality in C and will be usable in Python via a C-wrapper for rapid prototyping.

The goal is to provide researchers and astrodynamicists a fast, reliable and tested toolkit to rapidly prototype simulations
in Python, and use the safety of Rust for significantly faster repetitive simulation runs. To some extend, the ultimate goal of 
this library is to retire SPICE Toolkit and its analogs.

# License
The [LICENSE](https://github.com/ChristopherRabotin/nyx/blob/master/LICENSE) will be strictly enforced once (if?) this toolkit
reaches production-level quality.

# Features
Some of the following features may only be made available only through a commercial license.

- [ ] Multibody dynamics using SPICE Kernels
- [ ] Planetary and Solar eclipse and visibility computation
- [ ] Statistical Orbital Determination: Classical and Extended Kalman Filter
- [ ] Finite burns with fuel depletion (includes low thrust / ion propulsion)
- [ ] Monte Carlo simulations on different parameters
- [ ] Propagation with different stopping conditions
- [ ] Sub-Optimal Control of continuous thrust (Ruggerio, Naasz, Petropoulos)
- [ ] Link budget computations 
- [ ] Spacecraft attitude control and algorithms

# Who am I?
A happily-employed astrodynamics engineer with a heavy background in software. Nyx will rely on the fallbacks of 
[smd](https://github.com/ChristopherRabotin/smd), a library I wrote in Go while researching at the University
of Colorado at Boulder, and on my gained knowledge at Advanced Space ([check out the cool stuff we do!](http://advanced-space.com/)).
