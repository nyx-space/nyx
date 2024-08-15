# Orbit Determination of the Lunar Reconnaissance Orbiter

Spacecraft operations require high fidelity modeling of the orbital dynamics and high fidelity orbit determination.

In this example, you'll learn how to use an "as-flown" (_definitive_) SPICE BSP ephemeris file to simulate orbit determination measurements from ground stations. Then, you'll learn how to set up an orbit determination process in Nyx with high fidelity Moon dynamics and estimate the state of LRO. Finally, you'll learn how to compare two ephemerides in the Radial, In-track, Cross-track (RIC) frame.

To run this example, just execute:
```sh
RUST_LOG=info cargo run --example 04_lro_od --release
```

Building in `release` mode will make the computation significantly faster. Specifying `RUST_LOG=info` will allow you to see all of the information messages happening in ANISE and Nyx throughout the execution of the program.
