# GMAT_Script
Many of the features of this library are validated against [GMAT](https://software.nasa.gov/software/GSC-17177-1) ([wiki](http://gmatcentral.org/display/GW/GMAT+Wiki+Home)), version 2017a.
GMAT is validated on flown missions. It was also validated against other software, cf. [this PDF](./GMAT_V&V_ProcessAndResults.pdf).

# Propagators
The purpose of this test is solely to test the correct implementation of the propagator coefficients, error computation, and adaptive step size. The algorithms were taken from GMAT unless noted otherwise in the source code.
The tests in [`propagators.rs`](../propagators.rs) we executed in GMAT as well. Each script is in the subfolder [propagators](./propagators/) and is named after the propagator used.

## Results
The following table corresponds to the errors between the `nyx` output for the one day GEO propagation in two body dynamics and the GMAT output using the same step configuration.

Propagator  | x | y | z | vx | vy |  vz
--|---|---|---|---|---|--
RK4Fixed  | 1.4e-9 | 3.0e-8 | 4.0e-8 | 3.6e-11 | 2.4e-11 | 1.7e-11
Dormand45  | 2.1e-9 | 3.1e-7 | 4.5e-7 | 3.9e-10 | 2.6e-10 | 1.9e-10
Verner56  | 2.5e-9 | 8.6e-7 | 1.2e-6 | 1.1e-9 | 7.1e-10 | 5.2e-10
Dormand87  | 2.5e-9 | 9.6e-9 | 1.3e-8 | 1.2e-11 | 7.9e-12 | 5.6e-12
RK89  | 3.8e-10 | 7.7e-9 | 1.0e-8 | 9.4e-12 | 6.4e-12 | 4.4e-12
