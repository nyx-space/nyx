# GMAT_Script
Many of the features of this library are validated against [GMAT](https://software.nasa.gov/software/GSC-17177-1) ([wiki](http://gmatcentral.org/display/GW/GMAT+Wiki+Home)), version 2017a.
GMAT is validated on flown missions. It was also validated against other software, cf. [this PDF](./GMAT_V&V_ProcessAndResults.pdf).

# Propagators
The purpose of this test is solely to test the correct implementation of the propagator coefficients, error computation, and adaptive step size. The algorithms were taken from GMAT unless noted otherwise in the source code.
The tests in [`propagators.rs`](../propagators.rs) we executed in GMAT as well. Each script is in the subfolder [propagators](./propagators/) and is named after the propagator used.

The following table corresponds to the **errors** between the `nyx` output for the one day GEO propagation in two body dynamics and the GMAT output using the same step configuration.

Propagator  | x | y | z | vx | vy |  vz
--|---|---|---|---|---|--
RK4Fixed  | 1.4e-9 | 3.0e-8 | 4.0e-8 | 3.6e-11 | 2.4e-11 | 1.7e-11
Dormand45  | 2.1e-9 | 3.1e-7 | 4.5e-7 | 3.9e-10 | 2.6e-10 | 1.9e-10
Verner56  | 2.5e-9 | 8.6e-7 | 1.2e-6 | 1.1e-9 | 7.1e-10 | 5.2e-10
Dormand87  | 2.5e-9 | 9.6e-9 | 1.3e-8 | 1.2e-11 | 7.9e-12 | 5.6e-12
RK89  | 3.8e-10 | 7.7e-9 | 1.0e-8 | 9.4e-12 | 6.4e-12 | 4.4e-12

# Orbital state
This validation compares the computations of orbital elements in nyx with those in GMAT. Each scenario script is in the subfolder [state](./state/).
The validation is trivially the following: a spacecraft is defined in GMAT in the EarthMJ2000Eq frame by its Cartesian state. There is a "Propagate" segment of strictly zero seconds. The desired orbital state computations are exported to a report file.

The following table corresponds to the **errors** between the `nyx` computations and those of GMAT. Note that if `nyx` returned more significant digits than GMAT and all the GMAT digits matched those of nyx, the error is also marked as zero (this was only seen for the energy and orbital period where GMAT returns 13 digits and `nyx` 14).

## From a Cartesian state

Element / Scenario  | circular inclined  | circular equatorial  | elliptical |  
--|---|---|---|---|---|---|--
Earth.Energy  | 0.0 |   | 0.0
Earth.OrbitPeriod | 0.0 |   | 0.0
Earth.HX  | 7e-12  |   | 3e-12
Earth.HY  | 7e-12  |   | 0.0
Earth.HZ  | 0.0  |   | 0.0
Earth.SMA  | 0.0  |   | 0.0
Earth.ECC  |  0.0 |   | 1e-3<sup>(1)</sup>
EarthMJ2000Eq.INC  | 0.0 |   | 0.0
EarthMJ2000Eq.RAAN  | 0.0  |   | 0.0
EarthMJ2000Eq.AOP  | 0.0 |   | 5e-11
Earth.TA  | 0.0 |   | 4e-14
Earth.TLONG  | 0.0 |   | 1e-8<sup>(2)</sup>

### Footnotes
(1) This error _should_ be zero, but for transparency it's marked as 1e-3. The way this test was generated was by creating a given state using Keplerian orbital elements, and converting them to Cartesian in the GMAT spacecraft panel, and ensuring that `nyx` could convert the provided Cartesian state to the output state of GMAT. In this instance, `nyx` returns `0.15899999999999995` which is awfully close the `0.159` which I entered in GMAT in the first place.

(2) There is quite a large error in true longitude for this test. I am really not sure why given that `nyx` sums AoP, RAAN and TA to compute this, as per the definition. Summing these values leads _exactly_ to the value returned by `nyx`. I am very surprised that GMAT does not seem to use that same formula, I'll have to check why.

## From a Keplerian state

Element / Scenario  | circular inclined  | circular equatorial  | elliptical |  
--|---|---|---|---|---|---|--
Earth.Energy  |  |   |
Earth.OrbitPeriod |  |   |
Earth.HX  |   |   |
Earth.HY  |   |   |
Earth.HZ  |   |   |
Earth.SMA  |   |   |
Earth.ECC  |  |   |
EarthMJ2000Eq.INC  |  |   |
EarthMJ2000Eq.RAAN  |  |   |
EarthMJ2000Eq.AOP  | |   |
Earth.TA  | |   |
