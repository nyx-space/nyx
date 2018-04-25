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
Earth.Energy  | 0.0 | 0.0 | 0.0
Earth.OrbitPeriod | 0.0 | 0.0 | 0.0
Earth.HX  | 7e-12  | 0.0 | 3e-12
Earth.HY  | 7e-12  | 1e-16 | 0.0
Earth.HZ  | 0.0  | 0.0 | 0.0
Earth.SMA  | 0.0  | 0.0 | 0.0
Earth.ECC  |  0.0 | 1e-17 | 5e-17
EarthMJ2000Eq.INC  | 0.0 | 0.0  | 0.0
EarthMJ2000Eq.RAAN  | 0.0  | 0.0  | 0.0
EarthMJ2000Eq.AOP  | 0.0 | 3e-7 | 5e-11
Earth.TA  | 0.0 | 3e-7 | 4e-14
Earth.TLONG | 0.0 | 3e-1<sup>(3)</sup> | 1e-8<sup>(1)</sup>
Earth.EA | 0.0 | 0.0 | 0.0
Earth.MA | 0.0 | 0.0 | 0.0
Earth.RadApo | 0.0 | 0.0 | 0.0
Earth.RadPer | 0.0 | 0.0 | 0.0
Earth.SemilatusRectum | 2e-12 | 0.0 | 2e-12
### Footnotes
(1) There is quite a large error in true longitude for this test. I am really not sure why given that `nyx` sums AoP, RAAN and TA to compute this, as per the definition. Summing these values leads _exactly_ to the value returned by `nyx`. I am very surprised that GMAT does not seem to use that same formula, I'll have to check why.

(3) This error _should_ be zero, but for transparency it's marked as 1e-1. Using the same methodology as (1), the script was created using specific angles. These angles were all defined to within `1e-1`, and switching between a Keplerian and a Cartesian definition created rounding errors themselves within in the GMAT GUI. I admit being surprised that the TLONG returned by GMAT is precisely `159.6` since GMAT has a tendency to not round even floating-point computational errors.

## From a Keplerian state

Element / Scenario  | circular inclined  | circular equatorial  | elliptical |  
--|---|---|---|---|---|---|--
Earth.X  |  |  | 0.0
Earth.Y  |  |  | 4e-13
Earth.Z  |  |  | 3e-12
Earth.VX  |   |  | 1e-15
Earth.VY  |   |  | 2e-15
Earth.VZ  | |  | 1e-15
Earth.Energy  |  |  | 1e-14
Earth.OrbitPeriod |  |  | 5e-12
Earth.HX  |   |  | 0.0
Earth.HY  |   |  | 1e-11
Earth.HZ  |   |  | 2e-11
Earth.SMA  |  |  | 3e-12
Earth.ECC  |  |  | 3e-16
EarthMJ2000Eq.INC  |  |  | 1e-14
EarthMJ2000Eq.RAAN  |   |  | 0.0
EarthMJ2000Eq.AOP  |  |  | 1e-12
Earth.TA  |  |  | 1e-12
Earth.TLONG |  | | 1e-13
Earth.EA |  |  | 1e-12
Earth.MA |  |  | 1e-12
Earth.RadApo |  |  | 1e-12
Earth.RadPer |  |  | 7e-12
Earth.SemilatusRectum |  |  | 0.0
