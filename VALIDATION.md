# Validation
**WARNING:** This document is outdated. Please refer to the [MathSpec on nyxspace.com](https://nyxspace.com/MathSpec/) for the latest validation information.

The validation process ensures that the dynamics are correctly implemented.

Many of the features of this library are validated against [GMAT](https://software.nasa.gov/software/GSC-17177-1) ([wiki](http://gmatcentral.org/display/GW/GMAT+Wiki+Home)), version 2017a.
GMAT is validated on flown missions. It was also validated against other software, cf. [this PDF](./tests/GMAT_scripts/GMAT_V&V_ProcessAndResults.pdf).

- [Propagators](#propagators)
  * [With a fixed step](#fixed-step)
  * [Adaptive step](#adaptive-step)
- [Orbital state](#orbital-state)
  * [From a Cartesian state](#from-a-cartesian-state)
  * [From a Keplerian state](#from-a-keplerian-state)
  * [Frame transformations](#frame-transformations)
- [Multibody dynamics](#multibody-dynamics)
- [Propulsion](#propulsion)
- [SRP](#srp)
- [Orbit determination](#orbit-determination)

# Propagators
The purpose of this test is solely to test the correct implementation of the propagator coefficients, error computation, and adaptive step size. The algorithms were taken from GMAT unless noted otherwise in the source code.
The tests in [`propagators.rs`](./tests/propagators.rs) we executed in GMAT as well. Each script is in the subfolder [propagators](./tests/GMAT_scripts/propagators/) and is named after the propagator used. In addition to the GMAT tests, backward propagation is confirmed for two body dynamics using an adaptive step-size: all propagation tests are verified to return to their initial states with an error in position less than 1e-5 km and an error in velocity less than 1e-8 km/s.

The validation scenario is simply a two body propagation around Earth for one day.
The following tables corresponds to the **errors** between the `nyx` output for the one day LEO propagation in two body dynamics and the GMAT output using the same step configuration.

## Adaptive step

Step size configured as follows:
+ Minimum step: 0.1 seconds
+ Maximum step: 30.0 seconds
+ Accuracy: 1e-12
+ Max attempts: 50

*Note:* GMAT uses the secant method to determine the final step to take (as a fixed step). In this validation effort, we compute the exact number of seconds by which we overshot the expected elapsed time, and set the final step size to precisely that value. As such, the following table has the extra column `Δtime` which shows the difference in the propagation duration between GMAT and `nyx`. As explained, the latter will propagate for _exactly_ 86400 seconds (to machine precision), but GMAT doesn't, and this is mostly to blame for the (slight) differences in final states.

Propagator  | x | y | z | vx | vy | vz | Δtime
--|---|---|---|---|---|---|--
Dormand45  | 3.0e-07 | 9.2e-06 | 1.3e-05 | 1.1e-08 | 7.5e-09 | 5.3e-09 | 0.0
Verner56  | 6.4e-09 | 1.7e-07 | 2.5e-07 | 2.3e-10 | 1.5e-10 | 1.1e-10 | 1e-11
Dormand78  | 2.5e-10 | 9.7e-09 | 1.3e-08 | 1.2e-11 | 7.9e-12 | 5.6e-12 | 1e-11
RK89  | 3.9e-10 | 7.7e-09 | 1.0e-08 | 2.9e-12 | 6.4e-12 | 4.4e-12 | 0.0

## Fixed step

Step size set to 10.0 seconds.

Propagator  | x | y | z | vx | vy |  vz
--|---|---|---|---|---|--
RK4Fixed  | 1.4e-9 | 3.0e-8 | 4.0e-8 | 3.6e-11 | 2.4e-11 | 1.7e-11
Dormand45  | 2.9e-08 | 8.5e-07 | 1.2e-06 | 1.0e-09 | 7.0e-10 | 4.9e-10
Verner56  | 3.8e-10 | 3.9e-09 | 5.3e-09 | 4.8e-12 | 3.4e-12 | 2.1e-12
Dormand78  | 3.9e-10 | 1.2e-08 | 1.7e-08 | 1.5e-11 | 1.0e-11 | 7.0e-12
RK89  | 3.8e-11 | 2.3e-09 | 3.2e-09 | 2.9e-12 | 1.9e-12 | 1.5e-12

# Orbital state
This validation compares the computations of orbital elements in nyx with those in GMAT. Each scenario script is in the subfolder [state](./tests/GMAT_scripts/state/).
The validation is trivially the following: a spacecraft is defined in GMAT in the EarthMJ2000Eq frame by its Cartesian state. There is a "Propagate" segment of strictly zero seconds. The desired orbital state computations are exported to a report file.

The following table corresponds to the **absolute errors** between the `nyx` computations and those of GMAT. Note that if `nyx` returned more significant digits than GMAT and all the GMAT digits matched those of nyx, the error is also marked as zero (this was only seen for the energy and orbital period where GMAT returns 13 digits and `nyx` 14).

## From a Cartesian state

Element / Scenario  | circular inclined  | circular equatorial  | elliptical
--|---|---|--
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
Earth.TLONG | 0.0 | 3e-14 | **1e-8**<sup>(1)</sup>
Earth.EA | 0.0 | 0.0 | 0.0
Earth.MA | 0.0 | 0.0 | 0.0
Earth.RadApo | 0.0 | 0.0 | 0.0
Earth.RadPer | 0.0 | 0.0 | 0.0
Earth.SemilatusRectum | 2e-12 | 0.0 | 2e-12

### Footnotes
(1) There is quite a large error in true longitude for this test. I am really not sure why given that `nyx` sums AoP, RAAN and TA to compute this, as per the definition. Summing these values leads _exactly_ to the value returned by `nyx`. I am very surprised that GMAT does not seem to use that same formula, I'll have to check why.

## From a Keplerian state

Element / Scenario  | circular inclined  | circular equatorial  | elliptical
--|---|---|--
Earth.X  | 0.0 | 0.0 | 0.0
Earth.Y  | 0.0 | 1e-14 | 4e-13
Earth.Z  | 0.0 | **5e-5**<sup>(1)</sup> | 3e-12
Earth.VX  | 0.0 | 1e-15 | 1e-15
Earth.VY  | 0.0 | 1e-8<sup>(2)</sup> | 2e-15
Earth.VZ  | 0.0 | 0.0 | 1e-15
Earth.Energy  | 0.0 | 0.0 | 1e-14
Earth.OrbitPeriod | 0.0 | 1e-11 | 5e-12
Earth.HX  | 0.0 | **2e-4**<sup>(3)</sup> | 0.0
Earth.HY  | 0.0 | **1e-3**<sup>(3)</sup> | 1e-11
Earth.HZ  | 0.0 | **0.0**<sup>(3)</sup> | 2e-11
Earth.SMA  | 0.0 | 1e-11 | 3e-12
Earth.ECC  | 0.0 | 2e-16 | 3e-16
EarthMJ2000Eq.INC | 3e-15 | 0.0 | 1e-14
EarthMJ2000Eq.RAAN | 0.0 | 0.0 | 0.0
EarthMJ2000Eq.AOP | 0.0 | 3e-8 | 1e-12
Earth.TA  | 0.0 | 0.0 | 1e-12
Earth.TLONG | 9e-14 | 0.0 | 1e-13
Earth.EA | 0.0 | 3e-8 | 1e-12
Earth.MA | 0.0 | 3e-8 | 1e-12
Earth.RadApo | 0.0 | 2e-11 | 1e-12
Earth.RadPer | 0.0 | 0.0 | 7e-12
Earth.SemilatusRectum | 1e-12 | 0.0 | 0.0

### Footnotes
(1) This is one heck of an oddity! The error is quite large, although the algorithm used here is strictly a conversion of GMAT's C++ code to Rust. I _suspect_ this error is due to the orbit being circular equatorial with a value of Z of **2.5 meters** according to `nyx` and **3.0 meters** according to GMAT. As noted above, `nyx` often returns one to two additional digits than GMAT does in its report files, so maybe does GMAT suffer from multiplying several rounding errors, or could it be some shortcuts that Rust/`nyx` make?

(2) Similarly to (1), we get a large error in the velocity Z component. In this case, GMAT reports 5.9e-8 km/s and `nyx` computes 4.9e-8 km/s.

(3) Similarly to (1), we get a very significant error in the orbital momentum computation of both HX and HY. These components are small for the orbital momentum (both on the order of 1e-3 in GMAT and in `nyx`). I am not too concerned about these differences given that the orbital momentum component of the Z axis is exactly that returned by GMAT (all 16 digits are equal).

## Frame transformations
The frame transformation test, with methods to generate the test data, are in `src/celestia/cosm.rs` in the `test` module (starting with line `#[cfg(test)]`). Tests either compare with exact data and algorithm from [jplephem](https://github.com/brandon-rhodes/python-jplephem), SPICE or [JPL HORIZONS](https://ssd.jpl.nasa.gov/horizons.cgi), respectively for the "direct", "indirect" and "frame transformation" tests. The largest errors are found to be when comparing with SPICE, and these differences are similar to those (currently) found in [jplephem](https://github.com/brandon-rhodes/python-jplephem/issues/33).

Until version 0.0.19, the frame transformation _do not_ account for orientation transformation, only center of frame transformation.

# Multibody dynamics

*Note:* In the LLO scenarios, if propagation start on 2020 Jan 01 midnight TAI, the GMAT results and nyx results vary significantly. Nyx relies on the `hifitime` library for date time conversions, which has its own set of thorough validation. Hence, I have yet to find a good explanation for why the validation fails in 2020. In the following table, to avoid any misunderstandings, the start time of the one-day propagation is added.

The root mean squared errors in position and velocity are as follows:

The following scenario permutations are tested:
+ Halo orbit, LLO or LEO
+ RK8 10s time step, or RK89 with adaptive step size
+ Earth & Moon point masses, or Earth, Moon, Sun and Jupiter

Position errors are in **kilometers**, and velocity errors in **kilometers per second**.

Orbit | Adaptive/Fixed | Point masses | Prop. start date | x | y | z | vx | vy | vz | RSS position error | RSS velocity error
--|---|---|---|---|---|---|---|---|---|---|--
Halo | Fixed | Earth Moon  | 2020-01-01 | 3e-7  |  6e-7  |  1e-7  |  5e-12  | 9e-12  | 1e-12 | **6.71660e-7** | **1.08338e-11**
Halo | Adaptive | Earth Moon | 2002-02-07 | 5e-9  |  3e-7  |  1e-7  |  9e-13 |  6e-14 |  8e-16 | **3.24909e-7** | **9.18996e-13**
Halo | Fixed | Earth Moon Sun Jupiter | 2020-01-01 | 3e-7  |  6e-7  |  1e-7  |  5e-12 |  9e-12 |  1e-12 | **6.67462e-7** | **1.08000e-11**
Halo | Adaptive | Earth Moon Sun Jupiter | 2002-02-07 | 4e-9  |  2e-7 |   6e-8  |  5e-13 |  3e-1 |  0e0 | **1.79250e-7** | **5.08201e-13**
LLO | Adaptive | Earth Moon | 2002-02-07 | 2e-7  |  3e-7  |  8e-8  |  6e-13  |  2e-13 | 5e-14 | **3.70439e-7** | **6.61018e-13**
LLO | Adaptive | Earth Moon Sun Jupiter | 2002-02-07 | 4e-7  |  9e-7  |  2e-7  |  1e-12 |  4e-13 |  1e-13 |  **9.86872e-7** | **1.42994e-12**
LEO | Adaptive | Earth Moon Sun Jupiter | 2020-01-01 | 2e-8  |  2e-6  |  2e-6  |  2e-9  |  1e-9  |  9e-10 | **2.63186e-6** | **2.45434e-9**
LEO | Adaptive | Earth Sun Jupiter | 2020-01-01 | 3e-9  |  3e-7  |  4e-7  |  3e-10 |  2e-10 |  2e-10 | **4.80273e-7** | **4.47932e-10**

# Propulsion
## Finite burns
In both cases, we take a LEO spacecraft subjects to the point mass gravity of the Moon, the Sun and Jupiter. We set the dry mass to 1.0 ton/megagram/"metric ton" and the fuel mass to 756 kg. The spaceraft is equipped with a single thrusted whose thrust is 10 Newton and Isp 300 seconds. The finite burn is set to last 3000 seconds (50 minutes).

Mass depletion | Propagator | x | y | z | vx | vy | vz | RSS position error | RSS velocity error
--|---|---|---|---|---|---|---|---|--
Disabled | RK8 Fixed step | 5e-11 | 1e-10 | 3e-11 | 5e-14 | 3e-14 | 5e-14 | **1.35472e-10** | **7.54351e-14**
Enabled | RK8 Fixed step | 1e-11 | 9e-11 | 8e-12 | 4e-15 | 2e-14 | 2e-14 | **9.57349e-11** | **2.79494e-14**

# SRP
Solar radiation pressure using a cannonball model (i.e, spherical spacecraft approximation).

The difference between nyx and GMAT is **much** larger here than for other models: **1.83 km** in 24 hours! I assume that this is due to the ephemeris having a larger impact here than in multibody dynamics. **However**, the current frame transformations seem to be very good.

In fact, in the validation case used here, the spacecraft is always in full visibility of the Sun, and therefore the difference in Sun illumination does not matter. The propagation uses an RK89 and propagates for 24 hours.

However, note that nyx uses the IAU definition of an astronomical unit, whereas GMAT uses a pre-2012 definition of 1 AU (which is 19 meters shorter).

Note that further differences _may_ exist as GMAT uses a simpler algorithm for penumbra computations than nyx. The algorithm for nyx can be found in the [docs](./docs/).

# Orbit determination
Nyx has an easy-to-use UI to perform orbit determination scenario.

Usage example:

```
$ cargo run --release -- "data/od_validation/*"                                      
    Finished release [optimized] target(s) in 0.03s
     Running `target/release/nyx 'data/od_validation/*'`
 INFO  nyx > Loaded scenario `data/od_validation/*`
✔ 

Select the sequence to execute · leo
 INFO  nyx_space::celestia::cosm > Loaded 14 ephemerides in 0 seconds.
 INFO  nyx_space::celestia::cosm > Loaded frame iau uranus
 INFO  nyx_space::celestia::cosm > Loaded frame iau earth
 INFO  nyx_space::celestia::cosm > Loaded frame iau neptune
 INFO  nyx_space::celestia::cosm > Loaded frame iau sun
 INFO  nyx_space::celestia::cosm > Loaded frame iau venus
 INFO  nyx_space::celestia::cosm > Loaded frame iau jupiter
 INFO  nyx_space::celestia::cosm > Loaded frame iau moon
 INFO  nyx_space::celestia::cosm > Loaded frame iau saturn
 INFO  nyx_space::celestia::cosm > Loaded frame iau mars
 INFO  nyx_space::io::gravity    > data/Luna_jggrx_1500e_sha.tab.gz loaded with (degree, order) = (20, 20)
 INFO  nyx_space::io::gravity    > data/JGM3.cof.gz loaded with (degree, order) = (20, 20)
 INFO  nyx_space::io::gravity    > data/Luna_jggrx_1500e_sha.tab.gz loaded with (degree, order) = (20, 20)
 INFO  nyx_space::io::gravity    > data/JGM3.cof.gz loaded with (degree, order) = (20, 20)
 INFO  nyx_space::io::odp        > Saving truth to ./data/truth.csv
 INFO  nyx_space::io::odp        > Generating measurements over 28800 seconds (~ 0.333 days)
 INFO  nyx_space::io::odp        > Initial state: [399 (  0)] 2000-01-01T12:00:00 TAI   sma = 7712.186236 km    ecc = 0.001000  inc = 63.434003 deg     raan = 135.000000 deg   aop = 90.000000 deg     ta = 0.000000 deg       120 kg
 INFO  nyx_space::io::odp        > Final state:   [399 (  0)] 2000-01-01T20:00:00 TAI   sma = 7725.759915 km    ecc = 0.001975  inc = 63.456750 deg     raan = 134.242870 deg   aop = 103.095351 deg    ta = 82.948194 deg      120 kg (computed in 1.251 seconds)
 INFO  nyx_space::io::odp        > Generated 673 measurements in 0.031 seconds
 INFO  nyx_space::od::ui         > Navigation propagating for a total of 25360 seconds (~ 0.294 days)
 INFO  nyx_space::od::ui         > Processing 673 measurements with covariance mapping
 INFO  nyx_space::od::ui         > [ODProcess]   0% done (1 measurements processed)
 INFO  nyx_space::od::ui         > [ODProcess]  10% done (68 measurements processed)
 INFO  nyx_space::od::ui         > [ODProcess]  20% done (135 measurements processed)
 INFO  nyx_space::od::ui         > [ODProcess]  30% done (202 measurements processed)
 INFO  nyx_space::od::ui         > [ODProcess]  40% done (270 measurements processed)
 INFO  nyx_space::od::ui         > [ODProcess]  50% done (337 measurements processed)
 INFO  nyx_space::od::ui         > [ODProcess]  60% done (404 measurements processed)
 INFO  nyx_space::od::ui         > [ODProcess]  70% done (472 measurements processed)
 INFO  nyx_space::od::ui         > [ODProcess]  80% done (539 measurements processed)
 INFO  nyx_space::od::ui         > [ODProcess]  90% done (606 measurements processed)
 INFO  nyx_space::od::ui         > [ODProcess] 100% done (673 measurements processed)
 INFO  nyx_space::io::odp        > Saving output to ./data/estimates.csv
$ cargo run --release -- "data/od_validation/*"c
$ PLOT="LEO No SNC - No noise - 0.005 km state error" pipenv run python od_plotter.py
File loaded.
Delta avr. -2.0e+02, -2.1e+02, -1.5e+01, -5.9e-01, -8.2e-01, -1.7e-02
$
```
## Examples
In all examples, we are running a CKF _without_ any iteration. In the "perfect" cases, the same model is used to generate the measurements and execute the navigation. In the "non-perfect" models, the harmonics used in navigation are different, and we have enabled state noise compensation at the level shown in the plots. The subtitle of the plot shows the average error between the truth trajectory and the estimated trajectory.

Scenario | Truth harmonics (Earth & Moon) | Nav harmonics (Earth & Moon) | Position error (km)
--|---|---|--
NRHO | 20x20 | 10x10  | 0.5
LRO | 20x20 | 10x10  | 0.05
LEO | 20x20 | 10x10  | 0.005

### Near Rectilinear Halo Orbit (NRHO)
#### Perfect
![State deviation & covariance](./data/od_plots/nrho-perfect-fig1.png "State deviation & covariance (perfect)")
![Error with truth](./data/od_plots/nrho-perfect-fig2.png "Error with truth (perfect)")

#### Realistic
![State deviation & covariance](./data/od_plots/nrho-fig1.png "State deviation & covariance (perfect)")
![Error with truth](./data/od_plots/nrho-fig2.png "Error with truth (perfect)")

### Lunar Reconnaissance Orbiter
**Note:** The state errors we see are most likely due to integration and estimation happening in the EME2000 frame where the Earth is considered the primary body, even if the Moon is significantly closer.

#### Perfect
![State deviation & covariance](./data/od_plots/lro-perfect-fig1.png "State deviation & covariance (perfect)")
![Error with truth](./data/od_plots/lro-perfect-fig2.png "Error with truth (perfect)")

#### Realistic
![State deviation & covariance](./data/od_plots/lro-fig1.png "State deviation & covariance (perfect)")
![Error with truth](./data/od_plots/lro-fig2.png "Error with truth (perfect)")

### LEO
**Note:** I think the spikes in covariance correspond to ground station switching.

#### Perfect
![State deviation & covariance](./data/od_plots/leo-perfect-fig1.png "State deviation & covariance (perfect)")
![Error with truth](./data/od_plots/leo-perfect-fig2.png "Error with truth (perfect)")