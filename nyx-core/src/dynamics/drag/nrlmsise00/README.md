# NRLMSISE-00 empirical atmosphere model

This implementation is cloned from the [\`tobari\` crate](https://github.com/sksat/orts/tree/main/tobari),
developed by sksat. All credits for the Rust port of NRLMSISE-00 go to sksat and contributors of the \`tobari\` crate.

The implementation is a clean-room rewrite based on the following references:
- Picone, J.M. et al. (2002), "NRLMSISE-00 empirical model of the atmosphere:
  Statistical comparisons and scientific issues", J. Geophys. Res., 107(A12), 1468,
  doi:10.1029/2002JA009430
- Hedin, A.E. (1991), "Extension of the MSIS thermosphere model into the middle
  and lower atmosphere", J. Geophys. Res., 96(A2), 1159-1172.
- Hedin, A.E. (1987), "MSIS-86 thermospheric model",
  J. Geophys. Res., 92(A5), 4649-4662.

Coefficient values are from the official NRL distribution, treated as published data.
NRLMSISE-00 is believed to be in the public domain as a U.S. Government work
(17 U.S.C. § 105), though no explicit license was provided by NRL.

Note: MSIS is a registered trademark. This module uses the name "NRLMSISE-00"
for nominative fair use (identifying compatibility with the NRL model).

## Additional Nyx specific work

+ Additional validation against `pymsis` by comparing the density and atmopheric temperatures in different
orbital regimes, storm configurations, and times throughout one day.
+ Validation of the propagation against Orekit, GMAT, and others
+ Configuration flags for the NRLMSISE switches
+ Support for both mean local solar time and apparent local solar time

## Validation notes

One might assume that validating an atmospheric drag model across industry-standard flight dynamics engines should yield
sub-meter parity over multi-day propagations, much like vacuum kinematics. This is a false premise. Empirical models are
artifacts of the era and numerical constraints in which they were calibrated, and discrepancies between modern software
engines are rarely driven by the physics equations themselves. They are driven by environmental data pre-processing and
geometric coordinate frame conventions.

The Nyx implementation of the NRLMSISE-00 atmospheric model has been strictly validated against the original Fortran 77
empirical baseline and GMAT. When propagating a Low Earth Orbit in a static space weather environment over 3 days, Nyx
matches GMAT to within 65 meters in-track error when using a different gravitational parameter and an RK89.
This confirms that the RK89 numerical integration, geodetic altitude mapping, and core MSIS boundary evaluations are correct.

However, when dynamic space weather is introduced or when comparing against other open-source tools, single digit kilometer
deviations may arise. These are not integration defects, but conscious architectural decisions grounded in the model's
original scientific calibration.

### Local Solar Time: Mean versus Apparent

The original NRLMSISE-00 equations were calibrated to evaluate the diurnal density bulge using a linear proxy:
Mean Local Solar Time (calculated strictly as UTC seconds past midnight divided by 3600, plus geodetic longitude divided by 15).

Modern tools occasionally attempt to "upgrade" this by passing the true geometric Ephemeris Sun vector (Apparent LST).
This is an architectural error. Applying the true 13-minute Equation of Time offset shifts the diurnal bulge by over 3 degrees of longitude,
breaking the empirical fit of the model's spherical harmonics.

Furthermore, this has led to significant bugs in legacy software. Orekit 13.1, for example, generates a nearly 4-kilometer
along-track error against Nyx over 7 days. As acknowledged by its maintainers (Issues #1993 and #2003), Orekit attempted to
implement Apparent LST but introduced a severe coordinate frame defect—subtracting an inertial right ascension from a body-fixed
longitude—shifting the atmospheric bulge to a physically arbitrary location. Nyx enforces the required Mean LST by default to
maintain mathematical alignment with the original Naval Research Laboratory calibration.

### Space Weather Interpolation

When evaluating dynamic CSSI space weather files, Nyx will diverge from GMAT by 1 to 6 kilometers over a 3-day high-drag propagation:
this variance is an artifact of data interpolation.

The planetary $A_p$ and $K_p$ indices are discrete measurements representing maximum geomagnetic fluctuations over strict 3-hour bins.
The original scientists calibrated the MSIS model by feeding it these piecewise-constant step-functions.
GMAT, conversely, applies a global cubic spline to the space weather arrays to prevent sudden density discontinuities from rejecting
steps in its variable-step integrators.

Applying continuous splines to discrete step-functions acts as a low-pass filter, mathematically smoothing out the peak density responses
of severe geomagnetic storms before the MSIS equations ever evaluate them. My own internal testing proved that applying a linear smoothing
hack to Nyx's inputs immediately collapses the divergence with GMAT by up to 90%.
Nyx rejects this double-smoothing and natively feeds the raw 3-hourly step-functions exactly as the baseline requires.

### Operational Reality

Chasing sub-kilometer deterministic parity in a thermospheric environment is fundamentally an academic exercise. The NRLMSISE-00 model
inherently carries a 15% to 20% uncertainty due to the chaotic nature of solar weather.

In actual flight operations, an orbit determination filter cannot distinguish between a 0.5% spatial variation caused by cubic spline
interpolation, a missed solar flux prediction, or a fractional shift in a spacecraft's attitude. The estimator simply absorbs the
cumulative uncertainty entirely into the estimated drag coefficient ($C_d$). By ensuring the kinematics are pristine and the environmental
inputs trace strictly to the scientific baseline, Nyx provides a fully verified dynamics engine where the physical truth of the atmosphere
is handled exactly where it belongs: in the orbit determination filter.

