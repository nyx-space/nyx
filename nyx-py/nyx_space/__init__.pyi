# ruff: noqa
import typing

import numpy
from anise import astro, rotation, time

@typing.final
class DragData:
    area_m2: float
    coeff_drag: float

    def __init__(self, *args: typing.Any | None, **kwargs: typing.Any | None) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls, area_m2: typing.Any, coeff_drag: typing.Any=None) -> DragData:...

    @staticmethod
    def from_asn1(data: bytes) -> astro.DragData:
        """Decodes an ASN.1 DER encoded byte array into a DragData object."""

    def to_asn1(self) -> bytes:
        """Encodes this DragData object into an ASN.1 DER encoded byte array."""



@typing.final
class ExportCfg:
    """Configuration for exporting from Nyx to local disk."""

    def __init__(self, *args: typing.Any | None, **kwargs: typing.Any | None) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
Configuration for exporting from Nyx to local disk."""

    def __new__(cls, timestamped: typing.Any=False) -> ExportCfg:
        """Configuration for exporting from Nyx to local disk."""

    def __eq__(self, value: object) -> bool:
        """Return self==value."""

    def __ge__(self, value: typing.Any) -> bool:
        """Return self>=value."""

    def __gt__(self, value: typing.Any) -> bool:
        """Return self>value."""

    def __le__(self, value: typing.Any) -> bool:
        """Return self<=value."""

    def __lt__(self, value: typing.Any) -> bool:
        """Return self<value."""

    def __ne__(self, value: object) -> bool:
        """Return self!=value."""



@typing.final
class Frame:
    """A Frame uniquely defined by its ephemeris center and orientation. Refer to FrameDetail for frames combined with parameters."""
    ephemeris_id: int
    force_inertial: typing.Any
    frozen_epoch: time.Epoch
    orientation_id: int
    shape: astro.Ellipsoid

    def __init__(self, ephemeris_id: int, orientation_id: int, force_inertial: bool, mu_km3_s2: float | None, shape: astro.Ellipsoid | None, frozen_epoch: time.Epoch | None, *args: typing.Any | None, **kwargs: typing.Any | None) -> astro.Frame:
        """Initialize self.  See help(type(self)) for accurate signature.
A Frame uniquely defined by its ephemeris center and orientation. Refer to FrameDetail for frames combined with parameters."""

    def __new__(cls, ephemeris_id: int, orientation_id: int, mu_km3_s2: float | None=None, shape: astro.Ellipsoid | None=None, force_inertial: bool=False, frozen_epoch: time.Epoch | None=None) -> Frame:
        """A Frame uniquely defined by its ephemeris center and orientation. Refer to FrameDetail for frames combined with parameters."""

    def ephem_origin_id_match(self, other_id: int) -> bool:
        """Returns true if the ephemeris origin is equal to the provided ID"""

    def ephem_origin_match(self, other: astro.Frame) -> bool:
        """Returns true if the ephemeris origin is equal to the provided frame"""

    def flattening(self) -> float:
        """Returns the flattening ratio (unitless)"""

    @staticmethod
    def from_asn1(data: bytes) -> astro.Frame:
        """Decodes an ASN.1 DER encoded byte array into a Frame."""

    @staticmethod
    def from_frameuid(frameuid: typing.Any) -> astro.Frame:
        """Creates a Frame from a FrameUid"""

    def is_celestial(self) -> bool:
        """Returns whether this is a celestial frame"""

    def is_dynamic(self) -> bool:
        """Returns true if this is a dynamic frame, e.g. Mean/True of Date/Epoch"""

    def is_geodetic(self) -> bool:
        """Returns whether this is a geodetic frame"""

    def mean_equatorial_radius_km(self) -> float:
        """Returns the mean equatorial radius in km, if defined"""

    def mu_km3_s2(self) -> float:
        """Returns the gravitational parameters of this frame, if defined"""

    def orient_origin_id_match(self, other_id: int) -> bool:
        """Returns true if the orientation origin is equal to the provided ID"""

    def orient_origin_match(self, other: astro.Frame) -> bool:
        """Returns true if the orientation origin is equal to the provided frame"""

    def polar_radius_km(self) -> float:
        """Returns the polar radius in km, if defined"""

    def semi_major_radius_km(self) -> float:
        """Returns the semi major radius of the tri-axial ellipoid shape of this frame, if defined"""

    def strip(self) -> None:
        """Removes the graviational parameter and the shape information from this frame.
Use this to prevent astrodynamical computations."""

    def to_asn1(self) -> bytes:
        """Encodes this Frame into an ASN.1 DER encoded byte array."""

    def to_frameuid(self) -> astro.FrameUid:
        """Converts this Frame to a FrameUid"""

    def with_ephem(self, new_ephem_id: int) -> astro.Frame:
        """Returns a copy of this Frame whose ephemeris ID is set to the provided ID"""

    def with_mu_km3_s2(self, mu_km3_s2: float) -> astro.Frame:
        """Returns a copy of this frame with the graviational parameter set to the new value."""

    def with_orient(self, new_orient_id: int) -> astro.Frame:
        """Returns a copy of this Frame whose orientation ID is set to the provided ID"""

    def __eq__(self, value: object) -> bool:
        """Return self==value."""

    def __ge__(self, value: typing.Any) -> bool:
        """Return self>=value."""

    def __getnewargs__(self) -> tuple:
        """Allows for pickling the object"""

    def __gt__(self, value: typing.Any) -> bool:
        """Return self>value."""

    def __le__(self, value: typing.Any) -> bool:
        """Return self<=value."""

    def __lt__(self, value: typing.Any) -> bool:
        """Return self<value."""

    def __ne__(self, value: object) -> bool:
        """Return self!=value."""



@typing.final
class GuidanceMode:

    def __init__(self, *args: typing.Any | None, **kwargs: typing.Any | None) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> GuidanceMode:...

    def __int__(self) -> None:
        """int(self)"""

    Coast: GuidanceMode = ...
    Inhibit: GuidanceMode = ...
    Thrust: GuidanceMode = ...

@typing.final
class Mass:
    """Defines a spacecraft mass a the sum of the dry (structural) mass and the propellant mass, both in kilogram"""
    dry_mass_kg: float
    extra_mass_kg: float
    prop_mass_kg: float

    def __init__(self, *args: typing.Any | None, **kwargs: typing.Any | None) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
Defines a spacecraft mass a the sum of the dry (structural) mass and the propellant mass, both in kilogram"""

    def __new__(cls, dry_mass_kg: typing.Any, prop_mass_kg: typing.Any=None, extra_mass_kg: typing.Any=None) -> Mass:
        """Defines a spacecraft mass a the sum of the dry (structural) mass and the propellant mass, both in kilogram"""

    def abs(self) -> astro.Mass:
        """Returns a Mass structure that is guaranteed to be physically correct"""

    @staticmethod
    def from_asn1(data: bytes) -> astro.Mass:
        """Decodes an ASN.1 DER encoded byte array into a Mass object."""

    def is_valid(self) -> bool:
        """Returns true if all the masses are greater or equal to zero"""

    def to_asn1(self) -> bytes:
        """Encodes this Mass object into an ASN.1 DER encoded byte array."""

    def total_mass_kg(self) -> float:
        """Returns the total mass in kg"""



@typing.final
class Orbit:
    """Defines a Cartesian state in a given frame at a given epoch in a given time scale. Radius data is expressed in kilometers. Velocity data is expressed in kilometers per second.
Regardless of the constructor used, this struct stores all the state information in Cartesian coordinates as these are always non singular.

Unless noted otherwise, algorithms are from GMAT 2016a [StateConversionUtil.cpp](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/StateConversionUtil.cpp)."""
    epoch: time.Epoch
    frame: astro.Frame
    vx_km_s: float
    vy_km_s: float
    vz_km_s: float
    x_km: float
    y_km: float
    z_km: float

    def __init__(self, *args: typing.Any | None, **kwargs: typing.Any | None) -> astro.Orbit:
        """Initialize self.  See help(type(self)) for accurate signature.
Defines a Cartesian state in a given frame at a given epoch in a given time scale. Radius data is expressed in kilometers. Velocity data is expressed in kilometers per second.
Regardless of the constructor used, this struct stores all the state information in Cartesian coordinates as these are always non singular.

Unless noted otherwise, algorithms are from GMAT 2016a [StateConversionUtil.cpp](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/StateConversionUtil.cpp)."""

    def __new__(cls, *args: typing.Any | None) -> Orbit:
        """Defines a Cartesian state in a given frame at a given epoch in a given time scale. Radius data is expressed in kilometers. Velocity data is expressed in kilometers per second.
Regardless of the constructor used, this struct stores all the state information in Cartesian coordinates as these are always non singular.

Unless noted otherwise, algorithms are from GMAT 2016a [StateConversionUtil.cpp](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/StateConversionUtil.cpp)."""

    def abs_difference(self, other: astro.Orbit) -> tuple:
        """Returns the absolute position and velocity differences in km and km/s between this orbit and another.
Raises an error if the frames do not match (epochs do not need to match)."""

    def abs_pos_diff_km(self, other: astro.Orbit) -> float:
        """Returns the absolute position difference in kilometer between this orbit and another.
Raises an error if the frames do not match (epochs do not need to match)."""

    def abs_vel_diff_km_s(self, other: astro.Orbit) -> float:
        """Returns the absolute velocity difference in kilometer per second between this orbit and another.
Raises an error if the frames do not match (epochs do not need to match)."""

    def add_aop_deg(self, delta_aop_deg: float) -> astro.Orbit:
        """Returns a copy of the state with a provided AOP added to the current one"""

    def add_apoapsis_periapsis_km(self, delta_ra_km: float, delta_rp_km: float) -> astro.Orbit:
        """Returns a copy of this state with the provided apoasis and periapsis added to the current values"""

    def add_ecc(self, delta_ecc: float) -> astro.Orbit:
        """Returns a copy of the state with a provided ECC added to the current one"""

    def add_inc_deg(self, delta_inc_deg: float) -> astro.Orbit:
        """Returns a copy of the state with a provided INC added to the current one"""

    def add_raan_deg(self, delta_raan_deg: float) -> astro.Orbit:
        """Returns a copy of the state with a provided RAAN added to the current one"""

    def add_sma_km(self, delta_sma_km: float) -> astro.Orbit:
        """Returns a copy of the state with a provided SMA added to the current one"""

    def add_ta_deg(self, delta_ta_deg: float) -> astro.Orbit:
        """Returns a copy of the state with a provided TA added to the current one"""

    def altitude_km(self) -> float:
        """Returns the altitude in km"""

    def aol_deg(self) -> float:
        """Returns the argument of latitude in degrees

NOTE: If the orbit is near circular, the AoL will be computed from the true longitude
instead of relying on the ill-defined true anomaly."""

    def aop_brouwer_short_deg(self) -> float:
        """Returns the Brouwer-short mean Argument of Perigee in degrees."""

    def aop_deg(self) -> float:
        """Returns the argument of periapsis in degrees"""

    def apoapsis_altitude_km(self) -> float:
        """Returns the altitude of apoapsis (or apogee around Earth), in kilometers."""

    def apoapsis_km(self) -> float:
        """Returns the radius of apoapsis (or apogee around Earth), in kilometers."""

    def at_epoch(self, new_epoch: time.Epoch) -> astro.Orbit:
        """Adjusts the equinoctial mean longitude this orbit via the mean motion.

# Astrodynamics note
This is not a true propagation of the orbit. This is akin to a two body propagation ONLY without any other force models applied.
Use Nyx for high fidelity propagation. This implementation uses equinoctial elements and should be well behaved for circular equatorial orbits."""

    def c3_km2_s2(self) -> float:
        """Returns the $C_3$ of this orbit in km^2/s^2"""

    def cartesian_pos_vel(self) -> numpy.ndarray:
        """Returns this state as a Cartesian vector of size 6 in [km, km, km, km/s, km/s, km/s]

Note that the time is **not** returned in the vector."""

    def dcm3x3_from_rcn_to_inertial(self) -> rotation.DCM:
        """Builds the rotation matrix that rotates from this state's inertial frame to this state's RCN frame (radial, cross, normal)

# Frame warning
If the stattion is NOT in an inertial frame, then this computation is INVALID.

# Algorithm
1. Compute \\hat{r}, \\hat{h}, the unit vectors of the radius and orbital momentum.
2. Compute the cross product of these
3. Build the DCM with these unit vectors
4. Return the DCM structure"""

    def dcm3x3_from_ric_to_inertial(self) -> rotation.DCM:
        """Builds the rotation matrix that rotates from this state's inertial frame to this state's RIC frame

# Frame warning
If the state is NOT in an inertial frame, then this computation is INVALID.

# Algorithm
1. Build the c vector as the normalized orbital momentum vector
2. Build the i vector as the cross product of \\hat{r} and c
3. Build the RIC DCM as a 3x3 of the columns [\\hat{r}, \\hat{i}, \\hat{c}]
4. Return the DCM structure **without** accounting for the transport theorem."""

    def dcm3x3_from_topocentric_to_body_fixed(self) -> rotation.DCM:
        """Builds the rotation matrix that rotates from the topocentric frame (SEZ) into the body fixed frame of this state.

# Frame warning
If the state is NOT in a body fixed frame (i.e. ITRF93), then this computation is INVALID.

# Source
From the GMAT MathSpec, page 30 section 2.6.9 and from `Calculate_RFT` in `TopocentricAxes.cpp`, this returns the
rotation matrix from the topocentric frame (SEZ) to body fixed frame.
In the GMAT MathSpec notation, R_{IF} is the DCM from body fixed to inertial. Similarly, R{FT} is from topocentric
to body fixed."""

    def dcm3x3_from_vnc_to_inertial(self) -> rotation.DCM:
        """Builds the rotation matrix that rotates from this state's inertial frame to this state's VNC frame (velocity, normal, cross)

# Frame warning
If the stattion is NOT in an inertial frame, then this computation is INVALID.

# Algorithm
1. Compute \\hat{v}, \\hat{h}, the unit vectors of the radius and orbital momentum.
2. Compute the cross product of these
3. Build the DCM with these unit vectors
4. Return the DCM structure."""

    def dcm_from_rcn_to_inertial(self) -> rotation.DCM:
        """Builds the rotation matrix that rotates from this state's inertial frame to this state's RCN frame (radial, cross, normal)

# Frame warning
If the stattion is NOT in an inertial frame, then this computation is INVALID.

# Algorithm
1. Compute \\hat{r}, \\hat{h}, the unit vectors of the radius and orbital momentum.
2. Compute the cross product of these
3. Build the DCM with these unit vectors
4. Return the DCM structure with a 6x6 DCM with the time derivative of the VNC frame set.

# Note on the time derivative
If the pre or post states cannot be computed, then the time derivative of the DCM will _not_ be set.
Further note that most astrodynamics tools do *not* account for the time derivative in the RIC frame."""

    def dcm_from_ric_to_inertial(self) -> rotation.DCM:
        """Builds the rotation matrix that rotates from this state's inertial frame to this state's RIC frame

# Frame warning
If the state is NOT in an inertial frame, then this computation is INVALID.

# Algorithm
1. Compute the state data one millisecond before and one millisecond assuming two body dynamics
2. Compute the DCM for this state, and the pre and post states
3. Build the c vector as the normalized orbital momentum vector
4. Build the i vector as the cross product of \\hat{r} and c
5. Build the RIC DCM as a 3x3 of the columns [\\hat{r}, \\hat{i}, \\hat{c}], for the post, post, and current states
6. Compute the difference between the DCMs of the pre and post states, to build the DCM time derivative
7. Return the DCM structure with a 6x6 state DCM.

# Note on the time derivative
If the pre or post states cannot be computed, then the time derivative of the DCM will _not_ be set.
Further note that most astrodynamics tools do *not* account for the time derivative in the RIC frame."""

    def dcm_from_topocentric_to_body_fixed(self) -> rotation.DCM:
        """Builds the rotation matrix that rotates from the topocentric frame (SEZ) into the body fixed frame of this state.

# Frame warnings
+ If the state is NOT in a body fixed frame (i.e. ITRF93), then this computation is INVALID.
+ (Usually) no time derivative can be computed: the orbit is expected to be a body fixed frame where the `at_epoch` function will fail. Exceptions for Moon body fixed frames.

# UNUSED Arguments
+ `from`: ID of this new frame. Only used to set the "from" frame of the DCM. -- No longer used since 0.5.3

# Source
From the GMAT MathSpec, page 30 section 2.6.9 and from `Calculate_RFT` in `TopocentricAxes.cpp`, this returns the
rotation matrix from the topocentric frame (SEZ) to body fixed frame.
In the GMAT MathSpec notation, R_{IF} is the DCM from body fixed to inertial. Similarly, R{FT} is from topocentric
to body fixed."""

    def dcm_from_vnc_to_inertial(self) -> rotation.DCM:
        """Builds the rotation matrix that rotates from this state's inertial frame to this state's VNC frame (velocity, normal, cross)

# Frame warning
If the stattion is NOT in an inertial frame, then this computation is INVALID.

# Algorithm
1. Compute \\hat{v}, \\hat{h}, the unit vectors of the radius and orbital momentum.
2. Compute the cross product of these
3. Build the DCM with these unit vectors
4. Compute the difference between the DCMs of the pre and post states (+/- 1 ms), to build the DCM time derivative
4. Return the DCM structure with a 6x6 DCM with the time derivative of the VNC frame set.

# Note on the time derivative
If the pre or post states cannot be computed, then the time derivative of the DCM will _not_ be set.
Further note that most astrodynamics tools do *not* account for the time derivative in the RIC frame."""

    def dcm_to_inertial(self, local_frame: astro.LocalFrame) -> rotation.DCM:
        """Returns the DCM to rotate this orbit from the provided local frame to the inertial frame."""

    def declination_deg(self) -> float:
        """Returns the declination of this orbit in degrees"""

    def distance_to_km(self, other: astro.Orbit) -> float:
        """Returns the distance in kilometers between this state and another state, if both frame match (epoch does not need to match)."""

    def duration_to_radius(self, radius_km: float) -> time.Duration:
        """Calculates the duration to reach a specific radius in the orbit.

This function computes the time it will take for the orbiting body to reach
the given `radius_km` from its current position. The calculation assumes
two-body dynamics and considers the direction of motion.

# Assumptions & Limitations

- Assumes pure Keplerian motion.
- For elliptical orbits, if the radius is reachable at two points (ascending and descending parts
of the orbit), this function calculates the time to reach the radius corresponding to the
true anomaly in `[0, PI]` (typically the ascending part or up to apoapsis if starting past periapsis).
- For circular orbits, if the radius is within the apoapse and periapse, then a duration of zero is returned.
- For hyperbolic/parabolic orbits, the true anomaly at radius is also computed in `[0, PI]`. If this
point is in the past, the function returns an error, as it doesn't look for solutions on the
departing leg if `nu > PI` would be required (unless current TA is already > PI and target radius is further along).
The current implementation strictly uses the `acos` result, so `nu_rad_at_radius` is always `0 <= nu <= PI`.
This means it finds the time to reach the radius on the path from periapsis up to the point where true anomaly is PI."""

    def ea_deg(self) -> float:
        """Returns the eccentric anomaly in degrees

This is a conversion from GMAT's StateConversionUtil::TrueToEccentricAnomaly"""

    def ecc(self) -> float:
        """Returns the eccentricity (no unit)"""

    def ecc_brouwer_short(self) -> float:
        """Returns the Brouwer-short mean eccentricity."""

    def energy_km2_s2(self) -> float:
        """Returns the specific mechanical energy in km^2/s^2"""

    def eq_within(self, other: astro.Orbit, radial_tol_km: float, velocity_tol_km_s: float) -> bool:
        """Returns whether this orbit and another are equal within the specified radial and velocity absolute tolerances"""

    def equinoctial_a_km(self) -> float:
        """Returns the equinoctial semi-major axis (a) in km."""

    def equinoctial_elements(self) -> list[float]:
        """Returns the six equinoctial elements in order: sma (km), h, k, p, q, lambda (deg)"""

    def equinoctial_h(self) -> float:
        """Returns the equinoctial element h (ecc * sin(aop + raan))."""

    def equinoctial_k(self) -> float:
        """Returns the equinoctial element k (ecc * cos(aop + raan))."""

    def equinoctial_lambda_mean_deg(self) -> float:
        """Returns the equinoctial mean longitude (lambda = raan + aop + ma) in degrees."""

    def equinoctial_p(self) -> float:
        """Returns the equinoctial element p (sin(inc/2) * sin(raan))."""

    def equinoctial_q(self) -> float:
        """Returns the equinoctial element q (sin(inc/2) * cos(raan))."""

    def fpa_deg(self) -> float:
        """Returns the flight path angle in degrees"""

    @staticmethod
    def from_asn1(data: bytes) -> astro.Orbit:
        """Decodes an ASN.1 DER encoded byte array into a CartesianState (Orbit)."""

    @staticmethod
    def from_cartesian(x_km: float, y_km: float, z_km: float, vx_km_s: float, vy_km_s: float, vz_km_s: float, epoch: time.Epoch, frame: astro.Frame) -> astro.Orbit:
        """Creates a new Cartesian state in the provided frame at the provided Epoch.

**Units:** km, km, km, km/s, km/s, km/s"""

    @staticmethod
    def from_cartesian_npy(pos_vel: numpy.ndarray, epoch: time.Epoch, frame: astro.Frame) -> astro.Orbit:
        """Creates a new Cartesian state from a numpy array, an epoch, and a frame.

**Units:** km, km, km, km/s, km/s, km/s"""

    @staticmethod
    def from_equinoctial(sma_km: float, h: float, k: float, p: float, q: float, lambda_deg: float, epoch: time.Epoch, frame: astro.Frame) -> astro.Orbit:
        """Attempts to create a new Orbit from the Equinoctial orbital elements.

# Implementation notes
Note that this function computes the Keplerian elements from the equinoctial and then
calls the try_keplerian_mean_anomaly initializer."""

    @staticmethod
    def from_keplerian(sma_km: float, ecc: float, inc_deg: float, raan_deg: float, aop_deg: float, ta_deg: float, epoch: time.Epoch, frame: astro.Frame) -> astro.Orbit:
        """Creates a new Orbit around the provided Celestial or Geoid frame from the Keplerian orbital elements.

**Units:** km, none, degrees, degrees, degrees, degrees

NOTE: The state is defined in Cartesian coordinates as they are non-singular. This causes rounding
errors when creating a state from its Keplerian orbital elements (cf. the state tests).
One should expect these errors to be on the order of 1e-12."""

    @staticmethod
    def from_keplerian_altitude(sma_altitude_km: float, ecc: float, inc_deg: float, raan_deg: float, aop_deg: float, ta_deg: float, epoch: time.Epoch, frame: astro.Frame) -> astro.Orbit:
        """Creates a new Orbit from the provided semi-major axis altitude in kilometers"""

    @staticmethod
    def from_keplerian_apsis_altitude(apo_alt_km: float, peri_alt_km: float, inc_deg: float, raan_deg: float, aop_deg: float, ta_deg: float, epoch: time.Epoch, frame: astro.Frame) -> astro.Orbit:
        """Creates a new Orbit from the provided altitudes of apoapsis and periapsis, in kilometers"""

    @staticmethod
    def from_keplerian_apsis_radii(r_a_km: float, r_p_km: float, inc_deg: float, raan_deg: float, aop_deg: float, ta_deg: float, epoch: time.Epoch, frame: astro.Frame) -> astro.Orbit:
        """Attempts to create a new Orbit from the provided radii of apoapsis and periapsis, in kilometers"""

    @staticmethod
    def from_keplerian_mean_anomaly(sma_km: float, ecc: float, inc_deg: float, raan_deg: float, aop_deg: float, ma_deg: float, epoch: time.Epoch, frame: astro.Frame) -> astro.Orbit:
        """Initializes a new orbit from the Keplerian orbital elements using the mean anomaly instead of the true anomaly.

# Implementation notes
This function starts by converting the mean anomaly to true anomaly, and then it initializes the orbit
using the keplerian(..) method.
The conversion is from GMAT's MeanToTrueAnomaly function, transliterated originally by Claude and GPT4 with human adjustments."""

    @staticmethod
    def from_latlongalt(latitude_deg: float, longitude_deg: float, height_km: float, epoch: time.Epoch, frame: astro.Frame) -> astro.Orbit:
        """Creates a new Orbit from the latitude (φ), longitude (λ) and height (in km) with respect to the frame's ellipsoid, and with ZERO angular velocity in this frame.
Use this initializer for creating a fixed point in the ITRF93 frame for example: the correct angular velocity will be applied when transforming this to EME2000 for example.

Refer to [try_latlongalt_omega] if you need to build a fixed point with a non-zero angular velocity in the definition frame.

NOTE: This computation differs from the spherical coordinates because we consider the flattening of body.
Reference: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016"""

    @staticmethod
    def from_latlongalt_omega(latitude_deg: float, longitude_deg: float, height_km: float, angular_velocity_rad_s: numpy.ndarray, epoch: time.Epoch, frame: astro.Frame) -> astro.Orbit:
        """Creates a new Orbit from the latitude (φ), longitude (λ) and height (in km) with respect to the frame's ellipsoid given the angular velocity vector.
NOTE: Only specify the angular velocity if there's an EXTRA angular velocity of the lat/long/alt point in the provided frame.

Consider using the [Almanac]'s [angular_velocity_wrt_j2000_rad_s] function or [angular_velocity_rad_s] to retrieve the exact angular velocity vector between two orientations.
Example: build a lat/long/alt point referenced in the ITRF93 frame but by specifying the Frame as the EME2000 frame and providing the angular velocity between the ITRF93 and EME2000 frame at the desired time.

NOTE: This computation differs from the spherical coordinates because we consider the flattening of body.
Reference: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016"""

    def height_km(self) -> float:
        """Returns the geodetic height in km.

Reference: Vallado, 4th Ed., Algorithm 12 page 172."""

    def hmag(self) -> float:
        """Returns the norm of the orbital momentum"""

    def hx(self) -> float:
        """Returns the orbital momentum value on the X axis"""

    def hy(self) -> float:
        """Returns the orbital momentum value on the Y axis"""

    def hyperbolic_anomaly_deg(self) -> float:
        """Returns the hyperbolic anomaly in degrees between 0 and 360.0
Returns an error if the orbit is not hyperbolic."""

    def hz(self) -> float:
        """Returns the orbital momentum value on the Z axis"""

    def inc_brouwer_short_deg(self) -> float:
        """Returns the Brouwer-short mean inclination in degrees."""

    def inc_deg(self) -> float:
        """Returns the inclination in degrees"""

    def is_brouwer_short_valid(self) -> bool:
        """Returns whether this state satisfies the requirement to compute the Mean Brouwer Short orbital
element set.

This is a conversion from GMAT's StateConversionUtil::CartesianToBrouwerMeanShort.
The details are at the log level `info`.
NOTE: Mean Brouwer Short are only defined around Earth. However, `nyx` does *not* check the
main celestial body around which the state is defined (GMAT does perform this verification)."""

    def latitude_deg(self) -> float:
        """Returns the geodetic latitude (φ) in degrees. Value is between -180 and +180 degrees.

# Frame warning
This state MUST be in the body fixed frame (e.g. ITRF93) prior to calling this function, or the computation is **invalid**."""

    def latlongalt(self) -> tuple:
        """Returns the geodetic latitude, geodetic longitude, and geodetic height, respectively in degrees, degrees, and kilometers.

# Algorithm
This uses the Heikkinen procedure, which is not iterative. The results match Vallado and GMAT."""

    def light_time(self) -> time.Duration:
        """Returns the light time duration between this object and the origin of its reference frame."""

    def longitude_360_deg(self) -> float:
        """Returns the geodetic longitude (λ) in degrees. Value is between 0 and 360 degrees.

# Frame warning
This state MUST be in the body fixed frame (e.g. ITRF93) prior to calling this function, or the computation is **invalid**."""

    def longitude_deg(self) -> float:
        """Returns the geodetic longitude (λ) in degrees. Value is between -180 and 180 degrees.

# Frame warning
This state MUST be in the body fixed frame (e.g. ITRF93) prior to calling this function, or the computation is **invalid**."""

    def ltan_deg(self) -> float:
        """Returns the Longitude of the Ascending Node (LTAN), or an error of equatorial orbits"""

    def ma_brouwer_short_deg(self) -> float:
        """Returns the Brouwer-short mean Mean Anomaly in degrees."""

    def ma_deg(self) -> float:
        """Returns the mean anomaly in degrees

This is a conversion from GMAT's StateConversionUtil::TrueToMeanAnomaly"""

    def mean_motion_deg_s(self) -> float:
        """Returns the mean motion in degrees per seconds"""

    def periapsis_altitude_km(self) -> float:
        """Returns the altitude of periapsis (or perigee around Earth), in kilometers."""

    def periapsis_km(self) -> float:
        """Returns the radius of periapsis (or perigee around Earth), in kilometers."""

    def period(self) -> time.Duration:
        """Returns the period"""

    def raan_brouwer_short_deg(self) -> float:
        """Returns the Brouwer-short mean Right Ascension of the Ascending Node in degrees."""

    def raan_deg(self) -> float:
        """Returns the right ascension of the ascending node in degrees"""

    def radius_km(self) -> numpy.ndarray:
        """radius vector in km"""

    def rel_difference(self, other: astro.Orbit) -> tuple:
        """Returns the relative difference between this orbit and another for the position and velocity, respectively the first and second return values.
Both return values are UNITLESS because the relative difference is computed as the absolute difference divided by the rmag and vmag of this object.
Raises an error if the frames do not match, if the position is zero or the velocity is zero."""

    def rel_pos_diff(self, other: astro.Orbit) -> float:
        """Returns the relative position difference (unitless) between this orbit and another.
This is computed by dividing the absolute difference by the norm of this object's radius vector.
If the radius is zero, this function raises a math error.
Raises an error if the frames do not match or  (epochs do not need to match)."""

    def rel_vel_diff(self, other: astro.Orbit) -> float:
        """Returns the absolute velocity difference in kilometer per second between this orbit and another.
Raises an error if the frames do not match (epochs do not need to match)."""

    def ric_difference(self, other: astro.Orbit) -> astro.Orbit:
        """Returns a Cartesian state representing the RIC difference between self and other, in position and velocity (with transport theorem).
Refer to dcm_from_ric_to_inertial for details on the RIC frame.

# Algorithm
1. Compute the difference between `other` and `self`
2. Compute the RIC DCM of `self`
3. Rotate the difference into the RIC frame of `self`
4. Strip the astrodynamical information from the frame, enabling only computations from `CartesianState`"""

    def right_ascension_deg(self) -> float:
        """Returns the right ascension of this orbit in degrees"""

    def rmag_km(self) -> float:
        """Returns the magnitude of the radius vector in km"""

    def rms_radius_km(self, other: astro.Orbit) -> float:
        """Returns the root sum squared (RMS) radius difference between this state and another state, if both frames match (epoch does not need to match)"""

    def rms_velocity_km_s(self, other: astro.Orbit) -> float:
        """Returns the root sum squared (RMS) velocity difference between this state and another state, if both frames match (epoch does not need to match)"""

    def rss_radius_km(self, other: astro.Orbit) -> float:
        """Returns the root mean squared (RSS) radius difference between this state and another state, if both frames match (epoch does not need to match)"""

    def rss_velocity_km_s(self, other: astro.Orbit) -> float:
        """Returns the root mean squared (RSS) velocity difference between this state and another state, if both frames match (epoch does not need to match)"""

    def semi_minor_axis_km(self) -> float:
        """Returns the semi minor axis in km, includes code for a hyperbolic orbit"""

    def semi_parameter_km(self) -> float:
        """Returns the semi parameter (or semilatus rectum)"""

    def set_aop_deg(self, new_aop_deg: float) -> None:
        """Mutates this orbit to change the AOP"""

    def set_ecc(self, new_ecc: float) -> None:
        """Mutates this orbit to change the ECC"""

    def set_inc_deg(self, new_inc_deg: float) -> None:
        """Mutates this orbit to change the INC"""

    def set_raan_deg(self, new_raan_deg: float) -> None:
        """Mutates this orbit to change the RAAN"""

    def set_sma_km(self, new_sma_km: float) -> None:
        """Mutates this orbit to change the SMA"""

    def set_ta_deg(self, new_ta_deg: float) -> None:
        """Mutates this orbit to change the TA"""

    def sma_altitude_km(self) -> float:
        """Returns the SMA altitude in km"""

    def sma_brouwer_short_km(self) -> float:
        """Returns the Brouwer-short mean semi-major axis in km."""

    def sma_km(self) -> float:
        """Returns the semi-major axis in km"""

    def ta_deg(self) -> float:
        """Returns the true anomaly in degrees between 0 and 360.0

NOTE: This function will emit a warning stating that the TA should be avoided if in a very near circular orbit
Code from <https://github.com/ChristopherRabotin/GMAT/blob/80bde040e12946a61dae90d9fc3538f16df34190/src/gmatutil/util/StateConversionUtil.cpp#L6835>

LIMITATION: For an orbit whose true anomaly is (very nearly) 0.0 or 180.0, this function may return either 0.0 or 180.0 with a very small time increment.
This is due to the precision of the cosine calculation: if the arccosine calculation is out of bounds, the sign of the cosine of the true anomaly is used
to determine whether the true anomaly should be 0.0 or 180.0. **In other words**, there is an ambiguity in the computation in the true anomaly exactly at 180.0 and 0.0."""

    def ta_dot_deg_s(self) -> float:
        """Returns the time derivative of the true anomaly computed as the 360.0 degrees divided by the orbital period (in seconds)."""

    def tlong_deg(self) -> float:
        """Returns the true longitude in degrees"""

    def to_asn1(self) -> bytes:
        """Encodes this CartesianState (Orbit) into an ASN.1 DER encoded byte array."""

    def velocity_declination_deg(self) -> float:
        """Returns the velocity declination of this orbit in degrees"""

    def velocity_km_s(self) -> numpy.ndarray:
        """velocity vector in km/s"""

    def vinf_periapsis_km(self, turn_angle_degrees: float) -> float:
        """Returns the radius of periapse in kilometers for the provided turn angle of this hyperbolic orbit.
Returns an error if the orbit is not hyperbolic."""

    def vinf_turn_angle_deg(self, periapsis_km: float) -> float:
        """Returns the turn angle in degrees for the provided radius of periapse passage of this hyperbolic orbit
Returns an error if the orbit is not hyperbolic."""

    def vmag_km_s(self) -> float:
        """Returns the magnitude of the velocity vector in km/s"""

    def vnc_difference(self, other: astro.Orbit) -> astro.Orbit:
        """Returns a Cartesian state representing the VNC difference between self and other, in position and velocity (with transport theorem).
Refer to dcm_from_vnc_to_inertial for details on the VNC frame.

# Algorithm
1. Compute the difference between `other` and `self`
2. Compute the VNC DCM of `self`
3. Rotate the difference into the VNC frame of `self`
4. Strip the astrodynamical information from the frame, enabling only computations from `CartesianState`"""

    def with_aop_deg(self, new_aop_deg: float) -> astro.Orbit:
        """Returns a copy of the state with a new AOP"""

    def with_apoapsis_periapsis_km(self, new_ra_km: float, new_rp_km: float) -> astro.Orbit:
        """Returns a copy of this state with the provided apoasis and periapsis"""

    def with_ecc(self, new_ecc: float) -> astro.Orbit:
        """Returns a copy of the state with a new ECC"""

    def with_inc_deg(self, new_inc_deg: float) -> astro.Orbit:
        """Returns a copy of the state with a new INC"""

    def with_raan_deg(self, new_raan_deg: float) -> astro.Orbit:
        """Returns a copy of the state with a new RAAN"""

    def with_sma_km(self, new_sma_km: float) -> astro.Orbit:
        """Returns a copy of the state with a new SMA"""

    def with_ta_deg(self, new_ta_deg: float) -> astro.Orbit:
        """Returns a copy of the state with a new TA"""

    def __eq__(self, value: object) -> bool:
        """Return self==value."""

    def __ge__(self, value: typing.Any) -> bool:
        """Return self>=value."""

    def __getnewargs__(self) -> tuple:...

    def __gt__(self, value: typing.Any) -> bool:
        """Return self>value."""

    def __le__(self, value: typing.Any) -> bool:
        """Return self<=value."""

    def __lt__(self, value: typing.Any) -> bool:
        """Return self<value."""

    def __ne__(self, value: object) -> bool:
        """Return self!=value."""



@typing.final
class SRPData:
    area_m2: float
    coeff_reflectivity: float

    def __init__(self, *args: typing.Any | None, **kwargs: typing.Any | None) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls, area_m2: typing.Any, coeff_reflectivity: typing.Any=None) -> SRPData:...

    @staticmethod
    def from_asn1(data: bytes) -> astro.SRPData:
        """Decodes an ASN.1 DER encoded byte array into an SRPData object."""

    def to_asn1(self) -> bytes:
        """Encodes this SRPData object into an ASN.1 DER encoded byte array."""



@typing.final
class Spacecraft:
    """A spacecraft state, composed of its orbit, its masses (dry, prop, extra, all in kg), its SRP configuration, its drag configuration, its thruster configuration, and its guidance mode.

Optionally, the spacecraft state can also store the state transition matrix from the start of the propagation until the current time (i.e. trajectory STM, not step-size STM)."""
    drag: typing.Any
    mass: typing.Any
    orbit: typing.Any
    srp: typing.Any

    def __init__(self, *args: typing.Any | None, **kwargs: typing.Any | None) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
A spacecraft state, composed of its orbit, its masses (dry, prop, extra, all in kg), its SRP configuration, its drag configuration, its thruster configuration, and its guidance mode.

Optionally, the spacecraft state can also store the state transition matrix from the start of the propagation until the current time (i.e. trajectory STM, not step-size STM)."""

    def __new__(cls, orbit: typing.Any, mass: typing.Any=None, srp: typing.Any=None, drag: typing.Any=None, thruster: typing.Any=None, mode: typing.Any=None) -> Spacecraft:
        """A spacecraft state, composed of its orbit, its masses (dry, prop, extra, all in kg), its SRP configuration, its drag configuration, its thruster configuration, and its guidance mode.

Optionally, the spacecraft state can also store the state transition matrix from the start of the propagation until the current time (i.e. trajectory STM, not step-size STM)."""

    @staticmethod
    def from_asn1(data: bytes) -> astro.Mass:
        """Decodes an ASN.1 DER encoded byte array into a Mass object."""

    def rss(self, other: typing.Any) -> typing.Any:
        """Returns the root sum square error between this spacecraft and the other, in kilometers for the position, kilometers per second in velocity, and kilograms in prop"""

    def to_asn1(self) -> bytes:
        """Encodes this Mass object into an ASN.1 DER encoded byte array."""

    def __eq__(self, value: object) -> bool:
        """Return self==value."""

    def __ge__(self, value: typing.Any) -> bool:
        """Return self>=value."""

    def __gt__(self, value: typing.Any) -> bool:
        """Return self>value."""

    def __le__(self, value: typing.Any) -> bool:
        """Return self<=value."""

    def __lt__(self, value: typing.Any) -> bool:
        """Return self<value."""

    def __ne__(self, value: object) -> bool:
        """Return self!=value."""



@typing.final
class Thruster:
    """Defines a thruster with a maximum isp and a maximum thrust."""
    isp_s: typing.Any
    thrust_N: typing.Any

    def __init__(self, *args: typing.Any | None, **kwargs: typing.Any | None) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
Defines a thruster with a maximum isp and a maximum thrust."""

    def __new__(cls, thrust_N: typing.Any, isp_s: typing.Any) -> Thruster:
        """Defines a thruster with a maximum isp and a maximum thrust."""

    def exhaust_velocity_m_s(self) -> typing.Any:
        """Returns the exhaust velocity v_e in meters per second"""
