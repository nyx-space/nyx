from __future__ import annotations
from anise import Almanac
from anise import astro
from anise import time
import numpy
import nyx_space.od
import typing

@typing.final
class CN0:
    """Carrier power to noise density (C/N0) for stochastic modeling of Doppler observables.

    IMPORTANT: C/N0 governs the thermal noise of phase-locked loops (PLL) tracking
    the primary unmodulated carrier wave to measure frequency shift (velocity). It represents
    the total power of the carrier signal over the noise spectral density.

    Applying S/N0 to Doppler observables artificially inflates modeled velocity noise,
    as it fails to account for the unmodulated carrier power explicitly reserved for
    phase tracking."""

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Carrier power to noise density (C/N0) for stochastic modeling of Doppler observables.

        IMPORTANT: C/N0 governs the thermal noise of phase-locked loops (PLL) tracking
        the primary unmodulated carrier wave to measure frequency shift (velocity). It represents
        the total power of the carrier signal over the noise spectral density.

        Applying S/N0 to Doppler observables artificially inflates modeled velocity noise,
        as it fails to account for the unmodulated carrier power explicitly reserved for
        phase tracking."""

    def __new__(cls) -> CN0:
        """Carrier power to noise density (C/N0) for stochastic modeling of Doppler observables.

        IMPORTANT: C/N0 governs the thermal noise of phase-locked loops (PLL) tracking
        the primary unmodulated carrier wave to measure frequency shift (velocity). It represents
        the total power of the carrier signal over the noise spectral density.

        Applying S/N0 to Doppler observables artificially inflates modeled velocity noise,
        as it fails to account for the unmodulated carrier power explicitly reserved for
        phase tracking."""
    Average: type = ...
    ManualDbHz: type = ...
    Poor: type = ...
    Strong: type = ...

@typing.final
class Cadence:
    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> Cadence: ...
    @staticmethod
    def continuous() -> Cadence: ...
    @staticmethod
    def from_asn1(data: bytes) -> Cadence:
        """Decodes an ASN.1 DER encoded byte array into a Cadence object."""

    @staticmethod
    def intermittent(on: time.Duration, off: time.Duration) -> Cadence:
        """Set an intermittent cadence with specific on and off durations."""

    def to_asn1(self) -> bytes:
        """Encodes this Cadence object into an ASN.1 DER encoded byte array."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class CarrierFreq:
    """Carrier frequency helper enum, typical values."""

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Carrier frequency helper enum, typical values."""

    def __new__(cls) -> CarrierFreq:
        """Carrier frequency helper enum, typical values."""
    KaBand: type = ...
    ManualHz: type = ...
    SBand: type = ...
    XBand: type = ...

@typing.final
class ChipRate:
    """An enum helper with typical chip rates."""

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        An enum helper with typical chip rates."""

    def __new__(cls) -> ChipRate:
        """An enum helper with typical chip rates."""
    High: type = ...
    Low: type = ...
    Lowest: type = ...
    ManualHz: type = ...
    StandardT4B: type = ...
    VeryHigh: type = ...

@typing.final
class ExportCfg:
    """Configuration for exporting from Nyx to local disk."""

    def __init__(
        self,
        timestamped: bool,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Configuration for exporting from Nyx to local disk."""

    def __new__(cls, timestamped: typing.Optional[bool] = False) -> ExportCfg:
        """Configuration for exporting from Nyx to local disk."""

    def __eq__(self, value: typing.Any) -> bool:
        """Return self==value."""

    def __ge__(self, value: typing.Any) -> bool:
        """Return self>=value."""

    def __gt__(self, value: typing.Any) -> bool:
        """Return self>value."""

    def __le__(self, value: typing.Any) -> bool:
        """Return self<=value."""

    def __lt__(self, value: typing.Any) -> bool:
        """Return self<value."""

    def __ne__(self, value: typing.Any) -> bool:
        """Return self!=value."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class FrameUid:
    """A unique frame reference that only contains enough information to build the actual Frame object.
    It cannot be used for any computations, is it be used in any structure apart from error structures."""

    force_inertial: typing.Any
    frozen_epoch: time.Epoch

    def __init__(
        self,
        ephemeris_id: int,
        orientation_id: int,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        A unique frame reference that only contains enough information to build the actual Frame object.
        It cannot be used for any computations, is it be used in any structure apart from error structures."""

    def __new__(cls, ephemeris_id: int, orientation_id: int) -> FrameUid:
        """A unique frame reference that only contains enough information to build the actual Frame object.
        It cannot be used for any computations, is it be used in any structure apart from error structures."""

    @staticmethod
    def from_frame(frame: astro.Frame) -> astro.FrameUid:
        """Creates a FrameUid from a Frame"""

    def to_frame(self) -> astro.Frame:
        """Converts this FrameUid to a Frame"""

    def __eq__(self, value: typing.Any) -> bool:
        """Return self==value."""

    def __ge__(self, value: typing.Any) -> bool:
        """Return self>=value."""

    def __gt__(self, value: typing.Any) -> bool:
        """Return self>value."""

    def __le__(self, value: typing.Any) -> bool:
        """Return self<=value."""

    def __lt__(self, value: typing.Any) -> bool:
        """Return self<value."""

    def __ne__(self, value: typing.Any) -> bool:
        """Return self!=value."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class GaussMarkov:
    """A first order Gauss-Markov process for modeling biases as described in section 5.2.4 of the NASA Best Practices for Navigation Filters (D'Souza et al.).

    The process is defined by the following stochastic differential equation:

    \\dot{b(t)} = -1/τ * b(t) + w(t)

    Programmatically, it's calculated by sampling from b(t) ~ 𝓝(0, p_b(t)), where

    p_b(t) = exp((-2 / τ) * (t - t_0)) * p_b(t_0) + s(t - t_0)

    s(t - t_0) = ((q * τ) / 2) * (1 - exp((-2 / τ) * (t - t_0)))

    ## JPL DESCANSO Deep Space Network (DSN) Defaults

    - Range: 60 cm process noise over a 60 second average (tau, half life)
    - Doppler: 0.03 mm/s process noise over a 60 second average (tau, half life)"""

    constant: typing.Any
    init_sample: typing.Any
    prev_epoch: typing.Any
    process_noise: typing.Any
    tau: typing.Any

    def __init__(
        self,
        tau: time.Duration,
        process_noise: float,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        A first order Gauss-Markov process for modeling biases as described in section 5.2.4 of the NASA Best Practices for Navigation Filters (D'Souza et al.).

        The process is defined by the following stochastic differential equation:

        \\dot{b(t)} = -1/τ * b(t) + w(t)

        Programmatically, it's calculated by sampling from b(t) ~ 𝓝(0, p_b(t)), where

        p_b(t) = exp((-2 / τ) * (t - t_0)) * p_b(t_0) + s(t - t_0)

        s(t - t_0) = ((q * τ) / 2) * (1 - exp((-2 / τ) * (t - t_0)))

        ## JPL DESCANSO Deep Space Network (DSN) Defaults

        - Range: 60 cm process noise over a 60 second average (tau, half life)
        - Doppler: 0.03 mm/s process noise over a 60 second average (tau, half life)"""

    def __new__(cls, tau: time.Duration, process_noise: float) -> GaussMarkov:
        """A first order Gauss-Markov process for modeling biases as described in section 5.2.4 of the NASA Best Practices for Navigation Filters (D'Souza et al.).

        The process is defined by the following stochastic differential equation:

        \\dot{b(t)} = -1/τ * b(t) + w(t)

        Programmatically, it's calculated by sampling from b(t) ~ 𝓝(0, p_b(t)), where

        p_b(t) = exp((-2 / τ) * (t - t_0)) * p_b(t_0) + s(t - t_0)

        s(t - t_0) = ((q * τ) / 2) * (1 - exp((-2 / τ) * (t - t_0)))

        ## JPL DESCANSO Deep Space Network (DSN) Defaults

        - Range: 60 cm process noise over a 60 second average (tau, half life)
        - Doppler: 0.03 mm/s process noise over a 60 second average (tau, half life)"""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class GroundStation:
    """GroundStation defines a one-way or two-way ranging and doppler station. Set the integration time for two-way."""

    integration_time: typing.Any
    light_time_correction: typing.Any
    location: typing.Any
    measurement_types: typing.Any
    name: typing.Any
    stochastic_noises: typing.Any
    timestamp_noise_s: typing.Any

    def __init__(
        self,
        name: str,
        location: astro.Location,
        stochastic_noises: dict[MeasurementType, StochasticNoise],
        integration_time: time.Duration | None,
        light_time_correction: bool | None,
        timestamp_noise_s: StochasticNoise | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        GroundStation defines a one-way or two-way ranging and doppler station. Set the integration time for two-way."""

    def __new__(
        cls,
        name: str,
        location: astro.Location,
        stochastic_noises: dict[MeasurementType, StochasticNoise],
        integration_time: time.Duration | None = None,
        light_time_correction: bool | None = False,
        timestamp_noise_s: StochasticNoise | None = None,
    ) -> GroundStation:
        """GroundStation defines a one-way or two-way ranging and doppler station. Set the integration time for two-way."""

    def add_measurement_type(
        self, msr_type: MeasurementType, noise: StochasticNoise
    ) -> None:
        """Add a measurement type with stochastic noise."""

    def azimuth_elevation_of(
        self, rx: astro.Orbit, obstructing_body: astro.Frame | None, almanac: Almanac
    ) -> astro.AzElRange:
        """Computes the azimuth and elevation of the provided object seen from this ground station, both in degrees.
        This is a shortcut to almanac.azimuth_elevation_range_sez."""

    def clear_measurement_types(self) -> None:
        """Clear all measurement types"""

    def clear_stochastic_noises(self) -> None:
        """Clear stochastic noises"""

    @staticmethod
    def dump_many_yaml(stations: list[GroundStation], path: str) -> None:
        """Dump multiple GroundStations to a YAML file."""

    @staticmethod
    def dumps_many_yaml(stations: list[GroundStation]) -> str:
        """Dump multiple GroundStations to a YAML string."""

    @staticmethod
    def from_asn1(data: bytes) -> GroundStation:
        """Decodes an ASN.1 DER encoded byte array into a GroundStation object."""

    @staticmethod
    def from_yaml(yaml_str: str) -> GroundStation:
        """Load GroundStation from a YAML string."""

    def get_stochastic_noise(self, m_type: MeasurementType) -> StochasticNoise | None:
        """Get stochastic noise for a measurement type."""

    @staticmethod
    def load_many_yaml(path: str) -> list[GroundStation]:
        """Load multiple GroundStations from a YAML file."""

    @staticmethod
    def loads_many_yaml(yaml_str: str) -> list[GroundStation]:
        """Load multiple GroundStations from a YAML string."""

    def remove_measurement_type(self, msr_type: MeasurementType) -> bool:
        """Remove a measurement type."""

    def remove_stochastic_noise(
        self, m_type: MeasurementType
    ) -> StochasticNoise | None:
        """Remove stochastic noise for a measurement type."""

    def set_stochastic_noise(
        self, m_type: MeasurementType, noise: StochasticNoise
    ) -> None:
        """Set stochastic noise for a measurement type."""

    def to_asn1(self) -> bytes:
        """Encodes this GroundStation object into an ASN.1 DER encoded byte array."""

    def to_orbit(self, epoch: time.Epoch, almanac: Almanac) -> astro.Orbit:
        """Return this ground station as an orbit in its current frame"""

    def to_yaml(self) -> str: ...
    def __eq__(self, value: typing.Any) -> bool:
        """Return self==value."""

    def __ge__(self, value: typing.Any) -> bool:
        """Return self>=value."""

    def __gt__(self, value: typing.Any) -> bool:
        """Return self>value."""

    def __le__(self, value: typing.Any) -> bool:
        """Return self<=value."""

    def __lt__(self, value: typing.Any) -> bool:
        """Return self<value."""

    def __ne__(self, value: typing.Any) -> bool:
        """Return self!=value."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class GroundTrackingArcSim:
    """Simulated tracking architecture for a spacecraft."""

    configs: typing.Any
    devices: typing.Any

    def __init__(
        self,
        devices: dict[str, GroundStation],
        trajectory: PyTrajectory,
        configs: dict[str, TrkConfig],
        seed: int | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Simulated tracking architecture for a spacecraft."""

    def __new__(
        cls,
        devices: dict[str, GroundStation],
        trajectory: PyTrajectory,
        configs: dict[str, TrkConfig],
        seed: int | None = None,
    ) -> GroundTrackingArcSim:
        """Simulated tracking architecture for a spacecraft."""

    def build_schedule(self, almanac: Almanac) -> None:
        """Builds a schedule using the generate_schedule function, and set that schedule in this instance's configuration."""

    def generate_measurements(self, almanac: Almanac) -> TrackingDataArc:
        """Simulates operational tracking data across predefined tracking strands.

        This function strictly demands that a schedule already exists (stored in the `config` field).
        If a device is configured as a scheduler but lacks pre-computed
        strands, this function will raise an error rather than implicitly hallucinating a tracking
        pass. Call `generate_schedule` to build a schedule first.
        For each tracking device, the trajectory is sampled at the specific hardware rate,
        synthesizing measurements only when the spacecraft is visible.

        :raises ConfigError: If a scheduling configuration is present but the schedule was not built prior to execution."""

    def generate_schedule(self, almanac: Almanac) -> dict[str, TrkConfig]:
        """Builds the schedule provided the config.

        # Algorithm

        1. For each tracking device:
        2. Find when the vehicle elevation above ground station mask is greater or equal to zero, and use that as the first start of the first tracking arc for this station
        3. Find when the vehicle drops below the mask, after that initial epoch
        4. Repeat 2, 3 until the end of the trajectory
        5. Build each of these as "tracking strands" for this tracking device.
        6. Organize all of the built tracking strands chronologically.
        7. Iterate through all of the strands to adjust for tracker Greedy/Eager configuration.
        `Greedy` trackers will delay the start of subsequent station contacts, whereas `Eager` trackers will terminate
        current tracking strands prematurely to allow the next station to acquire.
        :raises AnalysisError: If underlying location dataset injection or visibility computation fails."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class Handoff:
    """Defines the handoff from a current ground station to the next one that is visible to prevent overlapping of measurements"""

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Defines the handoff from a current ground station to the next one that is visible to prevent overlapping of measurements"""

    def __new__(cls) -> Handoff:
        """Defines the handoff from a current ground station to the next one that is visible to prevent overlapping of measurements"""

    def __eq__(self, value: typing.Any) -> bool:
        """Return self==value."""

    def __ge__(self, value: typing.Any) -> bool:
        """Return self>=value."""

    def __gt__(self, value: typing.Any) -> bool:
        """Return self>value."""

    def __int__(self) -> None:
        """int(self)"""

    def __le__(self, value: typing.Any) -> bool:
        """Return self<=value."""

    def __lt__(self, value: typing.Any) -> bool:
        """Return self<value."""

    def __ne__(self, value: typing.Any) -> bool:
        """Return self!=value."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""
    Eager: Handoff = ...
    Greedy: Handoff = ...
    Overlap: Handoff = ...

@typing.final
class KalmanVariant:
    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> KalmanVariant: ...
    def __int__(self) -> None:
        """int(self)"""

    def __repr__(self) -> str:
        """Return repr(self)."""
    DeviationTracking: KalmanVariant = ...
    ReferenceUpdate: KalmanVariant = ...

@typing.final
class Location:
    """Location is defined by its latitude, longitude, height above the geoid, mean angular rotation of the geoid, and a frame UID.
    If the location includes a terrain mask, it will be used for obstruction checks when computing azimuth and elevation.
    **Note:** The mean Earth angular velocity is `0.004178079012116429` deg/s."""

    height_km: float
    latitude_deg: float
    longitude_deg: float
    terrain_mask: list
    terrain_mask_ignored: bool

    def __init__(
        self,
        latitude_deg: float,
        longitude_deg: float,
        height_km: float,
        frame: astro.FrameUid,
        terrain_mask: list,
        terrain_mask_ignored: bool,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Location is defined by its latitude, longitude, height above the geoid, mean angular rotation of the geoid, and a frame UID.
        If the location includes a terrain mask, it will be used for obstruction checks when computing azimuth and elevation.
        **Note:** The mean Earth angular velocity is `0.004178079012116429` deg/s."""

    def __new__(
        cls,
        latitude_deg: float,
        longitude_deg: float,
        height_km: float,
        frame: astro.FrameUid,
        terrain_mask: list,
        terrain_mask_ignored: bool,
    ) -> Location:
        """Location is defined by its latitude, longitude, height above the geoid, mean angular rotation of the geoid, and a frame UID.
        If the location includes a terrain mask, it will be used for obstruction checks when computing azimuth and elevation.
        **Note:** The mean Earth angular velocity is `0.004178079012116429` deg/s."""

    def elevation_mask_at_azimuth_deg(self, azimuth_deg: float) -> float:
        """Returns the elevation mask at the provided azimuth, does NOT account for whether the mask is ignored or not."""

    @staticmethod
    def from_dhall(repr: str) -> astro.Location:
        """Loads a Location from its Dhall representation"""

    def to_dhall(self) -> str:
        """Returns the Dhall representation of this Location"""

    def __eq__(self, value: typing.Any) -> bool:
        """Return self==value."""

    def __ge__(self, value: typing.Any) -> bool:
        """Return self>=value."""

    def __gt__(self, value: typing.Any) -> bool:
        """Return self>value."""

    def __le__(self, value: typing.Any) -> bool:
        """Return self<=value."""

    def __lt__(self, value: typing.Any) -> bool:
        """Return self<value."""

    def __ne__(self, value: typing.Any) -> bool:
        """Return self!=value."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class Measurement:
    """A type-agnostic simultaneous measurement storage structure. Allows storing any number of simultaneous measurement of a given taker.

    Note that two measurements are considered equal if the tracker and epoch match exactly, and if both have the same measurement types,
    and those measurements are equal to within 1e-10 (this allows for some leeway in TDM producers)."""

    def __init__(
        self,
        tracker: str,
        epoch: time.Epoch,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        A type-agnostic simultaneous measurement storage structure. Allows storing any number of simultaneous measurement of a given taker.

        Note that two measurements are considered equal if the tracker and epoch match exactly, and if both have the same measurement types,
        and those measurements are equal to within 1e-10 (this allows for some leeway in TDM producers)."""

    def __new__(cls, tracker: str, epoch: time.Epoch) -> Measurement:
        """A type-agnostic simultaneous measurement storage structure. Allows storing any number of simultaneous measurement of a given taker.

        Note that two measurements are considered equal if the tracker and epoch match exactly, and if both have the same measurement types,
        and those measurements are equal to within 1e-10 (this allows for some leeway in TDM producers)."""

    def correct(self, msr_type: MeasurementType, correction: float) -> None:
        """Correct the provided measurement type with the provided correction, if that measurement type is available"""

    def observation(self, msr_type: MeasurementType) -> float | None:
        """Returns the floating point value of this observation if this measurement contains the provided measurement type"""

    def push(self, msr_type: MeasurementType, msr_value: float) -> None:
        """Push a measurement type and value."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class MeasurementType:
    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> MeasurementType: ...
    def __int__(self) -> None:
        """int(self)"""

    def __repr__(self) -> str:
        """Return repr(self)."""
    Azimuth: MeasurementType = ...
    Doppler: MeasurementType = ...
    Elevation: MeasurementType = ...
    Range: MeasurementType = ...
    ReceiveFrequency: MeasurementType = ...
    TransmitFrequency: MeasurementType = ...
    TransmitFrequencyRate: MeasurementType = ...
    X: MeasurementType = ...
    Y: MeasurementType = ...
    Z: MeasurementType = ...

@typing.final
class PositionDevice:
    """Position device can be used to post-filter position measurements from GNSS/GPS devices.

    For GNSS devices, ensure to set the PositionDevice frame in the ITRF93 frame, the closest realization
    to WSG84."""

    frame: typing.Any
    name: typing.Any

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Position device can be used to post-filter position measurements from GNSS/GPS devices.

        For GNSS devices, ensure to set the PositionDevice frame in the ITRF93 frame, the closest realization
        to WSG84."""

    def __new__(cls, name: typing.Any, frame: typing.Any) -> PositionDevice:
        """Position device can be used to post-filter position measurements from GNSS/GPS devices.

        For GNSS devices, ensure to set the PositionDevice frame in the ITRF93 frame, the closest realization
        to WSG84."""

    def with_noise(
        self, msr_type: MeasurementType, noise: StochasticNoise
    ) -> PositionDevice:
        """Add a measurement type with stochastic noise."""

    def __eq__(self, value: typing.Any) -> bool:
        """Return self==value."""

    def __ge__(self, value: typing.Any) -> bool:
        """Return self>=value."""

    def __gt__(self, value: typing.Any) -> bool:
        """Return self>value."""

    def __le__(self, value: typing.Any) -> bool:
        """Return self<=value."""

    def __lt__(self, value: typing.Any) -> bool:
        """Return self<value."""

    def __ne__(self, value: typing.Any) -> bool:
        """Return self!=value."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class PositionResidual:
    epoch: typing.Any
    postfit: typing.Any
    prefit: typing.Any
    ratio: typing.Any
    rejected: typing.Any
    tracker: typing.Any

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> PositionResidual: ...
    def computed_obs(self, msr_type: MeasurementType) -> float | None:
        """Returns the computed/expected observation for this measurement type, if available"""

    def nis(self) -> float:
        """Returns the normalized innovation squared (NIS) as the norm squares of the whitened residual"""

    def real_obs(self, msr_type: MeasurementType) -> float | None:
        """Returns the real observation for this measurement type, if available"""

    def whitened_residual(self, msr_type: MeasurementType) -> float | None:
        """Returns the whitened residual for this measurement type, if available"""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class PositionTrackingArcSim:
    """Simulated tracking architecture for a spacecraft using position tracking devices."""

    configs: typing.Any
    devices: typing.Any

    def __init__(
        self,
        devices: dict[str, PositionDevice],
        trajectory: PyTrajectory,
        configs: dict[str, TrkConfig],
        seed: int | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Simulated tracking architecture for a spacecraft using position tracking devices."""

    def __new__(
        cls,
        devices: dict[str, PositionDevice],
        trajectory: PyTrajectory,
        configs: dict[str, TrkConfig],
        seed: int | None = None,
    ) -> PositionTrackingArcSim:
        """Simulated tracking architecture for a spacecraft using position tracking devices."""

    def generate_measurements(self, almanac: Almanac) -> TrackingDataArc:
        """Simulates operational tracking data across predefined tracking strands.

        :raises ConfigError: If a scheduling configuration is present but the schedule was not built prior to execution."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class ProcessNoise:
    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> ProcessNoise: ...
    @staticmethod
    def from_accel_m_s2(
        ax_m_s2: float,
        ay_m_s2: float,
        az_m_s2: float,
        disable_time: time.Duration,
        local_frame: astro.LocalFrame | None,
        x_decay_s: float | None,
        y_decay_s: float | None,
        z_decay_s: float | None,
    ) -> ProcessNoise:
        """Create process noise from acceleration standard deviations with optional exponential decay."""

    @staticmethod
    def from_velocity_m_s(
        vx_m_s: float,
        vy_m_s: float,
        vz_m_s: float,
        noise_duration: time.Duration,
        disable_time: time.Duration,
        local_frame: astro.LocalFrame | None,
    ) -> ProcessNoise:
        """Create process noise from velocity standard deviations."""

    def __eq__(self, value: typing.Any) -> bool:
        """Return self==value."""

    def __ge__(self, value: typing.Any) -> bool:
        """Return self>=value."""

    def __gt__(self, value: typing.Any) -> bool:
        """Return self>value."""

    def __le__(self, value: typing.Any) -> bool:
        """Return self<=value."""

    def __lt__(self, value: typing.Any) -> bool:
        """Return self<value."""

    def __ne__(self, value: typing.Any) -> bool:
        """Return self!=value."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class Residual:
    epoch: typing.Any
    postfit: typing.Any
    prefit: typing.Any
    ratio: typing.Any
    rejected: typing.Any
    tracker: typing.Any

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> Residual: ...
    def computed_obs(self, msr_type: MeasurementType) -> float | None:
        """Returns the computed/expected observation for this measurement type, if available"""

    def nis(self) -> float:
        """Returns the normalized innovation squared (NIS) as the norm squares of the whitened residual"""

    def real_obs(self, msr_type: MeasurementType) -> float | None:
        """Returns the real observation for this measurement type, if available"""

    def whitened_residual(self, msr_type: MeasurementType) -> float | None:
        """Returns the whitened residual for this measurement type, if available"""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class SN0:
    """Signal power to noise density (S/N0) for stochastic modeling of ranging observables.

    IMPORTANT: S/N0 governs the thermal noise of delay-locked loops (DLL) tracking
    the modulated ranging code or tone. Deep space architectures rely on phase modulation
    with a residual carrier. The total transmitted power is allocated fractionally among the
    main carrier wave, the telemetry subcarrier, and the ranging code, dictated by the modulation index.

    Because the power available for ranging is strictly a subset of the total carrier power,
    S/N0 <= C/N0. Applying C/N0 to ranging observables artificially suppresses the modeled thermal
    noise, yielding an overly optimistic covariance bound that ignores spacecraft power division."""

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Signal power to noise density (S/N0) for stochastic modeling of ranging observables.

        IMPORTANT: S/N0 governs the thermal noise of delay-locked loops (DLL) tracking
        the modulated ranging code or tone. Deep space architectures rely on phase modulation
        with a residual carrier. The total transmitted power is allocated fractionally among the
        main carrier wave, the telemetry subcarrier, and the ranging code, dictated by the modulation index.

        Because the power available for ranging is strictly a subset of the total carrier power,
        S/N0 <= C/N0. Applying C/N0 to ranging observables artificially suppresses the modeled thermal
        noise, yielding an overly optimistic covariance bound that ignores spacecraft power division."""

    def __new__(cls) -> SN0:
        """Signal power to noise density (S/N0) for stochastic modeling of ranging observables.

        IMPORTANT: S/N0 governs the thermal noise of delay-locked loops (DLL) tracking
        the modulated ranging code or tone. Deep space architectures rely on phase modulation
        with a residual carrier. The total transmitted power is allocated fractionally among the
        main carrier wave, the telemetry subcarrier, and the ranging code, dictated by the modulation index.

        Because the power available for ranging is strictly a subset of the total carrier power,
        S/N0 <= C/N0. Applying C/N0 to ranging observables artificially suppresses the modeled thermal
        noise, yielding an overly optimistic covariance bound that ignores spacecraft power division."""
    Average: type = ...
    ManualDbHz: type = ...
    Poor: type = ...
    Strong: type = ...

@typing.final
class Scheduler:
    """A scheduler allows building a scheduling of spaceraft tracking for a set of ground stations."""

    cadence: typing.Any
    handoff: typing.Any
    min_samples: typing.Any
    sample_alignment: typing.Any

    def __init__(
        self,
        handoff: Handoff,
        cadence: Cadence | None,
        min_samples: int,
        sample_alignment: time.Duration | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        A scheduler allows building a scheduling of spaceraft tracking for a set of ground stations."""

    def __new__(
        cls,
        handoff: typing.Optional[Handoff] = ...,
        cadence: Cadence | None = None,
        min_samples: typing.Optional[int] = 10,
        sample_alignment: time.Duration | None = None,
    ) -> Scheduler:
        """A scheduler allows building a scheduling of spaceraft tracking for a set of ground stations."""

    @staticmethod
    def from_asn1(data: bytes) -> Scheduler:
        """Decodes an ASN.1 DER encoded byte array into a Scheduler object."""

    def to_asn1(self) -> bytes:
        """Encodes this Scheduler object into an ASN.1 DER encoded byte array."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class SigmaRejection:
    """Reject measurements if the prefit is greater than the provided sigmas deviation from the measurement noise.

    # Important
    Some software, like ODTK, processes each measurement as a scalar. Nyx can process the measurements together.
    As such, if the prefit on range is bad, then the Doppler measurement with the same time stamp will also be rejected.
    This can lead to better convergence of the filter, and more appropriate results."""

    num_sigmas: typing.Any

    def __init__(
        self,
        num_sigmas: float,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Reject measurements if the prefit is greater than the provided sigmas deviation from the measurement noise.

        # Important
        Some software, like ODTK, processes each measurement as a scalar. Nyx can process the measurements together.
        As such, if the prefit on range is bad, then the Doppler measurement with the same time stamp will also be rejected.
        This can lead to better convergence of the filter, and more appropriate results."""

    def __new__(cls, num_sigmas: float) -> SigmaRejection:
        """Reject measurements if the prefit is greater than the provided sigmas deviation from the measurement noise.

        # Important
        Some software, like ODTK, processes each measurement as a scalar. Nyx can process the measurements together.
        As such, if the prefit on range is bad, then the Doppler measurement with the same time stamp will also be rejected.
        This can lead to better convergence of the filter, and more appropriate results."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class SpacecraftEstimate:
    covariance: numpy.ndarray
    nominal_state: typing.Any
    predicted: typing.Any
    state: typing.Any
    state_deviations: typing.Any

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> SpacecraftEstimate: ...
    @staticmethod
    def from_diag(nominal: Spacecraft, diag: numpy.ndarray) -> SpacecraftEstimate:
        """Initializes a new filter estimate from the nominal state (not dispersed) and the diagonal of the covariance"""

    @staticmethod
    def from_dispersions(
        nominal_state: Spacecraft,
        dispersions: list[StateDispersion],
        seed: int | None = None,
    ) -> SpacecraftEstimate:
        """Generates an initial Kalman filter state estimate dispersed from the nominal state using the provided standard deviation parameters.

        The resulting estimate will have a diagonal covariance matrix constructed from the variances of each parameter."""

    def to_random_variable(self) -> MvnSpacecraft:
        """Builds a multivariate random variable spacecraft from this estimate's nominal state and covariance, zero mean."""

    def within_3sigma(self) -> bool:
        """Returns whether this estimate is within three sigmas"""

    def within_sigma(self, sigma: float) -> bool:
        """Returns whether this estimate is within some bound
        The 68-95-99.7 rule is a good way to assess whether the filter is operating normally"""

    def __eq__(self, value: typing.Any) -> bool:
        """Return self==value."""

    def __ge__(self, value: typing.Any) -> bool:
        """Return self>=value."""

    def __gt__(self, value: typing.Any) -> bool:
        """Return self>value."""

    def __le__(self, value: typing.Any) -> bool:
        """Return self<=value."""

    def __lt__(self, value: typing.Any) -> bool:
        """Return self<value."""

    def __ne__(self, value: typing.Any) -> bool:
        """Return self!=value."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class SpacecraftODProcess:
    """Orbit determination process for a spacecraft."""

    sigma_rejection: typing.Any
    variant: typing.Any

    def __init__(
        self,
        prop: Propagator,
        kf_variant: KalmanVariant,
        devices: dict[str, GroundStation],
        sigma_reject: SigmaRejection | None,
        process_noise: ProcessNoise | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Orbit determination process for a spacecraft."""

    def __new__(
        cls,
        prop: Propagator,
        kf_variant: KalmanVariant,
        devices: dict[str, GroundStation],
        sigma_reject: SigmaRejection | None = ...,
        process_noise: ProcessNoise | None = None,
    ) -> SpacecraftODProcess:
        """Orbit determination process for a spacecraft."""

    def predict_for(
        self, initial_estimate: SpacecraftEstimate, duration: time.Duration
    ) -> SpacecraftODSolution:
        """Perform a time update. Continuously predicts the trajectory for the provided duration, with covariance mapping at each step."""

    def predict_until(
        self, initial_estimate: SpacecraftEstimate, end_epoch: time.Epoch
    ) -> SpacecraftODSolution:
        """Perform a time update. Continuously predicts the trajectory until the provided end epoch, with covariance mapping at each step."""

    def process_arc(
        self, initial_estimate: SpacecraftEstimate, arc: TrackingDataArc
    ) -> SpacecraftODSolution:
        """Process the provided tracking arc for this orbit determination process."""

@typing.final
class SpacecraftODSolution:
    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> SpacecraftODSolution: ...
    def accepted_residuals(self) -> list[Residual]: ...
    @staticmethod
    def from_parquet(
        path: str, devices: dict[str, GroundStation]
    ) -> SpacecraftODSolution:
        """Reconstruct an ODSolution from a parquet file."""

    def is_filter_run(self) -> bool: ...
    def is_nees_consistent(
        self, truth_traj: PyTrajectory, alpha: float | None = None
    ) -> bool:
        """Checks whether the filter estimates are statistically consistent
        by performing a Chi-squared test on the Normalized Estimation Error Squared (NEES).

        For each estimate, NEES is computed as:
        ```text
        error^T * P^-1 * error
        ```
        where `error` is the difference between the estimated state and the true state,
        and `P` is the estimated state covariance matrix.

        The sum of NEES values should fall within the confidence interval of a
        Chi-squared distribution with degrees of freedom `k = n * dim`, where `n`
        is the number of estimates and `dim` is the state dimension.

        Returns Ok(true) if the filter is consistent, Ok(false) if the filter
        is over-confident or under-confident, or an error if no estimates are available."""

    def is_nis_consistent(self, alpha: float | None = None) -> bool:
        """Checks whether the filter innovations are statistically consistent
        by performing a Chi-squared test on the Normalized Innovation Squared (NIS).

        For each accepted residual, NIS is computed as:
        ```text
        prefit^T * S_k^-1 * prefit
        ```

        The sum of NIS values should fall within the confidence interval of a
        Chi-squared distribution with degrees of freedom `k = n * m`, where `n`
        is the number of residuals and `m` is the measurement dimension.

        Returns Ok(true) if the filter is consistent, Ok(false) if the filter
        is over-confident or under-confident, or an error if no residuals are available."""

    def is_normal(self, alpha: float | None = None) -> bool:
        """Checks whether the whitened residuals of the accepted residuals pass a normality test at a given significance level `alpha`, default to 0.05.

        This uses a simplified KS-test threshold: D_alpha = c(α) / √n.
        For example, for α = 0.05, c(α) is approximately 1.36.
        α=0.05 means a 5% probability of Type I error (incorrectly rejecting the null hypothesis when it is true).
        This threshold is standard in many statistical tests to balance sensitivity and false positives

        The computation of the c(alpha) is from https://en.wikipedia.org/w/index.php?title=Kolmogorov%E2%80%93Smirnov_test&oldid=1280260701#Two-sample_Kolmogorov%E2%80%93Smirnov_test

        Returns Ok(true) if the residuals are consistent with a normal distribution,
        Ok(false) if not, or None if no residuals are available."""

    def is_smoother_run(self) -> bool: ...
    def ks_test_normality(self) -> float:
        """Computes the Kolmogorov–Smirnov statistic for the aggregated residual ratios of the accepted residuals.

        Returns Ok(ks_statistic) if residuals are available."""

    def nees_consistency(
        self, truth_traj: PyTrajectory, alpha: float | None = None
    ) -> NormalizedConsistency:
        """Checks whether the filter estimates are statistically consistent
        by performing a Chi-squared test on the Normalized Estimation Error Squared (NEES).

        For each estimate, NEES is computed as:
        ```text
        error^T * P^-1 * error
        ```
        where `error` is the difference between the estimated state and the true state,
        and `P` is the estimated state covariance matrix.

        The sum of NEES values should fall within the confidence interval of a
        Chi-squared distribution with degrees of freedom `k = n * dim`, where `n`
        is the number of estimates and `dim` is the state dimension.

        Returns Ok(true) if the filter is consistent, Ok(false) if the filter
        is over-confident or under-confident, or an error if no estimates are available."""

    def nis_consistency(self, alpha: float | None = None) -> NormalizedConsistency:
        """Checks whether the filter innovations are statistically consistent
        by performing a Chi-squared test on the Normalized Innovation Squared (NIS).

        For each accepted residual, NIS is computed as:
        ```text
        prefit^T * S_k^-1 * prefit
        ```

        The sum of NIS values should fall within the confidence interval of a
        Chi-squared distribution with degrees of freedom `k = n * m`, where `n`
        is the number of residuals and `m` is the measurement dimension.

        Returns Ok(true) if the filter is consistent, Ok(false) if the filter
        is over-confident or under-confident, or an error if no residuals are available."""

    def rejected_residuals(self) -> list[Residual]: ...
    def residual_ratio_within_threshold(self, threshold: float) -> float:
        """Computes the fraction of residual ratios that lie within ±threshold."""

    def rms_postfit_residuals(self) -> float:
        """Returns the root mean square of the postfit residuals"""

    def rms_prefit_residuals(self) -> float:
        """Returns the root mean square of the prefit residuals"""

    def rms_residual_ratios(self) -> float:
        """Returns the root mean square of the prefit residual ratios"""

    def smooth(self, almanac: Almanac) -> SpacecraftODSolution:
        """Smoothes this OD solution, returning a new OD solution and the filter-smoother consistency ratios, with updated **postfit** residuals, and where the ratio now represents the filter-smoother consistency ratio.

        Notes:
        1. Gains will be scrubbed because the smoother process does not recompute the gain.
        2. Prefit residuals, ratios, and measurement covariances are not updated, as these depend on the filtering process.
        3. Note: this function consumes the current OD solution to prevent reusing the wrong one.


        To assess whether the smoothing process improved the solution, compare the RMS of the postfit residuals from the filter and the smoother process.

        # Filter-Smoother consistency ratio

        The **filter-smoother consistency ratio** is used to evaluate the consistency between the state estimates produced by a filter (e.g., Kalman filter) and a smoother.
        This ratio is called "filter smoother consistency test" in the ODTK MathSpec.

        It is computed as follows:

        #### 1. Define the State Estimates
        **Filter state estimate**:
        $ \\hat{X}_{f,k} $
        This is the state estimate at time step $ k $ from the filter.

        **Smoother state estimate**:
        $ \\hat{X}_{s,k} $
        This is the state estimate at time step $ k $ from the smoother.

        #### 2. Define the Covariances

        **Filter covariance**:
        $ P_{f,k} $
        This is the covariance of the state estimate at time step $ k $ from the filter.

        **Smoother covariance**:
        $ P_{s,k} $
        This is the covariance of the state estimate at time step $ k $ from the smoother.

        #### 3. Compute the Differences

        **State difference**:
        $ \\Delta X_k = \\hat{X}_{s,k} - \\hat{X}_{f,k} $

        **Covariance difference**:
        $ \\Delta P_k = P_{s,k} - P_{f,k} $

        #### 4. Calculate the Consistency Ratio
        For each element $ i $ of the state vector, compute the ratio:

        $$
        R_{i,k} = \\frac{\\Delta X_{i,k}}{\\sqrt{\\Delta P_{i,k}}}
        $$

        #### 5. Evaluate Consistency
        - If $ |R_{i,k}| \\leq 3 $ for all $ i $ and $ k $, the filter-smoother consistency test is satisfied, indicating good consistency.
        - If $ |R_{i,k}| > 3 $ for any $ i $ or $ k $, the test fails, suggesting potential modeling inconsistencies or issues with the estimation process."""

    def to_ephemeris(self, object_id: str) -> astro.Ephemeris:
        """Export to an ANISE ephemeris, which can be converted to a CCSDS OEM"""

    def to_parquet(self, path: str, cfg: ExportCfg) -> str:
        """Export OD solutions, gains, ratios, residuals, sigmas, etc. to parquet"""

    def to_traj(self) -> Trajectory: ...
    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class SpacecraftPositionODProcess:
    """Orbit determination process for a spacecraft using position tracking devices."""

    sigma_rejection: typing.Any
    variant: typing.Any

    def __init__(
        self,
        prop: Propagator,
        kf_variant: KalmanVariant,
        devices: dict[str, PositionDevice],
        sigma_reject: SigmaRejection | None,
        process_noise: ProcessNoise | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Orbit determination process for a spacecraft using position tracking devices."""

    def __new__(
        cls,
        prop: Propagator,
        kf_variant: KalmanVariant,
        devices: dict[str, PositionDevice],
        sigma_reject: SigmaRejection | None = ...,
        process_noise: ProcessNoise | None = None,
    ) -> SpacecraftPositionODProcess:
        """Orbit determination process for a spacecraft using position tracking devices."""

    def predict_for(
        self, initial_estimate: SpacecraftEstimate, duration: time.Duration
    ) -> SpacecraftPositionODSolution:
        """Perform a time update. Continuously predicts the trajectory for the provided duration, with covariance mapping at each step."""

    def predict_until(
        self, initial_estimate: SpacecraftEstimate, end_epoch: time.Epoch
    ) -> SpacecraftPositionODSolution:
        """Perform a time update. Continuously predicts the trajectory until the provided end epoch, with covariance mapping at each step."""

    def process_arc(
        self, initial_estimate: SpacecraftEstimate, arc: TrackingDataArc
    ) -> SpacecraftPositionODSolution:
        """Process the provided tracking arc for this orbit determination process."""

@typing.final
class SpacecraftPositionODSolution:
    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> SpacecraftPositionODSolution: ...
    def accepted_residuals(self) -> list[PositionResidual]: ...
    @staticmethod
    def from_parquet(
        path: str, devices: dict[str, PositionDevice]
    ) -> SpacecraftPositionODSolution:
        """Reconstruct an ODSolution from a parquet file."""

    def is_filter_run(self) -> bool: ...
    def is_normal(self, alpha: float | None = None) -> bool:
        """Checks whether the whitened residuals of the accepted residuals pass a normality test at a given significance level `alpha`, default to 0.05."""

    def is_smoother_run(self) -> bool: ...
    def ks_test_normality(self) -> float:
        """Computes the Kolmogorov–Smirnov statistic for the aggregated residual ratios of the accepted residuals.

        Returns Ok(ks_statistic) if residuals are available."""

    def nees_consistency(
        self, truth_traj: PyTrajectory, alpha: float | None = None
    ) -> NormalizedConsistency:
        """Checks whether the filter estimates are statistically consistent
        by performing a Chi-squared test on the Normalized Estimation Error Squared (NEES)."""

    def nis_consistency(self, alpha: float | None = None) -> NormalizedConsistency:
        """Checks whether the filter innovations are statistically consistent
        by performing a Chi-squared test on the Normalized Innovation Squared (NIS)."""

    def rejected_residuals(self) -> list[PositionResidual]: ...
    def residual_ratio_within_threshold(self, threshold: float) -> float:
        """Computes the fraction of residual ratios that lie within ±threshold."""

    def rms_postfit_residuals(self) -> float:
        """Returns the root mean square of the postfit residuals"""

    def rms_prefit_residuals(self) -> float:
        """Returns the root mean square of the prefit residuals"""

    def rms_residual_ratios(self) -> float:
        """Returns the root mean square of the prefit residual ratios"""

    def smooth(self, almanac: Almanac) -> SpacecraftPositionODSolution:
        """Smoothes this OD solution."""

    def to_ephemeris(self, object_id: str) -> astro.Ephemeris:
        """Export to an ANISE ephemeris, which can be converted to a CCSDS OEM"""

    def to_parquet(self, path: str, cfg: ExportCfg) -> str:
        """Export OD solutions, gains, ratios, residuals, sigmas, etc. to parquet"""

    def to_traj(self) -> Trajectory: ...
    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class StochasticNoise:
    """Stochastic noise modeling used primarily for synthetic orbit determination measurements.

    This implementation distinguishes between the white noise model and the bias model. It also includes a constant offset."""

    bias: typing.Any
    white_noise: typing.Any

    def __init__(
        self,
        white_noise: WhiteNoise | None,
        bias: GaussMarkov | None,
        name: str | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Stochastic noise modeling used primarily for synthetic orbit determination measurements.

        This implementation distinguishes between the white noise model and the bias model. It also includes a constant offset."""

    def __new__(
        cls,
        white_noise: WhiteNoise | None = None,
        bias: GaussMarkov | None = None,
        name: str | None = None,
    ) -> StochasticNoise:
        """Stochastic noise modeling used primarily for synthetic orbit determination measurements.

        This implementation distinguishes between the white noise model and the bias model. It also includes a constant offset."""

    def covariance(self, epoch: time.Epoch) -> float:
        """Return the covariance of these stochastics at a given time."""

    @staticmethod
    def from_hardware_doppler_km_s(
        allan_deviation: float,
        integration_time: time.Duration,
        carrier: CarrierFreq,
        c_n0: CN0,
    ) -> StochasticNoise:
        """Constructs a hardware Doppler noise model."""

    @staticmethod
    def from_hardware_range_km(
        allan_deviation: float,
        integration_time: time.Duration,
        chip_rate: ChipRate,
        s_n0: SN0,
    ) -> StochasticNoise:
        """Constructs a high precision zero-mean range noise model (accounting for clock error and thermal error) from
        the Allan deviation of the clock, integration time, chip rate (depends on the ranging code), and
        signal-power-to-noise-density ratio (S/N₀).

        NOTE: The Allan Deviation should be provided given the integration time. For example, if the integration time
        is one second, the Allan Deviation should be the deviation over one second.

        IMPORTANT: These do NOT include atmospheric noises, which add up to ~10 cm one-sigma."""

    def simulate(
        self, path: str, runs: int | None, unit: str | None
    ) -> list[StochasticState]:
        """Executes a hardcoded 24-hour Monte Carlo simulation of the stochastic model, exporting the time history to a Parquet file.

        # Warning: Hardcoded Time Series & Diagnostic Data Gaps
        This method does *not* accept a user-defined tracking schedule or time series. It inherently evaluates the stochastic process
        over a strict 24-hour period, beginning at the exact system clock moment of method execution, utilizing a 1-minute step size.

        Furthermore, users will observe exactly 1,082 samples per simulation run, rather than the 1,441 samples expected from a
        continuous 24-hour 1-minute cadence. The simulation intentionally drops all epochs strictly greater than +6 hours and
        strictly less than +12 hours from the start time. This hardcoded artifact is designed to demonstrate variance bounds
        expansion in the absence of measurements (e.g., simulating a tracking dropout for a Gauss-Markov bias).

        # Algorithm
        1. Establish `start` as the system clock time at invocation.
        2. Construct an inclusive time series from `start` to `start + 24 hours` at 1-minute intervals.
        3. For each configured run, seed a PRNG (`Pcg64Mcg`) using system entropy.
        4. Evaluate the process covariance and sample the stochastic noise at each epoch.
        5. Discard all epochs inside the `(start + 6h, start + 12h)` open interval.
        6. Export the remaining 1,082 samples per run to an Apache Arrow RecordBatch and write to disk via Parquet.

        :param path: The filesystem path for the output Parquet file.
        :param runs: The number of Monte Carlo runs. Defaults to 25 if not provided.
        :param unit: An optional string appended to the Parquet column headers for plotting clarity.
        :raises Exception: If the underlying Apache Arrow RecordBatch fails to allocate or write to the specified filesystem path."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class StochasticState:
    dt_s: typing.Any
    run: typing.Any
    sample: typing.Any
    variance: typing.Any

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> StochasticState: ...
    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class Strand:
    """Stores a tracking strand with a start and end epoch"""

    end: typing.Any
    start: typing.Any

    def __init__(
        self,
        start: time.Epoch,
        end: time.Epoch,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Stores a tracking strand with a start and end epoch"""

    def __new__(cls, start: time.Epoch, end: time.Epoch) -> Strand:
        """Stores a tracking strand with a start and end epoch"""

    @staticmethod
    def from_asn1(data: bytes) -> Strand:
        """Decodes an ASN.1 DER encoded byte array into a Strand object."""

    def to_asn1(self) -> bytes:
        """Encodes this Strand object into an ASN.1 DER encoded byte array."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class TrackingDataArc:
    """Tracking data storing all of measurements as a B-Tree.
    It inherently does NOT support multiple concurrent measurements from several trackers.

    # Measurement Moduli, e.g. range modulus

    In the case of ranging, and possibly other data types, a code is used to measure the range to the spacecraft. The length of this code
    determines the ambiguity resolution, as per equation 9 in section 2.2.2.2 of the JPL DESCANSO, document 214, _Pseudo-Noise and Regenerative Ranging_.
    For example, using the JPL Range Code and a frequency range clock of 1 MHz, the range ambiguity is 75,660 km. In other words,
    as soon as the spacecraft is at a range of 75,660 + 1 km the JPL Range Code will report the vehicle to be at a range of 1 km.
    This is simply because the range code overlaps with itself, effectively loosing track of its own reference:
    it's due to the phase shift of the signal "lapping" the original signal length.

    ```text
    (Spacecraft)
    ^
    |    Actual Distance = 75,661 km
    |
    0 km                                         75,660 km (Wrap-Around)
    |-----------------------------------------------|
    When the "code length" is exceeded,
    measurements wrap back to 0.

    So effectively:
    Observed code range = Actual range (mod 75,660 km)
    75,661 km → 1 km

    ```

    Nyx can only resolve the range ambiguity if the tracking data specifies a modulus for this specific measurement type.
    For example, in the case of the JPL Range Code and a 1 MHz range clock, the ambiguity interval is 75,660 km.

    The measurement used in the Orbit Determination Process then becomes the following, where `//` represents the [Euclidian division](https://doc.rust-lang.org/std/primitive.f64.html#method.div_euclid).

    ```text
    k = computed_obs // ambiguity_interval
    real_obs = measured_obs + k * modulus
    ```

    Reference: JPL DESCANSO, document 214, _Pseudo-Noise and Regenerative Ranging_."""

    force_reject: typing.Any

    def __init__(
        self,
        measurements: list[Measurement],
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Tracking data storing all of measurements as a B-Tree.
        It inherently does NOT support multiple concurrent measurements from several trackers.

        # Measurement Moduli, e.g. range modulus

        In the case of ranging, and possibly other data types, a code is used to measure the range to the spacecraft. The length of this code
        determines the ambiguity resolution, as per equation 9 in section 2.2.2.2 of the JPL DESCANSO, document 214, _Pseudo-Noise and Regenerative Ranging_.
        For example, using the JPL Range Code and a frequency range clock of 1 MHz, the range ambiguity is 75,660 km. In other words,
        as soon as the spacecraft is at a range of 75,660 + 1 km the JPL Range Code will report the vehicle to be at a range of 1 km.
        This is simply because the range code overlaps with itself, effectively loosing track of its own reference:
        it's due to the phase shift of the signal "lapping" the original signal length.

        ```text
        (Spacecraft)
        ^
        |    Actual Distance = 75,661 km
        |
        0 km                                         75,660 km (Wrap-Around)
        |-----------------------------------------------|
        When the "code length" is exceeded,
        measurements wrap back to 0.

        So effectively:
        Observed code range = Actual range (mod 75,660 km)
        75,661 km → 1 km

        ```

        Nyx can only resolve the range ambiguity if the tracking data specifies a modulus for this specific measurement type.
        For example, in the case of the JPL Range Code and a 1 MHz range clock, the ambiguity interval is 75,660 km.

        The measurement used in the Orbit Determination Process then becomes the following, where `//` represents the [Euclidian division](https://doc.rust-lang.org/std/primitive.f64.html#method.div_euclid).

        ```text
        k = computed_obs // ambiguity_interval
        real_obs = measured_obs + k * modulus
        ```

        Reference: JPL DESCANSO, document 214, _Pseudo-Noise and Regenerative Ranging_."""

    def __new__(cls, measurements: list[Measurement]) -> TrackingDataArc:
        """Tracking data storing all of measurements as a B-Tree.
        It inherently does NOT support multiple concurrent measurements from several trackers.

        # Measurement Moduli, e.g. range modulus

        In the case of ranging, and possibly other data types, a code is used to measure the range to the spacecraft. The length of this code
        determines the ambiguity resolution, as per equation 9 in section 2.2.2.2 of the JPL DESCANSO, document 214, _Pseudo-Noise and Regenerative Ranging_.
        For example, using the JPL Range Code and a frequency range clock of 1 MHz, the range ambiguity is 75,660 km. In other words,
        as soon as the spacecraft is at a range of 75,660 + 1 km the JPL Range Code will report the vehicle to be at a range of 1 km.
        This is simply because the range code overlaps with itself, effectively loosing track of its own reference:
        it's due to the phase shift of the signal "lapping" the original signal length.

        ```text
        (Spacecraft)
        ^
        |    Actual Distance = 75,661 km
        |
        0 km                                         75,660 km (Wrap-Around)
        |-----------------------------------------------|
        When the "code length" is exceeded,
        measurements wrap back to 0.

        So effectively:
        Observed code range = Actual range (mod 75,660 km)
        75,661 km → 1 km

        ```

        Nyx can only resolve the range ambiguity if the tracking data specifies a modulus for this specific measurement type.
        For example, in the case of the JPL Range Code and a 1 MHz range clock, the ambiguity interval is 75,660 km.

        The measurement used in the Orbit Determination Process then becomes the following, where `//` represents the [Euclidian division](https://doc.rust-lang.org/std/primitive.f64.html#method.div_euclid).

        ```text
        k = computed_obs // ambiguity_interval
        real_obs = measured_obs + k * modulus
        ```

        Reference: JPL DESCANSO, document 214, _Pseudo-Noise and Regenerative Ranging_."""

    def apply_moduli(self) -> None:
        """Applies the moduli to each measurement, if defined."""

    def chunk(self, max_duration: time.Duration) -> list[TrackingDataArc]:
        """Splits a long tracking data arc into smaller chunks, each up to `max_duration` long."""

    def downsample(self, target_step: time.Duration) -> TrackingDataArc:
        """Downsamples the tracking data to a lower frequency using a simple moving average low-pass filter followed by decimation,
        returning new `TrackingDataArc` with downsampled measurements.

        It provides a computationally efficient approach to reduce the sampling rate while mitigating aliasing effects.

        # Algorithm

        1. A simple moving average filter is applied as a low-pass filter.
        2. Decimation is performed by selecting every Nth sample after filtering.

        # Advantages

        - Computationally efficient, suitable for large datasets common in spaceflight applications.
        - Provides basic anti-aliasing, crucial for preserving signal integrity in orbit determination and tracking.
        - Maintains phase information, important for accurate timing in spacecraft state estimation.

        # Limitations

        - The frequency response is not as sharp as more sophisticated filters (e.g., FIR, IIR).
        - May not provide optimal stopband attenuation for high-precision applications.

        ## Considerations for Spaceflight Applications

        - Suitable for initial data reduction in ground station tracking pipelines.
        - Adequate for many orbit determination and tracking tasks where computational speed is prioritized.
        - For high-precision applications (e.g., interplanetary navigation), consider using more advanced filtering techniques."""

    def duration(self) -> time.Duration | None:
        """Returns the duration this tracking arc"""

    def end_epoch(self) -> time.Epoch | None:
        """Returns the end epoch of this tracking arc"""

    def exclude_by_epoch(
        self, start: time.Epoch | None, end: time.Epoch | None
    ) -> TrackingDataArc:
        """Exclude measurements by epoch range."""

    def exclude_measurement_type(self, msr_type: MeasurementType) -> TrackingDataArc:
        """Exclude measurements by measurement type."""

    def exclude_tracker(self, tracker: str) -> TrackingDataArc:
        """Exclude measurements by tracker alias."""

    def filter_by_epoch(
        self, start: time.Epoch | None, end: time.Epoch | None
    ) -> TrackingDataArc:
        """Filter measurements by epoch range."""

    def filter_by_measurement_type(self, msr_type: MeasurementType) -> TrackingDataArc:
        """Filter measurements by measurement type."""

    def filter_by_offset(
        self, start: time.Duration | None, end: time.Duration | None
    ) -> TrackingDataArc:
        """Filter measurements by duration offset."""

    def filter_by_tracker(self, tracker: str) -> TrackingDataArc:
        """Filter measurements by tracker alias."""

    @staticmethod
    def from_ccsds_tdm(path: str, aliases: dict) -> nyx_space.od.TrackingDataArc:
        """Initializes a new Almanac from a file path to CCSDS OEM file, after converting to to SPICE SPK/BSP"""

    @staticmethod
    def from_parquet(path: str) -> TrackingDataArc:
        """Load TrackingDataArc from a parquet file."""

    def is_empty(self) -> bool:
        """Returns whether this arc has no measurements."""

    def len(self) -> int:
        """Returns the number of measurements in this data arc"""

    def min_duration_sep(self) -> time.Duration | None:
        """Returns the minimum duration between two subsequent measurements."""

    def resid_vs_ref_check(self) -> TrackingDataArc: ...
    def set_moduli(self, msr_type: MeasurementType, modulus: float) -> None:
        """Set (or overwrites) the modulus of the provided measurement type."""

    def sort(self) -> None:
        """Sort these measurements by epoch"""

    def start_epoch(self) -> time.Epoch | None:
        """Returns the start epoch of this tracking arc"""

    def to_parquet(self, path: str, cfg: ExportCfg) -> str:
        """Write tracking data arc to a parquet file."""

    def unique_aliases(self) -> list[str]: ...
    def unique_types(self) -> list[MeasurementType]: ...
    def write_ccsds_tdm(
        self, spacecraft_name: str, aliases: dict | None, path: str
    ) -> str:
        """Write tracking data in CCSDS TDM format."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class TrkConfig:
    """Stores a tracking configuration, there is one per tracking data simulator (e.g. one for ground station #1 and another for #2).
    By default, the tracking configuration is continuous and the tracking arc is from the beginning of the simulation to the end.
    In Python, any value that is set to None at initialization will use the default values: no scheduler, no strands, sampling at 1 min."""

    sampling: typing.Any
    scheduler: typing.Any
    strands: typing.Any

    def __init__(
        self,
        scheduler: Scheduler | None,
        sampling: time.Duration,
        strands: list[Strand] | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Stores a tracking configuration, there is one per tracking data simulator (e.g. one for ground station #1 and another for #2).
        By default, the tracking configuration is continuous and the tracking arc is from the beginning of the simulation to the end.
        In Python, any value that is set to None at initialization will use the default values: no scheduler, no strands, sampling at 1 min."""

    def __new__(
        cls,
        scheduler: Scheduler | None = None,
        sampling: typing.Optional[time.Duration] = ...,
        strands: list[Strand] | None = None,
    ) -> TrkConfig:
        """Stores a tracking configuration, there is one per tracking data simulator (e.g. one for ground station #1 and another for #2).
        By default, the tracking configuration is continuous and the tracking arc is from the beginning of the simulation to the end.
        In Python, any value that is set to None at initialization will use the default values: no scheduler, no strands, sampling at 1 min."""

    @staticmethod
    def from_asn1(data: bytes) -> TrkConfig:
        """Decodes an ASN.1 DER encoded byte array into a TrkConfig object."""

    def to_asn1(self) -> bytes:
        """Encodes this TrkConfig object into an ASN.1 DER encoded byte array."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class WhiteNoise:
    """White noise is an uncorrelated random variable."""

    mean: typing.Any
    sigma: typing.Any

    def __init__(
        self,
        mean: float,
        sigma: float,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        White noise is an uncorrelated random variable."""

    def __new__(cls, mean: float, sigma: float) -> WhiteNoise:
        """White noise is an uncorrelated random variable."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""
