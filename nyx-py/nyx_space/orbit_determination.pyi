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
    def continuous() -> typing.Any: ...
    @staticmethod
    def from_asn1(data: bytes) -> Cadence:
        """Decodes an ASN.1 DER encoded byte array into a Cadence object."""

    @staticmethod
    def intermittent(on: typing.Any, off: typing.Any) -> typing.Any: ...
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
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Configuration for exporting from Nyx to local disk."""

    def __new__(cls, timestamped: typing.Any = False) -> ExportCfg:
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
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
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

    def __new__(cls, tau: typing.Any, process_noise: typing.Any) -> GaussMarkov:
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

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        GroundStation defines a one-way or two-way ranging and doppler station. Set the integration time for two-way."""

    def __new__(
        cls,
        name: typing.Any,
        location: typing.Any,
        stochastic_noises: typing.Any,
        integration_time: typing.Any = None,
        light_time_correction: typing.Any = False,
        timestamp_noise_s: typing.Any = None,
    ) -> GroundStation:
        """GroundStation defines a one-way or two-way ranging and doppler station. Set the integration time for two-way."""

    def azimuth_elevation_of(
        self, rx: typing.Any, obstructing_body: typing.Any, almanac: typing.Any
    ) -> typing.Any:
        """Computes the azimuth and elevation of the provided object seen from this ground station, both in degrees.
        This is a shortcut to almanac.azimuth_elevation_range_sez."""

    @staticmethod
    def dump_many_yaml(stations: typing.Any, path: typing.Any) -> typing.Any: ...
    @staticmethod
    def dumps_many_yaml(stations: typing.Any) -> typing.Any: ...
    @staticmethod
    def from_asn1(data: bytes) -> GroundStation:
        """Decodes an ASN.1 DER encoded byte array into a GroundStation object."""

    @staticmethod
    def from_yaml(yaml_str: typing.Any) -> typing.Any: ...
    @staticmethod
    def load_many_yaml(path: typing.Any) -> typing.Any: ...
    @staticmethod
    def loads_many_yaml(yaml_str: typing.Any) -> typing.Any: ...
    def to_asn1(self) -> bytes:
        """Encodes this GroundStation object into an ASN.1 DER encoded byte array."""

    def to_orbit(self, epoch: typing.Any, almanac: typing.Any) -> typing.Any:
        """Return this ground station as an orbit in its current frame"""

    def to_yaml(self) -> typing.Any: ...
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
    configs: typing.Any
    devices: typing.Any

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(
        cls,
        devices: typing.Any,
        trajectory: typing.Any,
        configs: typing.Any,
        seed: typing.Any = None,
    ) -> GroundTrackingArcSim: ...
    def build_schedule(self, almanac: Almanac) -> typing.Any:
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
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        A type-agnostic simultaneous measurement storage structure. Allows storing any number of simultaneous measurement of a given taker.

        Note that two measurements are considered equal if the tracker and epoch match exactly, and if both have the same measurement types,
        and those measurements are equal to within 1e-10 (this allows for some leeway in TDM producers)."""

    def __new__(cls, tracker: typing.Any, epoch: typing.Any) -> Measurement:
        """A type-agnostic simultaneous measurement storage structure. Allows storing any number of simultaneous measurement of a given taker.

        Note that two measurements are considered equal if the tracker and epoch match exactly, and if both have the same measurement types,
        and those measurements are equal to within 1e-10 (this allows for some leeway in TDM producers)."""

    def correct(self, msr_type: typing.Any, correction: typing.Any) -> typing.Any:
        """Correct the provided measurement type with the provided correction, if that measurement type is available"""

    def observation(self, msr_type: typing.Any) -> typing.Any:
        """Returns the floating point value of this observation if this measurement contains the provided measurement type"""

    def push(self, msr_type: typing.Any, msr_value: typing.Any) -> typing.Any: ...
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
class ProcessNoise:
    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> ProcessNoise: ...
    @staticmethod
    def from_accel_m_s2(
        ax_m_s2: typing.Any,
        ay_m_s2: typing.Any,
        az_m_s2: typing.Any,
        disable_time: typing.Any,
        local_frame: typing.Any,
        x_decay_s: typing.Any,
        y_decay_s: typing.Any,
        z_decay_s: typing.Any,
    ) -> typing.Any: ...
    @staticmethod
    def from_velocity_m_s(
        vx_m_s: typing.Any,
        vy_m_s: typing.Any,
        vz_m_s: typing.Any,
        noise_duration: typing.Any,
        disable_time: typing.Any,
        local_frame: typing.Any,
    ) -> typing.Any: ...
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
    def computed_obs(self, msr_type: typing.Any) -> typing.Any:
        """Returns the computed/expected observation for this measurement type, if available"""

    def nis(self) -> typing.Any:
        """Returns the normalized innovation squared (NIS) as the norm squares of the whitened residual"""

    def real_obs(self, msr_type: typing.Any) -> typing.Any:
        """Returns the real observation for this measurement type, if available"""

    def whitened_residual(self, msr_type: typing.Any) -> typing.Any:
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
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        A scheduler allows building a scheduling of spaceraft tracking for a set of ground stations."""

    def __new__(
        cls,
        handoff: typing.Any = ...,
        cadence: typing.Any = None,
        min_samples: typing.Any = 10,
        sample_alignment: typing.Any = None,
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
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Reject measurements if the prefit is greater than the provided sigmas deviation from the measurement noise.

        # Important
        Some software, like ODTK, processes each measurement as a scalar. Nyx can process the measurements together.
        As such, if the prefit on range is bad, then the Doppler measurement with the same time stamp will also be rejected.
        This can lead to better convergence of the filter, and more appropriate results."""

    def __new__(cls, num_sigmas: typing.Any) -> SigmaRejection:
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
    def from_diag(nominal: typing.Any, diag: typing.Any) -> typing.Any:
        """Initializes a new filter estimate from the nominal state (not dispersed) and the diagonal of the covariance"""

    @staticmethod
    def from_dispersions(
        nominal_state: typing.Any, dispersions: typing.Any, seed: typing.Any = None
    ) -> typing.Any:
        """Generates an initial Kalman filter state estimate dispersed from the nominal state using the provided standard deviation parameters.

        The resulting estimate will have a diagonal covariance matrix constructed from the variances of each parameter."""

    def to_random_variable(self) -> typing.Any:
        """Builds a multivariate random variable spacecraft from this estimate's nominal state and covariance, zero mean."""

    def within_3sigma(self) -> typing.Any:
        """Returns whether this estimate is within three sigmas"""

    def within_sigma(self, sigma: typing.Any) -> typing.Any:
        """Returns whether this estimate is within some bound
        The 68-95-99.7 rule is a good way to assess whether the filter is operating normally"""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class SpacecraftODProcess:
    sigma_rejection: typing.Any
    variant: typing.Any

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(
        cls,
        prop: typing.Any,
        kf_variant: typing.Any,
        devices: typing.Any,
        sigma_reject: typing.Any = ...,
        process_noise: typing.Any = None,
    ) -> SpacecraftODProcess: ...
    def predict_for(
        self, initial_estimate: typing.Any, duration: typing.Any
    ) -> typing.Any:
        """Perform a time update. Continuously predicts the trajectory for the provided duration, with covariance mapping at each step."""

    def predict_until(
        self, initial_estimate: typing.Any, end_epoch: typing.Any
    ) -> typing.Any:
        """Perform a time update. Continuously predicts the trajectory until the provided end epoch, with covariance mapping at each step."""

    def process_arc(self, initial_estimate: typing.Any, arc: typing.Any) -> typing.Any:
        """Process the provided tracking arc for this orbit determination process."""

@typing.final
class SpacecraftODSolution:
    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> SpacecraftODSolution: ...
    def accepted_residuals(self) -> typing.Any: ...
    @staticmethod
    def from_parquet(path: typing.Any, devices: typing.Any) -> typing.Any: ...
    def is_filter_run(self) -> typing.Any: ...
    def is_nees_consistent(
        self, truth_traj: typing.Any, alpha: typing.Any = None
    ) -> typing.Any:
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

    def is_nis_consistent(self, alpha: typing.Any = None) -> typing.Any:
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

    def is_normal(self, alpha: typing.Any = None) -> typing.Any:
        """Checks whether the whitened residuals of the accepted residuals pass a normality test at a given significance level `alpha`, default to 0.05.

        This uses a simplified KS-test threshold: D_alpha = c(α) / √n.
        For example, for α = 0.05, c(α) is approximately 1.36.
        α=0.05 means a 5% probability of Type I error (incorrectly rejecting the null hypothesis when it is true).
        This threshold is standard in many statistical tests to balance sensitivity and false positives

        The computation of the c(alpha) is from https://en.wikipedia.org/w/index.php?title=Kolmogorov%E2%80%93Smirnov_test&oldid=1280260701#Two-sample_Kolmogorov%E2%80%93Smirnov_test

        Returns Ok(true) if the residuals are consistent with a normal distribution,
        Ok(false) if not, or None if no residuals are available."""

    def is_smoother_run(self) -> typing.Any: ...
    def ks_test_normality(self) -> typing.Any:
        """Computes the Kolmogorov–Smirnov statistic for the aggregated residual ratios of the accepted residuals.

        Returns Ok(ks_statistic) if residuals are available."""

    def rejected_residuals(self) -> typing.Any: ...
    def residual_ratio_within_threshold(self, threshold: typing.Any) -> typing.Any:
        """Computes the fraction of residual ratios that lie within ±threshold."""

    def rms_postfit_residuals(self) -> typing.Any:
        """Returns the root mean square of the postfit residuals"""

    def rms_prefit_residuals(self) -> typing.Any:
        """Returns the root mean square of the prefit residuals"""

    def rms_residual_ratios(self) -> typing.Any:
        """Returns the root mean square of the prefit residual ratios"""

    def smooth(self, almanac: typing.Any) -> typing.Any:
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

    def to_ephemeris(self, object_id: typing.Any) -> typing.Any:
        """Export to an ANISE ephemeris, which can be converted to a CCSDS OEM"""

    def to_parquet(self, path: typing.Any, cfg: typing.Any) -> typing.Any:
        """Export OD solutions, gains, ratios, residuals, sigmas, etc. to parquet"""

    def to_traj(self) -> typing.Any: ...
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
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Stochastic noise modeling used primarily for synthetic orbit determination measurements.

        This implementation distinguishes between the white noise model and the bias model. It also includes a constant offset."""

    def __new__(
        cls,
        white_noise: typing.Any = None,
        bias: typing.Any = None,
        name: typing.Any = None,
    ) -> StochasticNoise:
        """Stochastic noise modeling used primarily for synthetic orbit determination measurements.

        This implementation distinguishes between the white noise model and the bias model. It also includes a constant offset."""

    def covariance(self, epoch: typing.Any) -> typing.Any:
        """Return the covariance of these stochastics at a given time."""

    @staticmethod
    def from_hardware_doppler_km_s(
        allan_deviation: typing.Any,
        integration_time: typing.Any,
        carrier: typing.Any,
        c_n0: typing.Any,
    ) -> typing.Any: ...
    @staticmethod
    def from_hardware_range_km(
        allan_deviation: typing.Any,
        integration_time: typing.Any,
        chip_rate: typing.Any,
        s_n0: typing.Any,
    ) -> typing.Any:
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
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Stores a tracking strand with a start and end epoch"""

    def __new__(cls, start: typing.Any, end: typing.Any) -> Strand:
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
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
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

    def __new__(cls, measurements: typing.Any) -> TrackingDataArc:
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

    def apply_moduli(self) -> typing.Any:
        """Applies the moduli to each measurement, if defined."""

    def chunk(self, max_duration: typing.Any) -> TrackingDataArc:
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

    def duration(self) -> typing.Any:
        """Returns the duration this tracking arc"""

    def end_epoch(self) -> typing.Any:
        """Returns the end epoch of this tracking arc"""

    def exclude_by_epoch(self, start: typing.Any, end: typing.Any) -> nyx_space.orbit_determination.TrackingDataArc: ...
    def exclude_measurement_type(self, msr_type: typing.Any) -> nyx_space.orbit_determination.TrackingDataArc: ...
    def exclude_tracker(self, tracker: typing.Any) -> nyx_space.orbit_determination.TrackingDataArc: ...
    def filter_by_epoch(self, start: typing.Any, end: typing.Any) -> nyx_space.orbit_determination.TrackingDataArc: ...
    def filter_by_measurement_type(self, msr_type: typing.Any) -> nyx_space.orbit_determination.TrackingDataArc: ...
    def filter_by_offset(self, start: typing.Any, end: typing.Any) -> nyx_space.orbit_determination.TrackingDataArc: ...
    def filter_by_tracker(self, tracker: typing.Any) -> nyx_space.orbit_determination.TrackingDataArc: ...
    @staticmethod
    def from_ccsds_tdm(path: str, aliases: dict) -> nyx_space.orbit_determination.TrackingDataArc:
        """Initializes a new Almanac from a file path to CCSDS OEM file, after converting to to SPICE SPK/BSP"""

    def is_empty(self) -> typing.Any:
        """Returns whether this arc has no measurements."""

    def len(self) -> typing.Any:
        """Returns the number of measurements in this data arc"""

    def min_duration_sep(self) -> typing.Any:
        """Returns the minimum duration between two subsequent measurements."""

    def resid_vs_ref_check(self) -> typing.Any: ...
    def set_moduli(self, msr_type: typing.Any, modulus: typing.Any) -> typing.Any:
        """Set (or overwrites) the modulus of the provided measurement type."""

    def sort(self) -> typing.Any:
        """Sort these measurements by epoch"""

    def start_epoch(self) -> typing.Any:
        """Returns the start epoch of this tracking arc"""

    def unique_aliases(self) -> typing.Any: ...
    def unique_types(self) -> typing.Any: ...
    def write_ccsds_tdm(
        self, spacecraft_name: typing.Any, aliases: typing.Any, path: typing.Any
    ) -> typing.Any: ...
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
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Stores a tracking configuration, there is one per tracking data simulator (e.g. one for ground station #1 and another for #2).
        By default, the tracking configuration is continuous and the tracking arc is from the beginning of the simulation to the end.
        In Python, any value that is set to None at initialization will use the default values: no scheduler, no strands, sampling at 1 min."""

    def __new__(
        cls,
        scheduler: typing.Any = None,
        sampling: typing.Any = ...,
        strands: typing.Any = None,
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
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        White noise is an uncorrelated random variable."""

    def __new__(cls, mean: typing.Any, sigma: typing.Any) -> WhiteNoise:
        """White noise is an uncorrelated random variable."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

class PositionDevice:
    def __init__(
        self,
        name: str,
        stochastic_noises: dict[MeasurementType, StochasticNoise],
    ) -> None: ...
    @classmethod
    def from_yaml(cls, yaml_str: str) -> PositionDevice: ...
    def to_yaml(self) -> str: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...


class PositionTrackingArcSim:
    def __init__(
        self,
        devices: dict[str, PositionDevice],
        trajectory: PyTrajectory,
        configs: dict[str, TrkConfig],
        seed: int | None = None,
    ) -> None: ...
    def generate_measurements(self, almanac: Almanac) -> TrackingDataArc: ...
    @property
    def configs(self) -> dict[str, TrkConfig]: ...
    @property
    def devices(self) -> dict[str, PositionDevice]: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...


class SpacecraftPositionODProcess:
    def __init__(
        self,
        prop: Propagator,
        kf_variant: KalmanVariant,
        devices: dict[str, PositionDevice],
        sigma_reject: SigmaRejection | None = None,
        process_noise: ProcessNoise | None = None,
    ) -> None: ...
    def process_arc(
        self,
        initial_estimate: SpacecraftEstimate,
        arc: TrackingDataArc,
    ) -> SpacecraftPositionODSolution: ...
    def predict_until(
        self,
        initial_estimate: SpacecraftEstimate,
        end_epoch: Epoch,
    ) -> SpacecraftPositionODSolution: ...
    def predict_for(
        self,
        initial_estimate: SpacecraftEstimate,
        duration: Duration,
    ) -> SpacecraftPositionODSolution: ...


class SpacecraftPositionODSolution:
    def is_filter_run(self) -> bool: ...
    def is_smoother_run(self) -> bool: ...
    def to_traj(self) -> PyTrajectory: ...
    def to_parquet(self, path: str, cfg: ExportCfg) -> str: ...
    def smooth(self, almanac: Almanac) -> SpacecraftPositionODSolution: ...
    @classmethod
    def from_parquet(
        cls,
        path: str,
        devices: dict[str, PositionDevice],
    ) -> SpacecraftPositionODSolution: ...
    def __str__(self) -> str: ...
    def __repr__(self) -> str: ...
