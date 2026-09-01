from __future__ import annotations
from anise import Aberration
from anise import Almanac
from anise import analysis
from anise import astro
from anise import time
import typing

@typing.final
class AccelModels:
    """Acceleration models alter the orbital dynamics"""

    gravity_field: typing.Any
    point_masses: typing.Any
    solid_tides: typing.Any

    def __init__(
        self,
        point_masses: PointMasses | None,
        gravity_field: GravityFieldConfig | None,
        solid_tides: SolidTides | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Acceleration models alter the orbital dynamics"""

    def __new__(
        cls,
        point_masses: PointMasses | None = None,
        gravity_field: GravityFieldConfig | None = None,
        solid_tides: SolidTides | None = None,
    ) -> AccelModels:
        """Acceleration models alter the orbital dynamics"""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class AtmDensity:
    """Density in kg/m^3 and altitudes in kilometers"""

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Density in kg/m^3 and altitudes in kilometers"""

    def __new__(cls) -> AtmDensity:
        """Density in kg/m^3 and altitudes in kilometers"""

    @staticmethod
    def earth_exponential() -> AtmDensity:
        """Constructs a standard exponential drag model for Earth orbiters.

        Configured with nominal LEO reference parameters at $h_0 = 700\\text{ km}$:
        * $\\rho_0 = 3.614 \\times 10^{-13}\\text{ kg/m}^3$
        * $H = 88.667\\text{ km}$ ($88,667\\text{ m}$)"""
    Constant: type = ...
    Exponential: type = ...
    HarrisPriester: type = ...
    NRLMSISE00: type = ...
    StdAtm: type = ...

@typing.final
class Drag:
    """`Drag` implements all three drag models."""

    density: typing.Any
    estimate: typing.Any
    frame: typing.Any

    def __init__(
        self,
        density: AtmDensity,
        frame: astro.Frame,
        estimate: bool,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        `Drag` implements all three drag models."""

    def __new__(
        cls,
        density: AtmDensity,
        frame: astro.Frame,
        estimate: typing.Optional[bool] = True,
    ) -> Drag:
        """`Drag` implements all three drag models."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class Dynamics:
    """Dynamics defines the dynamical environment with a set of acceleration and force models"""

    def __init__(
        self,
        accel_models: AccelModels,
        force_models: ForceModels,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Dynamics defines the dynamical environment with a set of acceleration and force models"""

    def __new__(
        cls,
        accel_models: typing.Optional[AccelModels] = ...,
        force_models: typing.Optional[ForceModels] = ...,
    ) -> Dynamics:
        """Dynamics defines the dynamical environment with a set of acceleration and force models"""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

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
class ForceModels:
    """Force models alter the spacecraft dynamics (they need a mass)."""

    drag: typing.Any
    solar_pressure: typing.Any

    def __init__(
        self,
        solar_pressure: SolarPressure | None,
        drag: Drag | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Force models alter the spacecraft dynamics (they need a mass)."""

    def __new__(
        cls, solar_pressure: SolarPressure | None = None, drag: Drag | None = None
    ) -> ForceModels:
        """Force models alter the spacecraft dynamics (they need a mass)."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class GeomagneticMode:
    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> GeomagneticMode: ...
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
    ExtendedHistory57h: GeomagneticMode = ...
    Off: GeomagneticMode = ...
    StandardDailyAp: GeomagneticMode = ...

@typing.final
class GravityFieldConfig:
    """Configuration holder for gravity field.

    Data is first loaded as a SHADR, if that fails, Nyx will try to load it as a COF file."""

    degree: typing.Any
    filepath: typing.Any
    frame: typing.Any
    order: typing.Any

    def __init__(
        self,
        degree: int,
        order: int,
        filepath: str,
        frame: astro.FrameUid,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Configuration holder for gravity field.

        Data is first loaded as a SHADR, if that fails, Nyx will try to load it as a COF file."""

    def __new__(
        cls, degree: int, order: int, filepath: str, frame: astro.FrameUid
    ) -> GravityFieldConfig:
        """Configuration holder for gravity field.

        Data is first loaded as a SHADR, if that fails, Nyx will try to load it as a COF file."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class IntegratorMethod:
    """Enum of supported integration methods, all of which are part of the Runge Kutta family of ordinary differential equation (ODE) solvers.
    Nomenclature: X-Y means that this is an X order solver with a Y order error correction step."""

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Enum of supported integration methods, all of which are part of the Runge Kutta family of ordinary differential equation (ODE) solvers.
        Nomenclature: X-Y means that this is an X order solver with a Y order error correction step."""

    def __new__(cls) -> IntegratorMethod:
        """Enum of supported integration methods, all of which are part of the Runge Kutta family of ordinary differential equation (ODE) solvers.
        Nomenclature: X-Y means that this is an X order solver with a Y order error correction step."""

    def __int__(self) -> None:
        """int(self)"""

    def __repr__(self) -> str:
        """Return repr(self)."""
    CashKarp45: IntegratorMethod = ...
    DormandPrince45: IntegratorMethod = ...
    DormandPrince78: IntegratorMethod = ...
    RungeKutta4: IntegratorMethod = ...
    RungeKutta89: IntegratorMethod = ...
    Verner56: IntegratorMethod = ...

@typing.final
class IntegratorOptions:
    """Stores the integrator options, including the minimum and maximum step sizes, and the central body to perform the integration.

    Note that different step sizes and max errors are only used for adaptive
    methods. To use a fixed step integrator, initialize the options using `with_fixed_step`, and
    use whichever adaptive step integrator is desired.  For example, initializing an RK45 with
    fixed step options will lead to an RK4 being used instead of an RK45."""

    max_step: typing.Any
    min_step: typing.Any
    tolerance: typing.Any

    def __init__(
        self,
        min_step: time.Duration | None,
        max_step: time.Duration | None,
        tolerance: float | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Stores the integrator options, including the minimum and maximum step sizes, and the central body to perform the integration.

        Note that different step sizes and max errors are only used for adaptive
        methods. To use a fixed step integrator, initialize the options using `with_fixed_step`, and
        use whichever adaptive step integrator is desired.  For example, initializing an RK45 with
        fixed step options will lead to an RK4 being used instead of an RK45."""

    def __new__(
        cls,
        min_step: time.Duration | None = None,
        max_step: time.Duration | None = None,
        tolerance: float | None = None,
    ) -> IntegratorOptions:
        """Stores the integrator options, including the minimum and maximum step sizes, and the central body to perform the integration.

        Note that different step sizes and max errors are only used for adaptive
        methods. To use a fixed step integrator, initialize the options using `with_fixed_step`, and
        use whichever adaptive step integrator is desired.  For example, initializing an RK45 with
        fixed step options will lead to an RK4 being used instead of an RK45."""

    def info(self) -> str:
        """Returns a string with the information about these options"""

    def set_max_step(self, max_step: time.Duration) -> None:
        """Set the maximum step size and sets the initial step to that value if currently greater"""

    def set_min_step(self, min_step: time.Duration) -> None:
        """Set the minimum step size and sets the initial step to that value if currently smaller"""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class Msise00DailyWeather:
    """Target weather payload required by the NRLMSISE-00 density model."""

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Target weather payload required by the NRLMSISE-00 density model."""

    def __new__(cls) -> Msise00DailyWeather:
        """Target weather payload required by the NRLMSISE-00 density model."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class Nrlmsise00Flags:
    """Defines all of the available flags in the NRLMSISE00 model.
    NOTE By default, Nyx will use the mean local solar time computation. Set mean_lst=false
    to use the apparent local solar time."""

    annual_harmonics: typing.Any
    boundary_density_variations: typing.Any
    departures_from_diffusive_equilibrium: typing.Any
    diurnal_tides: typing.Any
    exospheric_temp_variations: typing.Any
    f107_solar_flux: typing.Any
    geomagnetic: typing.Any
    gradient_variations: typing.Any
    lower_boundary_temp_variations: typing.Any
    lower_mesosphere_temp_variations: typing.Any
    lower_thermosphere_temp_variations: typing.Any
    mean_lst: typing.Any
    semiannual_harmonics: typing.Any
    semidiurnal_tides: typing.Any
    terdiurnal_tides: typing.Any
    time_independent: typing.Any
    turbopause_scale_height_variations: typing.Any
    upper_stratosphere_temp_variations: typing.Any
    ut_and_longitude: typing.Any

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Defines all of the available flags in the NRLMSISE00 model.
        NOTE By default, Nyx will use the mean local solar time computation. Set mean_lst=false
        to use the apparent local solar time."""

    def __new__(
        cls,
        geomagnetic: typing.Any = None,
        f107_solar_flux: typing.Any = True,
        time_independent: typing.Any = True,
        annual_harmonics: typing.Any = True,
        semiannual_harmonics: typing.Any = True,
        diurnal_tides: typing.Any = True,
        semidiurnal_tides: typing.Any = True,
        terdiurnal_tides: typing.Any = True,
        ut_and_longitude: typing.Any = True,
        exospheric_temp_variations: typing.Any = True,
        lower_boundary_temp_variations: typing.Any = True,
        gradient_variations: typing.Any = True,
        departures_from_diffusive_equilibrium: typing.Any = True,
        lower_thermosphere_temp_variations: typing.Any = True,
        upper_stratosphere_temp_variations: typing.Any = True,
        boundary_density_variations: typing.Any = True,
        lower_mesosphere_temp_variations: typing.Any = True,
        turbopause_scale_height_variations: typing.Any = True,
        mean_lst: typing.Any = True,
    ) -> Nrlmsise00Flags:
        """Defines all of the available flags in the NRLMSISE00 model.
        NOTE By default, Nyx will use the mean local solar time computation. Set mean_lst=false
        to use the apparent local solar time."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class PointMasses:
    """PointMasses model"""

    celestial_objects: typing.Any
    correction: typing.Any

    def __init__(
        self,
        celestial_objects: list[int],
        correction: Aberration | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        PointMasses model"""

    def __new__(
        cls, celestial_objects: list[int], correction: Aberration | None = None
    ) -> PointMasses:
        """PointMasses model"""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class PropagationResult:
    state: typing.Any
    trajectory: typing.Any

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> PropagationResult: ...

@typing.final
class Propagator:
    """Numerical propagator for a spacecraft state."""

    dynamics: typing.Any
    method: typing.Any
    options: typing.Any

    def __init__(
        self,
        dynamics: Dynamics,
        almanac: Almanac,
        method: IntegratorMethod,
        options: IntegratorOptions,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Numerical propagator for a spacecraft state."""

    def __new__(
        cls,
        dynamics: Dynamics,
        almanac: Almanac,
        method: typing.Optional[IntegratorMethod] = ...,
        options: typing.Optional[IntegratorOptions] = ...,
    ) -> Propagator:
        """Numerical propagator for a spacecraft state."""

    def accel_km_s2(self, spacecraft: Spacecraft) -> list[float]:
        """Compute the instantaneous equations of motion for this spacecraft"""

    def for_duration(
        self,
        spacecraft: Spacecraft,
        duration: time.Duration,
        trajectory: typing.Optional[bool] = True,
    ) -> PropagationResult:
        """Propagates the initialization state for the desired duration, optionally not building the trajectory"""

    def many_for_duration(
        self,
        spacecraft: list[Spacecraft],
        duration: time.Duration,
        trajectory: typing.Optional[bool] = True,
    ) -> list[PropagationResult]:
        """Propagates the initialization state for the desired duration, optionally not building the trajectory"""

    def many_until_epoch(
        self,
        spacecraft: list[Spacecraft],
        epoch: time.Epoch,
        trajectory: typing.Optional[bool] = True,
    ) -> list[PropagationResult]:
        """Propagates the initialization state until the desired epoch, optionally not building the trajectory"""

    def many_until_event(
        self,
        spacecraft: list[Spacecraft],
        event: analysis.Event,
        max_duration: time.Duration,
        trigger: typing.Optional[int] = 1,
        event_frame: astro.Frame | None = None,
        trajectory: typing.Optional[bool] = True,
    ) -> list[PropagationResult]:
        """Propagates many states until event."""

    def until_epoch(
        self,
        spacecraft: Spacecraft,
        epoch: time.Epoch,
        trajectory: typing.Optional[bool] = True,
    ) -> PropagationResult:
        """Propagates the initialization state until the desired epoch, optionally not building the trajectory"""

    def until_event(
        self,
        spacecraft: Spacecraft,
        event: analysis.Event,
        max_duration: time.Duration,
        trigger: typing.Optional[int] = 1,
        event_frame: astro.Frame | None = None,
        trajectory: typing.Optional[bool] = True,
    ) -> PropagationResult:
        """Propagates the initialization state until the specified event has occurred `trigger` times, or until `max_duration` is reached.

        This method monitors the provided `event` during propagation. Once the event condition is met
        `trigger` number of times (e.g., set `trigger` to 1 for the first occurrence), the propagation stops
        at the end of that integration step.

        A root-finding algorithm (Brent's method) is then used to locate the exact time of the event
        within the final integration step. The returned state corresponds to this precise event time,
        interpolated from the trajectory.

        # Arguments

        * `max_duration` - The maximum duration to propagate if the event is not triggered the requested number of times.
        * `event` - The event definition (scalar expression and condition) to monitor.
        * `trigger` - The 1-based index of the event occurrence to stop at (e.g. 1 for the first crossing, 2 for the second).

        # Returns

        A tuple containing:
        1. The interpolated state exactly at the moment the $n$-th event occurred.
        2. The full trajectory recorded up to the end of the propagation step where the event occurred, unless explicitly ignored (but it is still built)

        # Errors

        * `PropagationError::NthEventError`: Returned if `max_duration` is reached before the event was triggered `trigger` times.
        * `PropagationError::TrajectoryEvent`: Returned if the interpolation of the event state fails.
        * `PropagationError::Analysis`: Returned if the event evaluation fails during the search."""

@typing.final
class PropagatorConfig:
    """Propagator config includes the method, options, and all dynamics"""

    def __init__(
        self,
        dynamics: Dynamics,
        method: IntegratorMethod,
        options: IntegratorOptions,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Propagator config includes the method, options, and all dynamics"""

    def __new__(
        cls, dynamics: Dynamics, method: IntegratorMethod, options: IntegratorOptions
    ) -> PropagatorConfig:
        """Propagator config includes the method, options, and all dynamics"""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class ShadowModel:
    correction: typing.Any
    light_source: typing.Any
    shadow_bodies: typing.Any

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> ShadowModel: ...

@typing.final
class SolarPressure:
    """Computation of solar radiation pressure is based on STK: <http://help.agi.com/stk/index.htm#gator/eq-solar.htm> ."""

    estimate: typing.Any
    phi: typing.Any
    shadow_model: typing.Any

    def __init__(
        self,
        shadow_bodies: list[astro.Frame],
        almanac: Almanac,
        flux_w_m2: float,
        estimate: bool,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Computation of solar radiation pressure is based on STK: <http://help.agi.com/stk/index.htm#gator/eq-solar.htm> ."""

    def __new__(
        cls,
        shadow_bodies: list[astro.Frame],
        almanac: Almanac,
        correction: typing.Any = None,
        flux_w_m2: typing.Optional[float] = ...,
        estimate: typing.Optional[bool] = True,
    ) -> SolarPressure:
        """Computation of solar radiation pressure is based on STK: <http://help.agi.com/stk/index.htm#gator/eq-solar.htm> ."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class SolidTides:
    """`SolidTides` implements the solid tide acceleration model.
    It accounts for the crust deformation due to the configured tidal perturbers.
    Formulas are based on IERS 2010 Conventions."""

    correction: typing.Any
    frame: typing.Any
    k2: typing.Any
    k3: typing.Any
    perturbers: typing.Any

    def __init__(
        self,
        frame: astro.Frame,
        k2: float,
        k3: float,
        perturbers: list[TidalPerturber],
        correction: Aberration,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        `SolidTides` implements the solid tide acceleration model.
        It accounts for the crust deformation due to the configured tidal perturbers.
        Formulas are based on IERS 2010 Conventions."""

    def __new__(
        cls,
        frame: astro.Frame,
        k2: float,
        k3: float,
        perturbers: list[TidalPerturber],
        correction: typing.Optional[Aberration] = None,
    ) -> SolidTides:
        """`SolidTides` implements the solid tide acceleration model.
        It accounts for the crust deformation due to the configured tidal perturbers.
        Formulas are based on IERS 2010 Conventions."""

    @staticmethod
    def earth_moon_system(
        earth_frame: astro.Frame,
        moon_frame: astro.Frame,
        almanac: Almanac,
        correction: typing.Optional[Aberration] = None,
    ) -> SolidTides:
        """Initializes solid tides with the Moon and the Sun, where the k3 is only computed for the Moon.
        Sets the k2 Love number to 0.3019 and the k3 Love number to 0.093"""

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
class SpaceWeatherData:
    """Stores SpaceWeather data as provided by [CelesTrak](https://celestrak.org/SpaceData/).
    Data may be provided either as original CSV or in a compressed (non-archived) gunzip (gz) format."""

    interpolate: typing.Any

    def __init__(
        self,
        path: str | None,
        fallback: StaticSpaceWeather | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Stores SpaceWeather data as provided by [CelesTrak](https://celestrak.org/SpaceData/).
        Data may be provided either as original CSV or in a compressed (non-archived) gunzip (gz) format."""

    def __new__(
        cls,
        path: str | None,
        fallback: StaticSpaceWeather | None,
        interpolate: typing.Any = False,
    ) -> SpaceWeatherData:
        """Stores SpaceWeather data as provided by [CelesTrak](https://celestrak.org/SpaceData/).
        Data may be provided either as original CSV or in a compressed (non-archived) gunzip (gz) format."""

    def build_ap_history(self, midnight: time.Epoch, bin_idx: int) -> list[float]:
        """Assembles the 7-element Ap array spanning current bin back 57 hours across 4 calendar days.

        Missing daily records or unforecasted bins are populated using the configured `SpaceWeatherFallback`."""

    def msise_weather(self, epoch: time.Epoch) -> Msise00DailyWeather:
        """Evaluates the space weather state at `epoch` and constructs the `Msise00DailyWeather` payload.

        Missing daily records or unforecasted fields are resolved using `SpaceWeatherFallback`."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class Spacecraft:
    """A spacecraft state, composed of its orbit, its masses (dry, prop, extra, all in kg), its SRP configuration, its drag configuration, its thruster configuration, and its guidance mode.

    Optionally, the spacecraft state can also store the state transition matrix from the start of the propagation until the current time (i.e. trajectory STM, not step-size STM)."""

    drag: typing.Any
    mass: typing.Any
    orbit: typing.Any
    srp: typing.Any

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        A spacecraft state, composed of its orbit, its masses (dry, prop, extra, all in kg), its SRP configuration, its drag configuration, its thruster configuration, and its guidance mode.

        Optionally, the spacecraft state can also store the state transition matrix from the start of the propagation until the current time (i.e. trajectory STM, not step-size STM)."""

    def __new__(
        cls,
        orbit: typing.Any,
        mass: typing.Any = None,
        srp: typing.Any = None,
        drag: typing.Any = None,
        thruster: typing.Any = None,
        mode: typing.Any = None,
    ) -> Spacecraft:
        """A spacecraft state, composed of its orbit, its masses (dry, prop, extra, all in kg), its SRP configuration, its drag configuration, its thruster configuration, and its guidance mode.

        Optionally, the spacecraft state can also store the state transition matrix from the start of the propagation until the current time (i.e. trajectory STM, not step-size STM)."""

    @staticmethod
    def from_asn1(data: bytes) -> astro.Mass:
        """Decodes an ASN.1 DER encoded byte array into a Mass object."""

    def rss(self, other: typing.Any) -> typing.Any:
        """Returns the root sum square error between this spacecraft and the other, in kilometers for the position, kilometers per second in velocity, and kilograms in prop"""

    def to_asn1(self) -> bytes:
        """Encodes this Mass object into an ASN.1 DER encoded byte array."""

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
class SpacecraftSequence:
    thruster_sets: typing.Any

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature."""

    def __new__(cls) -> SpacecraftSequence: ...
    @staticmethod
    def from_dhall(dhall_str: str) -> SpacecraftSequence:
        """Load SpacecraftSequence from Dhall."""

    @staticmethod
    def from_yaml(yaml_str: str) -> SpacecraftSequence:
        """Load SpacecraftSequence from YAML."""

    def propagate(
        self, state: Spacecraft, until_phase: str | None, almanac: Almanac
    ) -> list[str, str]:
        """Propagate the state through the sequence until a given phase."""

    def setup(self, almanac: Almanac) -> None:
        """Setup the sequence with the provided Almanac."""

    def thruster_set_insert(self, name: str, thruster: Thruster) -> None:
        """Insert a thruster with the given name into the thruster set."""

    def thruster_set_remove(self, name: str) -> None:
        """Remove a thruster with the given name from the thruster set."""

@typing.final
class StaticSpaceWeather:
    """Strategy for resolving missing predictive space weather parameters (F10.7, Ap, Kp)."""

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Strategy for resolving missing predictive space weather parameters (F10.7, Ap, Kp)."""

    def __new__(cls) -> StaticSpaceWeather:
        """Strategy for resolving missing predictive space weather parameters (F10.7, Ap, Kp)."""
    Custom: type = ...
    SolarAverage: type = ...
    SolarMaximum: type = ...
    SolarMinimum: type = ...

@typing.final
class Thruster:
    """Defines a thruster with a maximum isp and a maximum thrust."""

    isp_s: typing.Any
    thrust_N: typing.Any

    def __init__(
        self, *args: typing.Optional[typing.Any], **kwargs: typing.Optional[typing.Any]
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Defines a thruster with a maximum isp and a maximum thrust."""

    def __new__(cls, thrust_N: typing.Any, isp_s: typing.Any) -> Thruster:
        """Defines a thruster with a maximum isp and a maximum thrust."""

    def exhaust_velocity_m_s(self) -> typing.Any:
        """Returns the exhaust velocity v_e in meters per second"""

@typing.final
class TidalPerturber:
    """A celestial body raising tides on the central body."""

    compute_degree_3: typing.Any
    frame: typing.Any

    def __init__(
        self,
        frame: astro.Frame,
        compute_degree_3: bool,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        A celestial body raising tides on the central body."""

    def __new__(cls, frame: astro.Frame, compute_degree_3: bool) -> TidalPerturber:
        """A celestial body raising tides on the central body."""

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
class Trajectory:
    """Spacecraft Trajectory."""

    def __init__(
        self,
        path: str,
        template: Spacecraft | None,
        *args: typing.Optional[typing.Any],
        **kwargs: typing.Optional[typing.Any],
    ) -> None:
        """Initialize self.  See help(type(self)) for accurate signature.
        Spacecraft Trajectory."""

    def __new__(cls, path: str, template: Spacecraft | None) -> Trajectory:
        """Spacecraft Trajectory."""

    def append(self, many_spacecraft: list[Spacecraft]) -> None:
        """Append many spacecraft to this trajectory."""

    def at(self, epoch: time.Epoch) -> Spacecraft:
        """Evaluate the trajectory at this specific epoch."""

    def push(self, spacecraft: Spacecraft) -> None:
        """Add another state to this trajectory."""

    def ric_diff_to_parquet(self, other: Trajectory, path: str, cfg: ExportCfg) -> str:
        """Export the difference in RIC from of this trajectory compare to the "other" trajectory in parquet format.

        # Notes
        + The RIC frame accounts for the transport theorem by performing a finite differencing of the RIC frame."""

    def to_ephemeris(self, object_id: str, cfg: ExportCfg) -> astro.Ephemeris:
        """Export this spacecraft trajectory estimate to an ANISE Ephemeris"""

    def to_parquet(self, path: str, cfg: ExportCfg) -> str:
        """Write trajectory to a parquet file."""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""
