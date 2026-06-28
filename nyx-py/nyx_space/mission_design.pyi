from __future__ import annotations
import typing

@typing.final
class AccelModels:
    """Acceleration models alter the orbital dynamics"""

    gravity_field: typing.Any
    point_masses: typing.Any

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class AtmDensity:
    """Density in kg/m^3 and altitudes in meters, not kilometers!"""

    @staticmethod
    def earth_exponential() -> typing.Any: ...
    Constant: type = ...
    Exponential: type = ...
    StdAtm: type = ...

@typing.final
class Drag:
    """`Drag` implements all three drag models."""

    density: typing.Any
    estimate: typing.Any
    frame: typing.Any

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class Dynamics:
    """Dynamics defines the dynamical environment with a set of acceleration and force models"""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class ExportCfg:
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

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class GravityFieldConfig:
    """Configuration holder for gravity field.

    Data is first loaded as a SHADR, if that fails, Nyx will try to load it as a COF file."""

    degree: typing.Any
    filepath: typing.Any
    frame: typing.Any
    gunzipped: typing.Any
    order: typing.Any

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class IntegratorMethod:
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

    def info(self) -> typing.Any:
        """Returns a string with the information about these options"""

    def set_max_step(self, max_step: typing.Any) -> typing.Any:
        """Set the maximum step size and sets the initial step to that value if currently greater"""

    def set_min_step(self, min_step: typing.Any) -> typing.Any:
        """Set the minimum step size and sets the initial step to that value if currently smaller"""

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class PointMasses:
    """PointMasses model"""

    celestial_objects: typing.Any
    correction: typing.Any

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class Propagator:
    dynamics: typing.Any
    method: typing.Any
    options: typing.Any

    def for_duration(
        self,
        spacecraft: typing.Any,
        duration: typing.Any,
        trajectory: typing.Any = True,
    ) -> typing.Any:
        """Propagates the initialization state for the desired duration, optionally not building the trajectory"""

    def many_for_duration(
        self,
        spacecraft: typing.Any,
        duration: typing.Any,
        trajectory: typing.Any = True,
    ) -> typing.Any:
        """Propagates the initialization state for the desired duration, optionally not building the trajectory"""

    def many_until_epoch(
        self, spacecraft: typing.Any, epoch: typing.Any, trajectory: typing.Any = True
    ) -> typing.Any:
        """Propagates the initialization state until the desired epoch, optionally not building the trajectory"""

    def many_until_event(
        self,
        spacecraft: typing.Any,
        event: typing.Any,
        max_duration: typing.Any,
        trigger: typing.Any = 1,
        event_frame: typing.Any = None,
        trajectory: typing.Any = True,
    ) -> typing.Any: ...
    def until_epoch(
        self, spacecraft: typing.Any, epoch: typing.Any, trajectory: typing.Any = True
    ) -> typing.Any:
        """Propagates the initialization state until the desired epoch, optionally not building the trajectory"""

    def until_event(
        self,
        spacecraft: typing.Any,
        event: typing.Any,
        max_duration: typing.Any,
        trigger: typing.Any = 1,
        event_frame: typing.Any = None,
        trajectory: typing.Any = True,
    ) -> typing.Any:
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

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class ShadowModel:
    light_source: typing.Any
    shadow_bodies: typing.Any

@typing.final
class SolarPressure:
    """Computation of solar radiation pressure is based on STK: <http://help.agi.com/stk/index.htm#gator/eq-solar.htm> ."""

    estimate: typing.Any
    phi: typing.Any
    shadow_model: typing.Any

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class SpacecraftSequence:
    thruster_sets: typing.Any

    @staticmethod
    def from_dhall(dhall_str: typing.Any) -> typing.Any: ...
    @staticmethod
    def from_yaml(yaml_str: typing.Any) -> typing.Any: ...
    def propagate(
        self, state: typing.Any, until_phase: typing.Any, almanac: typing.Any
    ) -> typing.Any: ...
    def setup(self, almanac: typing.Any) -> typing.Any: ...
    def thruster_set_insert(
        self, name: typing.Any, thruster: typing.Any
    ) -> typing.Any: ...
    def thruster_set_remove(self, name: typing.Any) -> typing.Any: ...
