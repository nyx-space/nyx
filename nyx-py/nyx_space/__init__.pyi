from __future__ import annotations
from anise import astro
import typing

@typing.final
class DragData:
    area_m2: float
    coeff_drag: float

    @staticmethod
    def from_asn1(data: bytes) -> astro.DragData:
        """Decodes an ASN.1 DER encoded byte array into a DragData object."""

    def to_asn1(self) -> bytes:
        """Encodes this DragData object into an ASN.1 DER encoded byte array."""

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
class GuidanceMode:

    def __int__(self) -> None:
        """int(self)"""

    def __repr__(self) -> str:
        """Return repr(self)."""
    Coast: GuidanceMode = ...
    Inhibit: GuidanceMode = ...
    Thrust: GuidanceMode = ...

@typing.final
class Mass:
    """Defines a spacecraft mass a the sum of the dry (structural) mass and the propellant mass, both in kilogram"""
    dry_mass_kg: float
    extra_mass_kg: float
    prop_mass_kg: float

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

    def __repr__(self) -> str:
        """Return repr(self)."""

    def __str__(self) -> str:
        """Return str(self)."""

@typing.final
class SRPData:
    area_m2: float
    coeff_reflectivity: float

    @staticmethod
    def from_asn1(data: bytes) -> astro.SRPData:
        """Decodes an ASN.1 DER encoded byte array into an SRPData object."""

    def to_asn1(self) -> bytes:
        """Encodes this SRPData object into an ASN.1 DER encoded byte array."""

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
class Thruster:
    """Defines a thruster with a maximum isp and a maximum thrust."""
    isp_s: typing.Any
    thrust_N: typing.Any

    def exhaust_velocity_m_s(self) -> typing.Any:
        """Returns the exhaust velocity v_e in meters per second"""