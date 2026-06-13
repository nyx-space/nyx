from nyx_space import DragData, Mass, Spacecraft, SRPData
from nyx_space.anise import MetaAlmanac
from nyx_space.anise.analysis import OrbitalElement
from nyx_space.anise.astro import Orbit
from nyx_space.anise.constants import Frames
from nyx_space.anise.time import Epoch
from nyx_space.monte_carlo import MvnSpacecraft, StateDispersion, StateParameter


def test_def_sc():
    almanac = MetaAlmanac.latest()
    eme2k = almanac.frame_info(Frames.EME2000)
    # Define an orbit
    orbit = Orbit.from_keplerian(
        6800.0, 1e-4, 45.0, 60.0, 75.0, 90.0, Epoch("2020-02-29 01:02:03 TDB"), eme2k
    )

    sc = Spacecraft(orbit)
    sc_mass = Spacecraft(orbit, Mass(123.0))
    sc_mass_srp_drag = Spacecraft(
        orbit, Mass(123.0), SRPData(10.0, 1.2), DragData(10.0, 2.0)
    )

    # Test that we can serde to ASN1
    sc_mass2 = Spacecraft.from_asn1(sc_mass.to_asn1())
    assert sc_mass2.mass.dry_mass_kg == 123.0
    sc_mass_srp_drag2 = Spacecraft.from_asn1(sc_mass_srp_drag.to_asn1())
    assert sc_mass_srp_drag2.drag.coeff_drag == 2.0
    assert sc_mass_srp_drag2.srp.area_m2 == 10.0

    # Multivariate Normal Spacecraft test
    disp = [
        StateDispersion.zero_mean(
            StateParameter.Element(OrbitalElement.SemiMajorAxis), 15.0
        ),
        StateDispersion.zero_mean(StateParameter.Element(OrbitalElement.RAAN), 5.0),
    ]

    mvn = MvnSpacecraft(sc, disp)

    dispersed = mvn.sample(1000, 123)
    assert len(dispersed) == 1000


if __name__ == "__main__":
    test_def_sc()
