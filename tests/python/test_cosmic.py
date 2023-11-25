from nyx_space.cosmic import Cosm, Orbit, Spacecraft, SrpConfig, DragConfig
from nyx_space.time import Epoch, Unit, Duration
from nyx_space.monte_carlo import generate_orbits, generate_spacecraft, StateParameter
import pickle


def test_load_cosm():
    cosm = Cosm.de438()
    assert cosm.frames_get_names() == [
        "SSB J2000",
        "Sun J2000",
        "iau sun",
        "Mercury Barycenter J2000",
        "Venus Barycenter J2000",
        "iau venus",
        "Earth Barycenter J2000",
        "Earth J2000",
        "iau earth",
        "Moon J2000",
        "iau moon",
        "Mars Barycenter J2000",
        "iau mars",
        "Jupiter Barycenter J2000",
        "iau jupiter",
        "Saturn Barycenter J2000",
        "iau saturn",
        "Uranus Barycenter J2000",
        "iau uranus",
        "Neptune Barycenter J2000",
        "iau neptune",
        "Pluto Barycenter J2000",
    ]


def test_define_spacecraft():
    cosm = Cosm.de438()
    eme2k = cosm.frame("EME2000")

    assert eme2k.gm() == 398600.435392  # SPICE data, not GMAT data (398600.4415)
    assert eme2k.is_geoid()
    assert eme2k.equatorial_radius() == 6378.1363

    e = Epoch.system_now()

    orbit = Orbit.from_keplerian_altitude(
        400,
        ecc=1e-4,
        inc_deg=30.5,
        raan_deg=35.0,
        aop_deg=65.0,
        ta_deg=590,
        epoch=e,
        frame=eme2k,
    )

    assert orbit.period() == Unit.Second * 5553.623455582

    print(orbit)
    print(repr(orbit))

    # Define the SRP
    srp = SrpConfig(2.0)
    # Define the Drag
    drag = DragConfig(2.0)

    # Build a spacecraft
    sc = Spacecraft(
        orbit,
        dry_mass_kg=500.0,
        fuel_mass_kg=15.0,
        srp=srp,
        drag=drag,
    )

    print(sc)

    print(sc.orbit)

    # Check that we can pickle and dump the config of this spacecraft
    pkld = pickle.dumps(sc)
    print(sc.dumps())
    sc_loaded = Spacecraft.loads(sc.dumps())
    assert sc_loaded == sc
    unpkld = pickle.loads(pkld)
    assert unpkld == sc

    # NOTE: This does not return a pointer to the object, but a new object!
    # Therefore the equality _will_ fail!
    assert orbit == sc.orbit
    # But the printed information is identical
    assert f"{orbit}" == f"{sc.orbit}"
    # And the epochs match
    assert e.timedelta(sc.epoch) == Duration.zero()


def test_generate_states():
    """
    Tests that we can generate orbits from their state parameter deviations
    """
    # Build a demo orbit
    cosm = Cosm.de438()
    eme2k = cosm.frame("EME2000")

    assert eme2k.gm() == 398600.435392  # SPICE data, not GMAT data (398600.4415)
    assert eme2k.is_geoid()
    assert eme2k.equatorial_radius() == 6378.1363

    e = Epoch.system_now()

    orbit = Orbit.from_keplerian_altitude(
        400,
        ecc=1e-4,
        inc_deg=30.5,
        raan_deg=35.0,
        aop_deg=65.0,
        ta_deg=590,
        epoch=e,
        frame=eme2k,
    )

    # Generate a bunch of them in percentage of error
    orbits = generate_orbits(
        orbit,
        [
            # Note that we can create the state parameter from a string of its name
            (StateParameter("SMA"), 0.05),
            # Or directly as an enum (the preferred method)
            (StateParameter.Eccentricity, 0.1),
            (StateParameter.Inclination, 0.1),
        ],
        100,
        kind="prct",
    )
    assert len(orbits) >= 98

    # Define the SRP
    srp = SrpConfig(2.0)
    # Define the Drag
    drag = DragConfig(2.0)

    # Build a spacecraft
    sc = Spacecraft(
        orbit,
        dry_mass_kg=500.0,
        fuel_mass_kg=15.0,
        srp=srp,
        drag=drag,
    )

    # Generate a bunch of spacecraft with a specific seed, and using the absolute errors
    spacecraft = generate_spacecraft(
        sc,
        [
            (StateParameter.FuelMass, 1.9),
            (StateParameter.SMA, 1.0),
            (StateParameter.Eccentricity, 0.01),
            (StateParameter.Inclination, 3.5),
        ],
        100,
        kind="abs",
        seed=42,
    )
    assert len(spacecraft) >= 98


if __name__ == "__main__":
    test_define_spacecraft()
