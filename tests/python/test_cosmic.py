from nyx_space.cosmic import Cosm, Orbit, Spacecraft
from nyx_space.time import Epoch, Unit, Duration


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

    # Build a spacecraft
    sc = Spacecraft(
        orbit,
        dry_mass_kg=500.0,
        fuel_mass_kg=15.0,
        srp_area_m2=2.0,
        drag_area_m2=2.0,
        cr=1.8,
        cd=2.1,
    )

    print(sc)

    print(sc.orbit)

    # NOTE: This does not return a pointer to the object, but a new object!
    # Therefore the equality _will_ fail!
    assert orbit != sc.orbit
    # But the printed information is identical
    assert f"{orbit}" == f"{sc.orbit}"
    # And the epochs match
    assert e.timedelta(sc.epoch) == Duration.zero()
