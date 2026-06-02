from nyx_space.orbit_determination import (
    FrameUid,
    GroundStation,
    Location,
    Measurement,
    MeasurementType,
    TrackingDataArc,
)


def test_ground_station():
    """
    Demonstrate ground station config from a script, saving to a Yaml, and reloading
    """
    gs = GroundStation(
        "Paris, FR",
        Location(2.3522, 48.8566, 0.4, FrameUid(399, 399), [], True),
        [
            MeasurementType.Range,
            MeasurementType.Doppler,
            MeasurementType.Elevation,
            MeasurementType.Azimuth,
        ],
    )
    as_yaml = gs.to_yaml()
    gs_reloaded = GroundStation.from_yaml(as_yaml)
    assert str(gs) == str(gs_reloaded)
    assert gs == gs_reloaded
    print(gs)
    print(repr(gs))


def test_generate_tdm():
    """
    Demonstrate that we can generate TDMs
    """
