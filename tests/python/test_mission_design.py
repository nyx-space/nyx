import logging
import pickle
import sys
from pathlib import Path
from timeit import timeit

import pandas as pd
import yaml
from nyx_space.cosmic import Orbit, Spacecraft, SrpConfig
from nyx_space.mission_design import (
    Event,
    SpacecraftDynamics,
    StateParameter,
    TrajectoryLoader,
    propagate,
    two_body,
)
from nyx_space.monte_carlo import StateParameter, generate_orbits
from nyx_space.plots import plot_traj_errors
from nyx_space.time import Duration, Epoch, TimeSeries, Unit


def test_propagate():
    # Initialize logging
    FORMAT = "%(levelname)s %(name)s %(filename)s:%(lineno)d\t%(message)s"
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.INFO)

    # Base path
    root = Path(__file__).joinpath("../../../").resolve()

    config_path = root.joinpath("./data/tests/config/")

    sc = Spacecraft.load(str(config_path.joinpath("spacecraft.yaml")))
    # Check that we have loaded this correctly by checking the values from the YAML file
    assert sc.value_of(StateParameter.X) == -9042.862234
    assert sc.value_of(StateParameter.Y) == 18536.333069
    assert sc.value_of(StateParameter.Z) == 6999.957069
    assert sc.value_of(StateParameter.VX) == -3.288789
    assert sc.value_of(StateParameter.VY) == -2.226285
    assert sc.value_of(StateParameter.VZ) == 1.646738
    assert sc.value_of(StateParameter.Cr) == 1.0
    assert sc.value_of(StateParameter.Cd) == 2.2
    assert sc.value_of(StateParameter.FuelMass) == 50.0
    assert sc.epoch.timedelta(Epoch("2018-09-15T00:15:53.098 UTC")) == Duration.zero()
    # Check other stuff that is computed on request
    assert sc.value_of(StateParameter.SMA) == 21999.99774705774

    dynamics = SpacecraftDynamics.load_named(str(config_path.joinpath("dynamics.yaml")))
    # So far, there are no accessors on the dynamics, so we will check what they print out =(
    assert (
        f"{dynamics['lofi']}"
        == "Spacecraft dynamics (with guidance = false): No force models; Orbital dynamics: Point masses of Sun J2000, Earth J2000"
    )
    assert (
        f"{dynamics['hifi']}"
        == "Spacecraft dynamics (with guidance = false): SRP with φ = 1367 W/m^2 and eclipse light-source: Sun J2000, shadows casted by: Sun J2000, Moon J2000;  Orbital dynamics: Point masses of Sun J2000, Earth J2000, Moon J2000; IAU Earth gravity field 10x10 (order x degree)"
    )

    rslt, traj = propagate(sc, dynamics["lofi"], Unit.Day * 5.159)
    # Check that we propagated for the correct duration
    assert rslt.epoch.timedelta(sc.epoch) == Duration("5 days 3 h 48 min 57 s 600 ms")
    # We aren't in two body dynamics so the SMA has changed among other things.
    assert rslt.value_of(StateParameter.SMA) == 22000.024381883788
    # Grab a specific epoch from the trajectory
    sc_state = traj.at(Epoch("2018-09-16T00:16:53 TDB"))
    # Check that we got the right epoch
    assert sc_state.epoch.timedelta(Epoch("2018-09-16T00:16:53 TDB")) == Duration.zero()
    # Save the file to parquet
    traj.to_parquet("lofi.parquet")

    # We can also propagate with a different method
    rslt, traj = propagate(sc, dynamics["lofi"], Unit.Day * 5.159, method="Dormand78")
    assert rslt.epoch.timedelta(sc.epoch) == Duration("5 days 3 h 48 min 57 s 600 ms")

    # Let's now propagate the original spacecraft to its apoapsis, but let's search no more than 2 orbit periods
    rslt_apo, _ = propagate(
        sc,
        dynamics["lofi"],
        sc.orbit.period().unwrap() * 2,
        event=Event(StateParameter.Apoapsis, 0.0),
    )
    # Let's make sure that worked to within the default precision
    assert abs(rslt_apo.value_of(StateParameter.TrueAnomaly) - 180.0) < 0.001

    event = Event(StateParameter.Apoapsis, 0.0, value_precision=1e-6)

    # But we can also be even more precise
    rslt_apo, traj = propagate(
        sc,
        dynamics["lofi"],
        sc.orbit.period().unwrap() * 2,
        event=event,
    )
    assert abs(rslt_apo.orbit.ta_deg() - 180.0) <= 1e-6

    # Resample the trajectory at fixed step size of 25 seconds
    traj = traj.resample(Unit.Second * 25.0)

    # Rebuild the trajectory with specific epochs
    ts = TimeSeries(
        traj.first().epoch,
        traj.last().epoch,
        step=Duration("0.3 min 56 s 15 ns"),
        inclusive=True,
    )
    epochs = [epoch for epoch in ts]
    rebuilt_traj = traj.rebuild(epochs[1:-1])
    assert rebuilt_traj.first().epoch == epochs[1]
    assert rebuilt_traj.last().epoch == epochs[-2]

    # Export this trajectory with additional metadata and the events
    # Base path
    root = Path(__file__).joinpath("../../../").resolve()
    outpath = root.joinpath("output_data/")

    # Note: Python interface only supports strings for paths, not Path objects.
    traj.to_parquet(
        str(outpath.joinpath("./lofi_with_events.parquet")),
        metadata={"dynamics": str(dynamics["lofi"])},
        events=[event],
    )

    # Propagate until the lo-fi epoch to compare
    _, traj_hifi = propagate(
        sc,
        dynamics["hifi"],
        epoch=rslt_apo.epoch,
    )

    # Plot both trajectories in RIC frame
    if sys.platform != "win32":
        diff_path = str(outpath.joinpath("./lofi_hifi_ric_diff.parquet"))
        traj.ric_diff_to_parquet(traj_hifi, diff_path)

        traj_diff_ric = pd.read_parquet(diff_path)

        plot_traj_errors(
            traj_diff_ric,
            "LoFi vs HiFi",
            html_out=outpath.joinpath("./md_ric_lofi_hifi_diff.html"),
            show=False,
        )

        plot_traj_errors(
            traj_diff_ric,
            "LoFi vs HiFi",
            html_out=outpath.joinpath("./md_ric_lofi_hifi_diff_velocity.html"),
            vertical=True,
            velocity=True,
            show=False,
        )

    # Resample the trajectory at fixed step size of 25 seconds
    traj = traj.resample(Unit.Second * 25.0)

    # Also export this ground track
    cosm = Cosm.de438()
    iau_earth = cosm.frame("IAU Earth")
    traj.to_parquet("iau_earth_lofi.parquet", groundtrack=iau_earth)

    # Let's also search for this event in the trajectory
    for sc_at_event in traj.find(event):
        print(sc_at_event)
        assert abs(sc_at_event.value_of(StateParameter.TrueAnomaly) - 180.0) <= 1e-6


def test_build_spacecraft():
    """
    Tests that we can build a spacecraft from scratch without an input file
    """
    # Initialize logging
    FORMAT = "%(levelname)s %(name)s %(filename)s:%(lineno)d\t%(message)s"
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.INFO)

    cosm = Cosm.de438()
    eme2k = cosm.frame("EME2000")  # Earth Mean Equator J2000
    orbit = Orbit.from_try_keplerian_altitude(
        400.0, 0.01, 15.6, 45.0, 90.0, 75.0, Epoch.system_now(), eme2k
    )
    # Define the SRP
    srp = SrpConfig(2.0)  # 2.0 will be the area
    print(srp)
    assert srp.cr == 1.8  # Default value
    assert srp.area_m2 == 2.0
    sc = Spacecraft(orbit, 150.0, 15.0, srp)
    print(sc)

    assert sc.srp() == srp
    assert sc.drag().area_m2 == 0.0
    assert sc.drag().cd == 2.2  # Default value, but the area is zero, so it doesn't have any effect

    # Using this spacecraft as a template, let's load an OEM file, convert it to Parquet, and ensure we can load it back in.
    # The orbit data will be overwritten with data from the OEM file.
    root = Path(__file__).joinpath("../../../").resolve()

    config_path = root.joinpath("./data/tests/ccsds/oem/LEO_10s.oem")
    output_path = root.joinpath("./output_data/LEO_10s.parquet")
    # Convert from OEM to Parquet will happen on load
    traj = TrajectoryLoader(str(config_path), "oem", str(output_path), sc)
    print(traj)
    # Check that we can pickle the trajectory loader object
    traj_pkl = pickle.dumps(traj)
    traj_unpkl = pickle.loads(traj_pkl)
    assert traj_unpkl == traj
    # Check that we can convert this to a spacecraft trajectory
    traj_sc = traj.to_spacecraft_traj()
    traj_orbit = traj.to_orbit_traj()
    traj_orbit_dc = traj_sc.downcast()
    # Check that we can query it (will raise an exception if we can't, thereby failing the test)
    ts = TimeSeries(
        Epoch("2020-06-01T12:00:00.000000"),
        Epoch("2020-06-01T13:00:00.000000"),
        step=Unit.Minute * 17 + Unit.Second * 13.8,
        inclusive=True,
    )
    for epoch in ts:
        orbit = traj_orbit.at(epoch)
        dc_orbit = traj_orbit_dc.at(epoch)
        sc_orbit = traj_sc.at(epoch).orbit
        # Check params individually
        assert orbit.epoch == sc_orbit.epoch
        assert orbit.x_km == sc_orbit.x_km
        assert orbit.vx_km_s == sc_orbit.vx_km_s
        # Check the downcasted version
        assert dc_orbit.x_km == sc_orbit.x_km
        assert dc_orbit.vx_km_s == sc_orbit.vx_km_s
        assert dc_orbit == sc_orbit


def test_two_body():
    # Build a demo orbit
    cosm = Cosm.de438()
    eme2k = cosm.frame("EME2000")

    assert eme2k.gm() == 398600.435392  # SPICE data, not GMAT data (398600.4415)
    assert eme2k.is_geoid()
    assert eme2k.equatorial_radius() == 6378.1363

    e = Epoch.system_now()

    orbit = Orbit.from_try_keplerian_altitude(
        400,
        ecc=1e-4,
        inc_deg=30.5,
        raan_deg=35.0,
        aop_deg=65.0,
        ta_deg=590,
        epoch=e,
        frame=eme2k,
    )

    orbits = generate_orbits(
        orbit,
        [
            (StateParameter("SMA"), 0.05),
            (StateParameter.Eccentricity, 0.1),
            (StateParameter.Inclination, 0.1),
        ],
        1000,
        kind="prct",
    )

    # And propagate in parallel using a single duration
    proped_orbits = two_body(orbits, durations=[Unit.Day * 531.5])
    assert len(proped_orbits) >= len(orbits) - 2

    # And propagate in parallel using many epochs
    ts = TimeSeries(e, e + Unit.Day * 1000, step=Unit.Day * 1, inclusive=False)
    epochs = [e for e in ts]
    proped_orbits = two_body(orbits, new_epochs=epochs)
    # Allow up to two to fail
    assert len(proped_orbits) >= len(orbits) - 2

    timing = timeit(lambda: two_body(orbits, new_epochs=epochs), number=1)
    print(f"two body propagation of {len(orbits)} orbits in {timing} s")


def test_merge_traj():
    # Initialize logging
    FORMAT = "%(levelname)s %(name)s %(filename)s:%(lineno)d\t%(message)s"
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.INFO)

    # Base path
    root = Path(__file__).joinpath("../../../").resolve()

    config_path = root.joinpath("./data/tests/config/")

    sc1 = Spacecraft.load(str(config_path.joinpath("spacecraft.yaml")))

    dynamics = SpacecraftDynamics.load_named(str(config_path.joinpath("dynamics.yaml")))["lofi"]

    # Check loading from the YAML read from Python
    with open(config_path.joinpath("dynamics.yaml")) as fh:
        data = yaml.safe_load(fh)

    loaded = SpacecraftDynamics.loads(data["lofi"])
    assert loaded == dynamics

    sc2, traj1 = propagate(sc1, dynamics, Unit.Day * 5)
    # And propagate again
    sc3, traj2 = propagate(sc2, dynamics, Unit.Day * 5)
    # Add the trajectories
    traj = traj1 + traj2

    assert traj.last().epoch == sc3.epoch, f"{traj.last()} != {sc3}"

    # Convert into another frame and try to add them too.
    # We only check the epoch this time.
    traj1_moon = traj1.to_frame("Moon J2000")
    traj2_moon = traj2.to_frame("Moon J2000")

    traj_moon = traj1_moon + traj2_moon

    assert traj_moon.last().epoch == sc3.epoch
    print(traj_moon)


if __name__ == "__main__":
    test_propagate()
    test_merge_traj()
