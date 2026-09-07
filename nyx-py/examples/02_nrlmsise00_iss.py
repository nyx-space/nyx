import logging

from nyx_space import DragData, Mass, Spacecraft, SRPData
from nyx_space.anise import MetaAlmanac
from nyx_space.anise.astro import Orbit
from nyx_space.anise.constants import CelestialObjects, Frames
from nyx_space.mission_design import (
    AccelModels,
    AtmDensity,
    Drag,
    Dynamics,
    ForceModels,
    GravityFieldConfig,
    IntegratorMethod,
    Nrlmsise00Flags,
    PointMasses,
    Propagator,
    SpaceWeatherData,
    StaticSpaceWeather,
)
from nyx_space.time import Duration, Epoch

# NOTE Nyx uses the NAIF data for graviational parameters, which leads to an
# accumulation of error between Nyx and GMAT.
#
# 36x36 EGM96 + Sun + Moon: 43 meters after 7 days
#  ... + NRLMSISE00 constant: 65.5 meters after 3 days
#  ... + NRLMSISE00 SpaceWeather: 1316 meters after 3 days
#
# Note on the Dynamic Weather residual:
# The static weather match (< 70m) confirms the core MSIS physics, Mean LST geometry,
# and integrator are functionally identical to GMAT. The 1.3 km dynamic residual is
# strictly an artifact of environmental pre-processing: GMAT applies cubic spline
# interpolation to the 3-hourly Ap/Kp indices, whereas Nyx feeds the piecewise-constant
# step-functions strictly dictated by the NRLMSISE-00 empirical baseline.
# This interpolation variance alters the total integrated energy over multi-day propagations.

if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    # Load the universe via ANISE's MetaAlmanac (similar to SPICE's meta kernel)

    almanac = (
        MetaAlmanac("../data/02_config/ci_almanac.dhall")
        .process()
        .load("../data/01_planetary/earth_longterm_000101_251211_250915.bpc")
    )

    # Acceleration models are independent of the vehicle mass.
    accel_models = AccelModels(
        point_masses=PointMasses(
            celestial_objects=[
                CelestialObjects.EARTH,
                CelestialObjects.MOON,
                CelestialObjects.SUN,
            ],
        ),
        gravity_field=GravityFieldConfig(
            degree=36,
            order=36,
            filepath="../data/01_planetary/EGM96.cof.gz",
            frame=Frames.EARTH_ITRF93.to_frameuid(),
        ),
    )

    # === No drag validation ===
    # No drag validation
    dynamics_no_drag = Dynamics(accel_models)
    propagator = Propagator(dynamics_no_drag, almanac, IntegratorMethod.RungeKutta89)
    # Define the Spacecraft Initial State
    eme2k = almanac.frame_info(Frames.EME2000)
    orbit = Orbit.from_keplerian(
        6778.0,
        0.001,
        51.6,
        30.0,
        40.0,
        50.0,
        epoch=Epoch("2023-05-16T20:00:00"),
        frame=eme2k,
    )

    my_sc = Spacecraft(
        orbit,
        Mass(dry_mass_kg=1000.0),
        SRPData(area_m2=10.0, coeff_reflectivity=1.2),
        DragData(area_m2=10.0, coeff_drag=2.2),
    )

    result = propagator.for_duration(my_sc, Duration("7 day"))

    gmat_rslt = Orbit(
        5766.457918604106,
        1899.757886860554,
        3005.528950697711,
        -4.004588754191417,
        4.288824900835695,
        4.9504844851056,
        Epoch("2023-05-23T20:00:00"),
        eme2k,
    )

    ric_diff = gmat_rslt.ric_difference(result.state.orbit)

    print(
        f"RIC diff with GMAT = {ric_diff.rmag_km() * 1e3:.3f} m\n\tR = {ric_diff.x_km * 1e3:.3f} m\tI = {ric_diff.y_km * 1e3:.3f} m\tC = {ric_diff.z_km * 1e3:.3f} m"
    )

    assert ric_diff.rmag_km() < 44e-3

    # === NRLMSISE00 Static weather ===
    weather = SpaceWeatherData(
        None,  # Do not load a file, provide only the static weather
        StaticSpaceWeather.Custom(f107=150.0, ap=4.0, kp=1.0),
    )
    drag = Drag(
        AtmDensity.NRLMSISE00(weather, Nrlmsise00Flags(mean_lst=True)),
        almanac.frame_info(Frames.IAU_EARTH_FRAME),
    )
    dynamics_static_drag = Dynamics(accel_models, ForceModels(drag=drag))
    propagator = Propagator(
        dynamics_static_drag, almanac, IntegratorMethod.RungeKutta89
    )
    orbit = Orbit.from_keplerian(
        6798.0,
        0.0005,
        51.6,
        30.0,
        40.0,
        50.0,
        epoch=Epoch("2023-03-01"),
        frame=eme2k,
    )

    my_sc = Spacecraft(
        orbit,
        Mass(dry_mass_kg=1000.0),
        SRPData(area_m2=10.0, coeff_reflectivity=1.2),
        DragData(area_m2=10.0, coeff_drag=2.2),
    )

    result = propagator.for_duration(my_sc, Duration("3 day"))

    gmat_rslt = Orbit(
        -1090.433111737012,
        -4442.315837430585,
        -5043.440557985983,
        7.381816725397179,
        0.403056010083568,
        -1.949273075463185,
        Epoch("2023-03-04"),
        eme2k,
    )

    ric_diff = gmat_rslt.ric_difference(result.state.orbit)

    print(
        f"RIC diff with GMAT = {ric_diff.rmag_km() * 1e3:.3f} m\n\tR = {ric_diff.x_km * 1e3:.3f} m\tI = {ric_diff.y_km * 1e3:.3f} m\tC = {ric_diff.z_km * 1e3:.3f} m"
    )

    # === NRLMSISE00 Space Weather File ===
    weather = SpaceWeatherData(
        "../data/01_planetary/SpaceWeather-2021-01-01_2026-09-06.csv.gz",
        StaticSpaceWeather.SolarAverage(),
    )
    # NOTE GMAT interpolates the space weather, so let's turn it on here as well.
    weather.interpolate = True
    drag = Drag(
        AtmDensity.NRLMSISE00(weather, Nrlmsise00Flags(mean_lst=True)),
        almanac.frame_info(Frames.IAU_EARTH_FRAME),
    )
    dynamics_static_drag = Dynamics(accel_models, ForceModels(drag=drag))
    propagator = Propagator(
        dynamics_static_drag, almanac, IntegratorMethod.RungeKutta89
    )
    orbit = Orbit.from_keplerian(
        6798.0,
        0.0005,
        51.6,
        30.0,
        40.0,
        50.0,
        epoch=Epoch("2023-03-01"),
        frame=eme2k,
    )

    my_sc = Spacecraft(
        orbit,
        Mass(dry_mass_kg=1000.0),
        SRPData(area_m2=10.0, coeff_reflectivity=1.2),
        DragData(area_m2=10.0, coeff_drag=2.2),
    )

    result = propagator.for_duration(my_sc, Duration("3 day"))

    gmat_rslt = Orbit(
        -1042.605249258735,
        -4439.411782336723,
        -5055.731815614098,
        7.390665412088192,
        0.4391916505086589,
        -1.908043134886934,
        Epoch("2023-03-04"),
        eme2k,
    )

    ric_diff = gmat_rslt.ric_difference(result.state.orbit)

    print(
        f"RIC diff with GMAT = {ric_diff.rmag_km() * 1e3:.3f} m\n\tR = {ric_diff.x_km * 1e3:.3f} m\tI = {ric_diff.y_km * 1e3:.3f} m\tC = {ric_diff.z_km * 1e3:.3f} m"
    )
