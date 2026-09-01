import logging

from nyx_space import DragData, ExportCfg, Mass, Spacecraft, SRPData
from nyx_space.anise import Aberration, MetaAlmanac
from nyx_space.anise.analysis import Event
from nyx_space.anise.astro import DataType, Orbit
from nyx_space.anise.constants import CelestialObjects, Frames
from nyx_space.mission_design import (
    AccelModels,
    AtmDensity,
    Drag,
    Dynamics,
    ForceModels,
    GravityFieldConfig,
    IntegratorMethod,
    IntegratorOptions,
    Nrlmsise00Flags,
    PointMasses,
    Propagator,
    SolarPressure,
    SolidTides,
    SpaceWeatherData,
    StaticSpaceWeather,
)
from nyx_space.time import Duration, Epoch


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

    # Define Non-Gravitational Forces: force models depends on the vehicle mass.
    # Configure Solar Radiation Pressure targeting Earth and Moon frames, accounting for light-time aberration
    # because the direction of the photons depend on the position of the sun at time of emission.
    srp = SolarPressure(
        [Frames.EARTH_J2000, Frames.MOON_J2000], almanac, correction=Aberration("LT")
    )
    # Configure atmospheric drag using a standard exponential atmospheric model.
    # Load the Space Weather as provided by CelesTrak: https://celestrak.org/SpaceData/.
    # You should also provide a fallback for propagation that extends further from the
    # data set. This can be specified either as Solar Minimum, Solar Maximum, Solar Average,
    # or a custom number for the F10.7 (in SFU), Ap, and Kp values.
    # Data may be provided either as CSV or in non-archived gunzip (gz) format (decoded on the fly).
    weather = SpaceWeatherData(
        "../data/01_planetary/SpaceWeather-2021-01-01_2026-09-06.csv.gz",
        StaticSpaceWeather.SolarAverage(),
    )
    # NOTE While you must specify the flags for the NRLMSISE00 models, the default ones are typically what you want.
    # Call `help(Nrlmsise00Flag)` for details on available options.
    drag = Drag(
        AtmDensity.NRLMSISE00(weather, Nrlmsise00Flags()),
        almanac.frame_info(Frames.IAU_EARTH_FRAME),
    )

    # Combine accelerations and forces into a unified Dynamics model.
    # dynamics = Dynamics(accel_models, force_models=ForceModels(srp, drag))
    dynamics = Dynamics(accel_models)

    # Initialize the Propagator
    # The Propagator encapsulates the dynamics, the math (RungeKutta89), and the universe (almanac).
    # It remains stateless regarding the spacecraft, allowing it to be reused.
    propagator = Propagator(
        dynamics, almanac, IntegratorMethod.RungeKutta89, IntegratorOptions()
    )

    # Define the Spacecraft Initial State
    # We define the starting Keplerian orbital elements.
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

    # Define the physical characteristics of the spacecraft necessary to compute drag and SRP.
    # Surfaces are in square meters, masses are in kilograms.
    my_sc = Spacecraft(
        orbit,
        Mass(dry_mass_kg=1000.0),
        SRPData(area_m2=10.0, coeff_reflectivity=1.2),
        DragData(area_m2=10.0, coeff_drag=2.2),
    )

    # Propagate until a specific orbital event (periapsis).
    # This demonstrates the root-finding capabilities of Nyx to stop exactly at an event,
    # bounded by a maximum duration to prevent infinite searching.
    result = propagator.for_duration(my_sc, Duration("7 day"))

    print(result.state)

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
