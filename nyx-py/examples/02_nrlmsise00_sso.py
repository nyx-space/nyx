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

# NOTE Nyx uses the NAIF data for graviational parameters, which leads to an
# accumulation of error between Nyx and GMAT.
#
# 36x36 EGM96 + Sun + Moon: 43 meters after 7 days
#  ... + NRLMSISE00 constant:
#  ... + NRLMSISE00 SpaceWeather:

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
        6928.0,
        0.001,
        97.6,
        200.0,
        90.0,
        0.0,
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
        357.3206749596847,
        1138.437101108836,
        6817.742153372467,
        6.991584522403016,
        2.840019392253554,
        -0.835436591655589,
        Epoch("2023-03-04"),
        eme2k,
    )

    ric_diff = gmat_rslt.ric_difference(result.state.orbit)

    print(
        f"RIC diff with GMAT = {ric_diff.rmag_km() * 1e3:.3f} m\n\tR = {ric_diff.x_km * 1e3:.3f} m\tI = {ric_diff.y_km * 1e3:.3f} m\tC = {ric_diff.z_km * 1e3:.3f} m"
    )

    assert ric_diff.rmag_km()*1e3 < 78

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

    result = propagator.for_duration(my_sc, Duration("3 day"))

    gmat_rslt = Orbit(
        365.4022523482296,
        1141.713002498598,
        6816.733465998496,
        6.991100511804022,
        2.838453740184353,
        -0.8448775241656856,
        Epoch("2023-03-04"),
        eme2k,
    )

    ric_diff = gmat_rslt.ric_difference(result.state.orbit)

    print(
        f"RIC diff with GMAT = {ric_diff.rmag_km() * 1e3:.3f} m\n\tR = {ric_diff.x_km * 1e3:.3f} m\tI = {ric_diff.y_km * 1e3:.3f} m\tC = {ric_diff.z_km * 1e3:.3f} m"
    )
    assert ric_diff.rmag_km()*1e3 < 59
