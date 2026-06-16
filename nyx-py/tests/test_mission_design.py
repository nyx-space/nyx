from nyx_space import DragData, Mass, Spacecraft, SRPData
from nyx_space.anise import Aberration, MetaAlmanac
from nyx_space.anise.astro import FrameUid, Orbit
from nyx_space.anise.constants import CelestialObjects, Frames
from nyx_space.mission_design import (
    AtmDensity,
    AccelModels,
    ForceModels,
    SolarPressure,
    Drag,
    GravityFieldConfig,
    IntegratorMethod,
    IntegratorOptions,
    PointMasses,
    PropagatorConfig,
    PropagatorInstance,
)
from nyx_space.time import Duration, Epoch, Unit


def test_prop_cfg():
    """
    Demonstrate how to configure a propagator.
    """

    opts = IntegratorOptions()  # Default options is a variable step
    assert str(opts) == "min_step: 1 ms, max_step: 45 min, tol: 1e-12, attempts: 50"
    assert str(IntegratorOptions(Duration("1 s"), Unit.Second * 1)) == "fixed step: 1 s"

    # Build the acceleration models
    accel_models = AccelModels(
        PointMasses(
            [CelestialObjects.EARTH, CelestialObjects.MOON, CelestialObjects.SUN],
            Aberration("LT"),
        ),
        # Define the location of the gravity field, and then the frame in which it applies
        GravityFieldConfig(
            20,
            20,
            "../data/01_planetary/EGM2008_to2190_TideFree.gz",
            Frames.IAU_EARTH_FRAME.to_frameuid(),
        ),
    )

    almanac = MetaAlmanac.latest()

    srp = SolarPressure(
        [Frames.EARTH_J2000, Frames.MOON_J2000], almanac, estimate=False
    )
    drag = Drag(
        AtmDensity.earth_exponential(),
        almanac.frame_info(Frames.IAU_EARTH_FRAME),
        estimate=False,
    )
    print(f"srp = {srp}\ndrag = {drag}")

    cfg = PropagatorConfig(
        accel_models,
        ForceModels(srp, drag),
        IntegratorMethod.RungeKutta89,
        opts,
    )

    eme2k = almanac.frame_info(Frames.EME2000)
    # Define an orbit
    orbit = Orbit.from_keplerian(
        6800.0, 1e-4, 45.0, 60.0, 75.0, 90.0, Epoch("2030-12-29 01:02:03 TDB"), eme2k
    )

    my_sc = Spacecraft(orbit, Mass(123.0), SRPData(10.0, 1.2), DragData(10.0, 2.0))

    inst = PropagatorInstance(cfg, my_sc, almanac)

    result = inst.until_epoch(Epoch("2031-01-01 01:02:03 UTC"), trajectory=True)
    traj = result.trajectory
    breakpoint()


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    test_prop_cfg()
