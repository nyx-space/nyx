from nyx_space import DragData, Mass, Spacecraft, SRPData
from nyx_space.anise import MetaAlmanac
from nyx_space.anise.astro import Orbit
from nyx_space.anise.constants import Frames
from nyx_space.mission_design import (
    IntegratorMethod,
    IntegratorOptions,
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

    # TODO Export/config accel/force models

    cfg = PropagatorConfig(
        accel_models, force_models, IntegratorMethod.RungeKutta89, opts
    )

    almanac = MetaAlmanac.latest()
    eme2k = almanac.frame_info(Frames.EME2000)
    # Define an orbit
    orbit = Orbit.from_keplerian(
        6800.0, 1e-4, 45.0, 60.0, 75.0, 90.0, Epoch("2020-02-29 01:02:03 TDB"), eme2k
    )

    my_sc = Spacecraft(
        orbit, Mass(123.0), SRPData(10.0, 1.2), DragData(10.0, 2.0)
    )

    inst = PropagatorInstance(cfg, my_sc, almanac)

    final_sc, traj = inst.until_epoch(Epoch("2031-01-01"), trajectory=True)
