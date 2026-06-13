from nyx_space.mission_design import (
    PropagatorConfig,
    IntegratorMethod,
    AccelModels,
    ForceModels,
    PropagatorInstance,
)
from nyx_space.time import Duration, Unit
from nyx_space.anise import Almanac


def test_prop_cfg():
    """
    Demonstrate how to configure a propagator.
    """

    opts = IntegratorOptions()  # Default options is a variable step
    assert str(opts) == "min_step: 1 ms, max_step: 45 min, tol: 1e-12, attempts: 50"
    assert str(IntegratorOptions(Duration("1 s"), Unit.Second * 1)) == "fixed step: 1 s"

    # TODO Export/config accel/force models

    cfg = PropagatorConfig(
        accel_models, force_models, IntegratorMethod.RungeKutta89, options
    )

    inst = PropagatorInstance(cfg, my_spacecraft, almanac)

    final_sc, traj = inst.until_epoch(Epoch("2031-01-01"), trajectory=True)
