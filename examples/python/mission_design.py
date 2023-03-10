import logging
from pathlib import Path


from nyx_space.cosmic import Spacecraft
from nyx_space.mission_design import *
from nyx_space.time import Duration, Unit

if __name__ == '__main__':

    # Initialize logging
    FORMAT = '%(levelname)s %(name)s %(filename)s:%(lineno)d\t%(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.INFO)


    # Base path
    root = Path(__file__).joinpath("../../../").resolve()

    config_path = root.joinpath("./data/tests/config/")

    sc = Spacecraft.load_yaml(str(config_path.joinpath("spacecraft.yaml")))
    print(sc)

    dynamics = SpacecraftDynamics.load_named_yaml(str(config_path.joinpath("dynamics.yaml")))
    print(dynamics)

    rslt, traj = propagate(sc, dynamics['lofi'], Unit.Day * 5.159)
    print(rslt)
    print(traj)
    rslt, traj = propagate(sc, dynamics['lofi'], Duration("5.159 hours"))
    print(rslt)
    print(traj)

    rslt, traj = propagate(sc, dynamics['lofi'], event=Event(StateParameter.Apoapsis, 0.0))
    print(rslt)
    print(traj)