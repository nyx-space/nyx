import logging
from pathlib import Path

from nyx_space.orbit_determination import GroundStation, GroundTrackingArcSim, TrkConfig
from nyx_space.mission_design import DynamicTrajectory, SpacecraftDynamics, propagate
from nyx_space.cosmic import Spacecraft


def test_generate_msr():
    # Initialize logging
    FORMAT = "%(levelname)s %(name)s %(asctime)-15s %(filename)s:%(lineno)d %(message)s"
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.INFO)

    # Base path
    root = Path(__file__).joinpath("../../../").resolve()
    config_path = root.joinpath("./data/tests/config/")
    outpath = root.joinpath("output_data/")

    # Load the dynamics and spacecraft
    sc = Spacecraft.load_yaml(str(config_path.joinpath("spacecraft.yaml")))
    dynamics = SpacecraftDynamics.load_named_yaml(
        str(config_path.joinpath("dynamics.yaml"))
    )

    # An propagate for two periods (we only care about the trajectory)
    _, traj = propagate(sc, dynamics["hifi"], sc.orbit.period() * 2)
    # And save the trajectory
    traj_file = str(outpath.joinpath("./od_val_with_arc_truth_ephem.parquet"))
    traj.to_parquet(traj_file)

    # Now starts the measurement generation

    # Load the devices
    devices = GroundStation.load_many_yaml(
        str(root.joinpath("./data/tests/config/many_ground_stations.yaml"))
    )
    print(f"Loaded {devices}")
    for device in devices:
        # Check that we have loaded them correctly
        if device.name == "Demo ground station":
            assert device.latitude_deg == 2.3522
            assert device.longitude_deg == 48.8566
        else:
            # This is Canberra, inherited from the build-in ground station
            assert device.latitude_deg == -35.398333
            assert device.longitude_deg == 148.981944
        assert device.elevation_mask_deg == 5.0  # Both have the same mask

    # Load the tracking configuration as a dictionary
    trk_cfg = TrkConfig.load_named_yaml(str(config_path.joinpath("tracking_cfg.yaml")))

    # Load the trajectory
    traj = DynamicTrajectory(traj_file)
    print(f"Loaded {traj}")
    # Build the simulated tracking arc, setting the seed to zero
    arc_sim = GroundTrackingArcSim(devices, traj, trk_cfg, 0)
    # Generate the measurements
    path = arc_sim.generate_measurements(str(outpath.joinpath("./from_python.parquet")))
    print(f"Saved {arc_sim} to {path}")
