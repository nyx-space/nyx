import logging
from pathlib import Path

import nyx_space as nyx

if __name__ == "__main__":
    # Initialize logging
    FORMAT = '%(levelname)s %(name)s %(asctime)-15s %(filename)s:%(lineno)d %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.INFO)

    # Base path
    root = Path(__file__).joinpath("../../../").resolve()

    outpath = root.joinpath("output_data/")

    # Load the devices
    devices = nyx.orbit_determination.GroundStation.load_many_yaml(str(root.joinpath("./data/tests/config/many_ground_stations.yaml")))
    print(f"Loaded {devices}")
    for device in devices:
        # Check that we can set the latitude and longitude
        device.latitude_deg += 5.0
        device.longitude_deg -= 5.0
        print(f"{device.name} at lat = {device.latitude_deg} deg and long = {device.longitude_deg} deg")
    # Load the trajectory
    traj = nyx.mission_design.DynamicTrajectory(str(outpath.joinpath("./od_val_with_arc_truth_ephem.parquet")))
    print(f"Loaded {traj}")
    # Build the simulated tracking arc
    arc_sim = nyx.orbit_determination.GroundTrackingArcSim(devices, traj, 0)
    # Generate the measurements
    path = arc_sim.generate_measurements(str(outpath.joinpath("./from_python.parquet")))
    print(f"Saved {arc_sim} to {path}")