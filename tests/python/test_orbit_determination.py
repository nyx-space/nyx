import logging
from pathlib import Path
import numpy as np

from nyx_space.orbit_determination import (
    GroundStation,
    GroundTrackingArcSim,
    TrkConfig,
    OrbitEstimate,
    process_tracking_arc,
    DynamicTrackingArc,
    ExportCfg,
)
from nyx_space.mission_design import DynamicTrajectory, SpacecraftDynamics, propagate
from nyx_space.cosmic import Spacecraft
from nyx_space.time import Unit
from nyx_space.plots.od import (
    plot_covar,
    plot_state_deviation,
    plot_estimates,
    plot_residuals,
    plot_residual_histogram,
    plot_residual_qq,
    plot_qq,
)


def test_filter_arc():
    # Initialize logging
    FORMAT = "%(levelname)s %(name)s %(asctime)-15s %(filename)s:%(lineno)d %(message)s"
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.INFO)

    # Base path
    root = Path(__file__).joinpath("../../../").resolve()
    config_path = root.joinpath("./data/tests/config/")
    outpath = root.joinpath("output_data/")

    # Load the dynamics and spacecraft
    sc = Spacecraft.load(str(config_path.joinpath("spacecraft.yaml")))
    dynamics = SpacecraftDynamics.load_named(str(config_path.joinpath("dynamics.yaml")))

    # An propagate for two periods (we only care about the trajectory)
    _, traj = propagate(sc, dynamics["hifi"], sc.orbit.period() * 2)
    # And save the trajectory
    traj_file = str(outpath.joinpath("./od_val_with_arc_truth_ephem.parquet"))
    traj.to_parquet(traj_file)

    # Now starts the measurement generation

    # Load the devices
    devices = GroundStation.load_many(
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
    trk_cfg = TrkConfig.load_named(str(config_path.joinpath("tracking_cfg.yaml")))

    # Load the trajectory
    traj = DynamicTrajectory(traj_file)
    print(f"Loaded {traj}")

    # Set up the export -- We'll use the same config set up for both measurements and output of OD process
    cfg = ExportCfg(timestamp=True, metadata={"test key": "test value"})

    # Build the simulated tracking arc, setting the seed to zero
    arc_sim = GroundTrackingArcSim(devices, traj, trk_cfg, 0)
    # Generate the measurements
    path = arc_sim.generate_measurements(
        str(outpath.joinpath("./from_python.parquet")), cfg
    )
    print(f"Saved {arc_sim} to {path}")

    # Now let's filter this same data.
    # Load the tracking arc
    arc = DynamicTrackingArc(path)

    # Create the orbit estimate with the covariance diagonal (100 km on position and 1 km/s on velocity)
    orbit_est = OrbitEstimate(
        sc.orbit, covar=np.diag([100.0, 100.0, 100.0, 1.0, 1.0, 1.0])
    )

    msr_noise = [1e-3, 0, 0, 1e-6]  # TODO: Convert this to a numpy matrix
    # Switch from sequential to EKF after 100 measurements
    ekf_num_msr_trig = 100
    # Unless there is a 2 hour gap in the measurements, and then switch back to classical
    ekf_disable_time = Unit.Hour * 2

    rslt_path = process_tracking_arc(
        dynamics["hifi"],
        sc,
        orbit_est,
        msr_noise,
        arc,
        str(outpath.joinpath("./od_result.parquet")),
        cfg,
        ekf_num_msr_trig,
        ekf_disable_time,
    )

    print(f"Stored to {rslt_path}")

    # TODO: Add plotting and saving of the figures
    # TODO: Add figures to artifacts


def test_one_way_msr():
    """
    Test that we can perform a one-way measurement
    """

    # Base path
    root = Path(__file__).joinpath("../../../").resolve()
    config_path = root.joinpath("./data/tests/config/")

    # Load the dynamics and spacecraft
    sc = Spacecraft.load(str(config_path.joinpath("spacecraft.yaml")))
    dynamics = SpacecraftDynamics.load_named(str(config_path.joinpath("dynamics.yaml")))

    # Load the devices
    devices = GroundStation.load_many(
        str(config_path.joinpath("./many_ground_stations.yaml"))
    )
    print(f"Loaded {devices}")

    # One way measurement

    end_sc, _ = propagate(sc, dynamics["hifi"], sc.orbit.period() * 1.1)
    print(end_sc)

    range_km, doppler_km_s = devices[0].measure(end_sc.orbit)

    print(range_km, doppler_km_s)
    assert abs(range_km - 18097.562811514355) < 0.1
    assert abs(doppler_km_s - -0.2498238312640348) < 0.1


if __name__ == "__main__":
    test_one_way_msr()
