import logging
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from nyx_space.analysis import diff_traj_parquet
from nyx_space.cosmic import Cosm, Spacecraft
from nyx_space.mission_design import SpacecraftDynamics, TrajectoryLoader, propagate
from nyx_space.orbit_determination import (
    DynamicTrackingArc,
    ExportCfg,
    GroundStation,
    GroundTrackingArcSim,
    OrbitEstimate,
    TrkConfig,
    predictor,
    process_tracking_arc,
)
from nyx_space.plots.od import (
    plot_covar,
    plot_estimates,
    plot_residual_histogram,
    plot_residuals,
)
from nyx_space.plots.traj import plot_orbit_elements
from nyx_space.time import TimeSeries, Unit


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
    # Resample the trajectory at fixed step size
    traj = traj.resample(Unit.Second * 10.0)
    # And save the trajectory
    traj_file = str(outpath.joinpath("./python_ref_traj.parquet"))
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
    # Check that we can pickle and compare
    trk_cfg_demo = trk_cfg["Demo ground station"]
    print(trk_cfg_demo.dumps())
    unpkld = pickle.loads(pickle.dumps(trk_cfg_demo))
    assert unpkld == trk_cfg_demo

    # Load the trajectory
    traj = TrajectoryLoader(traj_file)
    print(f"Loaded {traj}")

    # Set up the export -- We'll use the same config set up for both measurements and output of OD process
    cfg = ExportCfg(timestamp=True, metadata={"test key": "test value"})

    # Build the simulated tracking arc, setting the seed to zero
    arc_sim = GroundTrackingArcSim(devices, traj, trk_cfg, 0)
    # Generate the measurements
    msr_path = arc_sim.generate_measurements(
        str(outpath.joinpath("./msr.parquet")), cfg
    )
    print(f"Saved {arc_sim} to {msr_path}")

    # Now let's filter this same data.
    # Load the tracking arc
    arc = DynamicTrackingArc(msr_path)

    # Create the orbit estimate with the covariance diagonal (100 km on position and 1 km/s on velocity)
    orbit_est = OrbitEstimate(
        sc.orbit, covar=np.diag([100.0, 100.0, 100.0, 1.0, 1.0, 1.0])
    )

    # Check loading from the YAML read from Python
    with open(config_path.joinpath("orbit_estimates.yaml")) as fh:
        data = yaml.safe_load(fh)

    loaded = OrbitEstimate.loads(data["example 1"])
    print(loaded)

    msr_noise = [1e-3, 0, 0, 1e-6]
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
        # predict_for=Unit.Hour * 12, # You can predict from the final estimate by uncommenting this line.
    )

    print(f"Stored to {rslt_path}")

    # Load the results
    oddf = pd.read_parquet(rslt_path)
    # Load the reference trajectory
    ref_traj = pd.read_parquet(traj_file)
    # Load the measurements
    msr_df = pd.read_parquet(msr_path)

    # There seems to be a bug when exporting the HTML on Github action windows
    # cf. https://github.com/nyx-space/nyx/actions/runs/5064848025/jobs/9092830624
    if sys.platform != "win32":
        # We'll plot the difference between the reference trajectory and the OD results, with the measurement bands overlaid.
        plot_estimates(
            oddf,
            "OD results from Python",
            cov_frame="RIC",
            ref_traj=ref_traj,
            msr_df=msr_df,
            html_out=str(outpath.joinpath("./od_estimate_plots.html")),
            show=False,
        )

        # Let's also plot the reference and the OD result's orbital elements
        plot_orbit_elements(
            [ref_traj, oddf],
            "Python OD result",
            ["Reference", "OD"],
            html_out=str(outpath.joinpath("./od_vs_ref_elements.html")),
            show=False,
        )

        # More often, the covariance is a better indicator
        plot_covar(
            oddf,
            "OD 1-sigma covar from Python",
            cov_sigma=1.0,
            msr_df=msr_df,
            html_out=str(outpath.joinpath("./od_covar_plots.html")),
            show=False,
        )

        # Now, we'll plot the prefit residuals
        plot_residuals(
            oddf,
            "OD residuals",
            msr_df=msr_df,
            html_out=str(outpath.joinpath("./od_residual_plots.html")),
            show=False,
        )
        # And the postfit histograms
        plot_residual_histogram(
            oddf,
            "OD residuals",
            kind="Postfit",
            html_out=str(outpath.joinpath("./od_residual_histograms.html")),
            show=False,
        )

        # Plot the errors between the nominal and estimation
        err_df = diff_traj_parquet(traj_file, rslt_path)
        plot_orbit_elements(
            err_df,
            "Error in orbital elements",
            html_out=str(outpath.joinpath("./od_vs_ref_error_elements.html")),
            show=False,
        )


def test_one_way_msr():
    """
    Test that we can manually perform a one-way measurement
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

    end_sc, traj = propagate(sc, dynamics["hifi"], sc.orbit.period() * 1.1)
    print(end_sc)
    print(traj)

    # Let's build a dataframe of the range, doppler, azimuth, and elevation as seen from a ground station that sees the spacecraft a bunch
    gs = devices[1]
    print(f"Using {gs}")
    cosm = Cosm.de438()
    data = {
        "epoch": [],
        "range (km)": [],
        "doppler (km/s)": [],
        "azimuth (deg)": [],
        "elevation (deg)": [],
    }
    # Start by building a time series
    ts = TimeSeries(
        traj.first().epoch, traj.last().epoch, step=Unit.Minute * 30, inclusive=True
    )
    # And iterate over it
    for epoch in ts:
        orbit = traj.at(epoch).orbit
        try:
            range_km, doppler_km_s = gs.measure(orbit)
        except:
            # Spacecraft is not visible then, nothing to store
            pass
        else:
            # Also grab the azimuth and elevation angles
            az_deg, el_deg = gs.compute_azimuth_elevation(orbit, cosm)
            # And push to the data dictionary
            data["epoch"] += [str(epoch)]
            data["azimuth (deg)"] += [az_deg]
            data["elevation (deg)"] += [el_deg]
            data["range (km)"] += [range_km]
            data["doppler (km/s)"] += [doppler_km_s]
    # And convert to a data frame
    df = pd.DataFrame(data, columns=data.keys())
    assert len(df) == 8
    print(df.describe())

    # Test values
    range_km, doppler_km_s = devices[0].measure(end_sc.orbit)

    print(range_km, doppler_km_s)
    assert abs(range_km - 18097.562811514355) < 0.1
    assert abs(doppler_km_s - -0.2498238312640348) < 0.1

    # Azimuth and elevation

    az_deg, el_deg = devices[0].compute_azimuth_elevation(end_sc.orbit, cosm)

    assert abs(az_deg - 128.66181520071825) < 1e-10
    assert abs(el_deg - 27.904687635388676) < 1e-10


def test_pure_prediction():
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

    # Set up the export -- We'll use the same config set up for both measurements and output of OD process
    cfg = ExportCfg(timestamp=True, metadata={"test key": "test value"})

    # We'll assume that we have a good estimate of the spacecraft's orbit before we predict it forward in time
    # Hence, create the orbit estimate with the covariance diagonal (100 m on position and 50 m/s on velocity)
    orbit_est = OrbitEstimate(
        sc.orbit,
        covar=np.diag([100.0e-3, 100.0e-3, 100.0e-3, 50.0e-3, 50.0e-3, 50.0e-3]),
    )

    rslt_path = predictor(
        dynamics["hifi"],
        sc,
        orbit_est,
        Unit.Second * 15.0,
        str(outpath.joinpath("./od_pred_result.parquet")),
        cfg,
        predict_for=Unit.Hour * 12,
    )

    print(f"Stored to {rslt_path}")

    # Load the prediction results
    oddf = pd.read_parquet(rslt_path)

    # There seems to be a bug when exporting the HTML on Github action windows
    # cf. https://github.com/nyx-space/nyx/actions/runs/5064848025/jobs/9092830624
    if sys.platform != "win32":
        # Let's plot the OD result's orbital elements
        plot_orbit_elements(
            oddf,
            "OD prediction result",
            html_out=str(outpath.joinpath("./od_pred_elements.html")),
            show=False,
        )

        # More often, the covariance is a better indicator
        plot_covar(
            oddf,
            "OD 1-sigma covar from Python",
            cov_sigma=1.0,
            cov_frame="Earth J2000",
            html_out=str(outpath.joinpath("./od_pred_covar_plots.html")),
            show=False,
        )


if __name__ == "__main__":
    test_filter_arc()
