use nyx_space::io::stations::StationSerde;
use nyx_space::io::tracking_data::DynamicTrackingArc;
use nyx_space::io::{ConfigRepr, Configurable};
use nyx_space::md::trajectory::ExportCfg;
use nyx_space::md::ui::*;
use nyx_space::od::msr::StdMeasurement;
use nyx_space::od::prelude::*;
use nyx_space::od::simulator::arc::TrackingArcSim;
use nyx_space::od::simulator::TrkConfig;
use std::collections::HashMap;
use std::env;
use std::path::PathBuf;
use std::str::FromStr;

#[test]
fn continuous_tracking() {
    // Test that continuous tracking
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    // Load cosm
    let cosm = Cosm::de438();

    // Dummy state
    let orbit = Orbit::keplerian_altitude(
        500.0,
        1e-3,
        30.0,
        45.0,
        75.0,
        23.4,
        Epoch::from_str("2023-02-22T19:18:17.16 UTC").unwrap(),
        cosm.frame("EME2000"),
    );

    // Generate a trajectory
    let (_, trajectory) = Propagator::default(OrbitalDynamics::two_body())
        .with(orbit)
        .for_duration_with_traj(1.5.days())
        .unwrap();

    println!("{trajectory}");

    // Save the trajectory to parquet
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "output_data",
        "tracking_truth_ephem.parquet",
    ]
    .iter()
    .collect();

    trajectory
        .to_parquet_with_cfg(
            path,
            ExportCfg {
                timestamp: true,
                ..Default::default()
            },
        )
        .unwrap();

    // Load the ground stations from the test data.
    let ground_station_yaml: PathBuf = [
        &env::var("CARGO_MANIFEST_DIR").unwrap(),
        "data",
        "tests",
        "config",
        "many_ground_stations.yaml",
    ]
    .iter()
    .collect();

    let stations_serde = StationSerde::load_many_yaml(ground_station_yaml).unwrap();
    let devices: Vec<GroundStation> = stations_serde
        .into_iter()
        .map(|station| GroundStation::from_config(station, cosm.clone()).unwrap())
        .collect();

    dbg!(&devices);

    // Load the tracking configuration from the test data.
    let trkconfg_yaml: PathBuf = [
        &env::var("CARGO_MANIFEST_DIR").unwrap(),
        "data",
        "tests",
        "config",
        "tracking_cfg.yaml",
    ]
    .iter()
    .collect();

    let configs: HashMap<String, TrkConfig> = TrkConfig::load_named_yaml(trkconfg_yaml).unwrap();

    dbg!(&configs);

    // Build the tracking arc simulation to generate a "standard measurement".
    let mut trk = TrackingArcSim::<Orbit, Orbit, StdMeasurement, _>::with_seed(
        devices, trajectory, configs, 12345,
    )
    .unwrap();

    let arc = trk.generate_measurements(cosm).unwrap();

    assert_eq!(arc.measurements.len(), 146);

    // And serialize to disk
    let path: PathBuf = [
        &env::var("CARGO_MANIFEST_DIR").unwrap(),
        "output_data",
        "simple_arc.parquet",
    ]
    .iter()
    .collect();

    let output_fn = arc.to_parquet_simple(path).unwrap();
    println!("[{}] {arc}", output_fn.to_string_lossy());

    // Now read this file back in.
    let dyn_arc = DynamicTrackingArc::from_parquet(output_fn).unwrap();
    // And convert to the same tracking arc as earlier
    let arc_concrete = dyn_arc.to_tracking_arc::<StdMeasurement>().unwrap();

    println!("{arc_concrete}");

    // Check that we've loaded all of the measurements
    assert_eq!(arc_concrete.measurements.len(), arc.measurements.len());
    // Check that we find the same device names too
    assert_eq!(arc_concrete.device_names(), arc.device_names());
    // Check that we've copied over the device configurations as well
    assert_eq!(arc_concrete.device_cfg, arc.device_cfg);
}
