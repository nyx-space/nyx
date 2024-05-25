use nyx_space::io::tracking_data::DynamicTrackingArc;
use nyx_space::io::ConfigRepr;
use nyx_space::md::prelude::*;
use nyx_space::md::trajectory::ExportCfg;
use nyx_space::od::msr::RangeDoppler;
use nyx_space::od::prelude::*;
use nyx_space::od::simulator::TrackingArcSim;
use nyx_space::od::simulator::TrkConfig;
use std::collections::BTreeMap;
use std::env;
use std::path::PathBuf;
use std::str::FromStr;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn continuous_tracking(almanac: Arc<Almanac>) {
    // Test that continuous tracking
    let _ = pretty_env_logger::try_init();

    // Dummy state
    let orbit = Orbit::try_keplerian_altitude(
        500.0,
        1e-3,
        30.0,
        45.0,
        75.0,
        23.4,
        Epoch::from_str("2023-02-22T19:18:17.16 UTC").unwrap(),
        almanac.frame_from_uid(EARTH_J2000).unwrap(),
    )
    .unwrap();

    // Generate a trajectory
    let (_, trajectory) = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()))
        .with(orbit.into(), almanac.clone())
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
    let ground_station_file: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "tests",
        "config",
        "many_ground_stations.yaml",
    ]
    .iter()
    .collect();

    let devices = GroundStation::load_many(ground_station_file).unwrap();

    // dbg!(&devices);

    // Load the tracking configuration from the test data.
    let trkconfg_yaml: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "tests",
        "config",
        "tracking_cfg.yaml",
    ]
    .iter()
    .collect();

    let configs: BTreeMap<String, TrkConfig> = TrkConfig::load_named(trkconfg_yaml).unwrap();

    dbg!(&configs);

    // Build the tracking arc simulation to generate a "standard measurement".
    let mut trk = TrackingArcSim::<Spacecraft, RangeDoppler, _>::with_seed(
        devices, trajectory, configs, 12345,
    )
    .unwrap();

    trk.build_schedule(almanac.clone()).unwrap();
    let arc = trk.generate_measurements(almanac).unwrap();

    // And serialize to disk
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
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
    let arc_concrete = dyn_arc.to_tracking_arc::<RangeDoppler>().unwrap();

    println!("{arc_concrete}");

    assert_eq!(arc.measurements.len(), 116);
    // Check that we've loaded all of the measurements
    assert_eq!(arc_concrete.measurements.len(), arc.measurements.len());
    // Check that we find the same device names too
    assert_eq!(arc_concrete.device_names(), arc.device_names());
    // Check that we've copied over the device configurations as well
    assert_eq!(arc_concrete.device_cfg, arc.device_cfg);
}
