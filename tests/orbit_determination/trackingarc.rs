use anise::constants::frames::{EARTH_J2000, IAU_EARTH_FRAME};
use nyx_space::io::ConfigRepr;
use nyx_space::md::prelude::*;
use nyx_space::od::prelude::*;
use nyx_space::od::simulator::TrackingArcSim;
use nyx_space::od::simulator::{Cadence, Strand, TrkConfig};
use rstest::*;
use std::collections::BTreeMap;
use std::env;
use std::path::PathBuf;
use std::str::FromStr;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[fixture]
fn traj(almanac: Arc<Almanac>) -> Traj<Spacecraft> {
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
        .with(Spacecraft::builder().orbit(orbit).build(), almanac)
        .for_duration_with_traj(3.days())
        .unwrap();

    println!("{trajectory}");

    trajectory
}

#[fixture]
fn devices() -> BTreeMap<String, GroundStation> {
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

    GroundStation::load_named(ground_station_file).unwrap()
}

#[rstest]
fn trk_simple(
    traj: Traj<Spacecraft>,
    devices: BTreeMap<String, GroundStation>,
    almanac: Arc<Almanac>,
) {
    // Path to output data
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "output_data",
        "tracking_truth_ephem.parquet",
    ]
    .iter()
    .collect();

    traj.to_parquet_simple(path.clone(), almanac.clone())
        .unwrap();

    traj.to_groundtrack_parquet(
        path.with_file_name("tracking_truth_ephem_groundtrack.parquet"),
        almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap(),
        None,
        None,
        almanac.clone(),
    )
    .unwrap();

    dbg!(&devices);

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
    let mut trk =
        TrackingArcSim::<Spacecraft, GroundStation>::with_seed(devices, traj, configs, 12345)
            .unwrap();

    // Test that building the schedule is deterministic
    let orig_sched = trk.generate_schedule(almanac.clone()).unwrap();
    for ii in 0..5 {
        let sched = trk.generate_schedule(almanac.clone()).unwrap();
        assert_eq!(
            sched, orig_sched,
            "{ii} was different:\n orig {orig_sched:?}\n sched {sched:?}"
        );
    }

    trk.build_schedule(almanac.clone()).unwrap();

    let arc = trk.generate_measurements(almanac).unwrap();

    // Regression
    assert_eq!(arc.measurements.len(), 197);

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
    let arc_concrete = TrackingDataArc::from_parquet(output_fn).unwrap();

    println!("{arc_concrete}");

    // Check that we've loaded all of the measurements
    assert_eq!(arc_concrete.measurements.len(), arc.measurements.len());
    assert_eq!(arc_concrete.unique(), arc.unique());
}

/// Tests that inclusion epochs work
#[rstest]
fn trkconfig_zero_inclusion(
    traj: Traj<Spacecraft>,
    devices: BTreeMap<String, GroundStation>,
    almanac: Arc<Almanac>,
) {
    // Build a tracking config that should always see this vehicle.
    let trkcfg_always = TrkConfig::builder()
        .strands(vec![Strand {
            start: traj.first().epoch(),
            end: traj.last().epoch(),
        }])
        .build();

    // Build the configs map, where we only have one of the two stations configured
    let mut configs = BTreeMap::new();
    configs.insert("Canberra".to_string(), trkcfg_always);

    let mut trk = TrackingArcSim::<Spacecraft, GroundStation>::new(devices, traj, configs).unwrap();

    trk.build_schedule(almanac.clone()).unwrap();

    let arc = trk.generate_measurements(almanac).unwrap();

    // Regression
    assert_eq!(arc.measurements.len(), 113);

    assert_eq!(
        arc.unique_aliases().len(),
        1,
        "only one device should have measurements"
    );
}

/// Test invalid tracking configurations
#[rstest]
fn trkconfig_invalid(traj: Traj<Spacecraft>, devices: BTreeMap<String, GroundStation>) {
    // Build a tracking config where the exclusion range is less than the sampling rate
    let trkcfg = TrkConfig::builder()
        .strands(vec![Strand {
            start: traj.first().epoch(),
            end: traj.first().epoch(),
        }])
        .build();

    // Build the configs map
    let mut configs = BTreeMap::new();
    for name in devices.keys() {
        configs.insert(name.clone(), trkcfg.clone());
    }

    assert!(TrackingArcSim::<Spacecraft, GroundStation>::new(devices, traj, configs).is_err());
}

/// Test a delayed start of the configuration
#[rstest]
fn trkconfig_delayed_start(
    traj: Traj<Spacecraft>,
    mut devices: BTreeMap<String, GroundStation>,
    almanac: Arc<Almanac>,
) {
    let trkcfg = TrkConfig::builder()
        .strands(vec![Strand {
            start: traj.first().epoch() + 2.hours(),
            end: traj.last().epoch(),
        }])
        .sampling(1.26.minutes())
        .build();

    devices.remove("Canberra").unwrap();

    // Build the configs map with a single ground station
    let mut configs = BTreeMap::new();
    configs.insert("Demo Ground Station".to_string(), trkcfg);

    let mut trk = TrackingArcSim::<Spacecraft, GroundStation>::new(devices, traj, configs).unwrap();

    trk.build_schedule(almanac.clone()).unwrap();

    let arc = trk.generate_measurements(almanac).unwrap();

    // Check the sampling of the arc.
    assert_eq!(
        arc.min_duration_sep().unwrap(),
        1.26.minutes(),
        "sampling invalid"
    );

    // Regression
    assert_eq!(arc.measurements.len(), 108);
}

/// Test different cadences and availabilities
#[rstest]
fn trkconfig_cadence(
    traj: Traj<Spacecraft>,
    devices: BTreeMap<String, GroundStation>,
    almanac: Arc<Almanac>,
) {
    // Build the configs map with a single ground station
    let mut configs = BTreeMap::new();

    configs.insert(
        "Demo Ground Station".to_string(),
        TrkConfig::builder()
            .scheduler(
                Scheduler::builder()
                    .cadence(Cadence::Intermittent {
                        on: 0.2.hours(),
                        off: 20.days(),
                    })
                    .build(),
            )
            .build(),
    );

    configs.insert(
        "Canberra".to_string(),
        TrkConfig::builder()
            .sampling(26.1.seconds())
            .scheduler(Scheduler::default())
            .build(),
    );

    let mut trk = TrackingArcSim::<Spacecraft, GroundStation>::new(devices, traj, configs).unwrap();

    trk.build_schedule(almanac.clone()).unwrap();

    let arc = trk.generate_measurements(almanac).unwrap();

    // Check the sampling of the arc is one minute: we don't have any overlap of availability and the default sampling is one minute.
    assert_eq!(
        arc.min_duration_sep().unwrap(),
        26.1.seconds(),
        "sampling should be the minimum of the two devices"
    );

    // Regression
    assert_eq!(arc.measurements.len(), 215);
}
