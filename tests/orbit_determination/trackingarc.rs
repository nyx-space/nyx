use nyx_space::io::tracking_data::DynamicTrackingArc;
use nyx_space::io::ConfigRepr;
use nyx_space::md::prelude::*;
use nyx_space::od::msr::RangeDoppler;
use nyx_space::od::prelude::*;
use nyx_space::od::simulator::TrackingArcSim;
use nyx_space::od::simulator::{Cadence, Strand, TrkConfig};
use rstest::*;
use std::collections::BTreeMap;
use std::env;
use std::path::PathBuf;
use std::str::FromStr;

#[fixture]
fn traj() -> Traj<Orbit> {
    let _ = pretty_env_logger::try_init();

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
        .for_duration_with_traj(3.days())
        .unwrap();

    println!("{trajectory}");

    trajectory
}

#[fixture]
fn devices() -> Vec<GroundStation> {
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

    GroundStation::load_many(ground_station_file).unwrap()
}

#[rstest]
fn trk_simple(traj: Traj<Orbit>, devices: Vec<GroundStation>) {
    // Load cosm
    let cosm = Cosm::de438();

    // Path to output data
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "output_data",
        "tracking_truth_ephem.parquet",
    ]
    .iter()
    .collect();

    traj.to_parquet_simple(path.clone()).unwrap();

    traj.to_groundtrack_parquet(
        path.with_file_name("tracking_truth_ephem_groundtrack.parquet"),
        cosm.frame("IAU Earth"),
        None,
        None,
        cosm.clone(),
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
        TrackingArcSim::<Orbit, RangeDoppler, _>::with_seed(devices, traj, configs, 12345).unwrap();

    // Test that building the schedule is deterministic
    let orig_sched = trk.generate_schedule(cosm.clone()).unwrap();
    for ii in 0..5 {
        let sched = trk.generate_schedule(cosm.clone()).unwrap();
        assert_eq!(
            sched, orig_sched,
            "{ii} was different:\n orig {orig_sched:?}\n sched {sched:?}"
        );
    }

    trk.build_schedule(cosm.clone()).unwrap();

    let arc = trk.generate_measurements(cosm).unwrap();

    // Test filtering by epoch
    let start_epoch = arc.measurements[0].1.epoch() + 1.minutes();
    for (_, msr) in arc.filter_by_epoch(start_epoch..).measurements {
        assert!(msr.epoch() >= start_epoch);
    }

    for (_, msr) in arc.filter_by_epoch(..=start_epoch).measurements {
        assert!(msr.epoch() <= start_epoch);
    }

    for (_, msr) in arc.filter_by_epoch(..start_epoch).measurements {
        assert!(msr.epoch() < start_epoch);
    }

    assert_eq!(
        arc.filter_by_epoch(start_epoch..start_epoch)
            .measurements
            .len(),
        0
    );

    // Test filtering by duration offset
    for (_, msr) in arc.filter_by_offset(1.minutes()..).measurements {
        assert!(msr.epoch() >= start_epoch);
    }

    for (_, msr) in arc.filter_by_offset(..=1.minutes()).measurements {
        assert!(msr.epoch() <= start_epoch);
    }

    for (_, msr) in arc.filter_by_offset(..1.minutes()).measurements {
        assert!(msr.epoch() < start_epoch);
    }

    assert_eq!(
        arc.filter_by_offset(1.minutes()..1.minutes())
            .measurements
            .len(),
        0
    );

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
    let dyn_arc = DynamicTrackingArc::from_parquet(output_fn).unwrap();
    // And convert to the same tracking arc as earlier
    let arc_concrete = dyn_arc.to_tracking_arc::<RangeDoppler>().unwrap();

    println!("{arc_concrete}");

    // Check that we've loaded all of the measurements
    assert_eq!(arc_concrete.measurements.len(), arc.measurements.len());
    // Check that we find the same device names too
    assert_eq!(arc_concrete.device_names(), arc.device_names());
    // Check that we've copied over the device configurations as well
    assert_eq!(arc_concrete.device_cfg, arc.device_cfg);
}

/// Tests that inclusion epochs work
#[rstest]
fn trkconfig_zero_inclusion(traj: Traj<Orbit>, devices: Vec<GroundStation>) {
    let cosm = Cosm::de438();

    // Build a tracking config that should always see this vehicle.
    let trkcfg_always = TrkConfig::builder()
        .strands(vec![Strand {
            start: traj.first().epoch(),
            end: traj.last().epoch(),
        }])
        .build();

    // Build the configs map, where we only have one of the two stations configured
    let mut configs = BTreeMap::new();
    configs.insert(devices[1].name.clone(), trkcfg_always);

    let mut trk = TrackingArcSim::<Orbit, RangeDoppler, _>::new(devices, traj, configs).unwrap();

    trk.build_schedule(cosm.clone()).unwrap();

    let arc = trk.generate_measurements(cosm).unwrap();

    // Regression
    assert_eq!(arc.measurements.len(), 113);

    assert_eq!(
        arc.device_names().len(),
        1,
        "only one device should have measurements"
    );
}

/// Test invalid tracking configurations
#[rstest]
fn trkconfig_invalid(traj: Traj<Orbit>, devices: Vec<GroundStation>) {
    // Build a tracking config where the exclusion range is less than the sampling rate
    let trkcfg = TrkConfig::builder()
        .strands(vec![Strand {
            start: traj.first().epoch(),
            end: traj.first().epoch(),
        }])
        .build();

    // Build the configs map
    let mut configs = BTreeMap::new();
    for device in &devices {
        configs.insert(device.name.clone(), trkcfg.clone());
    }

    assert!(TrackingArcSim::<Orbit, RangeDoppler, _>::new(devices, traj, configs).is_err());
}

/// Test a delayed start of the configuration
#[rstest]
fn trkconfig_delayed_start(traj: Traj<Orbit>, devices: Vec<GroundStation>) {
    let cosm = Cosm::de438();

    let trkcfg = TrkConfig::builder()
        .strands(vec![Strand {
            start: traj.first().epoch() + 2.hours(),
            end: traj.last().epoch(),
        }])
        .sampling(1.26.minutes())
        .build();

    // Build the configs map with a single ground station
    let mut configs = BTreeMap::new();

    configs.insert(devices[0].name.clone(), trkcfg);

    let mut trk =
        TrackingArcSim::<Orbit, RangeDoppler, _>::new(vec![devices[0].clone()], traj, configs)
            .unwrap();

    trk.build_schedule(cosm.clone()).unwrap();

    let arc = trk.generate_measurements(cosm).unwrap();

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
fn trkconfig_cadence(traj: Traj<Orbit>, devices: Vec<GroundStation>) {
    let cosm = Cosm::de438();

    // Build the configs map with a single ground station
    let mut configs = BTreeMap::new();

    configs.insert(
        devices[0].name.clone(),
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
        devices[1].name.clone(),
        TrkConfig::builder()
            .sampling(26.1.seconds())
            .scheduler(Scheduler::default())
            .build(),
    );

    let mut trk = TrackingArcSim::<Orbit, RangeDoppler, _>::new(devices, traj, configs).unwrap();

    trk.build_schedule(cosm.clone()).unwrap();

    let arc = trk.generate_measurements(cosm).unwrap();

    // Check the sampling of the arc is one minute: we don't have any overlap of availability and the default sampling is one minute.
    assert_eq!(
        arc.min_duration_sep().unwrap(),
        26.1.seconds(),
        "sampling should be the minimum of the two devices"
    );

    // Regression
    assert_eq!(arc.measurements.len(), 216);
}
