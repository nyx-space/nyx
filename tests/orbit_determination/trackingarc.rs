use nyx_space::io::tracking_data::DynamicTrackingArc;
use nyx_space::io::ConfigRepr;
use nyx_space::md::ui::*;
use nyx_space::od::msr::StdMeasurement;
use nyx_space::od::prelude::*;
use nyx_space::od::simulator::arc::TrackingArcSim;
use nyx_space::od::simulator::{Availability, EpochRanges, Schedule, TrkConfig};
use rstest::*;
use std::collections::HashMap;
use std::env;
use std::path::PathBuf;
use std::str::FromStr;

#[fixture]
fn traj() -> Traj<Orbit> {
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

    trajectory
}

#[fixture]
fn devices() -> Vec<GroundStation> {
    // Load the ground stations from the test data.
    let ground_station_file: PathBuf = [
        &env::var("CARGO_MANIFEST_DIR").unwrap(),
        "data",
        "tests",
        "config",
        "many_ground_stations.yaml",
    ]
    .iter()
    .collect();

    let devices = GroundStation::load_many(ground_station_file).unwrap();

    devices
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
        cosm.clone(),
    )
    .unwrap();

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

    let configs: HashMap<String, TrkConfig> = TrkConfig::load_named(trkconfg_yaml).unwrap();

    dbg!(&configs);

    // Build the tracking arc simulation to generate a "standard measurement".
    let mut trk = TrackingArcSim::<Orbit, StdMeasurement, _>::with_seed(
        devices,
        traj.clone(),
        configs,
        12345,
    )
    .unwrap();

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

/// Tests that exclusion epochs work
#[rstest]
fn trkconfig_zero_exclusion(traj: Traj<Orbit>, devices: Vec<GroundStation>) {
    let cosm = Cosm::de438();

    // Build a tracking config that should never see this vehicle.
    let trkcfg = TrkConfig {
        exclusion_epochs: Some(vec![EpochRanges {
            start: traj.first().epoch(),
            end: traj.last().epoch(),
        }]),
        ..Default::default()
    };
    // Build the configs map
    let mut configs = HashMap::new();
    for device in &devices {
        configs.insert(device.name.clone(), trkcfg.clone());
    }

    let mut trk = TrackingArcSim::<Orbit, StdMeasurement, _>::new(devices, traj, configs).unwrap();

    let arc = trk.generate_measurements(cosm).unwrap();

    assert_eq!(arc.measurements.len(), 0);
}

/// Tests that inclusion epochs work
#[rstest]
fn trkconfig_zero_inclusion(traj: Traj<Orbit>, devices: Vec<GroundStation>) {
    let cosm = Cosm::de438();

    // Build a tracking config that should always see this vehicle.
    let trkcfg_always = TrkConfig {
        inclusion_epochs: Some(vec![EpochRanges {
            start: traj.first().epoch(),
            end: traj.last().epoch(),
        }]),
        ..Default::default()
    };

    // And one that is never included
    let trkcfg_never = TrkConfig {
        exclusion_epochs: Some(vec![EpochRanges {
            start: traj.first().epoch(),
            end: traj.last().epoch(),
        }]),
        ..Default::default()
    };
    // Build the configs map
    let mut configs = HashMap::new();
    for (dno, device) in devices.iter().enumerate() {
        configs.insert(
            device.name.clone(),
            if dno == 0 {
                println!("{}", device.name);
                trkcfg_never.clone()
            } else {
                trkcfg_always.clone()
            },
        );
    }

    let mut trk = TrackingArcSim::<Orbit, StdMeasurement, _>::new(devices, traj, configs).unwrap();

    let arc = trk.generate_measurements(cosm).unwrap();

    // Regression
    assert_eq!(arc.measurements.len(), 79);

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
    let trkcfg = TrkConfig {
        exclusion_epochs: Some(vec![EpochRanges {
            start: traj.first().epoch(),
            end: traj.first().epoch(),
        }]),
        ..Default::default()
    };
    // Build the configs map
    let mut configs = HashMap::new();
    for device in &devices {
        configs.insert(device.name.clone(), trkcfg.clone());
    }

    assert!(TrackingArcSim::<Orbit, StdMeasurement, _>::new(devices, traj, configs).is_err());
}

/// Test a delayed start of the configuration
#[rstest]
fn trkconfig_delayed_start(traj: Traj<Orbit>, devices: Vec<GroundStation>) {
    let cosm = Cosm::de438();

    let trkcfg = TrkConfig {
        inclusion_epochs: Some(vec![EpochRanges {
            start: traj.first().epoch() + 2.hours(),
            end: traj.last().epoch(),
        }]),
        sampling: 1.26.minutes(),
        ..Default::default()
    };

    // Build the configs map with a single ground station
    let mut configs = HashMap::new();

    configs.insert(devices[0].name.clone(), trkcfg.clone());

    // Check that if if a device does not have an associated trkconfig, the tracking arc cannot be created.
    assert!(TrackingArcSim::<Orbit, StdMeasurement, _>::new(
        devices.clone(),
        traj.clone(),
        configs.clone()
    )
    .is_err());

    let mut trk =
        TrackingArcSim::<Orbit, StdMeasurement, _>::new(vec![devices[0].clone()], traj, configs)
            .unwrap();

    let arc = trk.generate_measurements(cosm).unwrap();

    // Check the sampling of the arc.
    assert_eq!(
        arc.min_duration_sep().unwrap(),
        1.26.minutes(),
        "sampling invalid"
    );

    // Regression
    assert_eq!(arc.measurements.len(), 53);
}

/// Test different cadences and availabilities
#[rstest]
fn trkconfig_cadence(traj: Traj<Orbit>, devices: Vec<GroundStation>) {
    let cosm = Cosm::de438();

    // Build the configs map with a single ground station
    let mut configs = HashMap::new();

    configs.insert(
        devices[0].name.clone(),
        TrkConfig {
            start: Availability::Visible,
            schedule: Schedule::Intermittent {
                on: 0.2.hours(),
                off: 20.days(),
            },
            ..Default::default()
        },
    );

    configs.insert(
        devices[1].name.clone(),
        TrkConfig {
            start: Availability::Epoch(traj.last().epoch() - 10.hours()),
            sampling: 26.1.seconds(),
            ..Default::default()
        },
    );

    let mut trk = TrackingArcSim::<Orbit, StdMeasurement, _>::new(devices, traj, configs).unwrap();

    let arc = trk.generate_measurements(cosm).unwrap();

    // Check the sampling of the arc is one minute: we don't have any overlap of availability and the default sampling is one minute.
    assert_eq!(
        arc.min_duration_sep().unwrap(),
        26.1.seconds(),
        "sampling should be the minimum of the two devices"
    );

    // Regression
    assert_eq!(arc.measurements.len(), 90);
}
