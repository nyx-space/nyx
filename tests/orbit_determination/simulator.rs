use nyx_space::dynamics::guidance::LocalFrame;
use nyx_space::io::ConfigRepr;
use nyx_space::md::prelude::*;
use nyx_space::md::trajectory::ExportCfg;
use nyx_space::od::prelude::*;
use nyx_space::od::simulator::TrackingArcSim;
use nyx_space::od::simulator::TrkConfig;
use std::collections::BTreeMap;
use std::collections::HashMap;
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

    let mut devices = BTreeMap::new();
    for gs in GroundStation::load_many(ground_station_file).unwrap() {
        devices.insert(gs.name.clone(), gs);
    }

    devices
}

#[fixture]
fn spacecraft(almanac: Arc<Almanac>) -> Spacecraft {
    // Dummy state
    let orbit = Orbit::try_keplerian_altitude(
        150_000.0,
        1e-2,
        30.0,
        45.0,
        75.0,
        23.4,
        Epoch::from_str("2023-02-22T19:18:17.16 UTC").unwrap(),
        almanac.frame_from_uid(EARTH_J2000).unwrap(),
    )
    .unwrap();

    orbit.into()
}

#[fixture]
fn trajectory(spacecraft: Spacecraft, almanac: Arc<Almanac>) -> Trajectory {
    // Generate a trajectory
    let (_, trajectory) = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()))
        .with(spacecraft, almanac.clone())
        .for_duration_with_traj(0.25 * spacecraft.orbit.period().unwrap())
        .unwrap();

    trajectory
}

#[fixture]
fn tracking_data(
    trajectory: Trajectory,
    devices: BTreeMap<String, GroundStation>,
    almanac: Arc<Almanac>,
) -> TrackingDataArc {
    // Test that continuous tracking
    let _ = pretty_env_logger::try_init();

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
            almanac.clone(),
        )
        .unwrap();

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
        TrackingArcSim::<Spacecraft, GroundStation>::with_seed(devices, trajectory, configs, 12345)
            .unwrap();

    trk.build_schedule(almanac.clone()).unwrap();
    trk.generate_measurements(almanac).unwrap()
}

#[rstest]
fn continuous_tracking_cov_test(tracking_data: TrackingDataArc) {
    let arc = tracking_data;

    let _ = pretty_env_logger::try_init();

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
    let arc_rtn = TrackingDataArc::from_parquet(output_fn).unwrap();

    println!("{arc_rtn}");

    assert_eq!(arc.measurements.len(), 96734);
    // Check that we've loaded all of the measurements
    assert_eq!(arc_rtn.measurements.len(), arc.measurements.len());
    assert_eq!(arc_rtn.unique(), arc.unique());

    // Serialize as TDM
    let path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "output_data", "simple_arc.tdm"]
        .iter()
        .collect();

    let mut aliases = HashMap::new();
    aliases.insert("Demo Ground Station".to_string(), "Fake GS".to_string());

    let output_fn = arc
        .clone()
        .to_tdm_file(
            "MySpacecraft".to_string(),
            Some(aliases.clone()),
            path,
            ExportCfg::default(),
        )
        .unwrap();

    // Read back from TDM
    let arc_tdm = TrackingDataArc::from_tdm(output_fn, None).unwrap();
    println!("{arc_tdm}");

    // Check everything but the source, since it'll be set when read from TDM.
    assert_eq!(arc_tdm.len(), arc.len());
    assert_eq!(arc_tdm.start_epoch(), arc.start_epoch());
    assert_eq!(arc_tdm.end_epoch(), arc.end_epoch());
    assert_eq!(arc_tdm.unique(), arc.unique());
    // Check that we have multiplied back and divided back correctly.
    assert_eq!(arc_tdm.measurements, arc.measurements);

    // Test the downsampling
    let tdm_failed_downsample = arc_tdm.clone().downsample(0.1.seconds());
    assert_eq!(
        tdm_failed_downsample.len(),
        arc_tdm.len(),
        "downsampling should have failed because it's upsampling"
    );

    let arc_downsample = arc_tdm.clone().downsample(10.seconds());
    println!("{arc_downsample}");
    assert_eq!(
        arc_downsample.len(),
        arc_tdm.len() / 10 + 1,
        "downsampling has wrong sample count"
    );

    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "output_data",
        "simple_arc_downsampled.parquet",
    ]
    .iter()
    .collect();

    arc_downsample
        .to_parquet(path, ExportCfg::default())
        .unwrap();
}

#[rstest]
fn od_with_modulus_cov_test(
    spacecraft: Spacecraft,
    tracking_data: TrackingDataArc,
    mut devices: BTreeMap<String, GroundStation>,
    trajectory: Trajectory,
    almanac: Arc<Almanac>,
) {
    let mut arc = tracking_data;

    // Assume JPL DSN Code is used, cf. DSN docs 214, section 2.2.2.2.
    let jpl_dsn_code_length_km = 75660.0;
    arc.set_moduli(MeasurementType::Range, jpl_dsn_code_length_km);
    arc.apply_moduli();

    // Increase the noise on the OD process
    // Set a bias instead of assuming a modulus.
    for device in devices.values_mut() {
        for (_, stochastics) in device.stochastic_noises.as_mut().unwrap().iter_mut() {
            *stochastics *= 2.0;
        }
    }

    let uncertainty = SpacecraftUncertainty::builder()
        .nominal(spacecraft)
        .frame(LocalFrame::RIC)
        .x_km(0.5)
        .y_km(0.5)
        .z_km(0.5)
        .vx_km_s(0.5e-3)
        .vy_km_s(0.5e-3)
        .vz_km_s(0.5e-3)
        .build();

    assert!((uncertainty.x_km - 0.5).abs() < f64::EPSILON);
    assert!((uncertainty.y_km - 0.5).abs() < f64::EPSILON);
    assert!((uncertainty.z_km - 0.5).abs() < f64::EPSILON);
    assert!((uncertainty.vx_km_s - 0.5e-3).abs() < f64::EPSILON);
    assert!((uncertainty.vy_km_s - 0.5e-3).abs() < f64::EPSILON);
    assert!((uncertainty.vz_km_s - 0.5e-3).abs() < f64::EPSILON);

    let estimate = uncertainty.to_estimate().unwrap();

    println!("{estimate}");

    let sigma_q = 1e-19_f64;
    let process_noise = ProcessNoise3D::from_velocity_km_s(
        &[sigma_q, sigma_q, sigma_q],
        20 * Unit::Minute,
        2 * Unit::Minute,
        Some(LocalFrame::RIC),
    );
    println!("{process_noise}");
    let kf = KF::new(estimate, process_noise);

    let setup = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let prop = setup.with(spacecraft.with_stm(), almanac.clone());

    let mut odp = SpacecraftODProcess::ekf(
        prop,
        kf,
        devices,
        EkfTrigger::new(10, Unit::Minute * 15),
        None,
        almanac.clone(),
    );

    let od_sol = odp.process_arc(&arc).unwrap();

    od_sol
        .to_parquet(
            "./output_data/od_with_modulus.parquet",
            ExportCfg::default(),
        )
        .unwrap();

    // Check the final error.
    let estimate = od_sol.estimates.last().unwrap();
    let rss_pos_km = trajectory
        .at(estimate.epoch())
        .unwrap()
        .orbit
        .rss_radius_km(&estimate.orbital_state())
        .unwrap();

    println!("rss_pos_km = {rss_pos_km}");

    let reject_count = od_sol.rejected_residuals().len();

    assert!(reject_count < 10, "wrong number of expected rejections");
}

#[rstest]
fn od_with_modulus_as_bias_cov_test(
    spacecraft: Spacecraft,
    mut tracking_data: TrackingDataArc,
    mut devices: BTreeMap<String, GroundStation>,
    trajectory: Trajectory,
    almanac: Arc<Almanac>,
) {
    // Assume JPL DSN Code is used, cf. DSN docs 214, section 2.2.2.2.
    let jpl_dsn_code_length_km = 75660.0;

    tracking_data.set_moduli(MeasurementType::Range, jpl_dsn_code_length_km);
    tracking_data.apply_moduli();
    // Forget there ever was a modulus!
    tracking_data.moduli = None;

    // Set a bias instead of assuming a modulus.
    for (name, device) in devices.clone() {
        let biased_device = device
            .with_msr_bias_constant(MeasurementType::Range, jpl_dsn_code_length_km)
            .unwrap();
        devices.insert(name, biased_device);
    }

    let uncertainty = SpacecraftUncertainty::builder()
        .nominal(spacecraft)
        .frame(LocalFrame::RIC)
        .x_km(0.5)
        .y_km(0.5)
        .z_km(0.5)
        .vx_km_s(0.5e-3)
        .vy_km_s(0.5e-3)
        .vz_km_s(0.5e-3)
        .build();

    let estimate = uncertainty.to_estimate().unwrap();

    let sigma_q = 1e-8_f64.powi(2);
    let process_noise =
        ProcessNoise3D::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);
    let kf = KF::new(estimate, process_noise);

    let setup = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let prop = setup.with(spacecraft.with_stm(), almanac.clone());

    let mut odp = SpacecraftODProcess::ekf(
        prop,
        kf,
        devices,
        EkfTrigger::new(10, Unit::Minute * 15),
        None,
        almanac,
    );

    let od_sol = odp.process_arc(&tracking_data).unwrap();

    od_sol
        .to_parquet(
            "./output_data/od_with_modulus.parquet",
            ExportCfg::default(),
        )
        .unwrap();

    // Check the final error.
    let estimate = od_sol.estimates.last().unwrap();
    let rss_pos_km = trajectory
        .at(estimate.epoch())
        .unwrap()
        .orbit
        .rss_radius_km(&estimate.orbital_state())
        .unwrap();

    assert!(
        rss_pos_km > 100_000.0,
        "expected bias to not correctly solve OD"
    )
}
