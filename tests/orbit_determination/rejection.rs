extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use anise::constants::frames::EARTH_J2000;
use nyx::cosmic::Orbit;
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::io::ConfigRepr;
use nyx::linalg::{SMatrix, SVector};
use nyx::od::prelude::*;
use nyx::propagators::{IntegratorOptions, Propagator};
use nyx_space::propagators::IntegratorMethod;
use std::collections::BTreeMap;
use std::env;
use std::path::PathBuf;

use anise::prelude::Almanac;
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[fixture]
fn sim_devices() -> BTreeMap<String, GroundStation> {
    let elevation_mask = 0.0;
    let dss65_madrid =
        GroundStation::dss65_madrid(elevation_mask, StochasticNoise::ZERO, StochasticNoise::ZERO);
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, StochasticNoise::ZERO, StochasticNoise::ZERO);
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        StochasticNoise::ZERO,
        StochasticNoise::ZERO,
    );

    let mut devices = BTreeMap::new();
    devices.insert("Madrid".to_string(), dss65_madrid);
    devices.insert("Canberra".to_string(), dss34_canberra);
    devices.insert("Goldstone".to_string(), dss13_goldstone);

    devices
}

/// Devices for processing the measurement, noise may not be zero.
#[fixture]
fn proc_devices() -> BTreeMap<String, GroundStation> {
    let elevation_mask = 0.0;
    let dss65_madrid =
        GroundStation::dss65_madrid(elevation_mask, StochasticNoise::MIN, StochasticNoise::MIN);
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, StochasticNoise::MIN, StochasticNoise::MIN);
    let dss13_goldstone =
        GroundStation::dss13_goldstone(elevation_mask, StochasticNoise::MIN, StochasticNoise::MIN);

    let mut devices = BTreeMap::new();
    devices.insert("Madrid".to_string(), dss65_madrid);
    devices.insert("Canberra".to_string(), dss34_canberra);
    devices.insert("Goldstone".to_string(), dss13_goldstone);

    devices
}

#[rstest]
fn od_rejection_test(
    almanac: Arc<Almanac>,
    sim_devices: BTreeMap<String, GroundStation>,
    proc_devices: BTreeMap<String, GroundStation>,
) {
    let _ = pretty_env_logger::try_init();

    // Load the tracking configurations
    let mut configs = BTreeMap::new();
    let trkconfig_yaml: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "03_tests",
        "config",
        "trk_cfg_od_val.yaml",
    ]
    .iter()
    .collect();

    let cfg = TrkConfig::load(trkconfig_yaml).unwrap();

    for name in sim_devices.keys() {
        configs.insert(name.clone(), cfg.clone());
    }

    let all_stations = sim_devices;

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // We're sharing both the propagator and the dynamics.
    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);

    let mut prop = setup.with(initial_state.into(), almanac.clone());
    let (_, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs.clone(), 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc_full = arc_sim.generate_measurements(almanac.clone()).unwrap();
    let initial_count = arc_full.len();
    println!("Full arc has {} measurements", initial_count);

    // Reject Madrid measurements
    let arc_rejected = arc_full.clone().reject_by_tracker("Madrid");
    // The number of measurements in the arc structure remains the same, but they are marked as rejected.
    assert_eq!(arc_rejected.len(), initial_count);

    let covar_radius_km = 1.0e-6;
    let covar_velocity_km_s = 1.0e-6;
    let init_covar = SMatrix::<f64, 9, 9>::from_diagonal(&SVector::<f64, 9>::from_iterator([
        covar_radius_km,
        covar_radius_km,
        covar_radius_km,
        covar_velocity_km_s,
        covar_velocity_km_s,
        covar_velocity_km_s,
        0.0,
        0.0,
        0.0,
    ]));

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(initial_state.into(), init_covar);

    let odp = SpacecraftKalmanOD::new(
        setup,
        KalmanVariant::ReferenceUpdate,
        None,
        proc_devices,
        almanac,
    );

    let od_sol_full = odp.process_arc(initial_estimate, &arc_full).unwrap();
    let accepted_full = od_sol_full.accepted_residuals().len();
    println!("Full arc accepted residuals: {}", accepted_full);

    let od_sol_rejected = odp.process_arc(initial_estimate, &arc_rejected).unwrap();
    let accepted_rejected = od_sol_rejected.accepted_residuals().len();
    println!("Rejected arc accepted residuals: {}", accepted_rejected);

    assert!(
        accepted_rejected < accepted_full,
        "Expected fewer accepted residuals after rejection"
    );

    // Check that we indeed rejected all Madrid measurements
    // We can check by iterating over residuals and checking if tracker "Madrid" is present in accepted residuals
    let madrid_accepted = od_sol_rejected
        .accepted_residuals()
        .iter()
        .filter(|r| r.tracker.as_ref().unwrap() == "Madrid")
        .count();
    assert_eq!(
        madrid_accepted, 0,
        "Expected 0 Madrid measurements to be accepted"
    );
}
