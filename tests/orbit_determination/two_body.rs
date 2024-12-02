extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use anise::constants::frames::IAU_EARTH_FRAME;
use nalgebra::U2;
use nyx::cosmic::Orbit;
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::sph_harmonics::Harmonics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::io::ConfigRepr;
use nyx::io::{gravity::*, ExportCfg};
use nyx::linalg::{SMatrix, SVector};
use nyx::od::prelude::*;
use nyx::propagators::{IntegratorOptions, Propagator};
use nyx::Spacecraft;
use nyx_space::propagators::IntegratorMethod;
use std::collections::BTreeMap;
use std::env;
use std::path::PathBuf;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[fixture]
fn sim_devices(almanac: Arc<Almanac>) -> BTreeMap<String, GroundStation> {
    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();
    let elevation_mask = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        StochasticNoise::ZERO,
        StochasticNoise::ZERO,
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        StochasticNoise::ZERO,
        StochasticNoise::ZERO,
        iau_earth,
    );
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        StochasticNoise::ZERO,
        StochasticNoise::ZERO,
        iau_earth,
    );

    let mut devices = BTreeMap::new();
    devices.insert("Madrid".to_string(), dss65_madrid);
    devices.insert("Canberra".to_string(), dss34_canberra);
    devices.insert("Goldstone".to_string(), dss13_goldstone);

    devices
}

/// Devices for processing the measurement, noise may not be zero.
#[fixture]
fn proc_devices(almanac: Arc<Almanac>) -> BTreeMap<String, GroundStation> {
    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();
    let elevation_mask = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        StochasticNoise::MIN,
        StochasticNoise::MIN,
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        StochasticNoise::MIN,
        StochasticNoise::MIN,
        iau_earth,
    );
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        StochasticNoise::MIN,
        StochasticNoise::MIN,
        iau_earth,
    );

    let mut devices = BTreeMap::new();
    devices.insert("Madrid".to_string(), dss65_madrid);
    devices.insert("Canberra".to_string(), dss34_canberra);
    devices.insert("Goldstone".to_string(), dss13_goldstone);

    devices
}

#[allow(clippy::identity_op)]
#[rstest]
fn od_tb_val_ekf_fixed_step_perfect_stations(
    almanac: Arc<Almanac>,
    sim_devices: BTreeMap<String, GroundStation>,
    proc_devices: BTreeMap<String, GroundStation>,
) {
    let _ = pretty_env_logger::try_init();

    // Define the ground stations.
    let ekf_num_meas = 100;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 5.0 * Unit::Second;

    // Load the tracking configurations
    let mut configs = BTreeMap::new();
    let trkconfig_yaml: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "tests",
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
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // We're sharing both the propagator and the dynamics.
    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);

    let mut prop = setup.with(initial_state.into(), almanac.clone());
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();
    println!("{}", final_truth);

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs.clone(), 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());
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
    println!("initial estimate:\n{}", initial_estimate);

    let kf = KF::no_snc(initial_estimate);

    let mut odp = ODProcess::<_, U2, _, _, _>::ekf(
        prop_est,
        kf,
        proc_devices,
        EkfTrigger::new(ekf_num_meas, ekf_disable_time),
        None,
        almanac,
    );

    odp.process_arc(&arc).unwrap();

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("Final estimate:\n{est}");
    assert!(
        est.state_deviation().norm() < 1e-12,
        "In perfect modeling, the state deviation should be near zero, got {:.3e}",
        est.state_deviation().norm()
    );
    for i in 0..6 {
        assert!(
            est.covar[(i, i)] >= 0.0,
            "covar diagonal element negative @ [{}, {}]",
            i,
            i
        );
    }
    for i in 0..6 {
        if i < 3 {
            assert!(
                est.covar[(i, i)] < covar_radius_km,
                "covar radius did not decrease"
            );
        } else {
            assert!(
                est.covar[(i, i)] < covar_velocity_km_s,
                "covar velocity did not decrease"
            );
        }
    }

    // Check the final estimate
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("{}\n\n{}\n{}", est.state_deviation(), est, final_truth);
    let delta = (est.state().orbit - final_truth.orbit).unwrap();
    println!(
        "RMAG error = {:.3} m\tVMAG error = {:.3} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 2e-16, "Position error should be zero");
    assert!(delta.vmag_km_s() < 2e-16, "Velocity error should be zero");
}

#[allow(clippy::identity_op)]
#[rstest]
fn od_tb_val_with_arc(
    almanac: Arc<Almanac>,
    sim_devices: BTreeMap<String, GroundStation>,
    proc_devices: BTreeMap<String, GroundStation>,
) {
    let _ = pretty_env_logger::try_init();

    // Define the ground stations.
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 5.0 * Unit::Second;

    let all_stations = sim_devices;

    let ekf_num_meas = 100;

    // Define the propagator information.
    let duration = 1 * Unit::Day;

    // Define the storages (channels for the states and a map for the measurements).

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // We're sharing both the propagator and the dynamics.
    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::new(
        orbital_dyn,
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step_s(10.0),
    );

    let mut prop = setup.with(initial_state.into(), almanac.clone());
    let (final_truth, traj) = prop.for_duration_with_traj(duration).unwrap();
    println!("{}", final_truth);

    // Save the trajectory to parquet
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "output_data",
        "od_val_with_arc_truth_ephem.parquet",
    ]
    .iter()
    .collect();
    traj.to_parquet_simple(path, almanac.clone()).unwrap();

    // Load the tracking configs
    let trkconfig_yaml: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "tests",
        "config",
        "trk_cfg_od_val_arc.yaml",
    ]
    .iter()
    .collect();

    let configs: BTreeMap<String, TrkConfig> = TrkConfig::load_named(trkconfig_yaml).unwrap();

    // Simulate tracking data of range and range rate
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs.clone(), 1).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // And serialize to disk
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "output_data",
        "two_body_od_val_arc.parquet",
    ]
    .iter()
    .collect();

    arc.to_parquet_simple(path).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());
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
    println!("initial estimate:\n{}", initial_estimate);

    let kf = KF::no_snc(initial_estimate);

    let mut odp = ODProcess::<_, U2, _, _, _>::ekf(
        prop_est,
        kf,
        proc_devices,
        EkfTrigger::new(ekf_num_meas, ekf_disable_time),
        Some(ResidRejectCrit::default()),
        almanac,
    );

    odp.process_arc(&arc).unwrap();

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("Final estimate:\n{est}");
    assert!(
        est.state_deviation().norm() < f64::EPSILON,
        "In perfect modeling, the state deviation should be near zero, got {:.3e}",
        est.state_deviation().norm()
    );
    for i in 0..6 {
        assert!(
            est.covar[(i, i)] >= 0.0,
            "covar diagonal element negative @ [{}, {}]",
            i,
            i
        );
    }
    for i in 0..6 {
        if i < 3 {
            assert!(
                est.covar[(i, i)] < covar_radius_km,
                "covar radius did not decrease"
            );
        } else {
            assert!(
                est.covar[(i, i)] < covar_velocity_km_s,
                "covar velocity did not decrease"
            );
        }
    }

    // Check the final estimate
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("{}\n\n{}\n{}", est.state_deviation(), est, final_truth);
    let delta = (est.state().orbit - final_truth.orbit).unwrap();
    println!(
        "RMAG error = {:.3} m\tVMAG error = {:.3} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 2e-16, "Position error should be zero");
    assert!(delta.vmag_km_s() < 2e-16, "Velocity error should be zero");
}

#[fixture]
fn cfg() -> TrkConfig {
    // Define the tracking configurations
    TrkConfig::builder()
        .sampling(10.seconds())
        .scheduler(Scheduler::builder().sample_alignment(10.seconds()).build())
        .build()
}

#[allow(clippy::identity_op)]
#[rstest]
fn od_tb_val_ckf_fixed_step_perfect_stations(
    almanac: Arc<Almanac>,
    sim_devices: BTreeMap<String, GroundStation>,
    proc_devices: BTreeMap<String, GroundStation>,
    cfg: TrkConfig,
) {
    /*
     * This tests that the state transition matrix computation is correct with two body dynamics.
     *
     * Specifically, the same dynamics are used for both the measurement generation and for the estimation.
     * However, only the estimation generation propagates the STM. When STM propagation is enabled, the code will compute
     * the dynamics using a hyperdual representation in 7 dimensions: 1 for the reals, 3 for the position partials,
     * 3 for the velocity partials.
     *
     * Hence, if the filter state estimation is any different from the truth data, then it means that the equations of
     * motion computed in hyperdual space differ from the ones computes in the reals.
     *
     * Thereby, this serves as a validation of the orbital dynamics implementation.
     **/
    let _ = pretty_env_logger::try_init();

    let mut configs = BTreeMap::new();
    for name in sim_devices.keys() {
        configs.insert(name.clone(), cfg.clone());
    }

    let all_stations = sim_devices;

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // Generate the truth data on one thread.
    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::two_body());

    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);
    let mut prop = setup.with(initial_state.into(), almanac.clone());

    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs.clone(), 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // And serialize to disk
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "output_data",
        "od_tb_val_ckf_fixed_step_perfect_stations.parquet",
    ]
    .iter()
    .collect();

    arc.to_parquet_simple(path).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let initial_state_est = Spacecraft::from(initial_state).with_stm();
    // Use the same setup as earlier
    let prop_est = setup.with(initial_state_est, almanac.clone());
    let covar_radius_km = 1.0e-3;
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

    // Define the initial orbit estimate
    let initial_estimate = KfEstimate::from_covar(initial_state_est, init_covar);

    let ckf = KF::no_snc(initial_estimate);

    let mut odp = ODProcess::<_, U2, _, _, _>::ckf(prop_est, ckf, proc_devices, None, almanac);

    odp.process_arc(&arc).unwrap();

    let path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "output_data", "tb_ckf.parquet"]
        .iter()
        .collect();

    odp.to_parquet(&arc, path, ExportCfg::default()).unwrap();

    // Check that there are no duplicates of epochs.
    let mut prev_epoch = odp.estimates[0].epoch();

    for est in odp.estimates.iter().skip(2) {
        let this_epoch = est.epoch();
        assert!(
            this_epoch > prev_epoch,
            "Estimates not continuously going forward: {this_epoch} <= {prev_epoch}"
        );
        prev_epoch = this_epoch;
    }

    for (no, est) in odp.estimates.iter().enumerate() {
        if no == 0 {
            // Skip the first estimate which is the initial estimate provided by user
            continue;
        }
        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0,
                "covar diagonal element negative @ [{}, {}] = {:e} @ {}",
                i,
                i,
                est.covar[(i, i)],
                est.epoch()
            );
        }
        assert!(
            est.state_deviation().norm() < 1e-12,
            "estimate error should be zero (perfect dynamics) ({:e})",
            est.state_deviation().norm()
        );
    }

    for res in odp.residuals.iter().flatten() {
        assert!(
            res.prefit.norm() < 1e-12,
            "prefit should be zero (perfect dynamics) ({:e})",
            res
        );
    }

    for res in odp.residuals.iter().flatten() {
        assert!(
            res.postfit.norm() < 1e-12,
            "postfit should be zero (perfect dynamics) ({:e})",
            res
        );
    }

    let est = odp.estimates.last().unwrap();
    println!("estimate error {:.2e}", est.state_deviation().norm());
    println!("estimate covariance {:.2e}", est.covar.diagonal().norm());

    assert!(
        est.state_deviation().norm() < 1e-12,
        "estimate error should be zero (perfect dynamics) ({:e})",
        est.state_deviation().norm()
    );

    assert!(
        est.covar.diagonal().norm() < 1e-6,
        "estimate covariance norm should be zero (perfect dynamics) ({:e})",
        est.covar.diagonal().norm()
    );

    let delta = (est.state().orbit - final_truth.orbit).unwrap();
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 2e-16, "Position error should be zero");
    assert!(delta.vmag_km_s() < 2e-16, "Velocity error should be zero");

    // Iterate
    odp.iterate_arc(
        &arc,
        IterationConf {
            smoother: SmoothingArc::TimeGap(10.0 * Unit::Second),
            ..Default::default()
        },
    )
    .unwrap();

    println!(
        "N-1 one iteration: \n{}",
        odp.estimates[odp.estimates.len() - 1]
    );

    println!(
        "Initial state after iteration: \n{:x}",
        odp.estimates[0].state()
    );

    // Check the final estimate
    let est = odp.estimates.last().unwrap();
    println!("estimate error {:.2e}", est.state_deviation().norm());
    println!("estimate covariance {:.2e}", est.covar.diagonal().norm());

    assert!(
        est.state_deviation().norm() < 1e-12,
        "estimate error should be zero (perfect dynamics) ({:e})",
        est.state_deviation().norm()
    );

    // Note we accept a larger covariance diagonal here because smoothing will increase the covariance
    assert!(
        est.covar.diagonal().norm() < 1e-4,
        "estimate covariance norm should be zero (perfect dynamics) ({:e})",
        est.covar.diagonal().norm()
    );

    let delta = (est.state().orbit - final_truth.orbit).unwrap();
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 1e-9, "More than 1 micrometer error");
    assert!(delta.vmag_km_s() < 1e-9, "More than 1 micrometer/s error");
}

#[allow(clippy::identity_op)]
#[rstest]
fn od_tb_val_az_el_ckf_fixed_step_perfect_stations(
    almanac: Arc<Almanac>,
    mut sim_devices: BTreeMap<String, GroundStation>,
    mut proc_devices: BTreeMap<String, GroundStation>,
    cfg: TrkConfig,
) {
    /*
     * This tests that the state transition matrix computation is correct with two body dynamics.
     *
     * Specifically, the same dynamics are used for both the measurement generation and for the estimation.
     * However, only the estimation generation propagates the STM. When STM propagation is enabled, the code will compute
     * the dynamics using a hyperdual representation in 7 dimensions: 1 for the reals, 3 for the position partials,
     * 3 for the velocity partials.
     *
     * Hence, if the filter state estimation is any different from the truth data, then it means that the equations of
     * motion computed in hyperdual space differ from the ones computes in the reals.
     *
     * Thereby, this serves as a validation of the orbital dynamics implementation.
     **/
    let _ = pretty_env_logger::try_init();

    for (_, dev) in sim_devices.iter_mut() {
        dev.measurement_types.swap_remove(&MeasurementType::Range);
        dev.measurement_types.swap_remove(&MeasurementType::Doppler);
        dev.measurement_types.insert(MeasurementType::Azimuth);
        dev.measurement_types.insert(MeasurementType::Elevation);
        dev.stochastic_noises
            .as_mut()
            .unwrap()
            .insert(MeasurementType::Azimuth, StochasticNoise::ZERO);
        dev.stochastic_noises
            .as_mut()
            .unwrap()
            .insert(MeasurementType::Elevation, StochasticNoise::ZERO);
    }

    for (_, dev) in proc_devices.iter_mut() {
        dev.measurement_types.swap_remove(&MeasurementType::Range);
        dev.measurement_types.swap_remove(&MeasurementType::Doppler);
        dev.measurement_types.insert(MeasurementType::Azimuth);
        dev.measurement_types.insert(MeasurementType::Elevation);
        dev.stochastic_noises
            .as_mut()
            .unwrap()
            .insert(MeasurementType::Azimuth, StochasticNoise::MIN);
        dev.stochastic_noises
            .as_mut()
            .unwrap()
            .insert(MeasurementType::Elevation, StochasticNoise::MIN);
    }

    let mut configs = BTreeMap::new();
    for name in sim_devices.keys() {
        configs.insert(name.clone(), cfg.clone());
    }

    let all_stations = sim_devices;

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // Generate the truth data on one thread.
    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::two_body());

    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);
    let mut prop = setup.with(initial_state.into(), almanac.clone());

    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs.clone(), 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // And serialize to disk
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "output_data",
        "od_tb_val_az_el_ckf_fixed_step_perfect_stations.parquet",
    ]
    .iter()
    .collect();

    arc.to_parquet_simple(path).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let initial_state_est = Spacecraft::from(initial_state).with_stm();
    // Use the same setup as earlier
    let prop_est = setup.with(initial_state_est, almanac.clone());
    let covar_radius_km = 1.0e-3;
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

    // Define the initial orbit estimate
    let initial_estimate = KfEstimate::from_covar(initial_state_est, init_covar);

    let ckf = KF::no_snc(initial_estimate);

    let mut odp = ODProcess::<_, U2, _, _, _>::ckf(prop_est, ckf, proc_devices, None, almanac);

    odp.process_arc(&arc).unwrap();

    let path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "output_data", "tb_ckf.parquet"]
        .iter()
        .collect();

    odp.to_parquet(&arc, path, ExportCfg::default()).unwrap();

    // Check that there are no duplicates of epochs.
    let mut prev_epoch = odp.estimates[0].epoch();

    for est in odp.estimates.iter().skip(2) {
        let this_epoch = est.epoch();
        assert!(
            this_epoch > prev_epoch,
            "Estimates not continuously going forward: {this_epoch} <= {prev_epoch}"
        );
        prev_epoch = this_epoch;
    }

    for (no, est) in odp.estimates.iter().enumerate() {
        if no == 0 {
            // Skip the first estimate which is the initial estimate provided by user
            continue;
        }
        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0,
                "covar diagonal element negative @ [{}, {}] = {:e} @ {}",
                i,
                i,
                est.covar[(i, i)],
                est.epoch()
            );
        }
        assert!(
            est.state_deviation().norm() < 1e-12,
            "estimate error should be zero (perfect dynamics) ({:e})",
            est.state_deviation().norm()
        );
    }

    for res in odp.residuals.iter().flatten() {
        assert!(
            res.prefit.norm() < 1e-12,
            "prefit should be zero (perfect dynamics) ({:e})",
            res
        );
    }

    for res in odp.residuals.iter().flatten() {
        assert!(
            res.postfit.norm() < 1e-12,
            "postfit should be zero (perfect dynamics) ({:e})",
            res
        );
    }

    let est = odp.estimates.last().unwrap();
    println!("estimate error {:.2e}", est.state_deviation().norm());
    println!("estimate covariance {:.2e}", est.covar.diagonal().norm());

    assert!(
        est.state_deviation().norm() < 1e-12,
        "estimate error should be zero (perfect dynamics) ({:e})",
        est.state_deviation().norm()
    );

    assert!(
        est.covar.diagonal().norm() < 1e-6,
        "estimate covariance norm should be zero (perfect dynamics) ({:e})",
        est.covar.diagonal().norm()
    );

    let delta = (est.state().orbit - final_truth.orbit).unwrap();
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 2e-16, "Position error should be zero");
    assert!(delta.vmag_km_s() < 2e-16, "Velocity error should be zero");

    // Iterate
    odp.iterate_arc(
        &arc,
        IterationConf {
            smoother: SmoothingArc::TimeGap(10.0 * Unit::Second),
            ..Default::default()
        },
    )
    .unwrap();

    println!(
        "N-1 one iteration: \n{}",
        odp.estimates[odp.estimates.len() - 1]
    );

    println!(
        "Initial state after iteration: \n{:x}",
        odp.estimates[0].state()
    );

    // Check the final estimate
    let est = odp.estimates.last().unwrap();
    println!("estimate error {:.2e}", est.state_deviation().norm());
    println!("estimate covariance {:.2e}", est.covar.diagonal().norm());

    assert!(
        est.state_deviation().norm() < 1e-12,
        "estimate error should be zero (perfect dynamics) ({:e})",
        est.state_deviation().norm()
    );

    // Note we accept a larger covariance diagonal here because smoothing will increase the covariance
    assert!(
        est.covar.diagonal().norm() < 1e-4,
        "estimate covariance norm should be zero (perfect dynamics) ({:e})",
        est.covar.diagonal().norm()
    );

    let delta = (est.state().orbit - final_truth.orbit).unwrap();
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 1e-9, "More than 1 micrometer error");
    assert!(delta.vmag_km_s() < 1e-9, "More than 1 micrometer/s error");
}

#[allow(clippy::identity_op)]
#[rstest]
fn od_tb_ckf_fixed_step_iteration_test(
    almanac: Arc<Almanac>,
    sim_devices: BTreeMap<String, GroundStation>,
    proc_devices: BTreeMap<String, GroundStation>,
) {
    let _ = pretty_env_logger::try_init();

    // Define the ground stations.
    let range_noise = 0.1; // in km (so 100 meters of error)
    let range_rate_noise = 0.001; // in km/s (or 1 meter per second of error)

    // Define the tracking configurations
    let cfg = TrkConfig::from_sample_rate(10.seconds());

    let mut configs = BTreeMap::new();
    for name in sim_devices.keys() {
        configs.insert(name.clone(), cfg.clone());
    }

    let all_stations = sim_devices;

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);

    let mut prop = setup.with(initial_state.into(), almanac.clone());
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs.clone(), 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());
    let covar_radius_km = 1.0e-3;
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
    // Define the initial estimate (x_hat): add 100 meters in X, remove 100 meters in Y and add 50 meters in Z
    let mut initial_state2 = initial_state;
    initial_state2.radius_km.x += 0.1;
    initial_state2.radius_km.y -= 0.1;
    initial_state2.radius_km.z += 0.05;
    let initial_estimate = KfEstimate::from_covar(initial_state2.into(), init_covar);

    let ckf = KF::no_snc(initial_estimate);

    let mut odp = ODProcess::<_, U2, _, _, _>::ckf(prop_est, ckf, proc_devices, None, almanac);

    odp.process_arc(&arc).unwrap();

    // Check the final estimate prior to iteration
    let delta = (odp.estimates.last().unwrap().state().orbit - final_truth.orbit).unwrap();
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(
        delta.rmag_km() < range_noise,
        "More than station level position error"
    );
    assert!(
        delta.vmag_km_s() < range_rate_noise,
        "More than station level velocity error"
    );

    // Iterate, and check that the initial state difference is lower
    odp.iterate_arc(
        &arc,
        IterationConf {
            smoother: SmoothingArc::TimeGap(10.0 * Unit::Second),
            ..Default::default()
        },
    )
    .unwrap();

    let dstate_no_iteration = (initial_state - initial_state2).unwrap();
    let dstate_iteration = (initial_state - odp.estimates[0].state().orbit).unwrap();

    println!("{}\n{}", initial_state2, odp.estimates[0].state());

    // Compute the order of magnitude of the errors, and check that iteration either decreases it or keeps it the same
    let err_it_oom = dstate_iteration.rmag_km().log10().floor() as i32;
    let err_no_it_oom = dstate_no_iteration.rmag_km().log10().floor() as i32;

    println!(
        "Difference in initial states radii without iterations: {} km (order of magnitude: {})",
        dstate_no_iteration.rmag_km(),
        err_no_it_oom
    );
    println!(
        "Difference in initial states radii with iterations: {} km (order of magnitude: {})",
        dstate_iteration.rmag_km(),
        err_it_oom
    );
    assert!(
        dstate_iteration.rmag_km() < dstate_no_iteration.rmag_km() || err_it_oom <= err_no_it_oom,
        "Iteration did not reduce initial error"
    );

    // Check the final estimate
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("{}\n\n{}\n{}", est.state_deviation(), est, final_truth);
    let delta = (est.state().orbit - final_truth.orbit).unwrap();
    println!(
        "RMAG error = {:.3} m\tVMAG error = {:.3} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 75e-3, "More than 75 meter error");
    assert!(delta.vmag_km_s() < 50e-6, "More than 50 mm/s error");
}

#[allow(clippy::identity_op)]
#[rstest]
fn od_tb_ckf_fixed_step_perfect_stations_snc_covar_map(
    almanac: Arc<Almanac>,
    sim_devices: BTreeMap<String, GroundStation>,
    proc_devices: BTreeMap<String, GroundStation>,
) {
    // Tests state noise compensation with covariance mapping
    let _ = pretty_env_logger::try_init();

    // Define the tracking configurations
    let mut configs = BTreeMap::new();
    let cfg = TrkConfig::from_sample_rate(10.seconds());
    for name in sim_devices.keys() {
        configs.insert(name.clone(), cfg.clone());
    }

    let all_stations = sim_devices;

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);

    let mut prop = setup.with(initial_state.into(), almanac.clone());
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs.clone(), 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());

    // Set up the filter
    let covar_radius_km = 1.0e-3;
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

    // Define the process noise to assume an unmodeled acceleration of 1e-3 km^2/s^2 on X, Y and Z in the ECI frame
    let sigma_q = 1e-8_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);

    let ckf = KF::new(initial_estimate, process_noise);

    let mut odp = ODProcess::<_, U2, _, _, _>::ckf(prop_est, ckf, proc_devices, None, almanac);

    odp.process_arc(&arc).unwrap();

    // Let's check that the covariance never falls below our sigma squared values
    for (no, est) in odp.estimates.iter().enumerate() {
        if no == 1 {
            println!("{}", est);
        }
        assert!(
            est.state_deviation().norm() < 1e-12,
            "estimate error should be zero (perfect dynamics) ({:e})",
            est.state_deviation().norm()
        );

        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0,
                "covar diagonal element negative @ [{}, {}]",
                i,
                i
            );
        }

        if est.predicted() {
            for i in 0..6 {
                assert!(
                    est.covar[(i, i)] >= sigma_q,
                    "covar diagonal less than SNC value @ {} = {:.3e}",
                    no,
                    est.covar[(i, i)]
                );
            }
        }
    }

    // Check the final estimate
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("{}\n\n{}\n{}", est.state_deviation(), est, final_truth);
    let delta = (est.state().orbit - final_truth.orbit).unwrap();
    println!(
        "RMAG error = {:.3} m\tVMAG error = {:.3} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 1e-3, "More than 1 meter error");
    assert!(delta.vmag_km_s() < 1e-6, "More than 1 mm/s error");
}

#[allow(clippy::identity_op)]
#[rstest]
fn od_tb_ckf_map_covar(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    // Define the propagator information.
    let duration = 2 * Unit::Day;
    let step_size = 10.0 * Unit::Second;

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.

    let setup = Propagator::new(
        SpacecraftDynamics::new(OrbitalDynamics::two_body()),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(step_size),
    );
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());
    let covar_radius_km = 1.0e-3;
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

    let initial_estimate = KfEstimate::from_covar(Spacecraft::from(initial_state), init_covar);

    let ckf = KF::no_snc(initial_estimate);

    let mut odp: SpacecraftODProcess =
        ODProcess::<_, U2, _, _, _>::ckf(prop_est, ckf, BTreeMap::new(), None, almanac);

    odp.predict_for(30.seconds(), duration).unwrap();

    // Check that the covariance inflated (we don't get the norm of the estimate because it's zero without any truth data)
    let estimates = odp.estimates;
    let est = &estimates[estimates.len() - 1];
    for i in 0..6 {
        assert!(
            est.covar[(i, i)] >= 0.0,
            "covar diagonal element negative @ [{}, {}]",
            i,
            i
        );
    }
    for i in 0..6 {
        if i < 3 {
            assert!(
                est.covar[(i, i)] > covar_radius_km,
                "covar radius did not increase"
            );
        } else {
            assert!(
                est.covar[(i, i)] > covar_velocity_km_s,
                "covar velocity did not increase"
            );
        }
    }
}

#[allow(clippy::identity_op)]
#[rstest]
fn od_tb_val_harmonics_ckf_fixed_step_perfect_cov_test(
    almanac: Arc<Almanac>,
    sim_devices: BTreeMap<String, GroundStation>,
    proc_devices: BTreeMap<String, GroundStation>,
) {
    // Tests state noise compensation with covariance mapping
    let _ = pretty_env_logger::try_init();

    // Define the tracking configurations
    let mut configs = BTreeMap::new();
    let cfg = TrkConfig::from_sample_rate(10.seconds());
    for name in sim_devices.keys() {
        configs.insert(name.clone(), cfg.clone());
    }

    let all_stations = sim_devices;

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let earth_sph_harm = HarmonicsMem::from_cof("data/JGM3.cof.gz", 70, 70, true).unwrap();
    let harmonics = Harmonics::from_stor(iau_earth, earth_sph_harm);
    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::from_model(harmonics));
    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);

    let mut prop = setup.with(initial_state.into(), almanac.clone());
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs.clone(), 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());

    // Set up the filter
    let covar_radius_km = 1.0e-3;
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

    let ckf = KF::no_snc(initial_estimate);

    let mut odp = ODProcess::<_, U2, _, _, _>::ckf(prop_est, ckf, proc_devices, None, almanac);

    odp.process_arc(&arc).unwrap();

    // Let's check that the covariance never falls below our sigma squared values
    for (no, est) in odp.estimates.iter().enumerate() {
        if no == 1 {
            println!("{}", est);
        }
        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0,
                "covar diagonal element negative @ [{}, {}]",
                i,
                i
            );
        }
        assert!(
            est.state_deviation().norm() < 1e-2,
            "estimate error should be good (perfect dynamics) ({:e})",
            est.state_deviation().norm()
        );
    }

    // Check the final estimate
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("{}\n\n{}\n{}", est.state_deviation(), est, final_truth);
    let delta = (est.state().orbit - final_truth.orbit).unwrap();
    println!(
        "RMAG error = {:.3} m\tVMAG error = {:.3} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 2e-16, "Position error should be zero");
    assert!(delta.vmag_km_s() < 2e-16, "Velocity error should be zero");
}

#[allow(clippy::identity_op)]
#[rstest]
fn od_tb_ckf_fixed_step_perfect_stations_several_snc_covar_map(
    almanac: Arc<Almanac>,
    sim_devices: BTreeMap<String, GroundStation>,
    proc_devices: BTreeMap<String, GroundStation>,
) {
    // Tests state noise compensation with covariance mapping
    let _ = pretty_env_logger::try_init();

    // Define the tracking configurations
    let mut configs = BTreeMap::new();
    let cfg = TrkConfig::from_sample_rate(10.seconds());
    for name in sim_devices.keys() {
        configs.insert(name.clone(), cfg.clone());
    }

    let all_stations = sim_devices;

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);
    let mut prop = setup.with(initial_state.into(), almanac.clone());
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs.clone(), 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());

    // Set up the filter
    let covar_radius_km = 1.0e-3;
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

    // Define the process noise to assume an unmodeled acceleration of 1e-3 km^2/s^2 on X, Y and Z in the ECI frame
    let sigma_q1 = 1e-7_f64.powi(2);
    let process_noise1 = SNC3::from_diagonal(2 * Unit::Day, &[sigma_q1, sigma_q1, sigma_q1]);

    let sigma_q2 = 1e-8_f64.powi(2);
    let sigma_q2_d = 3600.0;
    let mut process_noise2 = SNC3::with_decay(
        2 * Unit::Day,
        &[sigma_q2, sigma_q2, sigma_q2],
        &[sigma_q2_d, sigma_q2_d, sigma_q2_d],
    );
    process_noise2.start_time = Some(dt + 36_000.0); // Start the second process noise 10 hours into the tracking pass

    let ckf = KF::with_sncs(initial_estimate, vec![process_noise1, process_noise2]);

    let mut odp = ODProcess::<_, U2, _, _, _>::ckf(prop_est, ckf, proc_devices, None, almanac);

    odp.process_arc(&arc).unwrap();

    // Let's check that the covariance never falls below our sigma squared values
    for (no, est) in odp.estimates.iter().enumerate() {
        if no == 1 {
            println!("{}", est);
        }
        assert!(
            est.state_deviation().norm() < 1e-6,
            "estimate error should be zero (perfect dynamics) ({:e})",
            est.state_deviation().norm()
        );

        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0,
                "covar diagonal element negative @ [{}, {}]",
                i,
                i
            );
        }
    }

    // Check the final estimate
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("{}\n\n{}\n{}", est.state_deviation(), est, final_truth);
    let delta = (est.state().orbit - final_truth.orbit).unwrap();
    println!(
        "RMAG error = {:.3} m\tVMAG error = {:.3} m/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e3
    );

    assert!(delta.rmag_km() < 2e-16, "Position error should be zero");
    assert!(delta.vmag_km_s() < 2e-16, "Velocity error should be zero");
}
