extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MOON, SATURN_BARYCENTER, SUN};
use anise::constants::frames::IAU_EARTH_FRAME;
use nyx::cosmic::Orbit;
use nyx::dynamics::orbital::{OrbitalDynamics, PointMasses};
use nyx::dynamics::sph_harmonics::Harmonics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::io::gravity::*;
use nyx::linalg::{SMatrix, SVector};
use nyx::od::prelude::*;
use nyx::propagators::{IntegratorOptions, Propagator};
use nyx::utils::rss_orbit_errors;
use nyx::Spacecraft;
use nyx_space::io::ExportCfg;
use nyx_space::propagators::IntegratorMethod;
use std::collections::BTreeMap;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

/*
 * These tests check that if we start with a state deviation in the estimate, the filter will eventually converge back.
 * These tests do NOT check that the filter will converge if the initial state in the propagator has that state deviation.
 * The latter would require iteration and smoothing before playing with an EKF. This will be handled in a subsequent version.
**/

#[fixture]
fn devices(almanac: Arc<Almanac>) -> BTreeMap<String, GroundStation> {
    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();
    let elevation_mask = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        StochasticNoise::default_range_km(),
        StochasticNoise::default_doppler_km_s(),
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        StochasticNoise::default_range_km(),
        StochasticNoise::default_doppler_km_s(),
        iau_earth,
    );

    let mut devices = BTreeMap::new();
    devices.insert("Madrid".to_string(), dss65_madrid);
    devices.insert("Canberra".to_string(), dss34_canberra);
    devices
}

#[allow(clippy::identity_op)]
#[rstest]
fn xhat_dev_test_ekf_two_body(almanac: Arc<Almanac>, devices: BTreeMap<String, GroundStation>) {
    let _ = pretty_env_logger::try_init();

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    // Define the tracking configurations
    let mut configs = BTreeMap::new();
    configs.insert(
        "Madrid".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        "Canberra".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    // Define the propagator information.
    let prop_time = 0.01 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.radius_km.x += 5.0;
    initial_state_dev.radius_km.y -= 5.0;
    initial_state_dev.radius_km.z += 5.0;

    let (err_p, err_v) = rss_orbit_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\nDelta: {}",
        err_p * 1e3,
        err_v * 1e3,
        (initial_state - initial_state_dev).unwrap()
    );

    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);
    let (_, traj) = setup
        .with(Spacecraft::from(initial_state), almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Simulate tracking data
    println!("{traj}");
    let mut arc_sim = TrackingArcSim::with_seed(devices.clone(), traj.clone(), configs, 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(
        Spacecraft::from(initial_state_dev).with_stm(),
        almanac.clone(),
    );
    let covar_radius_km = 1.0e2;
    let covar_velocity_km_s = 1.0e1;
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
    let initial_estimate = KfEstimate::from_covar(initial_state_dev.into(), init_covar);
    println!("Initial estimate:\n{}", initial_estimate);

    let sigma_q = 1e-7_f64.powi(2);
    let process_noise =
        ProcessNoise3D::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);
    let kf = KF::new(initial_estimate, process_noise);

    let mut odp = SpacecraftODProcess::ckf(prop_est, kf, devices, None, almanac.clone());

    let od_sol = odp.process_arc(&arc).expect("OD process failed");
    let pre_smooth_first_est = od_sol.estimates[0];
    let pre_smooth_num_est = od_sol.estimates.len();

    let smoothed_od_sol = od_sol.clone().smooth(almanac).expect("OD smoothing failed");

    assert_eq!(
        pre_smooth_num_est,
        smoothed_od_sol.estimates.len(),
        "different number of estimates smoothed and not"
    );

    // Check the new initial estimate is better than at the start
    let smoothed_init_state = smoothed_od_sol.estimates[0].state().orbit;
    let (sm_err_p, sm_err_v) = rss_orbit_errors(&smoothed_init_state, &initial_state);
    println!(
        "New initial state dev: {:.3} m\t{:.3} m/s\n{}",
        sm_err_p * 1e3,
        sm_err_v * 1e3,
        (smoothed_init_state - initial_state_dev).unwrap()
    );

    assert!(
        sm_err_p < err_p,
        "initial position not improved by smoothing"
    );
    // We don't check the velocity because the initial error is zero, so the smoother will change the velocity for a better fit.

    // Check that the covariance deflated
    let est = &smoothed_od_sol.estimates.last().unwrap();
    println!("Estimate:\n{}", est);
    let final_truth_state = traj.at(est.epoch()).unwrap();
    println!("Truth:\n{}", final_truth_state);
    let (err_p, err_v) = rss_orbit_errors(&est.state().orbit, &final_truth_state.orbit);
    println!(
        "Delta state with truth (epoch match: {}): {:.3} m\t{:.3} m/s\n{}",
        final_truth_state.epoch() == est.epoch(),
        err_p * 1e3,
        err_v * 1e3,
        (final_truth_state.orbit - est.state().orbit).unwrap()
    );

    for i in 0..6 {
        if est.covar[(i, i)] < 0.0 {
            println!(
                "covar diagonal element negative @ [{}, {}] = {:.3e} -- issue #164",
                i,
                i,
                est.covar[(i, i)],
            );
        }
    }

    assert_eq!(
        final_truth_state.epoch(),
        est.epoch(),
        "time of final EST and TRUTH epochs differ"
    );
    let rmag_err = (final_truth_state.orbit - est.state().orbit)
        .unwrap()
        .rmag_km();
    assert!(
        rmag_err < sm_err_p,
        "final radius error ({:.3} m) should be better than initial state error ({:.3} m)",
        rmag_err * 1e3,
        sm_err_p * 1e3
    );

    let vmag_err = (final_truth_state.orbit - est.state().orbit)
        .unwrap()
        .vmag_km_s();
    assert!(
        vmag_err < sm_err_v,
        "final velocity error ({:.3} m) should be better than initial state error ({:.3} m)",
        vmag_err * 1e3,
        sm_err_v * 1e3
    );

    let post_smooth_first_est = smoothed_od_sol.estimates[0];

    let init_pos_rss = initial_state.rss_radius_km(&initial_state_dev).unwrap();
    let init_vel_rss = initial_state.rss_velocity_km_s(&initial_state_dev).unwrap();
    let zero_it_pos_rss = initial_state
        .rss_radius_km(&pre_smooth_first_est.state().orbit)
        .unwrap();
    let zero_it_vel_rss = initial_state
        .rss_velocity_km_s(&pre_smooth_first_est.state().orbit)
        .unwrap();

    let one_it_pos_rss = initial_state
        .rss_radius_km(&post_smooth_first_est.state().orbit)
        .unwrap();
    let one_it_vel_rss = initial_state
        .rss_velocity_km_s(&post_smooth_first_est.state().orbit)
        .unwrap();

    println!(
        "[pos] init: {}\tzero: {}\t one: {}",
        init_pos_rss, zero_it_pos_rss, one_it_pos_rss,
    );
    println!(
        "[vel] init: {}\tzero: {}\t one: {}",
        init_vel_rss, zero_it_vel_rss, one_it_vel_rss
    );

    println!(
        "RMS before smoothing: {}\t{}\t{}",
        od_sol.rms_prefit_residuals(),
        od_sol.rms_postfit_residuals(),
        od_sol.rms_residual_ratios()
    );
    println!(
        "RMS after smoothing: {}\t{}\t{}",
        smoothed_od_sol.rms_prefit_residuals(),
        smoothed_od_sol.rms_postfit_residuals(),
        smoothed_od_sol.rms_residual_ratios()
    );

    assert!(
        smoothed_od_sol.rms_postfit_residuals() < od_sol.rms_postfit_residuals(),
        "smoothing does not improve the residuals"
    );

    od_sol
        .to_parquet("od_sol_xhat.parquet", ExportCfg::default())
        .expect("could not export smoothed solutions");
    smoothed_od_sol
        .to_parquet("od_smooth_xhat.parquet", ExportCfg::default())
        .expect("could not export smoothed solutions");
}

#[allow(clippy::identity_op)]
#[rstest]
fn xhat_dev_test_ekf_multi_body(almanac: Arc<Almanac>, devices: BTreeMap<String, GroundStation>) {
    // We seed both propagators with the same initial state, but we let a large state deviation in the filter.
    // This does _not_ impact the prefits, but it impacts the state deviation and therefore the state estimate.
    // As such, it checks that the filter can return to a nominal state.
    let _ = pretty_env_logger::try_init();

    // Define the ground stations.
    let ekf_num_meas = 500;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 10.0 * Unit::Second;

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    // Define the tracking configurations
    let mut configs = BTreeMap::new();
    configs.insert(
        "Madrid".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        "Canberra".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.radius_km.x += 9.5;
    initial_state_dev.radius_km.y -= 9.5;
    initial_state_dev.radius_km.z += 9.5;

    let (err_p, err_v) = rss_orbit_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        (initial_state - initial_state_dev).unwrap()
    );

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);

    let (_, traj) = setup
        .with(initial_state.into(), almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(devices.clone(), traj.clone(), configs, 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());
    let covar_radius_km = 1.0e2;
    let covar_velocity_km_s = 1.0e1;
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
    let initial_estimate = KfEstimate::from_covar(initial_state_dev.into(), init_covar);
    println!("Initial estimate:\n{}", initial_estimate);

    let sigma_q = 1e-8_f64.powi(2);
    let process_noise =
        ProcessNoise3D::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);
    let kf = KF::new(initial_estimate, process_noise);

    let mut trig = EkfTrigger::new(ekf_num_meas, ekf_disable_time);
    trig.within_sigma = 3.0;

    let mut odp = SpacecraftODProcess::ekf(prop_est, kf, devices, trig, None, almanac);

    let od_sol = odp.process_arc(&arc).unwrap();
    // odp.iterate_arc(&arc, IterationConf::try_from(SmoothingArc::All).unwrap())
    //     .unwrap();

    // Check that the covariance deflated
    let est = &od_sol.estimates.last().unwrap();
    let final_truth_state = traj.at(est.epoch()).unwrap();

    println!("Estimate:\n{}", est);
    println!("Truth:\n{}", final_truth_state);

    // Some sanity checks to make sure that we have correctly indexed the estimates
    assert_eq!(est.epoch(), final_truth_state.epoch());

    let (err_p, err_v) = rss_orbit_errors(&est.state().orbit, &final_truth_state.orbit);

    // Some printing for debugging
    println!(
        "RSS error: estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        (final_truth_state.orbit - est.state().orbit).unwrap()
    );

    for i in 0..6 {
        if est.covar[(i, i)] < 0.0 {
            println!(
                "covar diagonal element negative @ [{}, {}] = {:.3e} -- issue #164",
                i,
                i,
                est.covar[(i, i)],
            );
        }
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

    assert_eq!(
        final_truth_state.epoch(),
        est.epoch(),
        "time of final EST and TRUTH epochs differ"
    );
    let rmag_err = (final_truth_state.orbit - est.state().orbit)
        .unwrap()
        .rmag_km();
    assert!(
        rmag_err < 0.1,
        "final radius error should be on 100 meter level (is instead {:.3} m)",
        rmag_err * 1e3
    );
}

#[allow(clippy::identity_op)]
#[rstest]
fn xhat_dev_test_ekf_harmonics(almanac: Arc<Almanac>, devices: BTreeMap<String, GroundStation>) {
    let _ = pretty_env_logger::try_init();

    // Define the ground stations.
    let ekf_num_meas = 5000;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 1 * Unit::Minute;

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    // Define the tracking configurations
    let mut configs = BTreeMap::new();
    configs.insert(
        "Madrid".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        "Canberra".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.radius_km.x += 9.5;
    initial_state_dev.radius_km.y -= 9.5;
    initial_state_dev.radius_km.z += 9.5;

    let (err_p, err_v) = rss_orbit_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        (initial_state - initial_state_dev).unwrap()
    );

    let hh_deg = 20;
    let hh_ord = 20;

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let earth_sph_harm = HarmonicsMem::from_cof("data/JGM3.cof.gz", hh_deg, hh_ord, true).unwrap();
    let harmonics = Harmonics::from_stor(iau_earth, earth_sph_harm);
    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::new(vec![
        harmonics,
        PointMasses::new(bodies),
    ]));

    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);

    let (_, traj) = setup
        .with(initial_state.into(), almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(devices.clone(), traj.clone(), configs, 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());
    let covar_radius_km = 1.0e2;
    let covar_velocity_km_s = 1.0e1;
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
    let initial_estimate = KfEstimate::from_covar(initial_state_dev.into(), init_covar);
    println!("Initial estimate:\n{}", initial_estimate);

    let sigma_q = 1e-7_f64.powi(2);
    let process_noise =
        ProcessNoise3D::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);
    let kf = KF::new(initial_estimate, process_noise);

    let mut trig = EkfTrigger::new(ekf_num_meas, ekf_disable_time);
    trig.within_sigma = 3.0;

    let mut odp = SpacecraftODProcess::ekf(prop_est, kf, devices, trig, None, almanac);

    let od_sol = odp.process_arc(&arc).unwrap();

    // Check that the covariance deflated
    let est = &od_sol.estimates.last().unwrap();
    let final_truth_state = traj.at(est.epoch()).unwrap();

    println!("Estimate:\n{}", est);
    println!("Truth:\n{}", final_truth_state);
    let (err_p, err_v) = rss_orbit_errors(&est.state().orbit, &final_truth_state.orbit);
    println!(
        "Delta state with truth (epoch match: {}): {:.3} m\t{:.3} m/s\n{}",
        final_truth_state.epoch() == est.epoch(),
        err_p * 1e3,
        err_v * 1e3,
        (final_truth_state.orbit - est.state().orbit).unwrap()
    );

    for i in 0..6 {
        if est.covar[(i, i)] < 0.0 {
            println!(
                "covar diagonal element negative @ [{}, {}] = {:.3e} -- issue #164",
                i,
                i,
                est.covar[(i, i)],
            );
        }
    }

    assert!(est.within_3sigma(), "Final estimate is not within 3 sigma!");

    assert_eq!(
        final_truth_state.epoch(),
        est.epoch(),
        "time of final EST and TRUTH epochs differ"
    );
    let rmag_err = (final_truth_state.orbit - est.state().orbit)
        .unwrap()
        .rmag_km();
    // XXX: Revisit this test
    assert!(
        rmag_err < 5e-1,
        "final radius error too large {:.3} m",
        rmag_err * 1e3
    );
}

#[allow(clippy::identity_op)]
#[rstest]
fn xhat_dev_test_ekf_realistic(almanac: Arc<Almanac>, devices: BTreeMap<String, GroundStation>) {
    let _ = pretty_env_logger::try_init();

    // Define the ground stations.
    let ekf_num_meas = 500;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 10.0 * Unit::Second;

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    // Define the tracking configurations
    let mut configs = BTreeMap::new();
    configs.insert(
        "Madrid".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        "Canberra".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.radius_km.x += 9.5;
    initial_state_dev.radius_km.y -= 9.5;
    initial_state_dev.radius_km.z += 9.5;

    println!(
        "Initial state dev:\n{}",
        (initial_state - initial_state_dev).unwrap()
    );

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER, SATURN_BARYCENTER];
    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let truth_setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);

    let (_, traj) = truth_setup
        .with(initial_state.into(), almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(devices.clone(), traj.clone(), configs, 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be _nearly_ perfect because we've removed SATURN_BARYCENTER from the estimated trajectory
    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let estimator = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let setup = Propagator::new(estimator, IntegratorMethod::RungeKutta4, opts);
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());
    let covar_radius_km = 1.0e2;
    let covar_velocity_km_s = 1.0e1;
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
    let initial_estimate = KfEstimate::from_covar(initial_state_dev.into(), init_covar);
    println!("Initial estimate:\n{}", initial_estimate);

    let kf = KF::no_snc(initial_estimate);

    let mut trig = EkfTrigger::new(ekf_num_meas, ekf_disable_time);
    trig.within_sigma = 3.0;

    let mut odp = SpacecraftODProcess::ekf(prop_est, kf, devices, trig, None, almanac);

    let od_sol = odp.process_arc(&arc).unwrap();

    // Check that the covariance deflated
    let est = &od_sol.estimates.last().unwrap();
    let final_truth_state = traj.at(est.epoch()).unwrap();

    println!("Estimate:\n{}", est);
    println!("Truth:\n{}", final_truth_state);
    println!(
        "Delta state with truth (epoch match: {}):\n{}",
        final_truth_state.epoch() == est.epoch(),
        (final_truth_state.orbit - est.state().orbit).unwrap()
    );

    for i in 0..6 {
        if est.covar[(i, i)] < 0.0 {
            println!(
                "covar diagonal element negative @ [{}, {}] = {:.3e}-- issue #164",
                i,
                i,
                est.covar[(i, i)],
            );
        }
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

    assert_eq!(
        final_truth_state.epoch(),
        est.epoch(),
        "time of final EST and TRUTH epochs differ"
    );
    let rmag_err = (final_truth_state.orbit - est.state().orbit)
        .unwrap()
        .rmag_km();
    assert!(
        rmag_err < 5e-1,
        "final radius error should be less than 500 m (is instead {:.3} m)",
        rmag_err * 1e3
    );
}

#[allow(clippy::identity_op)]
#[rstest]
fn xhat_dev_test_ckf_smoother_multi_body(
    almanac: Arc<Almanac>,
    devices: BTreeMap<String, GroundStation>,
) {
    let _ = pretty_env_logger::try_init();

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    // Define the tracking configurations
    let mut configs = BTreeMap::new();
    configs.insert(
        "Madrid".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        "Canberra".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.radius_km.x += 9.5;
    initial_state_dev.radius_km.y -= 9.5;
    initial_state_dev.radius_km.z += 9.5;

    let (err_p, err_v) = rss_orbit_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        (initial_state - initial_state_dev).unwrap()
    );

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);
    let (_, traj) = setup
        .with(initial_state.into(), almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(devices.clone(), traj.clone(), configs, 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());
    let covar_radius_km = 1.0e2;
    let covar_velocity_km_s = 1.0e1;
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
    let initial_estimate = KfEstimate::from_covar(initial_state_dev.into(), init_covar);
    println!("Initial estimate:\n{}", initial_estimate);

    let kf = KF::no_snc(initial_estimate);

    let mut odp = SpacecraftODProcess::ckf(prop_est, kf, devices, None, almanac.clone());

    let od_sol = odp.process_arc(&arc).unwrap();

    // Smoother
    let smoothed = od_sol.clone().smooth(almanac).unwrap();

    // Check that the estimates and smoothed estimates have the same epochs
    for (i, sm_est) in smoothed.estimates.iter().enumerate() {
        let this_epoch = sm_est.epoch();
        let est_epoch = od_sol.estimates[i].epoch();
        assert_eq!(
            this_epoch, est_epoch,
            "Smoothed estimate epoch different from ODP estimate epoch: {} != {}",
            this_epoch, est_epoch
        );
    }
    assert_eq!(
        od_sol.estimates.len(),
        smoothed.estimates.len(),
        "Different number of estimates and smoothed estimates"
    );

    // Check the first estimate, which should be better thanks to smoothing
    // Only the print the final estimate
    let est = &od_sol.estimates[0];
    let truth_state = traj.at(est.epoch()).unwrap().orbit;
    let smoothed_est = &smoothed.estimates[0];
    let (err_p, err_v) = rss_orbit_errors(&est.state().orbit, &truth_state);
    let (err_p_sm, err_v_sm) = rss_orbit_errors(&smoothed_est.state().orbit, &truth_state);

    println!("=== FIRST ===\nEstimate:\n{}", est);
    println!("Smoother estimate:\n{}", smoothed_est);
    println!("Truth:\n{}", truth_state);

    println!(
        "RSS error: estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        (truth_state - est.state().orbit).unwrap()
    );

    println!(
        "RSS error: smoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
        err_p_sm * 1e3,
        err_v_sm * 1e3,
        (truth_state - smoothed_est.state().orbit).unwrap()
    );

    let mut rss_pos_avr = 0.0;
    let mut rss_vel_avr = 0.0;
    let mut rss_pos_avr_sm = 0.0;
    let mut rss_vel_avr_sm = 0.0;
    let mut num_pos_ok = 0;
    let mut num_vel_ok = 0;

    // Test smoothed estimates
    for (i, est) in od_sol.estimates.iter().enumerate() {
        let smoothed_est = &smoothed.estimates[i];

        // Check that the covariance deflated
        let truth_state = traj.at(est.epoch()).unwrap().orbit;

        // Some sanity checks to make sure that we have correctly indexed the estimates
        assert_eq!(smoothed_est.epoch(), est.epoch());
        assert_eq!(est.epoch(), truth_state.epoch);

        let (err_p, err_v) = rss_orbit_errors(&est.state().orbit, &truth_state);
        let (err_p_sm, err_v_sm) = rss_orbit_errors(&smoothed_est.state().orbit, &truth_state);

        rss_pos_avr += err_p;
        rss_vel_avr += err_v;
        rss_pos_avr_sm += err_p_sm;
        rss_vel_avr_sm += err_v_sm;

        if err_p_sm <= err_p {
            num_pos_ok += 1;
        }

        if err_v_sm <= err_v {
            num_vel_ok += 1;
        }

        if i <= 1 {
            // Only the print the final estimate
            println!("Estimate:\n{}", est);
            println!("Smoother estimate:\n{}", smoothed_est);
            println!("Truth:\n{}", truth_state);

            println!(
                "RSS error: estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                err_p * 1e3,
                err_v * 1e3,
                (truth_state - est.state().orbit).unwrap()
            );

            println!(
                "RSS error: smoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                err_p_sm * 1e3,
                err_v_sm * 1e3,
                (truth_state - smoothed_est.state().orbit).unwrap()
            );
        }

        // The smoothed RSS errors should be better, or have the same order of magnitude or not significantly worse

        // Compute orders of magnitude
        let err_p_oom = err_p.log10().floor() as i32;
        let err_v_oom = err_v.log10().floor() as i32;
        let err_p_sm_oom = err_p_sm.log10().floor() as i32;
        let err_v_sm_oom = err_v_sm.log10().floor() as i32;

        if err_p_sm_oom - err_p_oom > 2 {
            println!(
                "[!!! POS !!!]RSS position error after smoothing not better @{} (#{}):\n\testimate vs truth: {:.3e} m\t{:.3e} m/s\n\tsmoothed estimate vs truth: {:.3e} m\t{:.3e} m/s",
                truth_state.epoch,
                i,
                err_p * 1e3,
                err_v * 1e3,
                err_p_sm * 1e3,
                err_v_sm * 1e3,
            );
        }

        if err_v_sm_oom - err_v_oom > 2 {
            println!(
                "[!!! VEL !!!] RSS velocity error after smoothing not better @{} (#{}):\n\testimate vs truth: {:.3e} m\t{:.3e} m/s\n{}\n\tsmoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                truth_state.epoch,
                i,
                err_p * 1e3,
                err_v * 1e3,
                (truth_state - est.state().orbit).unwrap(),
                err_p_sm * 1e3,
                err_v_sm * 1e3,
                (truth_state - smoothed_est.state().orbit).unwrap()
            );
        }

        for i in 0..6 {
            if est.covar[(i, i)] < 0.0 {
                println!(
                    "covar diagonal element negative @ [{}, {}] = {:.3e}: @{} (#{}) -- issue #164",
                    i,
                    i,
                    est.covar[(i, i)],
                    truth_state.epoch,
                    i,
                );
            }
        }
    }

    let cntf = od_sol.estimates.len() as f64;
    println!(
        "\nPos. better: {}/{}\tVel. better: {}/{}\nPre-smoothing  avr. RSS:\t{:.3e}\t{:.3e}\nPost-smoothing avr. RSS:\t{:.3e}\t{:.3e}\n",
        num_pos_ok,
        od_sol.estimates.len(),
        num_vel_ok,
        od_sol.estimates.len(),
        rss_pos_avr / cntf,
        rss_vel_avr / cntf,
        rss_pos_avr_sm / cntf,
        rss_vel_avr_sm / cntf,
    );

    // For the CKF, the average RSS errors are expected to be better.
    assert!(
        rss_pos_avr_sm / rss_pos_avr <= 2.0,
        "Average RSS position error worse by at least a factor of 2"
    );
    assert!(
        rss_vel_avr_sm / rss_vel_avr <= 2.0,
        "Average RSS velocity error worse by at least a factor of 2"
    );
}

#[allow(clippy::identity_op)]
#[rstest]
fn xhat_dev_test_ekf_snc_smoother_multi_body(
    almanac: Arc<Almanac>,
    devices: BTreeMap<String, GroundStation>,
) {
    let _ = pretty_env_logger::try_init();

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    // Define the tracking configurations
    let mut configs = BTreeMap::new();
    configs.insert(
        "Madrid".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        "Canberra".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.radius_km.x += 9.5;
    initial_state_dev.radius_km.y -= 9.5;
    initial_state_dev.radius_km.z += 9.5;

    let (err_p, err_v) = rss_orbit_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        (initial_state - initial_state_dev).unwrap()
    );

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);
    let (_, traj) = setup
        .with(initial_state.into(), almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(devices.clone(), traj.clone(), configs, 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());
    let covar_radius_km = 1.0e2;
    let covar_velocity_km_s = 1.0e1;
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

    // Define the ground stations.
    let ekf_num_meas = 100;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 1 * Unit::Hour;

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(initial_state_dev.into(), init_covar);
    println!("Initial estimate:\n{}", initial_estimate);

    let sigma_q = 1e-8_f64.powi(2);
    let process_noise =
        ProcessNoise3D::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);
    let kf = KF::new(initial_estimate, process_noise);

    let mut odp = SpacecraftODProcess::ekf(
        prop_est,
        kf,
        devices,
        EkfTrigger::new(ekf_num_meas, ekf_disable_time),
        None,
        almanac.clone(),
    );

    let od_sol = odp.process_arc(&arc).unwrap();

    // Smoother
    let smoothed = od_sol.clone().smooth(almanac.clone()).unwrap();

    let mut rss_pos_avr = 0.0;
    let mut rss_vel_avr = 0.0;
    let mut rss_pos_avr_sm = 0.0;
    let mut rss_vel_avr_sm = 0.0;
    let mut num_pos_ok = 0;
    let mut num_vel_ok = 0;

    // Test smoothed estimates
    // Skip the first 10 estimates which are surprisingly good in this case
    for offset in (1..od_sol.estimates.len() - 10).rev() {
        let smoothed_est = &smoothed.estimates[od_sol.estimates.len() - offset];

        // Check that the covariance deflated
        let est = &od_sol.estimates[od_sol.estimates.len() - offset];
        let truth_state = traj.at(est.epoch()).unwrap().orbit;

        // Some sanity checks to make sure that we have correctly indexed the estimates
        assert_eq!(smoothed_est.epoch(), est.epoch());
        assert_eq!(est.epoch(), truth_state.epoch);

        let (err_p, err_v) = rss_orbit_errors(&est.state().orbit, &truth_state);
        let (err_p_sm, err_v_sm) = rss_orbit_errors(&smoothed_est.state().orbit, &truth_state);

        rss_pos_avr += err_p;
        rss_vel_avr += err_v;
        rss_pos_avr_sm += err_p_sm;
        rss_vel_avr_sm += err_v_sm;

        if err_p_sm <= err_p {
            num_pos_ok += 1;
        }

        if err_v_sm <= err_v {
            num_vel_ok += 1;
        }

        if offset == 2 {
            // Only the print the final estimate
            println!("Estimate:\n{}", est);
            println!("Smoother estimate:\n{}", smoothed_est);
            println!("Truth:\n{}", truth_state);

            println!(
                "RSS error: estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                err_p * 1e3,
                err_v * 1e3,
                (truth_state - est.state().orbit).unwrap()
            );

            println!(
                "RSS error: smoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                err_p_sm * 1e3,
                err_v_sm * 1e3,
                (truth_state - smoothed_est.state().orbit).unwrap()
            );

            // Check that the covariance decreased for the final estimate
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

            let rmag_err_sm = (truth_state - smoothed_est.state().orbit)
                .unwrap()
                .rmag_km();
            let rmag_err = (truth_state - est.state().orbit).unwrap().rmag_km();
            assert!(
                rmag_err_sm < 0.150 || rmag_err_sm < rmag_err,
                "final radius error should be on ~ 150 meter level (is instead {:.3} m) OR the smoothed estimate should be better (but {:.3} m > {:.3} m)",
                rmag_err_sm * 1e3, rmag_err_sm * 1e3, rmag_err * 1e3
            );
        }

        // The smoothed RSS errors should be better, or have the same order of magnitude or not significantly worse

        // Compute orders of magnitude
        let err_p_oom = err_p.log10().floor() as i32;
        let err_v_oom = err_v.log10().floor() as i32;
        let err_p_sm_oom = err_p_sm.log10().floor() as i32;
        let err_v_sm_oom = err_v_sm.log10().floor() as i32;

        if err_p_sm_oom - err_p_oom > 2 {
            println!(
                "RSS position error after smoothing not better @{} (#{}):\n\testimate vs truth: {:.3e} m\t{:.3e} m/s\n{}\n\tsmoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                truth_state.epoch,
                od_sol.estimates.len() - offset,
                err_p * 1e3,
                err_v * 1e3,
                (truth_state - est.state().orbit).unwrap(),
                err_p_sm * 1e3,
                err_v_sm * 1e3,
                (truth_state - smoothed_est.state().orbit).unwrap()
            );
        }

        if err_v_sm_oom - err_v_oom > 2 {
            println!(
                "RSS velocity error after smoothing not better @{} (#{}):\n\testimate vs truth: {:.3e} m\t{:.3e} m/s\n{}\n\tsmoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                truth_state.epoch,
                od_sol.estimates.len() - offset,
                err_p * 1e3,
                err_v * 1e3,
                (truth_state - est.state().orbit).unwrap(),
                err_p_sm * 1e3,
                err_v_sm * 1e3,
                (truth_state - smoothed_est.state().orbit).unwrap()
            );
        }

        for i in 0..6 {
            if est.covar[(i, i)] < 0.0 {
                println!(
                    "covar diagonal element negative @ [{}, {}] = {:.3e}: @{} (#{}) -- issue #164",
                    i,
                    i,
                    est.covar[(i, i)],
                    truth_state.epoch,
                    od_sol.estimates.len() - offset,
                );
            }
        }
    }

    let cntf = od_sol.estimates.len() as f64;
    println!(
        "\nPos. better: {}/{}\tVel. better: {}/{}\nPre-smoothing  avr. RSS:\t{:.3e}\t{:.3e}\nPost-smoothing avr. RSS:\t{:.3e}\t{:.3e}\n",
        num_pos_ok,
        od_sol.estimates.len(),
        num_vel_ok,
        od_sol.estimates.len(),
        rss_pos_avr / cntf,
        rss_vel_avr / cntf,
        rss_pos_avr_sm / cntf,
        rss_vel_avr_sm / cntf,
    );

    // For the EKF with SNC, the estimate is already an order of magnitude better than the equivalent scenario with a CKF.
    // Hence, we expect the average RSS errors to be no worse than one order of magnitude away from what they initially were.
    assert!(
        rss_pos_avr_sm.log10().floor() - rss_pos_avr.log10().floor() < 2.0,
        "Average RSS position error more than two orders of magnitude worse"
    );
    assert!(
        rss_vel_avr_sm.log10().floor() - rss_vel_avr.log10().floor() < 2.0,
        "Average RSS velocity error more than two orders of magnitude worse"
    );
}

#[allow(clippy::identity_op)]
#[rstest]
fn xhat_dev_test_ckf_iteration_multi_body(
    almanac: Arc<Almanac>,
    devices: BTreeMap<String, GroundStation>,
) {
    let _ = pretty_env_logger::try_init();

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    // Define the tracking configurations
    let mut configs = BTreeMap::new();
    configs.insert(
        "Madrid".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        "Canberra".to_string(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.radius_km.x += 9.5;
    initial_state_dev.radius_km.y -= 9.5;
    initial_state_dev.radius_km.z += 9.5;

    let (err_p, err_v) = rss_orbit_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        (initial_state - initial_state_dev).unwrap()
    );

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let setup = Propagator::new(orbital_dyn, IntegratorMethod::RungeKutta4, opts);
    let (_, traj) = setup
        .with(initial_state.into(), almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(devices.clone(), traj.clone(), configs, 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());
    let covar_radius_km = 1.0e2;
    let covar_velocity_km_s = 1.0e1;
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
    let initial_estimate = KfEstimate::from_covar(initial_state_dev.into(), init_covar);
    println!("Initial estimate:\n{}", initial_estimate);

    let kf = KF::no_snc(initial_estimate);

    let mut odp = SpacecraftODProcess::ckf(prop_est, kf, devices, None, almanac.clone());

    let od_sol = odp.process_arc(&arc).unwrap();

    // Clone the initial estimates
    let pre_iteration_estimates = od_sol.estimates.clone();

    // Iterate
    let smoothed = od_sol.clone().smooth(almanac.clone()).unwrap();

    let mut rss_pos_avr = 0.0;
    let mut rss_vel_avr = 0.0;
    let mut rss_pos_avr_it = 0.0;
    let mut rss_vel_avr_it = 0.0;
    let mut num_pos_ok = 0;
    let mut num_vel_ok = 0;

    // Compare the initial estimates and the iterated estimates
    // Skip the first 10 estimates which are surprisingly good in this case
    for offset in (1..od_sol.estimates.len() - 10).rev() {
        let prior_est = &pre_iteration_estimates[od_sol.estimates.len() - offset];

        // Check that the covariance deflated
        let est = &smoothed.estimates[od_sol.estimates.len() - offset];
        let truth_state = traj.at(est.epoch()).unwrap().orbit;

        // Some sanity checks to make sure that we have correctly indexed the estimates
        assert_eq!(prior_est.epoch(), est.epoch());
        assert_eq!(est.epoch(), truth_state.epoch);

        let (err_p, err_v) = rss_orbit_errors(&prior_est.state().orbit, &truth_state);
        let (err_p_it, err_v_it) = rss_orbit_errors(&est.state().orbit, &truth_state);

        rss_pos_avr += err_p;
        rss_vel_avr += err_v;
        rss_pos_avr_it += err_p_it;
        rss_vel_avr_it += err_v_it;

        if err_p_it <= err_p {
            num_pos_ok += 1;
        }

        if err_v_it <= err_v {
            num_vel_ok += 1;
        }

        if offset == 2 {
            // Only the print the final estimate
            println!("Estimate:\n{}", prior_est);
            println!("Iterated estimate:\n{}", est);
            println!("Truth:\n{}", truth_state);

            println!(
                "RSS error: estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                err_p * 1e3,
                err_v * 1e3,
                (truth_state - prior_est.state().orbit).unwrap()
            );

            println!(
                "RSS error: iterated estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                err_p_it * 1e3,
                err_v_it * 1e3,
                (truth_state - est.state().orbit).unwrap()
            );

            // Check that the covariance decreased for the final estimate
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

            let rmag_err = (truth_state - est.state().orbit).unwrap().rmag_km();
            assert!(
                rmag_err < 1e-2,
                "final radius error should be on meter level (is instead {:.3} m)",
                rmag_err * 1e3
            );
        }

        // The smoothed RSS errors should be better, or have the same order of magnitude or not significantly worse

        // Compute orders of magnitude
        let err_p_oom = err_p.log10().floor() as i32;
        let err_v_oom = err_v.log10().floor() as i32;
        let err_p_it_oom = err_p_it.log10().floor() as i32;
        let err_v_it_oom = err_v_it.log10().floor() as i32;

        if err_p_it_oom - err_p_oom > 2 {
            println!(
                "RSS position error after iteration not better @{} (#{}):\n\testimate vs truth: {:.3e} m\t{:.3e} m/s\n{}\n\tsmoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                truth_state.epoch,
                od_sol.estimates.len() - offset,
                err_p * 1e3,
                err_v * 1e3,
                (truth_state - prior_est.state().orbit).unwrap(),
                err_p_it * 1e3,
                err_v_it * 1e3,
                (truth_state - est.state().orbit).unwrap()
            );
        }

        if err_v_it_oom - err_v_oom > 3 {
            println!(
                "RSS velocity error after smoothing not better @{} (#{}):\n\testimate vs truth: {:.3e} m\t{:.3e} m/s\n{}\n\tsmoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                truth_state.epoch,
                od_sol.estimates.len() - offset,
                err_p * 1e3,
                err_v * 1e3,
                (truth_state - prior_est.state().orbit).unwrap(),
                err_p_it * 1e3,
                err_v_it * 1e3,
                (truth_state - est.state().orbit).unwrap()
            );
        }

        for i in 0..6 {
            if est.covar[(i, i)] < 0.0 {
                println!(
                    "covar diagonal element negative @ [{}, {}] = {:.3e}: @{} (#{}) -- issue #164",
                    i,
                    i,
                    est.covar[(i, i)],
                    truth_state.epoch,
                    od_sol.estimates.len() - offset,
                );
            }
        }
    }

    let cntf = od_sol.estimates.len() as f64;
    println!(
        "\nPos. better: {}/{}\tVel. better: {}/{}\nPre-iteration  avr. RSS:\t{:.3e}\t{:.3e}\nPost-iteration avr. RSS:\t{:.3e}\t{:.3e}\n",
        num_pos_ok,
        od_sol.estimates.len(),
        num_vel_ok,
        od_sol.estimates.len(),
        rss_pos_avr / cntf,
        rss_vel_avr / cntf,
        rss_pos_avr_it / cntf,
        rss_vel_avr_it / cntf,
    );

    // For the CKF, the average RSS errors are expected to be better or on the same order of magnitude.
    assert!(
        rss_pos_avr_it.log10().floor() - rss_pos_avr.log10().floor() < 2.0,
        "Average RSS position error more than two orders of magnitude worse"
    );
    assert!(
        rss_vel_avr_it.log10().floor() - rss_vel_avr.log10().floor() < 2.0,
        "Average RSS velocity error more than two orders of magnitude worse"
    );
}
