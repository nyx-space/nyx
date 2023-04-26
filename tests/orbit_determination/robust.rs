extern crate csv;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use nyx::cosmic::{Bodies, Cosm, Orbit};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::io::formatter::{NavSolutionFormatter, StateFormatter};
use nyx::linalg::{Matrix2, Vector2};
use nyx::md::StateParameter;
use nyx::od::noise::GaussMarkov;
use nyx::od::prelude::*;
use nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use nyx::time::{Epoch, TimeUnits, Unit};
use nyx::utils::rss_orbit_errors;
use std::collections::HashMap;
use std::env;
use std::path::PathBuf;

/*
 * These tests check that if we start with a state deviation in the estimate, the filter will eventually converge back.
 * These tests do NOT check that the filter will converge if the initial state in the propagator has that state deviation.
 * The latter would require iteration and smoothing before playing with an EKF. This will be handled in a subsequent version.
**/

#[allow(clippy::identity_op)]
#[test]
fn od_robust_test_ekf_realistic() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();

    let iau_earth = cosm.frame("IAU Earth");
    // Define the ground stations.
    let ekf_num_meas = 300;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 3 * Unit::Minute;
    let elevation_mask = 0.0;

    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        GaussMarkov::high_precision_range_km(),
        GaussMarkov::high_precision_doppler_km_s(),
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        GaussMarkov::high_precision_range_km(),
        GaussMarkov::high_precision_doppler_km_s(),
        iau_earth,
    );

    // Define the tracking configurations
    let mut configs = HashMap::new();
    configs.insert(
        dss65_madrid.name.clone(),
        TrkConfig::from_sample_rate(60.seconds()),
    );
    configs.insert(
        dss34_canberra.name.clone(),
        TrkConfig::from_sample_rate(60.seconds()),
    );

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let all_stations = vec![dss65_madrid, dss34_canberra];

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // Disperse the initial state based on some orbital elements errors given from ULA Atlas 5 user guide, table 2.3.3-1 <https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf>
    // This assumes that the errors are ONE TENTH of the values given in the table. It assumes that the launch provider has provided an initial state vector, whose error is lower than the injection errors.
    // The initial covariance is computed based on the realized dispersions.
    let initial_estimate = KfEstimate::disperse_from_diag(
        initial_state,
        &[
            (StateParameter::Inclination, 0.0025),
            (StateParameter::RAAN, 0.022),
            (StateParameter::AoP, 0.02),
        ],
        Some(0),
    );
    println!("Initial estimate:\n{}", initial_estimate);

    let initial_state_dev = initial_estimate.nominal_state;
    let (init_rss_pos_km, init_rss_vel_km_s) = rss_orbit_errors(&initial_state, &initial_state_dev);

    println!("Truth initial state:\n{initial_state}\n{initial_state:x}");
    println!("Filter initial state:\n{initial_state_dev}\n{initial_state_dev:x}");
    println!(
        "Initial state dev:\t{:.3} km\t{:.3} km/s\n{}",
        init_rss_pos_km * 1e3,
        init_rss_vel_km_s * 1e3,
        initial_state - initial_state_dev
    );

    let bodies = vec![
        Bodies::Luna,
        Bodies::Sun,
        Bodies::JupiterBarycenter,
        Bodies::SaturnBarycenter,
    ];
    let orbital_dyn = OrbitalDynamics::point_masses(&bodies, cosm.clone());
    let truth_setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);
    let (_, traj) = truth_setup
        .with(initial_state)
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj.clone(), configs, 0).unwrap();
    arc_sim.disallow_overlap(); // Prevent overlapping measurements

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    // And serialize to disk
    let path: PathBuf = [
        &env::var("CARGO_MANIFEST_DIR").unwrap(),
        "output_data",
        "ekf_robust_msr.parquet",
    ]
    .iter()
    .collect();

    arc.to_parquet_simple(path).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be _nearly_ perfect because we've removed Saturn from the estimated trajectory
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let estimator = OrbitalDynamics::point_masses(&bodies, cosm.clone());
    let setup = Propagator::new::<RK4Fixed>(estimator, opts);
    let prop_est = setup.with(initial_state_dev.with_stm());

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    // Define the process noise to assume an unmodeled acceleration of 1e-3 km^2/s^2 on X, Y and Z in the ECI frame
    let sigma_q = 1e-8_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);

    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    let trig = EkfTrigger::new(ekf_num_meas, ekf_disable_time);

    let mut odp = ODProcess::ekf(prop_est, kf, trig, cosm.clone());

    // Let's filter and iterate on the first four hours of data
    let subset = arc.filter_by_offset(..3.hours());
    let remaining = arc.filter_by_offset(3.hours()..);

    odp.process_arc::<GroundStation>(&subset).unwrap();
    odp.iterate_arc::<GroundStation>(&subset, IterationConf::once())
        .unwrap();

    let (sm_rss_pos_km, sm_rss_vel_km_s) =
        rss_orbit_errors(&initial_state, &odp.estimates[0].state());

    println!(
        "Initial state error after smoothing:\t{:.3} m\t{:.3} m/s\n{}",
        sm_rss_pos_km * 1e3,
        sm_rss_vel_km_s * 1e3,
        initial_state - odp.estimates[0].state()
    );

    odp.process_arc::<GroundStation>(&remaining).unwrap();

    let estimate_fmtr = NavSolutionFormatter::default("robustness_test.csv".to_owned(), cosm);

    let mut wtr = csv::Writer::from_path("robustness_test.csv").unwrap();
    wtr.serialize(&estimate_fmtr.headers)
        .expect("could not write to output file");

    let mut err_wtr = csv::Writer::from_path("robustness_test_traj_err.csv").unwrap();
    err_wtr
        .serialize(vec![
            "epoch", "x_err", "y_err", "z_err", "vx_err", "vy_err", "vz_err",
        ])
        .expect("could not write to output file");

    for est in &odp.estimates {
        // Format the estimate
        wtr.serialize(estimate_fmtr.fmt(est))
            .expect("could not write to CSV");
        // Add the error data
        let truth_state = traj.at(est.epoch()).unwrap();
        let err = truth_state - est.state();
        err_wtr
            .serialize(vec![
                est.epoch().to_string(),
                format!("{}", err.x_km),
                format!("{}", err.y_km),
                format!("{}", err.z_km),
                format!("{}", err.vx_km_s),
                format!("{}", err.vy_km_s),
                format!("{}", err.vz_km_s),
            ])
            .expect("could not write to CSV");
    }

    let mut resid_wtr = csv::Writer::from_path("robustness_test_residuals.csv").unwrap();
    resid_wtr
        .serialize(vec![
            "epoch",
            "range_prefit",
            "doppler_prefit",
            "range_postfit",
            "doppler_postfit",
            "residual_ratio",
        ])
        .expect("could not write to output file");

    for res in &odp.residuals {
        resid_wtr
            .serialize(vec![
                res.dt.to_string(),
                format!("{}", res.prefit[0]),
                format!("{}", res.prefit[1]),
                format!("{}", res.postfit[0]),
                format!("{}", res.postfit[1]),
                format!("{}", res.ratio),
            ])
            .expect("could not write to CSV");
    }

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];
    let final_truth_state = traj.at(est.epoch()).unwrap();

    println!("Estimate:\n{}", est);
    println!("Truth:\n{}", final_truth_state);
    println!(
        "Delta state with truth (epoch match: {}):\n{}",
        final_truth_state.epoch == est.epoch(),
        final_truth_state - est.state()
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
        assert!(
            est.covar[(i, i)] < initial_estimate.covar[(i, i)],
            "covar[({i}, {i})] did not decrease"
        );
    }

    assert_eq!(
        final_truth_state.epoch,
        est.epoch(),
        "time of final EST and TRUTH epochs differ"
    );
    let delta = est.state() - final_truth_state;
    println!(
        "RMAG error = {:.6} m\tVMAG error = {:.6} m/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e3
    );

    assert!(
        delta.rmag_km() < 0.06,
        "Position error should be less than 50 meters"
    );
    assert!(
        delta.vmag_km_s() < 2e-4,
        "Velocity error should be on centimeter level"
    );
}

#[ignore]
#[allow(clippy::identity_op)]
#[test]
fn od_robust_ops_test() {
    // TODO: Enable this test after IOD (https://gitlab.com/nyx-space/nyx/-/issues/196), which will bring down the initial error and allow working with a converged state.
    // TODO: Update this test after #147
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let iau_earth = cosm.frame("IAU Earth");

    let elevation_mask = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        GaussMarkov::high_precision_range_km(),
        GaussMarkov::high_precision_doppler_km_s(),
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss13_goldstone(
        elevation_mask,
        GaussMarkov::high_precision_range_km(),
        GaussMarkov::high_precision_doppler_km_s(),
        iau_earth,
    );

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    // Define the tracking configurations
    let mut configs = HashMap::new();
    configs.insert(
        dss65_madrid.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        dss34_canberra.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let all_stations = vec![dss65_madrid, dss34_canberra];

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.9, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let orbital_dyn = OrbitalDynamics::point_masses(
        &[
            Bodies::Luna,
            Bodies::Sun,
            Bodies::JupiterBarycenter,
            Bodies::SaturnBarycenter,
        ],
        cosm.clone(),
    );
    let truth_setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);
    let (_, traj) = truth_setup
        .with(initial_state)
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Initialize the truth data output
    let mut initial_state_out = Some(initial_state);
    let truth_fmtr =
        StateFormatter::default("data/robust_test_ckf_truth.csv".to_string(), cosm.clone());
    let mut wtr =
        csv::Writer::from_path(truth_fmtr.filename.clone()).expect("could not create file");
    wtr.serialize(&truth_fmtr.headers)
        .expect("could not write headers");

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj.clone(), configs, 0).unwrap();
    arc_sim.disallow_overlap(); // Prevent overlapping measurements

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    for state in traj.every(10 * Unit::Second) {
        if let Some(first_state) = initial_state_out {
            wtr.serialize(&truth_fmtr.fmt(&first_state))
                .expect("could not format state");
            initial_state_out = None;
        }
        wtr.serialize(truth_fmtr.fmt(&state))
            .expect("could not format state");
    }

    let ekf_msr_trig = arc.measurements.len() / 10;

    println!(
        "Generated {} measurements in total (using {} for CKF)",
        arc.measurements.len(),
        ekf_msr_trig
    );

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.

    // Define the initial estimate
    let initial_estimate = KfEstimate::disperse_from_diag(
        initial_state,
        &[
            (StateParameter::X, 5.0),
            (StateParameter::Y, 5.0),
            (StateParameter::Z, 5.0),
        ],
        None,
    );

    let initial_state_dev = initial_estimate.nominal_state;

    let (err_p, err_v) = rss_orbit_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        initial_state - initial_state_dev
    );

    println!("Initial estimate:\n{}", initial_estimate);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    let sigma_q = 1e-7_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);
    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    let mut trig = EkfTrigger::new(ekf_msr_trig, 10.0 * Unit::Second);
    trig.within_sigma = 3.0;

    let orbital_dyn = OrbitalDynamics::point_masses(&[Bodies::Luna, Bodies::Sun], cosm.clone());
    let setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);
    let prop_est = setup.with(initial_state_dev.with_stm());

    let mut odp = ODProcess::ekf(prop_est, kf, trig, cosm.clone());

    odp.process_arc::<GroundStation>(&arc).unwrap();

    // Clone the initial estimate
    let pre_smooth_first_est = odp.estimates[0];
    // Output the pre-iteration estimates
    let fmtr = NavSolutionFormatter::default(
        "data/robust_test_ckf_pre_iteration.csv".to_string(),
        cosm.clone(),
    );
    let mut wtr = csv::Writer::from_path(fmtr.filename.clone()).expect("could not create file");
    wtr.serialize(&fmtr.headers)
        .expect("could not write headers");

    for est in &odp.estimates {
        wtr.serialize(fmtr.fmt(est))
            .expect("could not format state");
    }

    let est = &odp.estimates.last().unwrap();
    let final_truth_state = traj.at(est.epoch()).unwrap();

    println!("Pre-iteration Estimate:\n{}", est);
    println!("Pre-iteration Truth:\n{}", final_truth_state);

    let (err_p_no_it, err_v_no_it) = final_truth_state.rss(&est.state());

    println!(
        "Pre-iteration RSS error: estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
        err_p_no_it * 1e3,
        err_v_no_it * 1e3,
        final_truth_state - est.state()
    );

    // Iterate
    odp.iterate_arc::<GroundStation>(&arc, IterationConf::default())
        .unwrap();

    let fmtr =
        NavSolutionFormatter::default("data/robust_test_post_iteration.csv".to_string(), cosm);
    let mut wtr = csv::Writer::from_path(fmtr.filename.clone()).expect("could not create file");
    wtr.serialize(&fmtr.headers)
        .expect("could not write headers");

    for est in &odp.estimates {
        wtr.serialize(fmtr.fmt(est))
            .expect("could not format state");
    }

    let post_smooth_first_est = odp.estimates[0];

    let (init_pos_rss, init_vel_rss) = initial_state.rss(&initial_state_dev);
    let (zero_it_pos_rss, zero_it_vel_rss) = initial_state.rss(&pre_smooth_first_est.state());
    let (one_it_pos_rss, one_it_vel_rss) = initial_state.rss(&post_smooth_first_est.state());
    println!(
        "[pos] init: {}\tzero: {}\t one: {}",
        init_pos_rss, zero_it_pos_rss, one_it_pos_rss,
    );
    println!(
        "[vel] init: {}\tzero: {}\t one: {}",
        init_vel_rss, zero_it_vel_rss, one_it_vel_rss
    );

    // NOTE: From other tests, it's evident that there isn't enough visibility in this 0.9 eccentric orbit to change the initial estimate.
    assert!(
        one_it_pos_rss <= zero_it_pos_rss,
        "RSS position not better after iteration"
    );

    let est = &odp.estimates.last().unwrap();
    let final_truth_state = traj.at(est.epoch()).unwrap();

    println!("Estimate:\n{}", est);
    println!("Truth:\n{}", final_truth_state);

    let (err_p, err_v) = final_truth_state.rss(&est.state());

    println!(
        "RSS error: estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        final_truth_state - est.state()
    );

    // Check that the iteration leads to better results
    assert!(
        err_p <= err_p_no_it,
        "Position error not better after iteration: {:.3e} m < {:.3e} m",
        err_p * 1e3,
        err_p_no_it * 1e3
    );
    assert!(
        err_v <= err_v_no_it,
        "Velocity error not better after iteration {:.3e} m/s < {:.3e} m/s",
        err_v * 1e3,
        err_v_no_it * 1e3
    );

    // Reenable after #147
    let rmag_err = (final_truth_state - est.state()).rmag_km();
    assert!(
        rmag_err < 1e-2,
        "final radius error should be on meter level (is instead {:.3} m)",
        rmag_err * 1e3
    );
}
