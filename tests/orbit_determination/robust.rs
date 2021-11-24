extern crate csv;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use self::nyx::cosmic::{Bodies, Cosm, Orbit};
use self::nyx::dynamics::orbital::OrbitalDynamics;
use self::nyx::io::formatter::{NavSolutionFormatter, StateFormatter};
use self::nyx::linalg::{Matrix2, Matrix6, Vector2, Vector6};
use self::nyx::od::ui::*;
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use self::nyx::time::{Epoch, TimeUnit};
use self::nyx::utils::rss_orbit_errors;

/*
 * These tests check that if we start with a state deviation in the estimate, the filter will eventually converge back.
 * These tests do NOT check that the filter will converge if the initial state in the propagator has that state deviation.
 * The latter would require iteration and smoothing before playing with an EKF. This will be handled in a subsequent version.
**/

#[allow(clippy::identity_op)]
#[test]
fn od_robust_test_ekf_realistic() {
    // TODO: Add gravity field differences.
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();

    // Define the ground stations.
    let ekf_num_meas = 500;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 10.0 * TimeUnit::Second;
    let elevation_mask = 0.0;
    let range_noise = 1e-5;
    let range_rate_noise = 1e-7;
    let dss65_madrid =
        GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise, cosm.clone());
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise, cosm.clone());

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let all_stations = vec![dss65_madrid, dss34_canberra];

    // Define the propagator information.
    let prop_time = 1 * TimeUnit::Day;
    let step_size = 10.0 * TimeUnit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.x += 9.5;
    initial_state_dev.y -= 9.5;
    initial_state_dev.z += 9.5;

    println!("Initial state dev:\n{}", initial_state - initial_state_dev);

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

    for state in traj.every(10 * TimeUnit::Second) {
        for station in all_stations.iter() {
            let meas = station.measure(&state).unwrap();
            if meas.visible() {
                measurements.push(meas);
                break; // We know that only one station is in visibility at each time.
            }
        }
    }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be _nearly_ perfect because we've removed Saturn from the estimated trajectory
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let estimator = OrbitalDynamics::point_masses(&bodies, cosm);
    let setup = Propagator::new::<RK4Fixed>(estimator, opts);
    let prop_est = setup.with(initial_state.with_stm());
    let covar_radius = 1.0e2;
    let covar_velocity = 1.0e1;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius,
        covar_radius,
        covar_radius,
        covar_velocity,
        covar_velocity,
        covar_velocity,
    ));

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(initial_state_dev, init_covar);
    println!("Initial estimate:\n{}", initial_estimate);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    let kf = KF::no_snc(initial_estimate, measurement_noise);

    let mut trig = StdEkfTrigger::new(ekf_num_meas, ekf_disable_time);
    trig.within_sigma = 3.0;

    let mut odp = ODProcess::ekf(prop_est, kf, all_stations, trig);

    odp.process_measurements(&measurements).unwrap();
    odp.iterate(&measurements, IterationConf::default())
        .unwrap();
    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];
    let final_truth_state = traj.at(est.epoch()).unwrap();

    println!("Estimate:\n{}", est);
    println!("Truth:\n{}", final_truth_state);
    println!(
        "Delta state with truth (epoch match: {}):\n{}",
        final_truth_state.dt == est.epoch(),
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
        if i < 3 {
            assert!(
                est.covar[(i, i)] < covar_radius,
                "covar radius did not decrease"
            );
        } else {
            assert!(
                est.covar[(i, i)] < covar_velocity,
                "covar velocity did not decrease"
            );
        }
    }

    assert_eq!(
        final_truth_state.dt,
        est.epoch(),
        "time of final EST and TRUTH epochs differ"
    );
    let delta = est.state() - final_truth_state;
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag() * 1e3,
        delta.vmag() * 1e6
    );

    assert!(
        delta.rmag() < 1e-2,
        "Position error should be on meter level"
    );
    assert!(
        delta.vmag() < 1e-5,
        "Velocity error should be on millimeter level"
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

    let elevation_mask = 0.0;
    let range_noise = 0.0; // 1e-5
    let range_rate_noise = 0.0; // 1e-7
    let dss65_madrid =
        GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise, cosm.clone());
    let dss34_canberra =
        GroundStation::dss13_goldstone(elevation_mask, range_noise, range_rate_noise, cosm.clone());

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let all_stations = vec![dss65_madrid, dss34_canberra];
    // Bug identified in #147 means we need to redefine the stations here without noise
    let dss65_madrid = GroundStation::dss65_madrid(elevation_mask, 0.0, 0.0, cosm.clone());
    let dss34_canberra = GroundStation::dss13_goldstone(elevation_mask, 0.0, 0.0, cosm.clone());
    let all_stations_no_noise = vec![dss65_madrid, dss34_canberra];

    // Define the propagator information.
    let prop_time = 1 * TimeUnit::Day;
    let step_size = 10.0 * TimeUnit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.9, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.x += 0.0005;
    initial_state_dev.y -= 0.0005;
    initial_state_dev.z += 0.0005;

    let (err_p, err_v) = rss_orbit_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        initial_state - initial_state_dev
    );

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

    for state in traj.every(10 * TimeUnit::Second) {
        for station in all_stations.iter() {
            let meas = station.measure(&state).unwrap();
            if meas.visible() {
                // Always add it to the full list of measurements
                measurements.push(meas);
                break; // We know that only one station is in visibility at each time.
            }
        }

        if let Some(first_state) = initial_state_out {
            wtr.serialize(&truth_fmtr.fmt(&first_state))
                .expect("could not format state");
            initial_state_out = None;
        }
        wtr.serialize(truth_fmtr.fmt(&state))
            .expect("could not format state");
    }

    let ekf_msr_trig = measurements.len() / 1;

    println!(
        "Generated {} measurements in total (using {} for CKF)",
        measurements.len(),
        ekf_msr_trig
    );

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let orbital_dyn = OrbitalDynamics::point_masses(&[Bodies::Luna, Bodies::Sun], cosm.clone());
    let setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);
    let prop_est = setup.with(initial_state_dev.with_stm());
    let covar_radius = 1.0e2;
    let covar_velocity = 1.0e1;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius,
        covar_radius,
        covar_radius,
        covar_velocity,
        covar_velocity,
        covar_velocity,
    ));

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(initial_state_dev, init_covar);
    println!("Initial estimate:\n{}", initial_estimate);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    let sigma_q = 1e-7_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2 * TimeUnit::Minute, &[sigma_q, sigma_q, sigma_q]);
    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    let mut trig = StdEkfTrigger::new(ekf_msr_trig, 10.0 * TimeUnit::Second);
    trig.within_sigma = 3.0;

    let mut odp = ODProcess::ekf(prop_est, kf, all_stations_no_noise, trig);

    odp.process_measurements(&measurements).unwrap();

    // Clone the initial estimate
    let pre_smooth_first_est = odp.estimates[0].clone();
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
    odp.iterate(&measurements, IterationConf::default())
        .unwrap();

    let fmtr = NavSolutionFormatter::default(
        "data/robust_test_post_iteration.csv".to_string(),
        cosm.clone(),
    );
    let mut wtr = csv::Writer::from_path(fmtr.filename.clone()).expect("could not create file");
    wtr.serialize(&fmtr.headers)
        .expect("could not write headers");

    for est in &odp.estimates {
        wtr.serialize(fmtr.fmt(est))
            .expect("could not format state");
    }

    let post_smooth_first_est = odp.estimates[0].clone();

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
    let rmag_err = (final_truth_state - est.state()).rmag();
    assert!(
        rmag_err < 1e-2,
        "final radius error should be on meter level (is instead {:.3} m)",
        rmag_err * 1e3
    );
}
