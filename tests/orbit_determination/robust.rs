extern crate csv;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use nyx::cosmic::{Bodies, Cosm, Orbit};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::io::formatter::NavSolutionFormatter;
use nyx::io::ExportCfg;
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
fn od_robust_test_ekf_realistic_one_way() {
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
    let configs = HashMap::from([
        (
            dss65_madrid.name.clone(),
            TrkConfig::from_sample_rate(60.seconds()),
        ),
        (
            dss34_canberra.name.clone(),
            TrkConfig::from_sample_rate(60.seconds()),
        ),
    ]);

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
        "Initial state dev:\t{:.3} m\t{:.3} m/s\n{}",
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

    // Define the process noise to assume an unmodeled acceleration on X, Y and Z in the ECI frame
    let sigma_q = 5e-10_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);

    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    let trig = EkfTrigger::new(ekf_num_meas, ekf_disable_time);

    let mut odp = ODProcess::ekf(prop_est, kf, trig, None, cosm.clone());

    // Let's filter and iterate on the initial subset of the arc to refine the initial estimate
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
            "epoch",
            "x_err_km",
            "y_err_km",
            "z_err_km",
            "vx_err_km_s",
            "vy_err_km_s",
            "vz_err_km_s",
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

    for res in odp.residuals.iter().flatten() {
        resid_wtr
            .serialize(vec![
                res.epoch.to_string(),
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

#[allow(clippy::identity_op)]
#[test]
fn od_robust_test_ekf_realistic_two_way() {
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

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let mut dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        GaussMarkov::high_precision_range_km(),
        GaussMarkov::high_precision_doppler_km_s(),
        iau_earth,
    );
    // Set the integration time so as to generate two way measurements
    dss65_madrid.integration_time = Some(60.seconds());
    let mut dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        GaussMarkov::high_precision_range_km(),
        GaussMarkov::high_precision_doppler_km_s(),
        iau_earth,
    );
    dss34_canberra.integration_time = Some(60.seconds());

    // Define the tracking configurations
    let configs = HashMap::from([
        (
            dss65_madrid.name.clone(),
            TrkConfig {
                // Make sure to start the tracking one integration time after the start of the trajectory
                start: simulator::Availability::Epoch(dt + 60.seconds()),
                sampling: 60.seconds(),
                ..Default::default()
            },
        ),
        (
            dss34_canberra.name.clone(),
            TrkConfig {
                // Make sure to start the tracking one integration time after the start of the trajectory
                start: simulator::Availability::Epoch(dt + 60.seconds()),
                sampling: 60.seconds(),
                ..Default::default()
            },
        ),
    ]);

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let devices = vec![dss65_madrid, dss34_canberra];

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
        "Initial state dev:\t{:.3} m\t{:.3} m/s\n{}",
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
    let mut arc_sim = TrackingArcSim::with_seed(devices.clone(), traj.clone(), configs, 0).unwrap();
    arc_sim.disallow_overlap(); // Prevent overlapping measurements

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    // And serialize to disk
    let path: PathBuf = [
        &env::var("CARGO_MANIFEST_DIR").unwrap(),
        "output_data",
        "ekf_robust_two_way_msr.parquet",
    ]
    .iter()
    .collect();

    arc.to_parquet_simple(path.clone()).unwrap();

    println!("{arc}");

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be _nearly_ perfect because we've removed Saturn from the estimated trajectory
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let estimator = OrbitalDynamics::point_masses(&bodies, cosm.clone());
    let setup = Propagator::new::<RK4Fixed>(estimator, opts);
    let prop_est = setup.with(initial_state_dev.with_stm());

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    // Define the process noise to assume an unmodeled acceleration on X, Y and Z in the ECI frame
    let sigma_q = 5e-10_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);

    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    let trig = EkfTrigger::new(ekf_num_meas, ekf_disable_time);

    let mut odp = ODProcess::ekf(prop_est, kf, trig, None, cosm.clone());

    // TODO: Fix the deserialization of the measurements such that they also deserialize the integration time.
    // Without it, we're stuck having to rebuild them from scratch.
    // https://github.com/nyx-space/nyx/issues/140

    // Build the hashmap of devices from the vector using their names
    let mut devices_map = devices
        .into_iter()
        .map(|dev| (dev.name.clone(), dev))
        .collect::<HashMap<_, _>>();

    odp.process(
        &arc.measurements,
        &mut devices_map,
        arc.min_duration_sep().unwrap(),
    )
    .unwrap();

    // Export as Parquet
    odp.to_parquet(
        path.with_file_name("robustness_test_two_way.parquet"),
        ExportCfg::timestamped(),
    )
    .unwrap();

    let estimate_fmtr =
        NavSolutionFormatter::default("robustness_test_two_way.csv".to_owned(), cosm);

    let mut wtr = csv::Writer::from_path("robustness_test_two_way.csv").unwrap();
    wtr.serialize(&estimate_fmtr.headers)
        .expect("could not write to output file");

    let mut err_wtr = csv::Writer::from_path("robustness_test_traj_err_two_way.csv").unwrap();
    err_wtr
        .serialize(vec![
            "epoch",
            "x_err_km",
            "y_err_km",
            "z_err_km",
            "vx_err_km_s",
            "vy_err_km_s",
            "vz_err_km_s",
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

    let mut resid_wtr = csv::Writer::from_path("robustness_test_residuals_two_way.csv").unwrap();
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

    for res in odp.residuals.iter().flatten() {
        resid_wtr
            .serialize(vec![
                res.epoch.to_string(),
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
        delta.rmag_km() < 0.01,
        "Position error should be less than 10 meters"
    );
    assert!(
        delta.vmag_km_s() < 1e-5,
        "Velocity error should be on centimeter level"
    );
}
