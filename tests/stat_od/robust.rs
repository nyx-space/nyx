extern crate csv;
extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use self::hifitime::{Epoch, SECONDS_PER_DAY};
use self::na::{Matrix2, Matrix6, Vector2, Vector6};
use self::nyx::celestia::{bodies, Cosm, State};
use self::nyx::dynamics::orbital::{OrbitalDynamics, OrbitalDynamicsStm, PointMasses};
use self::nyx::dynamics::sph_harmonics::{Harmonics, HarmonicsDiff};
use self::nyx::io::gravity::*;
use self::nyx::od::ui::*;
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use self::nyx::utils::rss_state_errors;
use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};

#[test]
fn robust_test_ekf_two_body() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::thread;

    let cosm = Cosm::de438();

    // Define the ground stations.
    let ekf_num_meas = 100;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 3600.0;
    let elevation_mask = 0.0;
    let range_noise = 0.0;
    let range_rate_noise = 0.0;
    let dss65_madrid =
        GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise, &cosm);
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise, &cosm);

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let all_stations = vec![dss65_madrid, dss34_canberra];

    // Define the propagator information.
    let prop_time = SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx): (Sender<State>, Receiver<State>) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = State::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.x += 5.0;
    initial_state_dev.y -= 5.0;
    initial_state_dev.z += 5.0;

    let (err_p, err_v) = rss_state_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        initial_state - initial_state_dev
    );

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let mut dynamics = OrbitalDynamics::two_body(initial_state);
        let mut prop = Propagator::new::<RK4Fixed>(&mut dynamics, &opts);
        prop.tx_chan = Some(truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    let mut truth_states = Vec::with_capacity(10_000);
    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.recv() {
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state).unwrap();
            if meas.visible() {
                measurements.push(meas);
                break; // We know that only one station is in visibility at each time.
            }
        }
        truth_states.push(rx_state)
    }
    let final_truth_state = truth_states[truth_states.len() - 1];

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let opts_est = PropOpts::with_fixed_step(step_size);
    let mut tb_estimator = OrbitalDynamicsStm::two_body(initial_state_dev);
    let prop_est = Propagator::new::<RK4Fixed>(&mut tb_estimator, &opts_est);
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
    let process_noise = SNC3::from_diagonal(120.0, &[sigma_q, sigma_q, sigma_q]);
    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    let mut odp = ODProcess::ekf(
        prop_est,
        kf,
        all_stations,
        false,
        measurements.len(),
        StdEkfTrigger::new(ekf_num_meas, ekf_disable_time),
    );

    let rtn = odp.process_measurements(&measurements);
    assert!(rtn.is_none(), "ekf failed");

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("Estimate:\n{}", est);
    println!("Truth:\n{}", final_truth_state);
    let (err_p, err_v) = rss_state_errors(&est.state(), &final_truth_state);
    println!(
        "Delta state with truth (epoch match: {}): {:.3} m\t{:.3} m/s\n{}",
        final_truth_state.dt == est.epoch(),
        err_p * 1e3,
        err_v * 1e3,
        final_truth_state - est.state()
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
    let rmag_err = (final_truth_state - est.state()).rmag();
    assert!(
        rmag_err < 1e-2,
        "final radius error should be on meter level (is instead {:.3} m)",
        rmag_err * 1e3
    );

    assert_eq!(
        truth_states.len(),
        odp.estimates.len(),
        "different number of estimates"
    );
}

#[test]
fn robust_test_ekf_multi_body() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::thread;

    let cosm = Cosm::de438();

    // Define the ground stations.
    let ekf_num_meas = 500;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 10.0;
    let elevation_mask = 0.0;
    let range_noise = 1e-5;
    let range_rate_noise = 1e-7;
    let dss65_madrid =
        GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise, &cosm);
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise, &cosm);

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let all_stations = vec![dss65_madrid, dss34_canberra];

    // Define the propagator information.
    let prop_time = SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = State::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.x += 9.5;
    initial_state_dev.y -= 9.5;
    initial_state_dev.z += 9.5;

    let (err_p, err_v) = rss_state_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        initial_state - initial_state_dev
    );

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let cosm = Cosm::de438();
        let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
        let mut dynamics = OrbitalDynamics::point_masses(initial_state, bodies, &cosm);
        let mut prop = Propagator::new::<RK4Fixed>(&mut dynamics, &opts);
        prop.tx_chan = Some(truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    let mut truth_states = Vec::with_capacity(10_000);
    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.recv() {
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state).unwrap();
            if meas.visible() {
                measurements.push(meas);
                break; // We know that only one station is in visibility at each time.
            }
        }
        truth_states.push(rx_state)
    }

    // Note that we check the second to last state so we can compare with the smoother
    let offset = truth_states.len() - 10;

    let final_truth_state = truth_states[truth_states.len() - offset];

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let opts_est = PropOpts::with_fixed_step(step_size);
    let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
    let mut tb_estimator = OrbitalDynamicsStm::point_masses(initial_state, bodies, &cosm);
    let prop_est = Propagator::new::<RK4Fixed>(&mut tb_estimator, &opts_est);
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

    // let kf = KF::no_snc(initial_estimate, measurement_noise);
    let sigma_q = 1e-7_f64.powi(2);
    let process_noise = SNC3::from_diagonal(120.0, &[sigma_q, sigma_q, sigma_q]);
    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    let mut trig = StdEkfTrigger::new(ekf_num_meas, ekf_disable_time);
    trig.within_sigma = 3.0;

    let mut odp = ODProcess::ekf(prop_est, kf, all_stations, false, measurements.len(), trig);

    let rtn = odp.process_measurements(&measurements);
    assert!(rtn.is_none(), "ekf failed");

    // Smoother
    let smoothed_est = &odp.smooth().unwrap()[odp.estimates.len() - offset];

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - offset];

    println!("Estimate:\n{}", est);
    println!("Smoother estimate:\n{}", smoothed_est);
    println!("Truth:\n{}", final_truth_state);

    // Some sanity checks to make sure that we have correctly indexed the estimates
    assert_eq!(smoothed_est.epoch(), final_truth_state.dt);
    assert_eq!(est.epoch(), final_truth_state.dt);

    let (err_p, err_v) = rss_state_errors(&est.state(), &final_truth_state);
    let (err_p_sm, err_v_sm) = rss_state_errors(&smoothed_est.state(), &final_truth_state);

    // Some printing for debugging
    println!(
        "RSS error: estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        final_truth_state - est.state()
    );

    println!(
        "RSS error: smoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
        err_p_sm * 1e3,
        err_v_sm * 1e3,
        final_truth_state - smoothed_est.state()
    );

    // The smoothed RSS errors should be equal in the worst case

    assert!(
        err_p_sm <= err_p,
        "RSS position error after smoothing not better"
    );
    assert!(
        err_v_sm <= err_v,
        "RSS velocity error after smoothing not better"
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
    let rmag_err = (final_truth_state - est.state()).rmag();
    assert!(
        rmag_err < 1e-2,
        "final radius error should be on meter level (is instead {:.3} m)",
        rmag_err * 1e3
    );

    assert_eq!(
        truth_states.len(),
        odp.estimates.len() - 1,
        "different number of estimates"
    );
}

#[test]
#[ignore]
fn robust_test_ekf_harmonics() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::thread;

    let cosm = Cosm::de438();

    // Define the ground stations.
    let ekf_num_meas = 50000;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 3600.0;
    let elevation_mask = 0.0;
    let range_noise = 1e-5;
    let range_rate_noise = 1e-7;
    let dss65_madrid =
        GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise, &cosm);
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise, &cosm);

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let all_stations = vec![dss65_madrid, dss34_canberra];

    // Define the propagator information.
    let prop_time = SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let iau_earth = cosm.frame("IAU Earth");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = State::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.x += 0.00;
    initial_state_dev.y -= 0.00;
    initial_state_dev.z += 0.00;

    let (err_p, err_v) = rss_state_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        initial_state - initial_state_dev
    );

    let hh_deg = 20;
    let hh_ord = 20;

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let cosm = Cosm::de438();
        let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
        let mut dynamics = OrbitalDynamics::two_body(initial_state);
        let earth_sph_harm =
            HarmonicsMem::from_cof("data/JGM3.cof.gz", hh_deg, hh_ord, true).unwrap();
        let harmonics = Harmonics::from_stor(iau_earth, earth_sph_harm, &cosm);
        dynamics.add_model(Box::new(harmonics));
        dynamics.add_model(Box::new(PointMasses::new(eme2k, bodies, &cosm)));
        let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
        let mut dynamics = OrbitalDynamics::point_masses(initial_state, bodies, &cosm);
        let mut prop = Propagator::new::<RK4Fixed>(&mut dynamics, &opts);
        prop.tx_chan = Some(truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    let mut truth_states = Vec::with_capacity(10_000);
    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.recv() {
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state).unwrap();
            if meas.visible() {
                measurements.push(meas);
                break; // We know that only one station is in visibility at each time.
            }
        }
        truth_states.push(rx_state)
    }
    let final_truth_state = truth_states[truth_states.len() - 1];

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let opts_est = PropOpts::with_fixed_step(step_size);
    let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];

    let mut estimator = OrbitalDynamicsStm::two_body(initial_state);
    let earth_sph_harm = HarmonicsMem::from_cof("data/JGM3.cof.gz", hh_deg, hh_ord, true).unwrap();
    let harmonics = HarmonicsDiff::from_stor(iau_earth, earth_sph_harm, &cosm);
    estimator.add_model(Box::new(harmonics));
    estimator.add_model(Box::new(PointMasses::new(eme2k, bodies, &cosm)));

    let prop_est = Propagator::new::<RK4Fixed>(&mut estimator, &opts_est);
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

    // let kf = KF::no_snc(initial_estimate, measurement_noise);
    let sigma_q = 1e-6_f64.powi(2);
    let process_noise = SNC3::from_diagonal(120.0, &[sigma_q, sigma_q, sigma_q]);
    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    let mut trig = StdEkfTrigger::new(ekf_num_meas, ekf_disable_time);
    trig.within_sigma = 3.0;

    let mut odp = ODProcess::ekf(prop_est, kf, all_stations, false, measurements.len(), trig);

    let rtn = odp.process_measurements(&measurements);
    assert!(rtn.is_none(), "ekf failed");

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("Estimate:\n{}", est);
    println!("Truth:\n{}", final_truth_state);
    let (err_p, err_v) = rss_state_errors(&est.state(), &final_truth_state);
    println!(
        "Delta state with truth (epoch match: {}): {:.3} m\t{:.3} m/s\n{}",
        final_truth_state.dt == est.epoch(),
        err_p * 1e3,
        err_v * 1e3,
        final_truth_state - est.state()
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
    let rmag_err = (final_truth_state - est.state()).rmag();
    assert!(
        rmag_err < 1e-2,
        "final radius error should be on meter level (is instead {:.3} m)",
        rmag_err * 1e3
    );

    assert_eq!(
        truth_states.len(),
        odp.estimates.len() - 1,
        "different number of estimates"
    );
}

#[test]
#[ignore]
fn robust_test_ekf_realistic() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::thread;

    let cosm = Cosm::de438();

    // Define the ground stations.
    let ekf_num_meas = 500;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 10.0;
    let elevation_mask = 0.0;
    let range_noise = 1e-5;
    let range_rate_noise = 1e-7;
    let dss65_madrid =
        GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise, &cosm);
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise, &cosm);

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let all_stations = vec![dss65_madrid, dss34_canberra];

    // Define the propagator information.
    let prop_time = SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = State::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.x += 9.5;
    initial_state_dev.y -= 9.5;
    initial_state_dev.z += 9.5;

    println!("Initial state dev:\n{}", initial_state - initial_state_dev);

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let cosm = Cosm::de438();
        let bodies = vec![
            bodies::EARTH_MOON,
            bodies::SUN,
            bodies::JUPITER_BARYCENTER,
            bodies::SATURN_BARYCENTER,
        ];
        let mut dynamics = OrbitalDynamics::point_masses(initial_state, bodies, &cosm);
        let mut prop = Propagator::new::<RK4Fixed>(&mut dynamics, &opts);
        prop.tx_chan = Some(truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    let mut truth_states = Vec::with_capacity(10_000);
    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.recv() {
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state).unwrap();
            if meas.visible() {
                measurements.push(meas);
                break; // We know that only one station is in visibility at each time.
            }
        }
        truth_states.push(rx_state)
    }
    let final_truth_state = truth_states[truth_states.len() - 1];

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let opts_est = PropOpts::with_fixed_step(step_size);
    let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
    let mut tb_estimator = OrbitalDynamicsStm::point_masses(initial_state, bodies, &cosm);
    let prop_est = Propagator::new::<RK4Fixed>(&mut tb_estimator, &opts_est);
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

    let mut odp = ODProcess::ekf(prop_est, kf, all_stations, false, measurements.len(), trig);

    let rtn = odp.process_measurements(&measurements);
    assert!(rtn.is_none(), "ekf failed");

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("Estimate:\n{}", est);
    println!("Truth:\n{}", final_truth_state);
    println!(
        "Delta state with truth (epoch match: {}):\n{}",
        final_truth_state.dt == est.epoch(),
        final_truth_state - est.state()
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
    let rmag_err = (final_truth_state - est.state()).rmag();
    assert!(
        rmag_err < 1e-2,
        "final radius error should be on meter level (is instead {:.3} m)",
        rmag_err * 1e3
    );

    assert_eq!(
        truth_states.len(),
        odp.estimates.len() - 1,
        "different number of estimates"
    );
}
