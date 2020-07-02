extern crate csv;
extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use self::hifitime::{Epoch, SECONDS_PER_DAY};
use self::na::{Matrix2, Matrix3, Matrix6, Vector2, Vector3, Vector6};
use self::nyx::celestia::{bodies, Cosm, State};
use self::nyx::dynamics::orbital::{OrbitalDynamics, OrbitalDynamicsStm};
use self::nyx::dynamics::sph_harmonics::{Harmonics, HarmonicsDiff};
use self::nyx::io::gravity::*;
use self::nyx::od::ui::*;
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
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
    let ekf_num_meas = 400;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 10.0;
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
    let dt = Epoch::from_mjd_tai(21545.0);
    let initial_state = State::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.x += 9.5;
    initial_state_dev.y -= 9.5;
    initial_state_dev.z += 9.5;

    println!("Initial state dev:\n{}", initial_state - initial_state_dev);

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

    let kf = KF::no_snc(initial_estimate, measurement_noise);

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

    println!("Initial state dev:\n{}", initial_state - initial_state_dev);

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
    trig.when_converged = true;

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

#[test]
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
    trig.when_converged = true;

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
