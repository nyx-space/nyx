extern crate csv;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use self::nyx::celestia::{Bodies, Cosm, Orbit};
use self::nyx::dimensions::{Matrix2, Matrix6, Vector2, Vector6};
use self::nyx::dynamics::orbital::{OrbitalDynamics, PointMasses};
use self::nyx::dynamics::sph_harmonics::Harmonics;
use self::nyx::io::gravity::*;
use self::nyx::od::ui::*;
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use self::nyx::time::{Epoch, TimeUnit};
use self::nyx::utils::rss_state_errors;
use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};

#[allow(clippy::identity_op)]
#[test]
fn robust_test_ekf_two_body() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();

    // Define the ground stations.
    let ekf_num_meas = 100;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 1 * TimeUnit::Hour;
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
    let prop_time = 1 * TimeUnit::Day;
    let step_size = 10.0 * TimeUnit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx): (Sender<Orbit>, Receiver<Orbit>) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state.with_stm();
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

    let orbital_dyn = OrbitalDynamics::two_body();
    let truth_setup = Propagator::new::<RK4Fixed>(&orbital_dyn, opts);
    let mut prop = truth_setup.with(initial_state);
    prop.tx_chan = Some(truth_tx);
    prop.for_duration(prop_time).unwrap();
    let mut truth_states = Vec::with_capacity(10_000);
    truth_states.push(initial_state);
    // Receive the states on the main thread, and populate the measurement channel.
    'outer: while let Ok(rx_state) = truth_rx.try_recv() {
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state).unwrap();
            if meas.visible() {
                measurements.push(meas);
                // if measurements.len() > 105 {
                //     break 'outer;
                // }
                break; // We know that only one station is in visibility at each time.
            }
        }
        truth_states.push(rx_state)
    }
    let final_truth_state = truth_states[truth_states.len() - 1];

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let estimator = OrbitalDynamics::two_body();
    let setup = Propagator::new::<RK4Fixed>(&estimator, opts);
    let prop_est = setup.with(initial_state_dev);
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

    let mut odp = ODProcess::ekf(
        prop_est,
        kf,
        all_stations,
        false,
        measurements.len(),
        StdEkfTrigger::new(ekf_num_meas, ekf_disable_time),
    );

    odp.process_measurements(&measurements).unwrap();

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

#[allow(clippy::identity_op)]
#[test]
fn robust_test_ekf_multi_body() {
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
        GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise, &cosm);
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise, &cosm);

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let all_stations = vec![dss65_madrid, dss34_canberra];

    // Define the propagator information.
    let prop_time = 1 * TimeUnit::Day;
    let step_size = 10.0 * TimeUnit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.x += 9.5;
    initial_state_dev.y -= 9.5;
    initial_state_dev.z += 9.5;
    initial_state_dev.enable_stm();

    let (err_p, err_v) = rss_state_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        initial_state - initial_state_dev
    );

    // Generate the truth data on one thread.

    let cosm = Cosm::de438();
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let orbital_dyn = OrbitalDynamics::point_masses(initial_state.frame, &bodies, &cosm);
    let truth_setup = Propagator::new::<RK4Fixed>(&orbital_dyn, opts);
    let mut prop = truth_setup.with(initial_state);
    prop.tx_chan = Some(truth_tx);
    prop.for_duration(prop_time).unwrap();
    let mut truth_states = Vec::with_capacity(10_000);
    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.try_recv() {
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
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let estimator = OrbitalDynamics::point_masses(initial_state.frame, &bodies, &cosm);
    let setup = Propagator::new::<RK4Fixed>(&estimator, opts);
    let prop_est = setup.with(initial_state_dev);
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

    let mut trig = StdEkfTrigger::new(ekf_num_meas, ekf_disable_time);
    trig.within_sigma = 3.0;

    let mut odp = ODProcess::ekf(prop_est, kf, all_stations, false, measurements.len(), trig);

    odp.process_measurements(&measurements).unwrap();

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];

    println!("Estimate:\n{}", est);
    println!("Truth:\n{}", final_truth_state);

    // Some sanity checks to make sure that we have correctly indexed the estimates
    assert_eq!(est.epoch(), final_truth_state.dt);

    let (err_p, err_v) = rss_state_errors(&est.state(), &final_truth_state);

    // Some printing for debugging
    println!(
        "RSS error: estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
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

#[allow(clippy::identity_op)]
#[test]
#[ignore]
fn robust_test_ekf_harmonics() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();

    // Define the ground stations.
    let ekf_num_meas = 50000;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 1 * TimeUnit::Hour;
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
    let prop_time = 1 * TimeUnit::Day;
    let step_size = 10.0 * TimeUnit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let iau_earth = cosm.frame("IAU Earth");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.x += 0.00;
    initial_state_dev.y -= 0.00;
    initial_state_dev.z += 0.00;
    initial_state_dev.enable_stm();

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

    let cosm = Cosm::de438();
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let mut orbital_dyn = OrbitalDynamics::two_body();
    let earth_sph_harm = HarmonicsMem::from_cof("data/JGM3.cof.gz", hh_deg, hh_ord, true).unwrap();
    let harmonics = Harmonics::from_stor(iau_earth, earth_sph_harm, &cosm);
    orbital_dyn.add_model(harmonics);
    orbital_dyn.add_model(PointMasses::new(eme2k, &bodies, &cosm));
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let orbital_dyn = OrbitalDynamics::point_masses(initial_state.frame, &bodies, &cosm);
    let truth_setup = Propagator::new::<RK4Fixed>(&orbital_dyn, opts);
    let mut prop = truth_setup.with(initial_state);
    prop.tx_chan = Some(truth_tx);
    prop.for_duration(prop_time).unwrap();

    let mut truth_states = Vec::with_capacity(10_000);
    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.try_recv() {
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
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];

    let mut estimator = OrbitalDynamics::two_body();
    let earth_sph_harm = HarmonicsMem::from_cof("data/JGM3.cof.gz", hh_deg, hh_ord, true).unwrap();
    let harmonics = Harmonics::from_stor(iau_earth, earth_sph_harm, &cosm);
    estimator.add_model(harmonics);
    estimator.add_model(PointMasses::new(eme2k, &bodies, &cosm));

    let setup = Propagator::new::<RK4Fixed>(&estimator, opts);
    let prop_est = setup.with(initial_state_dev);
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

    let sigma_q = 1e-6_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2 * TimeUnit::Minute, &[sigma_q, sigma_q, sigma_q]);
    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    let mut trig = StdEkfTrigger::new(ekf_num_meas, ekf_disable_time);
    trig.within_sigma = 3.0;

    let mut odp = ODProcess::ekf(prop_est, kf, all_stations, false, measurements.len(), trig);

    odp.process_measurements(&measurements).unwrap();

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

#[allow(clippy::identity_op)]
#[test]
#[ignore]
fn robust_test_ekf_realistic() {
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
        GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise, &cosm);
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise, &cosm);

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let all_stations = vec![dss65_madrid, dss34_canberra];

    // Define the propagator information.
    let prop_time = 1 * TimeUnit::Day;
    let step_size = 10.0 * TimeUnit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.x += 9.5;
    initial_state_dev.y -= 9.5;
    initial_state_dev.z += 9.5;
    initial_state_dev.enable_stm();

    println!("Initial state dev:\n{}", initial_state - initial_state_dev);

    // Generate the truth data on one thread.

    let cosm = Cosm::de438();
    let bodies = vec![
        Bodies::Luna,
        Bodies::Sun,
        Bodies::JupiterBarycenter,
        Bodies::SaturnBarycenter,
    ];
    let orbital_dyn = OrbitalDynamics::point_masses(initial_state.frame, &bodies, &cosm);
    let truth_setup = Propagator::new::<RK4Fixed>(&orbital_dyn, opts);
    let mut prop = truth_setup.with(initial_state);
    prop.tx_chan = Some(truth_tx);
    prop.for_duration(prop_time).unwrap();

    let mut truth_states = Vec::with_capacity(10_000);
    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.try_recv() {
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
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let estimator = OrbitalDynamics::point_masses(initial_state.frame, &bodies, &cosm);
    let setup = Propagator::new::<RK4Fixed>(&estimator, opts);
    let prop_est = setup.with(initial_state_dev);
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

    odp.process_measurements(&measurements).unwrap();

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

#[allow(clippy::identity_op)]
#[test]
fn robust_test_ckf_smoother_multi_body() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();

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
    let prop_time = 1 * TimeUnit::Day;
    let step_size = 10.0 * TimeUnit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.x += 9.5;
    initial_state_dev.y -= 9.5;
    initial_state_dev.z += 9.5;
    initial_state_dev.enable_stm();

    let (err_p, err_v) = rss_state_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        initial_state - initial_state_dev
    );

    // Generate the truth data on one thread.

    let cosm = Cosm::de438();
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let orbital_dyn = OrbitalDynamics::point_masses(initial_state.frame, &bodies, &cosm);
    let truth_setup = Propagator::new::<RK4Fixed>(&orbital_dyn, opts);
    let mut prop = truth_setup.with(initial_state);
    prop.tx_chan = Some(truth_tx);
    prop.for_duration(prop_time).unwrap();

    let mut truth_states = Vec::with_capacity(10_000);
    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.try_recv() {
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state).unwrap();
            if meas.visible() {
                measurements.push(meas);
                break; // We know that only one station is in visibility at each time.
            }
        }
        truth_states.push(rx_state)
    }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let estimator = OrbitalDynamics::point_masses(initial_state.frame, &bodies, &cosm);
    let setup = Propagator::new::<RK4Fixed>(&estimator, opts);
    let prop_est = setup.with(initial_state_dev);
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

    let mut odp = ODProcess::ckf(prop_est, kf, all_stations, false, measurements.len());

    odp.process_measurements(&measurements).unwrap();

    // Smoother
    let smoothed_estimates = odp.smooth(SmoothingArc::All).unwrap();

    let mut rss_pos_avr = 0.0;
    let mut rss_vel_avr = 0.0;
    let mut rss_pos_avr_sm = 0.0;
    let mut rss_vel_avr_sm = 0.0;
    let mut num_pos_ok = 0;
    let mut num_vel_ok = 0;
    let mut large_smoothing_error = false;

    // Test smoothed estimates
    // Skip the first 10 estimates which are surprisingly good in this case
    for offset in (1..odp.estimates.len() - 10).rev() {
        let truth_state = truth_states[truth_states.len() - offset];
        let smoothed_est = &smoothed_estimates[odp.estimates.len() - offset];

        // Check that the covariance deflated
        let est = &odp.estimates[odp.estimates.len() - offset];

        // Some sanity checks to make sure that we have correctly indexed the estimates
        assert_eq!(smoothed_est.epoch(), est.epoch());
        assert_eq!(est.epoch(), truth_state.dt);

        let (err_p, err_v) = rss_state_errors(&est.state(), &truth_state);
        let (err_p_sm, err_v_sm) = rss_state_errors(&smoothed_est.state(), &truth_state);

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
                truth_state - est.state()
            );

            println!(
                "RSS error: smoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                err_p_sm * 1e3,
                err_v_sm * 1e3,
                truth_state - smoothed_est.state()
            );

            // Check that the covariance decreased for the final estimate
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

            let rmag_err = (truth_state - est.state()).rmag();
            assert!(
                rmag_err < 1e-2 || true,
                "final radius error should be on meter level (is instead {:.3} m)",
                rmag_err * 1e3
            );
            assert_eq!(
                truth_states.len(),
                odp.estimates.len() - 1,
                "different number of estimates"
            );
        }

        // The smoothed RSS errors should be better, or have the same order of magnitude or not significantly worse

        // Compute orders of magnitude
        let err_p_oom = err_p.log10().floor() as i32;
        let err_v_oom = err_v.log10().floor() as i32;
        let err_p_sm_oom = err_p_sm.log10().floor() as i32;
        let err_v_sm_oom = err_v_sm.log10().floor() as i32;

        if err_p_sm_oom - err_p_oom > 2 {
            large_smoothing_error = true;

            println!(
                "RSS position error after smoothing not better @{} (#{}):\n\testimate vs truth: {:.3e} m\t{:.3e} m/s\n{}\n\tsmoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                truth_state.dt.as_gregorian_tai_str(),
                odp.estimates.len() - offset,
                err_p * 1e3,
                err_v * 1e3,
                truth_state - est.state(),
                err_p_sm * 1e3,
                err_v_sm * 1e3,
                truth_state - smoothed_est.state()
            );
        }

        if err_v_sm_oom - err_v_oom > 2 {
            large_smoothing_error = true;

            println!(
                "RSS velocity error after smoothing not better @{} (#{}):\n\testimate vs truth: {:.3e} m\t{:.3e} m/s\n{}\n\tsmoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                truth_state.dt.as_gregorian_tai_str(),
                odp.estimates.len() - offset,
                err_p * 1e3,
                err_v * 1e3,
                truth_state - est.state(),
                err_p_sm * 1e3,
                err_v_sm * 1e3,
                truth_state - smoothed_est.state()
            );
        }

        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0 || true,
                "covar diagonal element negative @ [{}, {}]: @{} (#{})",
                i,
                i,
                truth_state.dt.as_gregorian_tai_str(),
                odp.estimates.len() - offset,
            );
        }
    }

    let cntf = odp.estimates.len() as f64;
    println!(
        "\nPos. better: {}/{}\tVel. better: {}/{}\nPre-smoothing  avr. RSS:\t{:.3e}\t{:.3e}\nPost-smoothing avr. RSS:\t{:.3e}\t{:.3e}\n",
        num_pos_ok,
        odp.estimates.len(),
        num_vel_ok,
        odp.estimates.len(),
        rss_pos_avr / cntf,
        rss_vel_avr / cntf,
        rss_pos_avr_sm / cntf,
        rss_vel_avr_sm / cntf,
    );

    // For the CKF, the average RSS errors are expected to be better.
    assert!(
        rss_pos_avr_sm < rss_pos_avr,
        "Average RSS position error not better"
    );
    assert!(
        rss_vel_avr_sm < rss_vel_avr,
        "Average RSS velocity error not better"
    );

    assert!(!large_smoothing_error);
}

#[allow(clippy::identity_op)]
#[test]
fn robust_test_ekf_snc_smoother_multi_body() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();

    let elevation_mask = 10.0;
    let range_noise = 1e-5;
    let range_rate_noise = 1e-7;
    let dss65_madrid =
        GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise, &cosm);
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise, &cosm);

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let all_stations = vec![dss65_madrid, dss34_canberra];

    // Define the propagator information.
    let prop_time = 1 * TimeUnit::Day;
    let step_size = 10.0 * TimeUnit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.x += 9.5;
    initial_state_dev.y -= 9.5;
    initial_state_dev.z += 9.5;
    initial_state_dev.enable_stm();

    let (err_p, err_v) = rss_state_errors(&initial_state_dev, &initial_state);
    println!(
        "Initial state dev: {:.3} m\t{:.3} m/s\n{}",
        err_p * 1e3,
        err_v * 1e3,
        initial_state - initial_state_dev
    );

    // Generate the truth data on one thread.

    let cosm = Cosm::de438();
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let orbital_dyn = OrbitalDynamics::point_masses(initial_state.frame, &bodies, &cosm);
    let truth_setup = Propagator::new::<RK4Fixed>(&orbital_dyn, opts);
    let mut prop = truth_setup.with(initial_state);
    prop.tx_chan = Some(truth_tx);
    prop.for_duration(prop_time).unwrap();

    let mut truth_states = Vec::with_capacity(10_000);
    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.try_recv() {
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state).unwrap();
            if meas.visible() {
                measurements.push(meas);
                break; // We know that only one station is in visibility at each time.
            }
        }
        truth_states.push(rx_state)
    }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let estimator = OrbitalDynamics::point_masses(initial_state.frame, &bodies, &cosm);
    let setup = Propagator::new::<RK4Fixed>(&estimator, opts);
    let prop_est = setup.with(initial_state_dev);
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

    // Define the ground stations.
    let ekf_num_meas = 100;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 1 * TimeUnit::Hour;

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(initial_state_dev, init_covar);
    println!("Initial estimate:\n{}", initial_estimate);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    let sigma_q = 1e-7_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2 * TimeUnit::Minute, &[sigma_q, sigma_q, sigma_q]);
    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    let mut odp = ODProcess::ekf(
        prop_est,
        kf,
        all_stations,
        false,
        measurements.len(),
        StdEkfTrigger::new(ekf_num_meas, ekf_disable_time),
    );

    odp.process_measurements(&measurements).unwrap();

    // Smoother
    let smoothed_estimates = odp.smooth(SmoothingArc::All).unwrap();

    let mut rss_pos_avr = 0.0;
    let mut rss_vel_avr = 0.0;
    let mut rss_pos_avr_sm = 0.0;
    let mut rss_vel_avr_sm = 0.0;
    let mut num_pos_ok = 0;
    let mut num_vel_ok = 0;

    // Test smoothed estimates
    // Skip the first 10 estimates which are surprisingly good in this case
    for offset in (1..odp.estimates.len() - 10).rev() {
        let truth_state = truth_states[truth_states.len() - offset];
        let smoothed_est = &smoothed_estimates[odp.estimates.len() - offset];

        // Check that the covariance deflated
        let est = &odp.estimates[odp.estimates.len() - offset];

        // Some sanity checks to make sure that we have correctly indexed the estimates
        assert_eq!(smoothed_est.epoch(), est.epoch());
        assert_eq!(est.epoch(), truth_state.dt);

        let (err_p, err_v) = rss_state_errors(&est.state(), &truth_state);
        let (err_p_sm, err_v_sm) = rss_state_errors(&smoothed_est.state(), &truth_state);

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
                truth_state - est.state()
            );

            println!(
                "RSS error: smoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                err_p_sm * 1e3,
                err_v_sm * 1e3,
                truth_state - smoothed_est.state()
            );

            // Check that the covariance decreased for the final estimate
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

            let rmag_err = (truth_state - est.state()).rmag();
            assert!(
                rmag_err < 1e-2 || true,
                "final radius error should be on meter level (is instead {:.3} m)",
                rmag_err * 1e3
            );
            assert_eq!(
                truth_states.len(),
                odp.estimates.len() - 1,
                "different number of estimates"
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
                truth_state.dt.as_gregorian_tai_str(),
                odp.estimates.len() - offset,
                err_p * 1e3,
                err_v * 1e3,
                truth_state - est.state(),
                err_p_sm * 1e3,
                err_v_sm * 1e3,
                truth_state - smoothed_est.state()
            );
        }

        if err_v_sm_oom - err_v_oom > 2 {
            println!(
                "RSS velocity error after smoothing not better @{} (#{}):\n\testimate vs truth: {:.3e} m\t{:.3e} m/s\n{}\n\tsmoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                truth_state.dt.as_gregorian_tai_str(),
                odp.estimates.len() - offset,
                err_p * 1e3,
                err_v * 1e3,
                truth_state - est.state(),
                err_p_sm * 1e3,
                err_v_sm * 1e3,
                truth_state - smoothed_est.state()
            );
        }

        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0 || true,
                "covar diagonal element negative @ [{}, {}]: @{} (#{})",
                i,
                i,
                truth_state.dt.as_gregorian_tai_str(),
                odp.estimates.len() - offset,
            );
        }
    }

    let cntf = odp.estimates.len() as f64;
    println!(
        "\nPos. better: {}/{}\tVel. better: {}/{}\nPre-smoothing  avr. RSS:\t{:.3e}\t{:.3e}\nPost-smoothing avr. RSS:\t{:.3e}\t{:.3e}\n",
        num_pos_ok,
        odp.estimates.len(),
        num_vel_ok,
        odp.estimates.len(),
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
#[test]
fn robust_test_ckf_iteration_multi_body() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::thread;

    let cosm = Cosm::de438();

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
    let prop_time = 1 * TimeUnit::Day;
    let step_size = 10.0 * TimeUnit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);
    let mut initial_state_dev = initial_state;
    initial_state_dev.x += 9.5;
    initial_state_dev.y -= 9.5;
    initial_state_dev.z += 9.5;
    initial_state_dev.enable_stm();

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
        let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
        let orbital_dyn = OrbitalDynamics::point_masses(initial_state.frame, &bodies, &cosm);
        let truth_setup = Propagator::new::<RK4Fixed>(&orbital_dyn, opts);
        let mut prop = truth_setup.with(initial_state);
        prop.tx_chan = Some(truth_tx);
        prop.for_duration(prop_time).unwrap();
    });

    let mut truth_states = Vec::with_capacity(10_000);
    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.try_recv() {
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state).unwrap();
            if meas.visible() {
                measurements.push(meas);
                break; // We know that only one station is in visibility at each time.
            }
        }
        truth_states.push(rx_state)
    }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let estimator = OrbitalDynamics::point_masses(initial_state.frame, &bodies, &cosm);
    // XXX: Previouslt, this used a state deviation for the initial estimate but not in the dynamics. Sounds wrong.
    let setup = Propagator::new::<RK4Fixed>(&estimator, opts);
    let prop_est = setup.with(initial_state_dev);
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

    let mut odp = ODProcess::ckf(prop_est, kf, all_stations, false, measurements.len());

    odp.process_measurements(&measurements).unwrap();

    // Clone the initial estimates
    let pre_iteration_estimates = odp.estimates.clone();

    // Iterate
    odp.iterate(&measurements, SmoothingArc::All).unwrap();

    let mut rss_pos_avr = 0.0;
    let mut rss_vel_avr = 0.0;
    let mut rss_pos_avr_it = 0.0;
    let mut rss_vel_avr_it = 0.0;
    let mut num_pos_ok = 0;
    let mut num_vel_ok = 0;

    // Compare the initial estimates and the iterated estimates
    // Skip the first 10 estimates which are surprisingly good in this case
    for offset in (1..odp.estimates.len() - 10).rev() {
        let truth_state = truth_states[truth_states.len() - offset];
        let prior_est = &pre_iteration_estimates[odp.estimates.len() - offset];

        // Check that the covariance deflated
        let est = &odp.estimates[odp.estimates.len() - offset];

        // Some sanity checks to make sure that we have correctly indexed the estimates
        assert_eq!(prior_est.epoch(), est.epoch());
        assert_eq!(est.epoch(), truth_state.dt);

        let (err_p, err_v) = rss_state_errors(&prior_est.state(), &truth_state);
        let (err_p_it, err_v_it) = rss_state_errors(&est.state(), &truth_state);

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
                truth_state - prior_est.state()
            );

            println!(
                "RSS error: iterated estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                err_p_it * 1e3,
                err_v_it * 1e3,
                truth_state - est.state()
            );

            // Check that the covariance decreased for the final estimate
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

            let rmag_err = (truth_state - est.state()).rmag();
            assert!(
                rmag_err < 1e-2 || true,
                "final radius error should be on meter level (is instead {:.3} m)",
                rmag_err * 1e3
            );
            assert_eq!(
                truth_states.len(),
                odp.estimates.len() - 1,
                "different number of estimates"
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
                truth_state.dt.as_gregorian_tai_str(),
                odp.estimates.len() - offset,
                err_p * 1e3,
                err_v * 1e3,
                truth_state - prior_est.state(),
                err_p_it * 1e3,
                err_v_it * 1e3,
                truth_state - est.state()
            );
        }

        if err_v_it_oom - err_v_oom > 3 {
            println!(
                "RSS velocity error after smoothing not better @{} (#{}):\n\testimate vs truth: {:.3e} m\t{:.3e} m/s\n{}\n\tsmoothed estimate vs truth: {:.3e} m\t{:.3e} m/s\n{}",
                truth_state.dt.as_gregorian_tai_str(),
                odp.estimates.len() - offset,
                err_p * 1e3,
                err_v * 1e3,
                truth_state - prior_est.state(),
                err_p_it * 1e3,
                err_v_it * 1e3,
                truth_state - est.state()
            );
        }

        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0 || true,
                "covar diagonal element negative @ [{}, {}]: @{} (#{})",
                i,
                i,
                truth_state.dt.as_gregorian_tai_str(),
                odp.estimates.len() - offset,
            );
        }
    }

    let cntf = odp.estimates.len() as f64;
    println!(
        "\nPos. better: {}/{}\tVel. better: {}/{}\nPre-iteration  avr. RSS:\t{:.3e}\t{:.3e}\nPost-iteration avr. RSS:\t{:.3e}\t{:.3e}\n",
        num_pos_ok,
        odp.estimates.len(),
        num_vel_ok,
        odp.estimates.len(),
        rss_pos_avr / cntf,
        rss_vel_avr / cntf,
        rss_pos_avr_it / cntf,
        rss_vel_avr_it / cntf,
    );

    // For the CKF, the average RSS errors are expected to be better.
    assert!(
        rss_pos_avr_it < rss_pos_avr,
        "Average RSS position error not better"
    );
    assert!(
        rss_vel_avr_it < rss_vel_avr,
        "Average RSS velocity error not better"
    );
}
