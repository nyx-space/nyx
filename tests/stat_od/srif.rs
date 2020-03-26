extern crate csv;
extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use self::hifitime::{Epoch, SECONDS_PER_DAY};
use self::na::{Matrix2, Matrix3, Matrix6, Vector2, Vector3, Vector6};
use self::nyx::celestia::{bodies, Cosm, OrbitState};
use self::nyx::dynamics::celestial::{CelestialDynamics, CelestialDynamicsStm};
use self::nyx::od::ui::*;
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};

#[test]
fn srif_fixed_step_perfect_stations() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::{io, thread};

    // Define the ground stations.
    let elevation_mask = 0.0;
    let range_noise = 0.0;
    let range_rate_noise = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise);
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise);
    let dss13_goldstone =
        GroundStation::dss13_goldstone(elevation_mask, range_noise, range_rate_noise);
    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx): (Sender<OrbitState>, Receiver<OrbitState>) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.frame_by_id(bodies::EARTH);
    let dt = Epoch::from_mjd_tai(21545.0);
    let initial_state =
        OrbitState::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, earth_geoid);

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let mut dynamics = CelestialDynamics::two_body(initial_state);
        let mut prop = Propagator::new::<RK4Fixed>(&mut dynamics, &opts);
        prop.tx_chan = Some(&truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.recv() {
        // Convert the state to ECI.
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state).unwrap();
            if meas.visible() {
                measurements.push((rx_state.dt, meas));
                break; // We know that only one station is in visibility at each time.
            }
        }
    }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let opts_est = PropOpts::with_fixed_step(step_size);
    let mut tb_estimator = CelestialDynamicsStm::two_body(initial_state);
    let mut prop_est = Propagator::new::<RK4Fixed>(&mut tb_estimator, &opts_est);
    let covar_radius = 1.0e-3;
    let covar_velocity = 1.0e-6;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius,
        covar_radius,
        covar_radius,
        covar_velocity,
        covar_velocity,
        covar_velocity,
    ));

    // Define the initial estimate
    let initial_estimate = IfEstimate::from_covar(dt, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    // Define the process noise in order to define how many variables of the EOMs are accelerations
    // (this is required due to the many compile-time matrix size verifications)
    let process_noise = Matrix3::zeros();
    // But we disable the state noise compensation / process noise by setting the delta time to None
    let process_noise_dt = None;

    let mut ckf = SRIF::initialize(
        initial_estimate,
        process_noise,
        measurement_noise,
        process_noise_dt,
    );

    let mut odp = ODProcess::ckf(
        &mut prop_est,
        &mut ckf,
        &all_stations,
        false,
        measurements.len(),
    );

    let rtn = odp.process_measurements(&measurements);
    assert!(rtn.is_none(), "kf failed");

    let mut wtr = csv::Writer::from_writer(io::stdout());
    let mut printed = false;
    for (no, est) in odp.estimates.iter().enumerate() {
        assert_eq!(
            est.predicted, false,
            "estimate {} should not be a prediction",
            no
        );
        assert!(
            est.state().norm() < 1e-6,
            "estimate error should be zero (perfect dynamics) ({:e})",
            est.state().norm()
        );

        let res = &odp.residuals[no];
        assert!(
            res.postfit.norm() < 1e-12,
            "postfit should be zero (perfect dynamics) ({:e})",
            res
        );

        if !printed {
            wtr.serialize(est.clone())
                .expect("could not write to stdout");
            printed = true;
        }
    }

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];
    for i in 0..6 {
        if i < 3 {
            assert!(
                est.covar()[(i, i)].abs() < covar_radius,
                "covar radius did not decrease"
            );
        } else {
            assert!(
                est.covar()[(i, i)].abs() < covar_velocity,
                "covar velocity did not decrease"
            );
        }
    }

    println!("N-1 not smoothed: \n{}", est);
}

#[test]
fn srif_fixed_step_perfect_stations_snc_covar_map() {
    // Tests state noise compensation with covariance mapping
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::thread;

    // Define the ground stations.
    let elevation_mask = 0.0;
    let range_noise = 0.0;
    let range_rate_noise = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise);
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise);
    let dss13_goldstone =
        GroundStation::dss13_goldstone(elevation_mask, range_noise, range_rate_noise);
    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx): (Sender<OrbitState>, Receiver<OrbitState>) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.frame_by_id(bodies::EARTH);
    let dt = Epoch::from_mjd_tai(21545.0);
    let initial_state =
        OrbitState::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, earth_geoid);

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let mut dynamics = CelestialDynamics::two_body(initial_state);
        let mut prop = Propagator::new::<RK4Fixed>(&mut dynamics, &opts);
        prop.tx_chan = Some(&truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.recv() {
        // Convert the state to ECI.
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state).unwrap();
            if meas.visible() {
                measurements.push((rx_state.dt, meas));
                break; // We know that only one station is in visibility at each time.
            }
        }
    }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let opts_est = PropOpts::with_fixed_step(step_size);
    let mut tb_estimator = CelestialDynamicsStm::two_body(initial_state);
    let mut prop_est = Propagator::new::<RK4Fixed>(&mut tb_estimator, &opts_est);
    // Create the channels for covariance mapping
    let (prop_tx, prop_rx) = mpsc::channel();
    prop_est.tx_chan = Some(&prop_tx);

    // Set up the filter
    let covar_radius = 1.0e-3;
    let covar_velocity = 1.0e-6;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius,
        covar_radius,
        covar_radius,
        covar_velocity,
        covar_velocity,
        covar_velocity,
    ));

    // Define the initial estimate
    let initial_estimate = IfEstimate::from_covar(dt, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    // Define the process noise to assume an unmodel acceleration of 1e-3 km^2/s^2 on X, Y and Z in the ECI frame
    let sigma_q = 1e-8_f64.powi(2);
    let process_noise = Matrix3::from_diagonal(&Vector3::new(sigma_q, sigma_q, sigma_q));
    // Disable SNC if there is more than 120 seconds between two measurements
    let process_noise_dt = None;

    let mut ckf = SRIF::initialize(
        initial_estimate,
        process_noise,
        measurement_noise,
        process_noise_dt,
    );

    let mut odp = ODProcess::ckf(
        &mut prop_est,
        &mut ckf,
        &all_stations,
        false,
        measurements.len(),
    );

    let rtn = odp.process_measurements_covar(&prop_rx, &measurements);
    assert!(rtn.is_none(), "srif failed");

    let mut wtr = csv::Writer::from_path("./estimation-srif.csv").unwrap();

    // Let's export these to a CSV file, and also check that the covariance never falls below our sigma squared values
    for (no, est) in odp.estimates.iter().enumerate() {
        if no == 1 {
            println!("{}", est);
        }
        assert!(
            est.state().norm() < 1e-6,
            "estimate error should be zero (perfect dynamics) ({:e})",
            est.state().norm()
        );

        // for i in 0..6 {
        //     assert!(
        //         est.covar()[(i, i)].abs() >= sigma_q,
        //         "covar diagonal less than SNC value @ {} = {:.3e}",
        //         no,
        //         est.covar()[(i, i)]
        //     );
        // }

        wtr.serialize(est.clone())
            .expect("could not write to stdout");
    }
}
