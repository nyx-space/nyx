extern crate csv;
extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use self::hifitime::{Epoch, SECONDS_PER_DAY};
use self::na::{Matrix2, Matrix3, Matrix6, Vector2, Vector6};
use self::nyx::celestia::{bodies, Cosm, Geoid, State};
use self::nyx::dynamics::celestial::{CelestialDynamics, CelestialDynamicsStm};
use self::nyx::od::ui::*;
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};

#[test]
fn multi_body_ckf_perfect_stations() {
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
    let (truth_tx, truth_rx): (Sender<State<Geoid>>, Receiver<State<Geoid>>) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(bodies::EARTH);
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state =
        State::from_keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, earth_geoid);

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let cosm = Cosm::from_xb("./de438s");
        let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
        let mut dynamics = CelestialDynamics::new(initial_state, bodies, &cosm);
        let mut prop = Propagator::new::<RK4Fixed>(&mut dynamics, &opts);
        prop.tx_chan = Some(&truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.recv() {
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state);
            if meas.visible() {
                measurements.push((rx_state.dt, meas));
                break;
            }
        }
    }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let opts_est = PropOpts::with_fixed_step(step_size);
    let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
    let mut estimator = CelestialDynamicsStm::new(initial_state, bodies, &cosm);
    let mut prop_est = Propagator::new::<RK4Fixed>(&mut estimator, &opts_est);
    let covar_radius = 1.0e-3_f64.powi(2);
    let covar_velocity = 1.0e-6_f64.powi(2);
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius,
        covar_radius,
        covar_radius,
        covar_velocity,
        covar_velocity,
        covar_velocity,
    ));

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(dt, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise =
        Matrix2::from_diagonal(&Vector2::new(15e-3_f64.powi(2), 1e-5_f64.powi(2)));

    // Define the process noise in order to define how many variables of the EOMs are accelerations
    // (this is required due to the many compile-time matrix size verifications)
    let process_noise = Matrix3::zeros();
    // But we disable the state noise compensation / process noise by setting the delta time to None
    let process_noise_dt = None;

    let mut ckf = KF::initialize(
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
    let mut last_est = None;
    for (no, est) in odp.estimates.iter().enumerate() {
        assert_eq!(
            est.predicted, false,
            "estimate {} should not be a prediction",
            no
        );
        assert!(
            est.state.norm() < 1e-6,
            "estimate error should be zero (perfect dynamics) ({:e})",
            est.state.norm()
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

        last_est = Some(est);
    }

    // NOTE: We do not check whether the covariance has deflated because it is possible that it inflates before deflating.
    // The filter in multibody dynamics has been validated against JPL tools using a proprietary scenario.
    let est = last_est.unwrap();
    assert!(est.state.norm() < 1e-8);
    assert!(est.covar.norm() < 1e-5);
}

#[test]
fn multi_body_ckf_covar_map() {
    // For this test, we're only enabling one station so we can check that the covariance inflates between visibility passes.
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::{io, thread};

    // Define the ground stations.
    let elevation_mask = 0.0;
    let range_noise = 0.0;
    let range_rate_noise = 0.0;
    let dss13_goldstone =
        GroundStation::dss13_goldstone(elevation_mask, range_noise, range_rate_noise);
    let all_stations = vec![dss13_goldstone];

    // Define the propagator information.
    let prop_time = SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx): (Sender<State<Geoid>>, Receiver<State<Geoid>>) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(bodies::EARTH);
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state =
        State::from_keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, earth_geoid);

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let cosm = Cosm::from_xb("./de438s");
        let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
        let mut dynamics = CelestialDynamics::new(initial_state, bodies, &cosm);
        let mut prop = Propagator::new::<RK4Fixed>(&mut dynamics, &opts);
        prop.tx_chan = Some(&truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.recv() {
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state);
            if meas.visible() {
                measurements.push((rx_state.dt, meas));
                break;
            }
        }
    }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let opts_est = PropOpts::with_fixed_step(step_size);
    let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
    let mut estimator = CelestialDynamicsStm::new(initial_state, bodies, &cosm);

    let (pest_tx, pest_rx): (
        Sender<(State<Geoid>, Matrix6<f64>)>,
        Receiver<(State<Geoid>, Matrix6<f64>)>,
    ) = mpsc::channel();

    let mut prop_est = Propagator::new::<RK4Fixed>(&mut estimator, &opts_est);
    prop_est.tx_chan = Some(&pest_tx);
    let covar_radius = 1.0e-3_f64.powi(2);
    let covar_velocity = 1.0e-6_f64.powi(2);
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius,
        covar_radius,
        covar_radius,
        covar_velocity,
        covar_velocity,
        covar_velocity,
    ));

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(dt, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise =
        Matrix2::from_diagonal(&Vector2::new(15e-3_f64.powi(2), 1e-5_f64.powi(2)));

    // Define the process noise in order to define how many variables of the EOMs are accelerations
    // (this is required due to the many compile-time matrix size verifications)
    let process_noise = Matrix3::zeros();
    // But we disable the state noise compensation / process noise by setting the delta time to None
    let process_noise_dt = None;

    let mut ckf = KF::initialize(
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

    let rtn = odp.process_measurements_covar(&pest_rx, &measurements);

    assert!(rtn.is_none(), "kf failed");

    let mut last_est = None;
    let mut num_pred = 0_u32;
    for est in odp.estimates.iter() {
        if est.predicted {
            num_pred += 1;
        } else {
            // Only check that the covariance is low IF this isn't a predicted estimate
            assert!(
                est.state.norm() < 1e-6,
                "estimate error should be zero (perfect dynamics) ({:e})",
                est.state.norm()
            );
        }
        last_est = Some(est);
    }

    // Note that we check the residuals separately from the estimates because we have many predicted estimates which do not have any associated residuals.
    for res in odp.residuals.iter() {
        assert!(
            res.postfit.norm() < 1e-12,
            "postfit should be zero (perfect dynamics) ({:e})",
            res
        );
    }

    assert!(num_pred > 0, "no predicted estimates");

    // NOTE: We do not check whether the covariance has deflated because it is possible that it inflates before deflating.
    // The filter in multibody dynamics has been validated against JPL Monte using a proprietary scenario.
    let est = last_est.unwrap();

    let mut wtr = csv::Writer::from_writer(io::stdout());
    wtr.serialize(est.clone())
        .expect("could not write to stdout");

    assert!(est.state.norm() < 1e-8);
    assert!(est.covar.norm() < 1e-4);
}
