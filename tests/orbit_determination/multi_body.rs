extern crate csv;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use self::nyx::celestia::{Bodies, Cosm, Orbit};
use self::nyx::dimensions::{Matrix2, Matrix6, Vector2, Vector6};
use self::nyx::dynamics::orbital::OrbitalDynamics;
use self::nyx::od::ui::*;
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use self::nyx::time::{Epoch, TimeUnit};
use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};

#[allow(clippy::identity_op)]
#[test]
fn multi_body_ckf_perfect_stations() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::io;

    let cosm = Cosm::de438();

    // Define the ground stations.
    let elevation_mask = 0.0;
    let range_noise = 0.0;
    let range_rate_noise = 0.0;
    let dss65_madrid =
        GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise, cosm.clone());
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise, cosm.clone());
    let dss13_goldstone =
        GroundStation::dss13_goldstone(elevation_mask, range_noise, range_rate_noise, cosm.clone());
    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

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

    // Generate the truth data on one thread.

    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let orbital_dyn = OrbitalDynamics::point_masses(initial_state.frame, &bodies, cosm);
    let setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);
    let mut prop = setup.with(initial_state);
    prop.tx_chan = Some(truth_tx);
    prop.for_duration(prop_time).unwrap();

    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.try_recv() {
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state).unwrap();
            if meas.visible() {
                measurements.push(meas);
                break;
            }
        }
    }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let mut init_est_state = initial_state;
    init_est_state.enable_stm();
    let prop_est = setup.with(init_est_state);
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
    let initial_estimate = KfEstimate::from_covar(initial_state, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise =
        Matrix2::from_diagonal(&Vector2::new(15e-3_f64.powi(2), 1e-5_f64.powi(2)));

    let ckf = KF::no_snc(initial_estimate, measurement_noise);

    let mut odp = ODProcess::ckf(prop_est, ckf, all_stations, false, measurements.len());

    odp.process_measurements(&measurements).unwrap();

    let mut wtr = csv::Writer::from_writer(io::stdout());
    let mut printed = false;
    let mut last_est = None;
    for (no, est) in odp.estimates.iter().enumerate() {
        if no == 0 {
            // Skip the first estimate which is the initial estimate provided by user
            continue;
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
            est.state_deviation().norm() < 1e-3,
            "estimate error should be very good (perfect dynamics) ({:e})",
            est.state_deviation().norm()
        );

        if !printed {
            wtr.serialize(est.clone())
                .expect("could not write to stdout");
            printed = true;
        }

        last_est = Some(est);
    }

    for res in &odp.residuals {
        assert!(
            res.postfit.norm() < 1e-5,
            "postfit should be zero (perfect dynamics) ({:e})",
            res
        );
    }

    // NOTE: We do not check whether the covariance has deflated because it is possible that it inflates before deflating.
    // The filter in multibody dynamics has been validated against JPL tools using a proprietary scenario.
    let est = last_est.unwrap();
    assert!(est.state_deviation().norm() < 5e-5);
    assert!(est.covar.norm() < 1e-5);
}

#[allow(clippy::identity_op)]
#[test]
fn multi_body_ckf_covar_map() {
    // For this test, we're only enabling one station so we can check that the covariance inflates between visibility passes.
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::io;

    let cosm = Cosm::de438();

    // Define the ground stations.
    let elevation_mask = 0.0;
    let range_noise = 0.0;
    let range_rate_noise = 0.0;
    let dss13_goldstone =
        GroundStation::dss13_goldstone(elevation_mask, range_noise, range_rate_noise, cosm.clone());
    let all_stations = vec![dss13_goldstone];

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

    // Generate the truth data on one thread.

    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let orbital_dyn = OrbitalDynamics::point_masses(initial_state.frame, &bodies, cosm);
    let setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);
    let mut prop = setup.with(initial_state);
    prop.tx_chan = Some(truth_tx);
    prop.for_duration(prop_time).unwrap();

    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_state) = truth_rx.try_recv() {
        for station in all_stations.iter() {
            let meas = station.measure(&rx_state).unwrap();
            if meas.visible() {
                measurements.push(meas);
                break;
            }
        }
    }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let mut init_est_state = initial_state;
    init_est_state.enable_stm();

    let prop_est = setup.with(init_est_state);
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
    let initial_estimate = KfEstimate::from_covar(initial_state, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise =
        Matrix2::from_diagonal(&Vector2::new(15e-3_f64.powi(2), 1e-5_f64.powi(2)));

    let ckf = KF::no_snc(initial_estimate, measurement_noise);

    let mut odp = ODProcess::ckf(prop_est, ckf, all_stations, false, measurements.len());

    odp.process_measurements(&measurements).unwrap();

    let mut last_est = None;
    let mut num_pred = 0_u32;
    for est in odp.estimates.iter() {
        if est.predicted {
            num_pred += 1;
        } else {
            // Only check that the covariance is low IF this isn't a predicted estimate
            assert!(
                est.state_deviation().norm() < 1e-3,
                "estimate error should be zero (perfect dynamics) ({:e})",
                est.state_deviation().norm()
            );
        }
        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0,
                "covar diagonal element negative @ [{}, {}]",
                i,
                i
            );
        }
        last_est = Some(est);
    }

    // Note that we check the residuals separately from the estimates because we have many predicted estimates which do not have any associated residuals.
    for res in odp.residuals.iter() {
        assert!(
            res.postfit.norm() < 1e-5,
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

    println!("{:.2e}", est.state_deviation().norm());
    println!("{:.2e}", est.covar.norm());

    assert!(est.state_deviation().norm() < 5e-5);
    assert!(est.covar.norm() < 1e-4);
}
