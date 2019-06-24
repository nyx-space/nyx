extern crate csv;
extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::hifitime::datetime::*;
use self::hifitime::julian::*;
use self::hifitime::SECONDS_PER_DAY;
use self::na::{Matrix2, Matrix6, Vector2, Vector6, VectorN, U36, U42, U6};
use self::nyx::celestia::{Cosm, Geoid, State};
use self::nyx::dynamics::celestial::{TwoBody, TwoBodyWithStm};
use self::nyx::dynamics::Dynamics;
use self::nyx::io::cosmo::Cosmographia;
use self::nyx::od::kalman::{Estimate, KF};
use self::nyx::od::ranging::GroundStation;
use self::nyx::od::Measurement;
use self::nyx::propagators::error_ctrl::{LargestError, RSSStepPV};
use self::nyx::propagators::{PropOpts, Propagator, RK89};
use std::f64::EPSILON;
use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};
use std::thread;

fn main() {
    // Define the output file with all the estimates and covariance
    let mut est_csv = csv::Writer::from_path("estimation.csv").expect("could not open estimation.csv for writing");

    // Define the ground stations.
    let num_meas_for_ekf = 15;
    let elevation_mask = 0.0;
    let range_noise = 0.0;
    let range_rate_noise = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise);
    let dss34_canberra = GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise);
    let dss13_goldstone = GroundStation::dss13_goldstone(elevation_mask, range_noise, range_rate_noise);
    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = 10.0 * SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = PropOpts::with_fixed_step(step_size, RSSStepPV {});

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx): (Sender<(f64, Vector6<f64>)>, Receiver<(f64, Vector6<f64>)>) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let dt = ModifiedJulian::from_instant(Datetime::new(2018, 2, 27, 0, 0, 0, 0).expect("ugh?").into_instant());
    let cosm = Cosm::from_xb("../../de438s");
    let earth_geoid = cosm.geoid_from_id(3).unwrap();
    let initial_state = State::from_keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, earth_geoid.clone());

    // Initialize the XYZV exporter
    let mut outfile = Cosmographia::from_path("truth.xyzv".to_owned());

    // Generate the truth data on one thread.
    let geoid = earth_geoid.clone();
    thread::spawn(move || {
        let mut dynamics = TwoBody::from_state_vec(initial_state.to_cartesian_vec(), geoid);
        let mut prop = Propagator::new::<RK89>(&mut dynamics, &opts.clone());
        prop.tx_chan = Some(&truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    // Receive the states on the main thread, and populate the measurement channel.
    loop {
        match truth_rx.recv() {
            Ok((t, state_vec)) => {
                let this_dt =
                    ModifiedJulian::from_instant(dt.into_instant() + Instant::from_precise_seconds(t, Era::Present).duration());
                let rx_state = State::<Geoid>::from_cartesian_vec(&state_vec, this_dt, earth_geoid.clone());
                // Export state
                outfile.append(rx_state);
                // Check visibility
                for station in all_stations.iter() {
                    let meas = station.measure(rx_state, this_dt.into_instant());
                    if meas.visible() {
                        // XXX: Instant does not implement Eq, only PartialEq, so can't use it as an index =(
                        measurements.push((t, meas));
                        break; // We know that only one station is in visibility at each time.
                    }
                }
            }
            Err(_) => {
                break;
            }
        }
    }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let opts_est = PropOpts::with_fixed_step(step_size, LargestError {});
    let (est_tx, est_rx): (Sender<(f64, VectorN<f64, U42>)>, Receiver<(f64, VectorN<f64, U42>)>) = mpsc::channel();
    // WARNING: The dynamics must be created after the channel is created due to lifetime issues.
    let mut tb_estimator = TwoBodyWithStm::from_state(initial_state);
    let mut prop_est = Propagator::new::<RK89>(&mut tb_estimator, &opts_est);
    prop_est.tx_chan = Some(&est_tx);

    // And propagate
    prop_est.until_time_elapsed(prop_time);

    let covar_radius = 1.0e-6;
    let covar_velocity = 1.0e-6;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius,
        covar_radius,
        covar_radius,
        covar_velocity,
        covar_velocity,
        covar_velocity,
    ));

    let initial_estimate = Estimate {
        state: Vector6::zeros(),
        covar: init_covar,
        stm: prop_est.dynamics.stm.clone(),
        predicted: false,
    };

    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));
    let mut kf = KF::initialize(initial_estimate, measurement_noise);

    println!("Will process {} measurements", measurements.len());

    let mut prev_dt = dt;
    let mut meas_no = 0;
    loop {
        match est_rx.recv() {
            Ok((t, full_state_vec)) => {
                let state_vec = full_state_vec.fixed_rows::<U6>(0).into_owned();
                let stm_as_vec = full_state_vec.fixed_rows::<U36>(6).into_owned();
                let stm = Matrix6::from_row_slice(&stm_as_vec.as_slice());

                // Update the time
                let this_dt = ModifiedJulian::from_instant(
                    prev_dt.into_instant() + Instant::from_precise_seconds(step_size, Era::Present).duration(),
                );
                if prev_dt == this_dt {
                    panic!("already handled that time: {}", prev_dt);
                }
                prev_dt = this_dt;

                // Start by setting the next STM
                kf.update_stm(stm.clone());
                // Check to see if we have a measurement at this time
                let (meas_time, real_meas) = measurements[meas_no];

                if t == meas_time {
                    // We've got a measurement here, so let's get ready to move onto the next measurement
                    meas_no += 1;
                    // Get the computed observation
                    let rx_state = State::<Geoid>::from_cartesian_vec(&state_vec, this_dt, earth_geoid.clone());
                    let mut still_empty = true;
                    for station in all_stations.iter() {
                        let computed_meas = station.measure(rx_state, this_dt.into_instant());
                        if computed_meas.visible() {
                            kf.update_h_tilde(*computed_meas.sensitivity());
                            let mut latest_est = kf
                                .measurement_update(*real_meas.observation(), *computed_meas.observation())
                                .expect("wut?");
                            still_empty = false;
                            assert_eq!(latest_est.predicted, false, "estimate should not be a prediction");
                            assert!(
                                latest_est.state.norm() < EPSILON,
                                "estimate error should be zero (perfect dynamics)"
                            );
                            if kf.ekf {
                                // It's an EKF, so let's update the state in the dynamics.
                                // let now = tb_estimator.time(); // Needed because we can't do a mutable borrow while doing an immutable one too.
                                let new_state = state_vec + latest_est.state;
                                prop_est.dynamics.set_state(
                                    t,
                                    &VectorN::<f64, U42>::from_iterator(new_state.iter().chain(stm_as_vec.iter()).cloned()),
                                );
                                // tb_estimator.two_body_dyn.set_state(t, &new_state);
                            }
                            // We want to show the 3 sigma covariance, so le'ts multiply the covariance by 3
                            latest_est.covar *= 3.0;
                            // Let's export this estimation to the CSV file
                            est_csv.serialize(latest_est).expect("could not write to stdout");
                            break; // We know that only one station is in visibility at each time.
                        }
                    }
                    if still_empty {
                        // We're doing perfect everything, so we should always be in visibility if there is a measurement
                        panic!("T+{} : not in visibility", this_dt);
                    }
                    // If we've reached the last measurement, let's break this loop.
                    if meas_no == measurements.len() {
                        break;
                    }
                } else {
                    // No measurement, so let's do a time update only.
                    let mut latest_est = kf.time_update().expect("time update failed");
                    // We want to show the 3 sigma covariance, so le'ts multiply the covariance by 3
                    latest_est.covar *= 3.0;
                    // Let's export this estimation to the CSV file
                    est_csv.serialize(latest_est).expect("could not write to stdout");
                }
                if meas_no > num_meas_for_ekf && !kf.ekf {
                    println!("switched to EKF");
                    kf.ekf = true;
                }
            }
            Err(_) => {
                break;
            }
        }
    }
    println!("DONE");
}
