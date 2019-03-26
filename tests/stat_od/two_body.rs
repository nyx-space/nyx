extern crate csv;
extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::hifitime::instant::*;
use self::hifitime::julian::*;
use self::hifitime::SECONDS_PER_DAY;
use self::na::{Matrix2, Matrix2x6, Matrix6, Vector2, Vector6, U6};
use self::nyx::celestia::{State, EARTH, ECI};
use self::nyx::dynamics::celestial::{TwoBody, TwoBodyWithDualStm, TwoBodyWithStm};
use self::nyx::dynamics::Dynamics;
use self::nyx::od::kalman::{Estimate, FilterError, KF};
use self::nyx::od::ranging::GroundStation;
use self::nyx::od::Measurement;
use self::nyx::propagators::error_ctrl::{LargestError, RSSStepPV};
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use std::f64::EPSILON;
use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};

macro_rules! f64_nil {
    ($x:expr, $msg:expr) => {
        assert!($x.abs() < EPSILON, $msg)
    };
}

#[test]
fn empty_estimate() {
    let empty = Estimate::<U6>::empty();
    f64_nil!(empty.state.norm(), "expected state norm to be nil");
    f64_nil!(empty.covar.norm(), "expected covar norm to be nil");
    f64_nil!(empty.stm.norm(), "expected STM norm to be nil");
    assert_eq!(empty.predicted, true, "expected predicted to be true");
}

#[test]
fn csv_serialize_empty_estimate() {
    use std::io;
    let empty = Estimate::<U6>::empty();
    let mut wtr = csv::Writer::from_writer(io::stdout());
    wtr.serialize(empty).expect("could not write to stdout");
}

#[test]
fn filter_errors() {
    let initial_estimate = Estimate::<U6>::empty();
    let measurement_noise = Matrix2::zeros();
    let real_obs = Vector2::zeros();
    let computed_obs = Vector2::zeros();
    let sensitivity = Matrix2x6::zeros();
    let stm = Matrix6::zeros();

    let mut ckf = KF::initialize(initial_estimate, measurement_noise);
    match ckf.time_update() {
        Ok(_) => panic!("expected the time update to fail"),
        Err(e) => {
            assert_eq!(e, FilterError::StateTransitionMatrixNotUpdated);
        }
    }
    match ckf.measurement_update(real_obs, computed_obs) {
        Ok(_) => panic!("expected the measurement update to fail"),
        Err(e) => {
            assert_eq!(e, FilterError::StateTransitionMatrixNotUpdated);
        }
    }
    ckf.update_stm(stm);
    match ckf.measurement_update(real_obs, computed_obs) {
        Ok(_) => panic!("expected the measurement update to fail"),
        Err(e) => {
            assert_eq!(e, FilterError::SensitivityNotUpdated);
        }
    }
    ckf.update_h_tilde(sensitivity);
    match ckf.measurement_update(real_obs, computed_obs) {
        Ok(_) => panic!("expected the measurement update to fail"),
        Err(e) => {
            assert_eq!(e, FilterError::GainIsSingular);
        }
    }
}

#[test]
fn ckf_fixed_step_perfect_stations_std() {
    use std::{io, thread};

    // Define the ground stations.
    let elevation_mask = 0.0;
    let range_noise = 0.0;
    let range_rate_noise = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise);
    let dss34_canberra = GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise);
    let dss13_goldstone = GroundStation::dss13_goldstone(elevation_mask, range_noise, range_rate_noise);
    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = PropOpts::with_fixed_step(step_size, RSSStepPV {});

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx): (Sender<(f64, Vector6<f64>)>, Receiver<(f64, Vector6<f64>)>) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let dt = ModifiedJulian { days: 21545.0 };
    let initial_state = State::from_keplerian_eci(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt);

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let mut dyn = TwoBody::from_state_vec::<EARTH>(initial_state.to_cartesian_vec());
        let mut prop = Propagator::new::<RK4Fixed>(&mut dyn, &opts.clone());
        prop.tx_chan = Some(&truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    // Receive the states on the main thread, and populate the measurement channel.
    let mut prev_dt = dt.into_instant();
    loop {
        match truth_rx.recv() {
            Ok((t, state_vec)) => {
                // Convert the state to ECI.
                let this_dt =
                    ModifiedJulian::from_instant(dt.into_instant() + Instant::from_precise_seconds(t, Era::Present).duration());
                let rx_state = State::from_cartesian_vec::<EARTH, ModifiedJulian>(&state_vec, this_dt, ECI {});
                for station in all_stations.iter() {
                    let meas = station.measure(rx_state, this_dt.into_instant());
                    if meas.visible() {
                        // XXX: Instant does not implement Eq, only PartialEq, so can't use it as an index =(
                        let delta = (prev_dt - this_dt.into_instant().duration()).duration().as_secs();
                        measurements.push((delta, meas));
                        prev_dt = this_dt.into_instant();
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
    let mut tb_estimator = TwoBodyWithStm::from_state::<EARTH, ECI>(initial_state);
    let mut prop_est = Propagator::new::<RK4Fixed>(&mut tb_estimator, &opts_est);
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
    let mut ckf = KF::initialize(initial_estimate, measurement_noise);

    println!("Will process {} measurements", measurements.len());

    let mut prev_dt = dt;
    let mut printed = false;

    let mut wtr = csv::Writer::from_writer(io::stdout());

    for (duration, real_meas) in measurements.iter() {
        // Propagate the dynamics to the measurement, and then start the filter.
        let delta_time = (*duration) as f64;
        prop_est.until_time_elapsed(delta_time);
        // Update the STM of the KF
        ckf.update_stm(prop_est.dynamics.stm.clone());
        // Get the computed observation
        let this_dt = ModifiedJulian::from_instant(
            prev_dt.into_instant() + Instant::from_precise_seconds(delta_time, Era::Present).duration(),
        );
        if prev_dt == this_dt {
            panic!("already handled that time: {}", prev_dt);
        }
        prev_dt = this_dt;
        let rx_state =
            State::from_cartesian_vec::<EARTH, ModifiedJulian>(&prop_est.dynamics.two_body_dyn.state(), this_dt, ECI {});
        let mut still_empty = true;
        for station in all_stations.iter() {
            let computed_meas = station.measure(rx_state, this_dt.into_instant());
            if computed_meas.visible() {
                ckf.update_h_tilde(*computed_meas.sensitivity());
                let latest_est = ckf
                    .measurement_update(*real_meas.observation(), *computed_meas.observation())
                    .expect("wut?");
                still_empty = false;
                assert_eq!(latest_est.predicted, false, "estimate should not be a prediction");
                assert!(
                    latest_est.state.norm() < EPSILON,
                    "estimate error should be zero (perfect dynamics)"
                );
                if !printed {
                    wtr.serialize(latest_est).expect("could not write to stdout");
                    printed = true;
                }
                break; // We know that only one station is in visibility at each time.
            }
        }
        if still_empty {
            // We're doing perfect everything, so we should always be in visibility
            panic!("T+{} : not in visibility", this_dt);
        }
    }
}

#[test]
fn ekf_fixed_step_perfect_stations() {
    use std::thread;

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
    let prop_time = SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = PropOpts::with_fixed_step(step_size, RSSStepPV {});

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx): (Sender<(f64, Vector6<f64>)>, Receiver<(f64, Vector6<f64>)>) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let dt = ModifiedJulian { days: 21545.0 };
    let initial_state = State::from_keplerian_eci(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt);

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let mut dyn = TwoBody::from_state_vec::<EARTH>(initial_state.to_cartesian_vec());
        let mut prop = Propagator::new::<RK4Fixed>(&mut dyn, &opts.clone());
        prop.tx_chan = Some(&truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    // Receive the states on the main thread, and populate the measurement channel.
    let mut prev_dt = dt.into_instant();
    loop {
        match truth_rx.recv() {
            Ok((t, state_vec)) => {
                // Convert the state to ECI.
                let this_dt =
                    ModifiedJulian::from_instant(dt.into_instant() + Instant::from_precise_seconds(t, Era::Present).duration());
                let rx_state = State::from_cartesian_vec::<EARTH, ModifiedJulian>(&state_vec, this_dt, ECI {});
                for station in all_stations.iter() {
                    let meas = station.measure(rx_state, this_dt.into_instant());
                    if meas.visible() {
                        // XXX: Instant does not implement Eq, only PartialEq, so can't use it as an index =(
                        let delta = (prev_dt - this_dt.into_instant().duration()).duration().as_secs();
                        measurements.push((delta, meas));
                        prev_dt = this_dt.into_instant();
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
    let mut tb_estimator = TwoBodyWithStm::from_state::<EARTH, ECI>(initial_state);
    let mut prop_est = Propagator::new::<RK4Fixed>(&mut tb_estimator, &opts_est);
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
    for (meas_no, (duration, real_meas)) in measurements.iter().enumerate() {
        // Propagate the dynamics to the measurement, and then start the filter.
        let delta_time = (*duration) as f64;
        prop_est.until_time_elapsed(delta_time);
        if meas_no > num_meas_for_ekf && !kf.ekf {
            println!("switched to EKF");
            kf.ekf = true;
        }
        // Update the STM of the KF
        kf.update_stm(prop_est.dynamics.stm.clone());
        // Get the computed observation
        let this_dt = ModifiedJulian::from_instant(
            prev_dt.into_instant() + Instant::from_precise_seconds(delta_time, Era::Present).duration(),
        );
        if prev_dt == this_dt {
            panic!("already handled that time: {}", prev_dt);
        }
        prev_dt = this_dt;
        let rx_state =
            State::from_cartesian_vec::<EARTH, ModifiedJulian>(&prop_est.dynamics.two_body_dyn.state(), this_dt, ECI {});
        let mut still_empty = true;
        for station in all_stations.iter() {
            let computed_meas = station.measure(rx_state, this_dt.into_instant());
            if computed_meas.visible() {
                kf.update_h_tilde(*computed_meas.sensitivity());
                let latest_est = kf
                    .measurement_update(*real_meas.observation(), *computed_meas.observation())
                    .expect("wut?");
                still_empty = false;
                assert_eq!(latest_est.predicted, false, "estimate should not be a prediction");
                assert!(
                    latest_est.state.norm() < EPSILON,
                    "estimate error should be zero (perfect dynamics)"
                );
                // It's an EKF, so let's update the state in the dynamics.
                let now = prop_est.time(); // Needed because we can't do a mutable borrow while doing an immutable one too.
                let new_state = prop_est.dynamics.two_body_dyn.state() + latest_est.state;
                prop_est.dynamics.two_body_dyn.set_state(now, &new_state);
                break; // We know that only one station is in visibility at each time.
            }
        }
        if still_empty {
            // We're doing perfect everything, so we should always be in visibility
            panic!("T+{} : not in visibility", this_dt);
        }
    }
}

/*
// #[test]
fn ckf_fixed_step_perfect_stations_dual() {
    use std::{io, thread};

    // Define the ground stations.
    let elevation_mask = 0.0;
    let range_noise = 0.0;
    let range_rate_noise = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise);
    let dss34_canberra = GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise);
    let dss13_goldstone = GroundStation::dss13_goldstone(elevation_mask, range_noise, range_rate_noise);
    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = PropOpts::with_fixed_step(step_size, RSSStepPV {});

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx): (Sender<(f64, Vector6<f64>)>, Receiver<(f64, Vector6<f64>)>) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000); // Assume that we won't get more than 10k measurements.

    // Define state information.
    let dt = ModifiedJulian { days: 21545.0 };
    let initial_state = State::from_keplerian_eci(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt);

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let mut dyn = TwoBody::from_state_vec::<EARTH>(initial_state.to_cartesian_vec());
        let mut prop = Propagator::new::<RK4Fixed>(dyn, &opts.clone());
        dyn.tx_chan = Some(&truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    // Receive the states on the main thread, and populate the measurement channel.
    let mut prev_dt = dt.into_instant();
    loop {
        match truth_rx.recv() {
            Ok((t, state_vec)) => {
                // Convert the state to ECI.
                let this_dt =
                    ModifiedJulian::from_instant(dt.into_instant() + Instant::from_precise_seconds(t, Era::Present).duration());
                let rx_state = State::from_cartesian_vec::<EARTH, ModifiedJulian>(&state_vec, this_dt, ECI {});
                for station in all_stations.iter() {
                    let meas = station.measure(rx_state, this_dt.into_instant());
                    if meas.visible() {
                        // XXX: Instant does not implement Eq, only PartialEq, so can't use it as an index =(
                        let delta = (prev_dt - this_dt.into_instant().duration()).duration().as_secs();
                        measurements.push((delta, meas));
                        prev_dt = this_dt.into_instant();
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
    let mut tb_estimator = TwoBodyWithDualStm::from_state::<EARTH, ECI>(initial_state);
    let mut prop_est = Propagator::new::<RK4Fixed>(tb_estimator, &opts_est);
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
        stm: tb_estimator.stm.clone(),
        predicted: false,
    };

    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));
    let mut ckf = KF::initialize(initial_estimate, measurement_noise);

    println!("Will process {} measurements", measurements.len());

    let mut prev_dt = dt;
    let mut printed = false;

    let mut wtr = csv::Writer::from_writer(io::stdout());

    for (duration, real_meas) in measurements.iter() {
        // Propagate the dynamics to the measurement, and then start the filter.
        let delta_time = (*duration) as f64;
        prop_est.until_time_elapsed(delta_time);
        // Update the STM of the KF
        ckf.update_stm(tb_estimator.stm.clone());
        // Get the computed observation
        let this_dt = ModifiedJulian::from_instant(
            prev_dt.into_instant() + Instant::from_precise_seconds(delta_time, Era::Present).duration(),
        );
        if prev_dt == this_dt {
            panic!("already handled that time: {}", prev_dt);
        }
        prev_dt = this_dt;
        let rx_state = State::from_cartesian_vec::<EARTH, ModifiedJulian>(&tb_estimator.pos_vel, this_dt, ECI {});
        let mut still_empty = true;
        for station in all_stations.iter() {
            let computed_meas = station.measure(rx_state, this_dt.into_instant());
            if computed_meas.visible() {
                ckf.update_h_tilde(*computed_meas.sensitivity());
                let latest_est = ckf
                    .measurement_update(*real_meas.observation(), *computed_meas.observation())
                    .expect("wut?");
                still_empty = false;
                assert_eq!(latest_est.predicted, false, "estimate should not be a prediction");
                assert!(
                    latest_est.state.norm() < 1e-9,
                    "estimate error should be zero (perfect dynamics) ({:e})",
                    latest_est.state.norm()
                );
                if !printed {
                    wtr.serialize(latest_est).expect("could not write to stdout");
                    printed = true;
                }
                break; // We know that only one station is in visibility at each time.
            }
        }
        if still_empty {
            // We're doing perfect everything, so we should always be in visibility
            panic!("T+{} : not in visibility", this_dt);
        }
    }
}
*/
