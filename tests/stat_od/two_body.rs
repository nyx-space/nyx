extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::hifitime::instant::*;
use self::hifitime::julian::*;
use self::hifitime::SECONDS_PER_DAY;
use self::na::{Matrix2, Matrix6, U42, U6, Vector2, Vector6};
use self::nyx::celestia::{State, EARTH, ECI};
use self::nyx::dynamics::celestial::{TwoBody, TwoBodyWithStm};
use self::nyx::dynamics::Dynamics;
use self::nyx::od::kalman::{Estimate, FilterError, KF};
use self::nyx::od::ranging::{GroundStation, StdMeasurement};
use self::nyx::od::{Linearization, Measurement};
use self::nyx::propagators::{error_ctrl, Options, Propagator, RK4Fixed};
use std::collections::BTreeMap;
use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};

#[test]
fn fixed_step_perfect_stations() {
    use std::thread;
    // Tests that we can generate measurements on one side and get a proper estimate on the other.

    // Define the ground stations.
    let dss65_madrid = GroundStation::from_noise_values("Madrid", 0.0, 40.427222, 4.250556, 0.834939, 0.0, 0.0);
    let dss34_canberra = GroundStation::from_noise_values("Canberra", 0.0, -35.398333, 148.981944, 0.691750, 0.0, 0.0);
    let dss13_goldstone = GroundStation::from_noise_values("Goldstone", 0.0, 35.247164, 243.205, 1.07114904, 0.0, 0.0);
    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = Options::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx): (Sender<(f64, Vector6<f64>)>, Receiver<(f64, Vector6<f64>)>) = mpsc::channel();
    let mut measurements = BTreeMap::new();

    // Define state information.
    let dt = ModifiedJulian { days: 21545.0 };
    let initial_state = State::from_cartesian_eci(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0, dt);

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let mut prop = Propagator::new::<RK4Fixed>(&opts.clone());
        let mut dyn = TwoBody::from_state_vec::<EARTH>(initial_state.to_cartesian_vec());
        dyn.tx_chan = Some(&truth_tx);
        let (final_t, final_state_vec) = prop.until_time_elapsed(prop_time, &mut dyn, error_ctrl::rss_step_pos_vel);
        dyn.set_state(final_t, &final_state_vec);
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
                        measurements.insert(delta, meas);
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
    let mut prop_est = Propagator::new::<RK4Fixed>(&opts);
    let mut tb_estimator = TwoBodyWithStm::from_state::<EARTH, ECI>(initial_state);
    let covar_radius = 1.0;
    let covar_velocity = 1.0e-3;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius,
        covar_radius,
        covar_radius,
        covar_velocity,
        covar_velocity,
        covar_velocity,
    ));

    let initial_estimate = Estimate {
        state: tb_estimator.two_body_dyn.state(),
        covar: init_covar,
        stm: tb_estimator.stm.clone(),
        predicted: false,
    };

    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-3, 1e-3));
    let mut ckf = KF::initialize(initial_estimate, measurement_noise);

    println!("Will process {} measurements", measurements.keys().len());

    for (meas_no, (duration, real_meas)) in measurements.iter().enumerate() {
        // Propagate the dynamics to the measurement, and then start the filter.
        let delta_time = (*duration) as f64;
        print!("propagating for {:?}", delta_time);
        let (final_t, final_state) =
            prop_est.until_time_elapsed(delta_time, &mut tb_estimator, error_ctrl::largest_error::<U42>);
        tb_estimator.set_state(final_t, &final_state);
        println!("...done",);
        // Update the STM of the KF
        ckf.update_stm(tb_estimator.stm.clone());
        // Get the computed observation
        let this_dt = ModifiedJulian::from_instant(
            dt.into_instant() + Instant::from_precise_seconds(delta_time, Era::Present).duration(),
        );
        let rx_state = State::from_cartesian_vec::<EARTH, ModifiedJulian>(&tb_estimator.two_body_dyn.state(), this_dt, ECI {});
        let mut latest_est = Estimate::<U6>::empty();
        let mut still_empty = true;
        for station in all_stations.iter() {
            let computed_meas = station.measure(rx_state, this_dt.into_instant());
            if computed_meas.visible() {
                // We've got a visible measurement, so let's do a KF measurement update and stop searching for measurements.
                ckf.update_h_tilde(*computed_meas.sensitivity());
                latest_est = ckf
                    .measurement_update(*real_meas.observation(), *computed_meas.observation())
                    .expect("wut?");
                still_empty = false;
                break; // We know that only one station is in visibility at each time.
            }
        }
        if still_empty {
            // Do a time update instead
            latest_est = ckf.time_update().expect("huh?");
        }
        println!(
            "=== #{} PREDICTED: {} ===\nState {}Covariance {}\n=====================",
            meas_no, latest_est.predicted, latest_est.state, latest_est.covar
        );
    }
}
