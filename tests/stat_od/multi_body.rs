extern crate csv;
extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::hifitime::{Epoch, SECONDS_PER_DAY};
use self::na::{Matrix2, Matrix2x6, Matrix6, Vector2, Vector6, U6};
use self::nyx::celestia::{Cosm, Geoid, State};
use self::nyx::dynamics::celestial::{CelestialDynamics, CelestialDynamicsStm};
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
fn muli_body_ckf_perfect_stations() {
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
    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(3);
    let dt = Epoch::from_mjd_tai(21545.0);
    let initial_state = State::from_keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, earth_geoid);

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
        let mut dynamics = CelestialDynamics::new(initial_state, bodies, &cosm);
        // let mut dynamics = CelestialDynamics::two_body(initial_state);
        let mut prop = Propagator::new::<RK4Fixed>(&mut dynamics, &opts);
        prop.tx_chan = Some(&truth_tx);
        prop.until_time_elapsed(prop_time);
    });

    // Receive the states on the main thread, and populate the measurement channel.
    let mut prev_t = 0.0;
    while let Ok((t, state_vec)) = truth_rx.recv() {
        // Convert the state to ECI.
        let this_dt = Epoch::from_mjd_tai(dt.as_mjd_tai_days() + t / SECONDS_PER_DAY);
        let rx_state = State::<Geoid>::from_cartesian_vec(&state_vec, this_dt, earth_geoid);
        for station in all_stations.iter() {
            let meas = station.measure(rx_state, this_dt);
            if meas.visible() {
                let delta = t - prev_t;
                measurements.push((delta, meas));
                prev_t = t;
                break; // We know that only one station is in visibility at each time.
            }
        }
    }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let opts_est = PropOpts::with_fixed_step(step_size, LargestError {});
    let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
    let mut tb_estimator = CelestialDynamicsStm::new(initial_state, bodies, &cosm);
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
        stm: prop_est.dynamics.stm,
        predicted: false,
    };

    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));
    let mut ckf = KF::initialize(initial_estimate, measurement_noise);

    println!("Will process {} measurements", measurements.len());

    let mut this_dt = dt;
    let mut printed = false;

    let mut wtr = csv::Writer::from_writer(io::stdout());

    for (duration, real_meas) in measurements.iter() {
        // Propagate the dynamics to the measurement, and then start the filter.
        let delta_time = (*duration) as f64;
        prop_est.until_time_elapsed(delta_time);
        // Update the STM of the KF
        ckf.update_stm(prop_est.dynamics.stm);
        // Get the computed observation
        assert!(delta_time > 0.0, "repeated time");
        this_dt.mut_add_secs(delta_time);
        let rx_state = State::<Geoid>::from_cartesian_vec(&prop_est.dynamics.state.to_cartesian_vec(), this_dt, earth_geoid);
        let mut still_empty = true;
        for station in all_stations.iter() {
            let computed_meas = station.measure(rx_state, this_dt);
            if computed_meas.visible() {
                ckf.update_h_tilde(computed_meas.sensitivity());
                let latest_est = ckf
                    .measurement_update(real_meas.observation(), computed_meas.observation())
                    .expect("wut?");
                still_empty = false;
                assert_eq!(latest_est.predicted, false, "estimate should not be a prediction");
                assert!(
                    latest_est.state.norm() < EPSILON,
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
            panic!("T {} : not in visibility", this_dt.as_mjd_tai_days());
        }
    }
}
