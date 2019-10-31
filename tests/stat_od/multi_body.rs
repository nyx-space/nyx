extern crate csv;
extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::hifitime::{Epoch, SECONDS_PER_DAY};
use self::na::{Matrix2, Matrix6, Vector2, Vector6};
use self::nyx::celestia::{bodies, Cosm, Geoid, State};
use self::nyx::dynamics::celestial::{CelestialDynamics, CelestialDynamicsStm};
use self::nyx::od::kalman::{Estimate, KF};
use self::nyx::od::ranging::GroundStation;
use self::nyx::od::Measurement;
use self::nyx::propagators::error_ctrl::{RSSStepPV, RSSStepPVStm};
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};

#[test]
fn multi_body_ckf_perfect_stations() {
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
            let meas = station.measure(rx_state, rx_state.dt);
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

    let initial_estimate = Estimate::from_covar(dt, init_covar);

    let measurement_noise =
        Matrix2::from_diagonal(&Vector2::new(15e-3_f64.powi(2), 1e-5_f64.powi(2)));
    let mut ckf = KF::initialize(initial_estimate, measurement_noise);

    println!("Will process {} measurements", measurements.len());

    let mut printed = false;

    let mut wtr = csv::Writer::from_writer(io::stdout());

    let mut last_est = None;
    let mut prev_dt = dt;
    for (abs_dt, real_meas) in measurements.iter() {
        // Propagate the dynamics to the measurement, and then start the filter.
        let delta_time = abs_dt.as_tai_seconds() - prev_dt.as_tai_seconds();
        let (new_t, _) = prop_est.until_time_elapsed(delta_time);
        prev_dt = *abs_dt;
        // Update the STM of the KF
        ckf.update_stm(prop_est.dynamics.stm);
        // Get the computed observation
        assert!(delta_time > 0.0, "repeated time");
        let mut this_dt = dt;
        this_dt.mut_add_secs(new_t);
        let rx_state = State::<Geoid>::from_cartesian_vec(
            &prop_est.dynamics.state.to_cartesian_vec(),
            this_dt,
            earth_geoid,
        );
        let mut still_empty = true;
        for station in all_stations.iter() {
            let computed_meas = station.measure(rx_state, this_dt);
            if computed_meas.visible() {
                ckf.update_h_tilde(computed_meas.sensitivity());
                let latest_est = ckf
                    .measurement_update(
                        this_dt,
                        real_meas.observation(),
                        computed_meas.observation(),
                    )
                    .expect("wut?");
                still_empty = false;
                assert_eq!(
                    latest_est.predicted, false,
                    "estimate should not be a prediction"
                );
                assert!(
                    latest_est.state.norm() < 1e-6,
                    "estimate error should be zero (perfect dynamics) ({:e})",
                    latest_est.state.norm()
                );

                if !printed {
                    wtr.serialize(latest_est.clone())
                        .expect("could not write to stdout");
                    printed = true;
                }

                last_est = Some(latest_est);

                break; // We know that only one station is in visibility at each time.
            }
        }
        if still_empty {
            // We're doing perfect everything, so we should always be in visibility
            panic!("T {} : not in visibility", this_dt.as_mjd_tai_days());
        }
    }

    // NOTE: We do not check whether the covariance has deflated because it is possible that it inflates before deflating.
    // The filter in multibody dynamics has been validated against JPL tools using a proprietary scenario.
    let est = last_est.unwrap();
    assert!(est.state.norm() < 1e-8);
    assert!(est.covar.norm() < 1e-5);
}
