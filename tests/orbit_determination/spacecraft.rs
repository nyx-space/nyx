extern crate csv;
extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use self::hifitime::{Epoch, SECONDS_PER_DAY};
use self::na::{Matrix2, Matrix6, Vector2, Vector6};
use self::nyx::celestia::{bodies, Cosm, Orbit};
use self::nyx::dynamics::orbital::{OrbitalDynamics, OrbitalDynamicsStm};
use self::nyx::dynamics::spacecraft::{SolarPressure, Spacecraft};
use self::nyx::dynamics::Dynamics;
use self::nyx::io::formatter::NavSolutionFormatter;
use self::nyx::od::ui::*;
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use std::sync::mpsc;

#[test]
fn sc_ckf_perfect_stations() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::{io, thread};

    let cosm = Cosm::de438();

    // Define the ground stations.
    let elevation_mask = 0.0;
    let range_noise = 0.0;
    let range_rate_noise = 0.0;
    let dss65_madrid =
        GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise, &cosm);
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise, &cosm);
    let dss13_goldstone =
        GroundStation::dss13_goldstone(elevation_mask, range_noise, range_rate_noise, &cosm);
    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let sc_dry_mass = 100.0; // in kg
    let sc_area = 5.0; // m^2

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let cosm = Cosm::de438();
        let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
        let orbital_dyn = OrbitalDynamics::point_masses(initial_state, &bodies, &cosm);
        let mut dynamics = Spacecraft::new(orbital_dyn, sc_dry_mass);
        dynamics.add_model(Box::new(SolarPressure::default(
            sc_area,
            vec![cosm.frame("EME2000")],
            &cosm,
        )));
        let mut prop = Propagator::new::<RK4Fixed>(&mut dynamics, opts);
        prop.tx_chan = Some(truth_tx);
        prop.until_time_elapsed(prop_time).unwrap();
    });

    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_sc_state) = truth_rx.recv() {
        for station in all_stations.iter() {
            let rx_state = rx_sc_state.orbit;
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
    let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
    let orbital_dyn = OrbitalDynamicsStm::point_masses(initial_state, &bodies, &cosm);
    let mut estimator = Spacecraft::with_stm(orbital_dyn, sc_dry_mass);
    let init_sc_state = estimator.state();
    estimator.add_model(Box::new(SolarPressure::default(
        sc_area,
        vec![cosm.frame("EME2000")],
        &cosm,
    )));
    let prop_est = Propagator::new::<RK4Fixed>(&mut estimator, opts);
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
    let initial_estimate = KfEstimate::from_covar(init_sc_state, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise =
        Matrix2::from_diagonal(&Vector2::new(15e-3_f64.powi(2), 1e-5_f64.powi(2)));

    let ckf = KF::no_snc(initial_estimate, measurement_noise);

    let mut odp = ODProcess::ckf(prop_est, ckf, all_stations, false, measurements.len());

    odp.process_measurements(&measurements).unwrap();

    // Initialize the formatter
    let estimate_fmtr = NavSolutionFormatter::default("sc_ckf.csv".to_owned(), &cosm);

    let mut wtr = csv::Writer::from_writer(io::stdout());
    wtr.serialize(&estimate_fmtr.headers)
        .expect("could not write to stdout");
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
            est.state_deviation().norm() < 1e-6,
            "estimate error should be zero (perfect dynamics) ({:e})",
            est.state_deviation().norm()
        );

        if !printed {
            // Format the estimate
            wtr.serialize(estimate_fmtr.fmt(est))
                .expect("could not write to stdout");
            printed = true;
        }

        last_est = Some(est);
    }

    for res in &odp.residuals {
        assert!(
            res.postfit.norm() < 1e-9,
            "postfit should be zero (perfect dynamics) ({:e})",
            res
        );
    }

    // NOTE: We do not check whether the covariance has deflated because it is possible that it inflates before deflating.
    // The filter in multibody dynamics has been validated against JPL tools using a proprietary scenario.
    let est = last_est.unwrap();
    assert!(est.state_deviation().norm() < 1e-8);
    assert!(est.covar.norm() < 1e-5);
}
