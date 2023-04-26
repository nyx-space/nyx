extern crate csv;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use nyx::cosmic::{Cosm, Orbit};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::sph_harmonics::Harmonics;
use nyx::io::formatter::NavSolutionFormatter;
use nyx::io::gravity::*;
use nyx::io::ConfigRepr;
use nyx::linalg::{Matrix2, Matrix6, Vector2, Vector6};
use nyx::od::noise::GaussMarkov;
use nyx::od::prelude::*;
use nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use std::collections::HashMap;
use std::env;
use std::path::PathBuf;

#[allow(clippy::identity_op)]
#[test]
fn od_tb_val_ekf_fixed_step_perfect_stations() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let iau_earth = cosm.frame("IAU Earth");

    // Define the ground stations.
    let ekf_num_meas = 100;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 5.0 * Unit::Second;
    let elevation_mask = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );

    // Define the tracking configurations
    let mut configs = HashMap::new();
    configs.insert(
        dss65_madrid.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        dss34_canberra.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        dss13_goldstone.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // We're sharing both the propagator and the dynamics.
    let orbital_dyn = OrbitalDynamics::two_body();
    let setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);

    let mut prop = setup.with(initial_state);
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs, 0).unwrap();
    arc_sim.disallow_overlap(); // Prevent overlapping measurements

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    println!("{}", final_truth);

    // let mut rng = thread_rng();
    // // Receive the states on the main thread, and populate the measurement channel.
    // while let Ok(rx_state) = truth_rx.try_recv() {
    //     for station in all_stations.iter_mut() {
    //         if let Some(meas) = station.measure(&rx_state, &mut rng, cosm.clone()) {
    //             measurements.push(meas);
    //             break; // We know that only one station is in visibility at each time.
    //         }
    //     }
    // }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(initial_state.with_stm());
    let covar_radius_km = 1.0e-6;
    let covar_velocity_km_s = 1.0e-6;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius_km,
        covar_radius_km,
        covar_radius_km,
        covar_velocity_km_s,
        covar_velocity_km_s,
        covar_velocity_km_s,
    ));

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(initial_state, init_covar);
    println!("initial estimate:\n{}", initial_estimate);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    let kf = KF::no_snc(initial_estimate, measurement_noise);

    let mut odp = ODProcess::ekf(
        prop_est,
        kf,
        EkfTrigger::new(ekf_num_meas, ekf_disable_time),
        cosm,
    );

    odp.process_arc::<GroundStation>(&arc).unwrap();

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("Final estimate:\n{est}");
    assert!(
        est.state_deviation().norm() < 1e-12,
        "In perfect modeling, the state deviation should be near zero"
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
                est.covar[(i, i)] < covar_radius_km,
                "covar radius did not decrease"
            );
        } else {
            assert!(
                est.covar[(i, i)] < covar_velocity_km_s,
                "covar velocity did not decrease"
            );
        }
    }

    // Check the final estimate
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("{}\n\n{}\n{}", est.state_deviation(), est, final_truth);
    let delta = est.state() - final_truth;
    println!(
        "RMAG error = {:.3} m\tVMAG error = {:.3} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 2e-16, "Position error should be zero");
    assert!(delta.vmag_km_s() < 2e-16, "Velocity error should be zero");
}

#[allow(clippy::identity_op)]
#[test]
fn od_tb_val_with_arc() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let iau_earth = cosm.frame("IAU Earth");

    // Define the ground stations.
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 5.0 * Unit::Second;
    let elevation_mask = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    let ekf_num_meas = 100;

    // Define the propagator information.
    let duration = 1 * Unit::Day;

    // Define the storages (channels for the states and a map for the measurements).

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // We're sharing both the propagator and the dynamics.
    let orbital_dyn = OrbitalDynamics::two_body();
    let setup = Propagator::new::<RK4Fixed>(orbital_dyn, PropOpts::with_fixed_step_s(10.0));

    let mut prop = setup.with(initial_state);
    let (final_truth, traj) = prop.for_duration_with_traj(duration).unwrap();
    println!("{}", final_truth);

    // Save the trajectory to parquet
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "output_data",
        "od_val_with_arc_truth_ephem.parquet",
    ]
    .iter()
    .collect();
    traj.to_parquet_simple(path).unwrap();

    // Load the tracking configs
    let trkconfig_yaml: PathBuf = [
        &env::var("CARGO_MANIFEST_DIR").unwrap(),
        "data",
        "tests",
        "config",
        "trk_cfg_od_val_arc.yaml",
    ]
    .iter()
    .collect();

    let configs: HashMap<String, TrkConfig> = TrkConfig::load_named(trkconfig_yaml).unwrap();

    // Simulate tracking data of range and range rate
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs, 1).unwrap();

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    // And serialize to disk
    let path: PathBuf = [
        &env::var("CARGO_MANIFEST_DIR").unwrap(),
        "output_data",
        "two_body_od_val_arc.parquet",
    ]
    .iter()
    .collect();

    arc.to_parquet_simple(path).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(initial_state.with_stm());
    let covar_radius_km = 1.0e-6;
    let covar_velocity_km_s = 1.0e-6;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius_km,
        covar_radius_km,
        covar_radius_km,
        covar_velocity_km_s,
        covar_velocity_km_s,
        covar_velocity_km_s,
    ));

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(initial_state, init_covar);
    println!("initial estimate:\n{}", initial_estimate);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    let kf = KF::no_snc(initial_estimate, measurement_noise);

    let mut odp = ODProcess::ekf(
        prop_est,
        kf,
        EkfTrigger::new(ekf_num_meas, ekf_disable_time),
        cosm,
    );

    odp.process_arc::<GroundStation>(&arc).unwrap();

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("Final estimate:\n{est}");
    assert!(
        dbg!(est.state_deviation().norm()) < f64::EPSILON,
        "In perfect modeling, the state deviation should be near zero"
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
                est.covar[(i, i)] < covar_radius_km,
                "covar radius did not decrease"
            );
        } else {
            assert!(
                est.covar[(i, i)] < covar_velocity_km_s,
                "covar velocity did not decrease"
            );
        }
    }

    // Check the final estimate
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("{}\n\n{}\n{}", est.state_deviation(), est, final_truth);
    let delta = est.state() - final_truth;
    println!(
        "RMAG error = {:.3} m\tVMAG error = {:.3} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 2e-16, "Position error should be zero");
    assert!(delta.vmag_km_s() < 2e-16, "Velocity error should be zero");
}

#[allow(clippy::identity_op)]
#[test]
fn od_tb_val_ckf_fixed_step_perfect_stations() {
    /*
     * This tests that the state transition matrix computation is correct with two body dynamics.
     *
     * Specifically, the same dynamics are used for both the measurement generation and for the estimation.
     * However, only the estimation generation propagates the STM. When STM propagation is enabled, the code will compute
     * the dynamics using a hyperdual representation in 7 dimensions: 1 for the reals, 3 for the position partials,
     * 3 for the velocity partials.
     *
     * Hence, if the filter state estimation is any different from the truth data, then it means that the equations of
     * motion computed in hyperdual space differ from the ones computes in the reals.
     *
     * Thereby, this serves as a validation of the orbital dynamics implementation.
     **/
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::io;

    let cosm = Cosm::de438();
    let iau_earth = cosm.frame("IAU Earth");

    // Define the ground stations.
    let elevation_mask = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );

    // Define the tracking configurations
    let mut configs = HashMap::new();
    configs.insert(
        dss65_madrid.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        dss34_canberra.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        dss13_goldstone.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // Generate the truth data on one thread.
    let orbital_dyn = OrbitalDynamics::two_body();

    let setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);
    let mut prop = setup.with(initial_state);

    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs, 0).unwrap();
    arc_sim.disallow_overlap(); // Prevent overlapping measurements

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let initial_state_est = initial_state.with_stm();
    // Use the same setup as earlier
    let prop_est = setup.with(initial_state_est);
    let covar_radius_km = 1.0e-3;
    let covar_velocity_km_s = 1.0e-6;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius_km,
        covar_radius_km,
        covar_radius_km,
        covar_velocity_km_s,
        covar_velocity_km_s,
        covar_velocity_km_s,
    ));

    // Define the initial orbit estimate
    let initial_estimate = KfEstimate::from_covar(initial_state_est, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    let ckf = KF::no_snc(initial_estimate, measurement_noise);

    let mut odp = ODProcess::ckf(prop_est, ckf, cosm.clone());

    odp.process_arc::<GroundStation>(&arc).unwrap();

    // Initialize the formatter
    let estimate_fmtr = NavSolutionFormatter::default("tb_ckf.csv".to_owned(), cosm);

    let mut wtr = csv::Writer::from_writer(io::stdout());
    wtr.serialize(&estimate_fmtr.headers)
        .expect("could not write to stdout");
    // Check that we have as many estimates as steps taken by the propagator.
    // Note that this test cannot work when using a variable step propagator in that same setup.
    // We're adding +1 because the propagation time is inclusive on both ends.
    // We add another +1 because we have a measurement on the initial estimate of the propagator, so we propagate by zero seconds but still process that measurement.
    let expected_num_estimates = (prop_time.to_seconds() / step_size.to_seconds()) as usize + 2;

    // Check that there are no duplicates of epochs.
    let mut prev_epoch = odp.estimates[0].epoch();

    for est in odp.estimates.iter().skip(2) {
        let this_epoch = est.epoch();
        assert!(
            this_epoch > prev_epoch,
            "Estimates not continuously going forward: {this_epoch} <= {prev_epoch}"
        );
        prev_epoch = this_epoch;
    }

    assert_eq!(
        odp.estimates.len(),
        expected_num_estimates,
        "Different number of estimates received"
    );

    let mut printed = false;
    for (no, est) in odp.estimates.iter().enumerate() {
        if no == 0 {
            // Skip the first estimate which is the initial estimate provided by user
            continue;
        }
        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0,
                "covar diagonal element negative @ [{}, {}] = {:e} @ {}",
                i,
                i,
                est.covar[(i, i)],
                est.epoch()
            );
        }
        assert!(
            est.state_deviation().norm() < 1e-12,
            "estimate error should be zero (perfect dynamics) ({:e})",
            est.state_deviation().norm()
        );

        if !printed {
            // Format the estimate
            wtr.serialize(estimate_fmtr.fmt(est))
                .expect("could not write to stdout");
            printed = true;
        }
    }

    for res in &odp.residuals {
        assert!(
            res.postfit.norm() < 1e-12,
            "postfit should be zero (perfect dynamics) ({:e})",
            res
        );
    }

    let est = odp.estimates.last().unwrap();
    println!("estimate error {:.2e}", est.state_deviation().norm());
    println!("estimate covariance {:.2e}", est.covar.diagonal().norm());

    assert!(
        est.state_deviation().norm() < 1e-12,
        "estimate error should be zero (perfect dynamics) ({:e})",
        est.state_deviation().norm()
    );

    assert!(
        est.covar.diagonal().norm() < 1e-6,
        "estimate covariance norm should be zero (perfect dynamics) ({:e})",
        est.covar.diagonal().norm()
    );

    let delta = est.state() - final_truth;
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 2e-16, "Position error should be zero");
    assert!(delta.vmag_km_s() < 2e-16, "Velocity error should be zero");

    // Iterate
    odp.iterate_arc::<GroundStation>(
        &arc,
        IterationConf {
            smoother: SmoothingArc::TimeGap(10.0 * Unit::Second),
            ..Default::default()
        },
    )
    .unwrap();

    println!(
        "N-1 one iteration: \n{}",
        odp.estimates[odp.estimates.len() - 1]
    );

    println!(
        "Initial state after iteration: \n{:x}",
        odp.estimates[0].state()
    );

    // Check the final estimate
    let est = odp.estimates.last().unwrap();
    println!("estimate error {:.2e}", est.state_deviation().norm());
    println!("estimate covariance {:.2e}", est.covar.diagonal().norm());

    assert!(
        est.state_deviation().norm() < 1e-12,
        "estimate error should be zero (perfect dynamics) ({:e})",
        est.state_deviation().norm()
    );

    // Note we accept a larger covariance diagonal here because smoothing will increase the covariance
    assert!(
        est.covar.diagonal().norm() < 1e-4,
        "estimate covariance norm should be zero (perfect dynamics) ({:e})",
        est.covar.diagonal().norm()
    );

    let delta = est.state() - final_truth;
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 1e-9, "More than 1 micrometer error");
    assert!(delta.vmag_km_s() < 1e-9, "More than 1 micrometer/s error");
}

#[allow(clippy::identity_op)]
#[test]
fn od_tb_ckf_fixed_step_iteration_test() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let iau_earth = cosm.frame("IAU Earth");

    // Define the ground stations.
    let elevation_mask = 0.0;
    let range_noise = 0.1; // in km (so 100 meters of error)
    let range_rate_noise = 0.001; // in km/s (or 1 meter per second of error)
    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );

    // Define the tracking configurations
    let mut configs = HashMap::new();
    configs.insert(
        dss65_madrid.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        dss34_canberra.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        dss13_goldstone.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let orbital_dyn = OrbitalDynamics::two_body();
    let setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);

    let mut prop = setup.with(initial_state);
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs, 0).unwrap();
    arc_sim.disallow_overlap(); // Prevent overlapping measurements

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    // Check that we have the same number of measurements as before the behavior change.
    assert_eq!(arc.measurements.len(), 7954);

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(initial_state.with_stm());
    let covar_radius_km = 1.0e-3;
    let covar_velocity_km_s = 1.0e-6;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius_km,
        covar_radius_km,
        covar_radius_km,
        covar_velocity_km_s,
        covar_velocity_km_s,
        covar_velocity_km_s,
    ));

    // Define the initial estimate (x_hat): add 100 meters in X, remove 100 meters in Y and add 50 meters in Z
    let mut initial_state2 = initial_state;
    initial_state2.x_km += 0.1;
    initial_state2.y_km -= 0.1;
    initial_state2.z_km += 0.05;
    let initial_estimate = KfEstimate::from_covar(initial_state2, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    let ckf = KF::no_snc(initial_estimate, measurement_noise);

    let mut odp = ODProcess::ckf(prop_est, ckf, cosm);

    odp.process_arc::<GroundStation>(&arc).unwrap();

    // Check the final estimate prior to iteration
    let delta = odp.estimates.last().unwrap().state() - final_truth;
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(
        delta.rmag_km() < range_noise,
        "More than station level position error"
    );
    assert!(
        delta.vmag_km_s() < range_rate_noise,
        "More than station level velocity error"
    );

    // Iterate, and check that the initial state difference is lower
    odp.iterate_arc::<GroundStation>(
        &arc,
        IterationConf {
            smoother: SmoothingArc::TimeGap(10.0 * Unit::Second),
            ..Default::default()
        },
    )
    .unwrap();

    let dstate_no_iteration = initial_state - initial_state2;
    let dstate_iteration = initial_state - odp.estimates[0].state();

    println!("{}\n{}", initial_state2, odp.estimates[0].state());

    // Compute the order of magnitude of the errors, and check that iteration either decreases it or keeps it the same
    let err_it_oom = dstate_iteration.rmag_km().log10().floor() as i32;
    let err_no_it_oom = dstate_no_iteration.rmag_km().log10().floor() as i32;

    println!(
        "Difference in initial states radii without iterations: {} km (order of magnitude: {})",
        dstate_no_iteration.rmag_km(),
        err_no_it_oom
    );
    println!(
        "Difference in initial states radii with iterations: {} km (order of magnitude: {})",
        dstate_iteration.rmag_km(),
        err_it_oom
    );
    assert!(
        dstate_iteration.rmag_km() < dstate_no_iteration.rmag_km() || err_it_oom <= err_no_it_oom,
        "Iteration did not reduce initial error"
    );

    // Check the final estimate
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("{}\n\n{}\n{}", est.state_deviation(), est, final_truth);
    let delta = est.state() - final_truth;
    println!(
        "RMAG error = {:.3} m\tVMAG error = {:.3} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 75e-3, "More than 75 meter error");
    assert!(delta.vmag_km_s() < 50e-6, "More than 50 mm/s error");
}

#[allow(clippy::identity_op)]
#[test]
fn od_tb_ckf_fixed_step_perfect_stations_snc_covar_map() {
    // Tests state noise compensation with covariance mapping
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let iau_earth = cosm.frame("IAU Earth");

    // Define the ground stations.
    let elevation_mask = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );

    // Define the tracking configurations
    let mut configs = HashMap::new();
    configs.insert(
        dss65_madrid.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        dss34_canberra.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        dss13_goldstone.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let orbital_dyn = OrbitalDynamics::two_body();
    let setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);

    let mut prop = setup.with(initial_state);
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs, 0).unwrap();
    arc_sim.disallow_overlap(); // Prevent overlapping measurements

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(initial_state.with_stm());

    // Set up the filter
    let covar_radius_km = 1.0e-3;
    let covar_velocity_km_s = 1.0e-6;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius_km,
        covar_radius_km,
        covar_radius_km,
        covar_velocity_km_s,
        covar_velocity_km_s,
        covar_velocity_km_s,
    ));

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(initial_state, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    // Define the process noise to assume an unmodeled acceleration of 1e-3 km^2/s^2 on X, Y and Z in the ECI frame
    let sigma_q = 1e-8_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);

    let ckf = KF::new(initial_estimate, process_noise, measurement_noise);

    let mut odp = ODProcess::ckf(prop_est, ckf, cosm);

    odp.process_arc::<GroundStation>(&arc).unwrap();

    let mut wtr = csv::Writer::from_path("./estimation.csv").unwrap();

    // Let's export these to a CSV file, and also check that the covariance never falls below our sigma squared values
    for (no, est) in odp.estimates.iter().enumerate() {
        if no == 1 {
            println!("{}", est);
        }
        assert!(
            est.state_deviation().norm() < 1e-12,
            "estimate error should be zero (perfect dynamics) ({:e})",
            est.state_deviation().norm()
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
            assert!(
                est.covar[(i, i)] >= sigma_q,
                "covar diagonal less than SNC value @ {} = {:.3e}",
                no,
                est.covar[(i, i)]
            );
        }

        wtr.serialize(*est).expect("could not write to stdout");
    }

    // Check the final estimate
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("{}\n\n{}\n{}", est.state_deviation(), est, final_truth);
    let delta = est.state() - final_truth;
    println!(
        "RMAG error = {:.3} m\tVMAG error = {:.3} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 1e-3, "More than 1 meter error");
    assert!(delta.vmag_km_s() < 1e-6, "More than 1 mm/s error");
}

#[allow(clippy::identity_op)]
#[test]
fn od_tb_ckf_map_covar() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();

    // Define the propagator information.
    let prop_time = 2 * Unit::Day;
    let step_size = 10.0 * Unit::Second;

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.

    let setup = Propagator::new::<RK4Fixed>(
        OrbitalDynamics::two_body(),
        PropOpts::with_fixed_step(step_size),
    );
    let prop_est = setup.with(initial_state.with_stm());
    let covar_radius_km = 1.0e-3;
    let covar_velocity_km_s = 1.0e-6;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius_km,
        covar_radius_km,
        covar_radius_km,
        covar_velocity_km_s,
        covar_velocity_km_s,
        covar_velocity_km_s,
    ));

    let initial_estimate = KfEstimate::from_covar(initial_state, init_covar);

    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));
    let ckf = KF::no_snc(initial_estimate, measurement_noise);

    let mut odp: ODProcess<
        OrbitalDynamics,
        nyx::propagators::RSSCartesianStep,
        RangeDoppler,
        CkfTrigger,
        nalgebra::Const<3>,
        Orbit,
        KF<Orbit, nalgebra::Const<3>, nalgebra::Const<2>>,
    > = ODProcess::ckf(prop_est, ckf, cosm);

    odp.map_covar(dt + prop_time).unwrap();

    // Check that the covariance inflated (we don't get the norm of the estimate because it's zero without any truth data)
    let estimates = odp.estimates;
    let est = &estimates[estimates.len() - 1];
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
                est.covar[(i, i)] > covar_radius_km,
                "covar radius did not increase"
            );
        } else {
            assert!(
                est.covar[(i, i)] > covar_velocity_km_s,
                "covar velocity did not increase"
            );
        }
    }
}

#[allow(clippy::identity_op)]
#[test]
fn od_tb_val_harmonics_ckf_fixed_step_perfect() {
    // Tests state noise compensation with covariance mapping
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let iau_earth = cosm.frame("IAU Earth");

    // Define the ground stations.
    let elevation_mask = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );

    // Define the tracking configurations
    let mut configs = HashMap::new();
    configs.insert(
        dss65_madrid.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        dss34_canberra.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        dss13_goldstone.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let iau_earth = cosm.frame("IAU Earth");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let earth_sph_harm = HarmonicsMem::from_cof("data/JGM3.cof.gz", 70, 70, true).unwrap();
    let harmonics = Harmonics::from_stor(iau_earth, earth_sph_harm, cosm.clone());
    let orbital_dyn = OrbitalDynamics::from_model(harmonics);
    let setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);

    let mut prop = setup.with(initial_state);
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs, 0).unwrap();
    arc_sim.disallow_overlap(); // Prevent overlapping measurements

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(initial_state.with_stm());

    // Set up the filter
    let covar_radius_km = 1.0e-3;
    let covar_velocity_km_s = 1.0e-6;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius_km,
        covar_radius_km,
        covar_radius_km,
        covar_velocity_km_s,
        covar_velocity_km_s,
        covar_velocity_km_s,
    ));

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(initial_state, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    let ckf = KF::no_snc(initial_estimate, measurement_noise);

    let mut odp = ODProcess::ckf(prop_est, ckf, cosm);

    odp.process_arc::<GroundStation>(&arc).unwrap();
    let mut wtr = csv::Writer::from_path("./estimation.csv").unwrap();

    // Let's export these to a CSV file, and also check that the covariance never falls below our sigma squared values
    for (no, est) in odp.estimates.iter().enumerate() {
        if no == 1 {
            println!("{}", est);
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
            est.state_deviation().norm() < 1e-2,
            "estimate error should be good (perfect dynamics) ({:e})",
            est.state_deviation().norm()
        );

        wtr.serialize(*est).expect("could not write to stdout");
    }

    // Check the final estimate
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("{}\n\n{}\n{}", est.state_deviation(), est, final_truth);
    let delta = est.state() - final_truth;
    println!(
        "RMAG error = {:.3} m\tVMAG error = {:.3} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 2e-16, "Position error should be zero");
    assert!(delta.vmag_km_s() < 2e-16, "Velocity error should be zero");
}

#[allow(clippy::identity_op)]
#[test]
fn od_tb_ckf_fixed_step_perfect_stations_several_snc_covar_map() {
    // Tests state noise compensation with covariance mapping
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let iau_earth = cosm.frame("IAU Earth");

    // Define the ground stations.
    let elevation_mask = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );

    // Define the tracking configurations
    let mut configs = HashMap::new();
    configs.insert(
        dss65_madrid.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        dss34_canberra.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );
    configs.insert(
        dss13_goldstone.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let orbital_dyn = OrbitalDynamics::two_body();
    let setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);
    let mut prop = setup.with(initial_state);
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs, 0).unwrap();
    arc_sim.disallow_overlap(); // Prevent overlapping measurements

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(initial_state.with_stm());

    // Set up the filter
    let covar_radius_km = 1.0e-3;
    let covar_velocity_km_s = 1.0e-6;
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius_km,
        covar_radius_km,
        covar_radius_km,
        covar_velocity_km_s,
        covar_velocity_km_s,
        covar_velocity_km_s,
    ));

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(initial_state, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    // Define the process noise to assume an unmodeled acceleration of 1e-3 km^2/s^2 on X, Y and Z in the ECI frame
    let sigma_q1 = 1e-7_f64.powi(2);
    let process_noise1 = SNC3::from_diagonal(2 * Unit::Day, &[sigma_q1, sigma_q1, sigma_q1]);

    let sigma_q2 = 1e-8_f64.powi(2);
    let sigma_q2_d = 3600.0;
    let mut process_noise2 = SNC3::with_decay(
        2 * Unit::Day,
        &[sigma_q2, sigma_q2, sigma_q2],
        &[sigma_q2_d, sigma_q2_d, sigma_q2_d],
    );
    process_noise2.start_time = Some(dt + 36_000.0); // Start the second process noise 10 hours into the tracking pass

    let ckf = KF::with_sncs(
        initial_estimate,
        vec![process_noise1, process_noise2],
        measurement_noise,
    );

    let mut odp = ODProcess::ckf(prop_est, ckf, cosm);

    odp.process_arc::<GroundStation>(&arc).unwrap();

    let mut wtr = csv::Writer::from_path("./estimation.csv").unwrap();

    // Let's export these to a CSV file, and also check that the covariance never falls below our sigma squared values
    for (no, est) in odp.estimates.iter().enumerate() {
        if no == 1 {
            println!("{}", est);
        }
        assert!(
            est.state_deviation().norm() < 1e-6,
            "estimate error should be zero (perfect dynamics) ({:e})",
            est.state_deviation().norm()
        );

        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0,
                "covar diagonal element negative @ [{}, {}]",
                i,
                i
            );
        }

        wtr.serialize(*est).expect("could not write to stdout");
    }

    // Check the final estimate
    let est = &odp.estimates[odp.estimates.len() - 1];
    println!("{}\n\n{}\n{}", est.state_deviation(), est, final_truth);
    let delta = est.state() - final_truth;
    println!(
        "RMAG error = {:.3} m\tVMAG error = {:.3} m/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e3
    );

    assert!(delta.rmag_km() < 2e-16, "Position error should be zero");
    assert!(delta.vmag_km_s() < 2e-16, "Velocity error should be zero");
}
