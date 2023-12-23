extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use nyx::cosmic::{Bodies, Cosm, Orbit};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::io::ExportCfg;
use nyx::linalg::{Matrix2, Vector2};
use nyx::md::StateParameter;
use nyx::od::noise::GaussMarkov;
use nyx::od::prelude::*;
use nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use nyx::time::{Epoch, TimeUnits, Unit};
use nyx::utils::rss_orbit_errors;
use polars::prelude::*;
use std::collections::HashMap;
use std::env;
use std::fs::File;
use std::path::PathBuf;

/*
 * These tests check that if we start with a state deviation in the estimate, the filter will eventually converge back.
 * These tests do NOT check that the filter will converge if the initial state in the propagator has that state deviation.
 * The latter would require iteration and smoothing before playing with an EKF. This will be handled in a subsequent version.
**/

#[allow(clippy::identity_op)]
#[test]
fn od_robust_test_ekf_realistic_one_way() {
    let _ = pretty_env_logger::try_init();

    let cosm = Cosm::de438();

    let iau_earth = cosm.frame("IAU Earth");
    // Define the ground stations.
    let ekf_num_meas = 300;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 3 * Unit::Minute;
    let elevation_mask = 0.0;

    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        GaussMarkov::high_precision_range_km(),
        GaussMarkov::high_precision_doppler_km_s(),
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        GaussMarkov::high_precision_range_km(),
        GaussMarkov::high_precision_doppler_km_s(),
        iau_earth,
    );

    // Define the tracking configurations
    let configs = HashMap::from([
        (
            dss65_madrid.name.clone(),
            TrkConfig::from_sample_rate(60.seconds()),
        ),
        (
            dss34_canberra.name.clone(),
            TrkConfig::from_sample_rate(60.seconds()),
        ),
    ]);

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let all_stations = vec![dss65_madrid, dss34_canberra];

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // Disperse the initial state based on some orbital elements errors given from ULA Atlas 5 user guide, table 2.3.3-1 <https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf>
    // This assumes that the errors are ONE TENTH of the values given in the table. It assumes that the launch provider has provided an initial state vector, whose error is lower than the injection errors.
    // The initial covariance is computed based on the realized dispersions.
    let initial_estimate = KfEstimate::disperse_from_diag(
        initial_state,
        &[
            (StateParameter::Inclination, 0.0025),
            (StateParameter::RAAN, 0.022),
            (StateParameter::AoP, 0.02),
        ],
        Some(0),
    );
    println!("Initial estimate:\n{}", initial_estimate);

    let initial_state_dev = initial_estimate.nominal_state;
    let (init_rss_pos_km, init_rss_vel_km_s) = rss_orbit_errors(&initial_state, &initial_state_dev);

    println!("Truth initial state:\n{initial_state}\n{initial_state:x}");
    println!("Filter initial state:\n{initial_state_dev}\n{initial_state_dev:x}");
    println!(
        "Initial state dev:\t{:.3} m\t{:.3} m/s\n{}",
        init_rss_pos_km * 1e3,
        init_rss_vel_km_s * 1e3,
        initial_state - initial_state_dev
    );

    let bodies = vec![
        Bodies::Luna,
        Bodies::Sun,
        Bodies::JupiterBarycenter,
        Bodies::SaturnBarycenter,
    ];
    let orbital_dyn = OrbitalDynamics::point_masses(&bodies, cosm.clone());
    let truth_setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);
    let (_, traj) = truth_setup
        .with(initial_state)
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj.clone(), configs, 0).unwrap();
    arc_sim.build_schedule(cosm.clone()).unwrap();

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    // And serialize to disk
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "output_data",
        "ekf_robust_msr.parquet",
    ]
    .iter()
    .collect();

    arc.to_parquet_simple(&path).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be _nearly_ perfect because we've removed Saturn from the estimated trajectory
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let estimator = OrbitalDynamics::point_masses(&bodies, cosm.clone());
    let setup = Propagator::new::<RK4Fixed>(estimator, opts);
    let prop_est = setup.with(initial_state_dev.with_stm());

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    // Define the process noise to assume an unmodeled acceleration on X, Y and Z in the ECI frame
    let sigma_q = 5e-10_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);

    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    let trig = EkfTrigger::new(ekf_num_meas, ekf_disable_time);

    let mut odp = ODProcess::ekf(prop_est, kf, trig, None, cosm);

    // Let's filter and iterate on the initial subset of the arc to refine the initial estimate
    let subset = arc.filter_by_offset(..3.hours());
    let remaining = arc.filter_by_offset(3.hours()..);

    odp.process_arc::<GroundStation>(&subset).unwrap();
    odp.iterate_arc::<GroundStation>(&subset, IterationConf::once())
        .unwrap();

    let (sm_rss_pos_km, sm_rss_vel_km_s) =
        rss_orbit_errors(&initial_state, &odp.estimates[0].state());

    println!(
        "Initial state error after smoothing:\t{:.3} m\t{:.3} m/s\n{}",
        sm_rss_pos_km * 1e3,
        sm_rss_vel_km_s * 1e3,
        initial_state - odp.estimates[0].state()
    );

    odp.process_arc::<GroundStation>(&remaining).unwrap();

    odp.to_parquet(
        path.with_file_name("robustness_test_one_way.parquet"),
        ExportCfg::timestamped(),
    )
    .unwrap();

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];
    let final_truth_state = traj.at(est.epoch()).unwrap();

    println!("Estimate:\n{}", est);
    println!("Truth:\n{}", final_truth_state);
    println!(
        "Delta state with truth (epoch match: {}):\n{}",
        final_truth_state.epoch == est.epoch(),
        final_truth_state - est.state()
    );

    for i in 0..6 {
        if est.covar[(i, i)] < 0.0 {
            println!(
                "covar diagonal element negative @ [{}, {}] = {:.3e}-- issue #164",
                i,
                i,
                est.covar[(i, i)],
            );
        }
    }
    for i in 0..6 {
        assert!(
            est.covar[(i, i)] < initial_estimate.covar[(i, i)],
            "covar[({i}, {i})] did not decrease"
        );
    }

    assert_eq!(
        final_truth_state.epoch,
        est.epoch(),
        "time of final EST and TRUTH epochs differ"
    );
    let delta = est.state() - final_truth_state;
    println!(
        "RMAG error = {:.6} m\tVMAG error = {:.6} m/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e3
    );

    assert!(
        delta.rmag_km() < 0.06,
        "Position error should be less than 50 meters"
    );
    assert!(
        delta.vmag_km_s() < 2e-4,
        "Velocity error should be on centimeter level"
    );
}

#[allow(clippy::identity_op)]
#[test]
fn od_robust_test_ekf_realistic_two_way() {
    let _ = pretty_env_logger::try_init();

    let cosm = Cosm::de438();

    let iau_earth = cosm.frame("IAU Earth");
    // Define the ground stations.
    let ekf_num_meas = 300;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 3 * Unit::Minute;
    let elevation_mask = 0.0;

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let mut dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        GaussMarkov::high_precision_range_km(),
        GaussMarkov::high_precision_doppler_km_s(),
        iau_earth,
    );
    // Set the integration time so as to generate two way measurements
    dss65_madrid.integration_time = Some(60.seconds());
    let mut dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        GaussMarkov::high_precision_range_km(),
        GaussMarkov::high_precision_doppler_km_s(),
        iau_earth,
    );
    dss34_canberra.integration_time = Some(60.seconds());

    // Define the tracking configurations
    let configs = HashMap::from([
        (dss65_madrid.name.clone(), TrkConfig::default()),
        (dss34_canberra.name.clone(), TrkConfig::default()),
    ]);

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let devices = vec![dss65_madrid, dss34_canberra];

    // Disperse the initial state based on some orbital elements errors given from ULA Atlas 5 user guide, table 2.3.3-1 <https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf>
    // This assumes that the errors are ONE TENTH of the values given in the table. It assumes that the launch provider has provided an initial state vector, whose error is lower than the injection errors.
    // The initial covariance is computed based on the realized dispersions.
    let initial_estimate = KfEstimate::disperse_from_diag(
        initial_state,
        &[
            (StateParameter::Inclination, 0.0025),
            (StateParameter::RAAN, 0.022),
            (StateParameter::AoP, 0.02),
        ],
        Some(0),
    );
    println!("Initial estimate:\n{}", initial_estimate);

    let initial_state_dev = initial_estimate.nominal_state;
    let (init_rss_pos_km, init_rss_vel_km_s) = rss_orbit_errors(&initial_state, &initial_state_dev);

    println!("Truth initial state:\n{initial_state}\n{initial_state:x}");
    println!("Filter initial state:\n{initial_state_dev}\n{initial_state_dev:x}");
    println!(
        "Initial state dev:\t{:.3} m\t{:.3} m/s\n{}",
        init_rss_pos_km * 1e3,
        init_rss_vel_km_s * 1e3,
        initial_state - initial_state_dev
    );

    let bodies = vec![
        Bodies::Luna,
        Bodies::Sun,
        Bodies::JupiterBarycenter,
        Bodies::SaturnBarycenter,
    ];
    let orbital_dyn = OrbitalDynamics::point_masses(&bodies, cosm.clone());
    let truth_setup = Propagator::default(orbital_dyn);
    let (_, traj) = truth_setup
        .with(initial_state)
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(devices.clone(), traj.clone(), configs, 0).unwrap();
    arc_sim.build_schedule(cosm.clone()).unwrap();

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    // And serialize to disk
    let path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "output_data"].iter().collect();

    traj.to_parquet_simple(path.join("ekf_robust_two_way_traj.parquet"))
        .unwrap();
    arc.to_parquet_simple(path.join("ekf_robust_two_way_msr.parquet"))
        .unwrap();

    println!("{arc}");

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be _nearly_ perfect because we've removed Saturn from the estimated trajectory
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let estimator = OrbitalDynamics::point_masses(&bodies, cosm.clone());
    let setup = Propagator::default(estimator);
    let prop_est = setup.with(initial_state_dev.with_stm());

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    // Define the process noise to assume an unmodeled acceleration on X, Y and Z in the ECI frame
    let sigma_q = 5e-10_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);

    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    let trig = EkfTrigger::new(ekf_num_meas, ekf_disable_time);

    let mut odp = ODProcess::ekf(prop_est, kf, trig, None, cosm);

    // TODO: Fix the deserialization of the measurements such that they also deserialize the integration time.
    // Without it, we're stuck having to rebuild them from scratch.
    // https://github.com/nyx-space/nyx/issues/140

    // Build the hashmap of devices from the vector using their names
    let mut devices_map = devices
        .into_iter()
        .map(|dev| (dev.name.clone(), dev))
        .collect::<HashMap<_, _>>();

    // Check that exporting an empty set returns an error.
    assert!(odp
        .to_parquet(
            path.with_file_name("robustness_test_two_way.parquet"),
            ExportCfg::timestamped(),
        )
        .is_err());

    odp.process(
        &arc.measurements,
        &mut devices_map,
        arc.min_duration_sep().unwrap(),
    )
    .unwrap();

    let mut num_residual_none = 0;
    let mut num_residual_some = 0;
    odp.residuals.iter().for_each(|opt_v| match opt_v {
        Some(_) => num_residual_some += 1,
        None => {
            num_residual_none += 1;
        }
    });

    // Export as Parquet
    let timestamped_path = odp
        .to_parquet(
            path.with_file_name("robustness_test_two_way.parquet"),
            ExportCfg::timestamped(),
        )
        .unwrap();

    // Read in the Parquet file and assert proper data was written.

    // let df = LazyFrame::scan_parquet(timestamped_path, Default::default()).unwrap();
    let df = ParquetReader::new(File::open(timestamped_path).unwrap())
        .finish()
        .unwrap();

    // Note: this also checks that the columns that match the given measurement kind exist.
    let df_residuals = df
        .columns([
            "Prefit residual: Range (km)",
            "Prefit residual: Doppler (km/s)",
            "Postfit residual: Range (km)",
            "Postfit residual: Doppler (km/s)",
            "Residual ratio",
        ])
        .unwrap();

    for series in df_residuals.iter() {
        assert_eq!(series.len(), odp.estimates.len());
        let mut num_none = 0;
        let mut num_some = 0;
        series
            .f64()
            .unwrap()
            .into_iter()
            .for_each(|opt_v| match opt_v {
                Some(_) => num_some += 1,
                None => {
                    num_none += 1;
                }
            });

        assert_eq!(num_none, num_residual_none);
        assert_eq!(num_some, num_residual_some);
    }

    // Check that the position and velocity estimates are present, along with the epochs
    assert!(df
        .columns([
            "Epoch:Gregorian UTC",
            "Epoch:Gregorian TAI",
            "Epoch:TAI (s)",
            "x (km)",
            "y (km)",
            "z (km)",
            "vx (km/s)",
            "vy (km/s)",
            "vz (km/s)",
        ])
        .is_ok());

    // Check that the covariance in the integration frame is present
    assert!(df
        .columns([
            "Covariance XX (Earth J2000)",
            "Covariance XY (Earth J2000)",
            "Covariance XZ (Earth J2000)",
            "Covariance XVx (Earth J2000)",
            "Covariance XVy (Earth J2000)",
            "Covariance XVz (Earth J2000)",
            "Covariance YY (Earth J2000)",
            "Covariance YZ (Earth J2000)",
            "Covariance YVx (Earth J2000)",
            "Covariance YVy (Earth J2000)",
            "Covariance YVz (Earth J2000)",
            "Covariance ZZ (Earth J2000)",
            "Covariance ZVx (Earth J2000)",
            "Covariance ZVy (Earth J2000)",
            "Covariance ZVz (Earth J2000)",
            "Covariance VxVx (Earth J2000)",
            "Covariance VxVy (Earth J2000)",
            "Covariance VxVz (Earth J2000)",
            "Covariance VyVy (Earth J2000)",
            "Covariance VyVz (Earth J2000)",
            "Covariance VzVz (Earth J2000)",
        ])
        .is_ok());

    // Check that the covariance in the RIC frame is present
    assert!(df
        .columns([
            "Covariance XX (RIC)",
            "Covariance XY (RIC)",
            "Covariance XZ (RIC)",
            "Covariance XVx (RIC)",
            "Covariance XVy (RIC)",
            "Covariance XVz (RIC)",
            "Covariance YY (RIC)",
            "Covariance YZ (RIC)",
            "Covariance YVx (RIC)",
            "Covariance YVy (RIC)",
            "Covariance YVz (RIC)",
            "Covariance ZZ (RIC)",
            "Covariance ZVx (RIC)",
            "Covariance ZVy (RIC)",
            "Covariance ZVz (RIC)",
            "Covariance VxVx (RIC)",
            "Covariance VxVy (RIC)",
            "Covariance VxVz (RIC)",
            "Covariance VyVy (RIC)",
            "Covariance VyVz (RIC)",
            "Covariance VzVz (RIC)",
        ])
        .is_ok());

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];
    let final_truth_state = traj.at(est.epoch()).unwrap();

    println!("Estimate:\n{}", est);
    println!("Truth:\n{}", final_truth_state);
    println!(
        "Delta state with truth (epoch match: {}):\n{}",
        final_truth_state.epoch == est.epoch(),
        final_truth_state - est.state()
    );

    for i in 0..6 {
        if est.covar[(i, i)] < 0.0 {
            println!(
                "covar diagonal element negative @ [{}, {}] = {:.3e}-- issue #164",
                i,
                i,
                est.covar[(i, i)],
            );
        }
    }
    for i in 0..6 {
        assert!(
            est.covar[(i, i)] < initial_estimate.covar[(i, i)],
            "covar[({i}, {i})] did not decrease"
        );
    }

    assert_eq!(
        final_truth_state.epoch,
        est.epoch(),
        "time of final EST and TRUTH epochs differ"
    );
    let delta = est.state() - final_truth_state;
    println!(
        "RMAG error = {:.6} m\tVMAG error = {:.6} m/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e3
    );

    assert!(
        delta.rmag_km() < 0.2,
        "Position error should be less than 200 meters (down from >3 km)"
    );
    assert!(
        delta.vmag_km_s() < 1e-4,
        "Velocity error should be on decimeter level"
    );
}
