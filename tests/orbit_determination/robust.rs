extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MOON, SATURN_BARYCENTER, SUN};
use anise::constants::frames::IAU_EARTH_FRAME;
use nyx::cosmic::Orbit;
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::io::ExportCfg;
use nyx::md::StateParameter;
use nyx::od::prelude::*;
use nyx::propagators::{IntegratorOptions, Propagator};
use nyx::time::{Epoch, TimeUnits, Unit};
use nyx::utils::rss_orbit_errors;
use nyx::Spacecraft;
use nyx_space::mc::StateDispersion;
use nyx_space::propagators::IntegratorMethod;
use polars::prelude::*;
use std::collections::BTreeMap;
use std::env;
use std::fs::File;
use std::path::PathBuf;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

/*
 * These tests check that if we start with a state deviation in the estimate, the filter will eventually converge back.
 * These tests do NOT check that the filter will converge if the initial state in the propagator has that state deviation.
 * The latter would require iteration and smoothing before playing with an EKF. This will be handled in a subsequent version.
**/

#[allow(clippy::identity_op)]
#[rstest]
fn od_robust_test_ekf_realistic_one_way_cov_test(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();
    // Define the ground stations.
    let ekf_num_meas = 300;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 3 * Unit::Minute;
    let elevation_mask = 0.0;

    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        StochasticNoise::default_range_km(),
        StochasticNoise::default_doppler_km_s(),
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        StochasticNoise::default_range_km(),
        StochasticNoise::default_doppler_km_s(),
        iau_earth,
    );

    // Define the tracking configurations
    let configs = BTreeMap::from([
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
    let mut devices = BTreeMap::new();
    devices.insert("Madrid".to_string(), dss65_madrid);
    devices.insert("Canberra".to_string(), dss34_canberra);

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Spacecraft::from(Orbit::keplerian(
        22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k,
    ));

    // Disperse the initial state based on some orbital elements errors given from ULA Atlas 5 user guide, table 2.3.3-1 <https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf>
    // This assumes that the errors are ONE TENTH of the values given in the table. It assumes that the launch provider has provided an initial state vector, whose error is lower than the injection errors.
    // The initial covariance is computed based on the realized dispersions.
    let initial_estimate = KfEstimate::disperse_from_diag(
        initial_state,
        vec![
            StateDispersion::zero_mean(StateParameter::Inclination, 0.0025),
            StateDispersion::zero_mean(StateParameter::RAAN, 0.022),
            StateDispersion::zero_mean(StateParameter::AoP, 0.02),
        ],
        Some(0),
    )
    .unwrap();

    println!("Initial estimate:\n{}", initial_estimate);

    let initial_state_dev = initial_estimate.nominal_state;
    let (init_rss_pos_km, init_rss_vel_km_s) =
        rss_orbit_errors(&initial_state.orbit, &initial_state_dev.orbit);

    println!("Truth initial state:\n{initial_state}\n{initial_state:x}");
    println!("Filter initial state:\n{initial_state_dev}\n{initial_state_dev:x}");
    println!(
        "Initial state dev:\t{:.3} m\t{:.3} m/s\n{}",
        init_rss_pos_km * 1e3,
        init_rss_vel_km_s * 1e3,
        (initial_state.orbit - initial_state_dev.orbit).unwrap()
    );

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER, SATURN_BARYCENTER];
    let orbital_dyn = OrbitalDynamics::point_masses(bodies);
    let truth_setup = Propagator::new(
        SpacecraftDynamics::new(orbital_dyn),
        IntegratorMethod::RungeKutta4,
        opts,
    );
    let (_, traj) = truth_setup
        .with(initial_state, almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(devices.clone(), traj.clone(), configs, 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // And serialize to disk
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "output_data",
        "ekf_robust_msr.parquet",
    ]
    .iter()
    .collect();

    arc.to_parquet_simple(&path).unwrap();

    println!("{arc}\n{arc:?}");
    // Reload
    let reloaded = TrackingDataArc::from_parquet(&path).unwrap();
    assert_eq!(reloaded, arc);

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be _nearly_ perfect because we've removed SATURN_BARYCENTER from the estimated trajectory
    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let estimator = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let setup = Propagator::new(estimator, IntegratorMethod::RungeKutta4, opts);
    let prop_est = setup.with(initial_state_dev.with_stm(), almanac.clone());

    // Define the process noise to assume an unmodeled acceleration on X, Y and Z in the ECI frame
    let sigma_q = 5e-10_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);

    let kf = KF::new(initial_estimate, process_noise);

    let trig = EkfTrigger::new(ekf_num_meas, ekf_disable_time);

    let mut odp = SpacecraftODProcess::ekf(prop_est, kf, devices, trig, None, almanac);

    // Let's filter and iterate on the initial subset of the arc to refine the initial estimate
    let subset = arc.clone().filter_by_offset(..3.hours());
    let remaining = arc.filter_by_offset(3.hours()..);

    odp.process_arc(&subset).unwrap();
    odp.iterate_arc(&subset, IterationConf::once()).unwrap();

    let (sm_rss_pos_km, sm_rss_vel_km_s) =
        rss_orbit_errors(&initial_state.orbit, &odp.estimates[0].state().orbit);

    println!(
        "Initial state error after smoothing:\t{:.3} m\t{:.3} m/s\n{}",
        sm_rss_pos_km * 1e3,
        sm_rss_vel_km_s * 1e3,
        (initial_state.orbit - odp.estimates[0].state().orbit).unwrap()
    );

    odp.process_arc(&remaining).unwrap();

    odp.to_parquet(
        &remaining,
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
        final_truth_state.epoch() == est.epoch(),
        (final_truth_state.orbit - est.state().orbit).unwrap()
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

    assert_eq!(
        final_truth_state.epoch(),
        est.epoch(),
        "time of final EST and TRUTH epochs differ"
    );
    let delta = (est.state().orbit - final_truth_state.orbit).unwrap();
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
#[rstest]
fn od_robust_test_ekf_realistic_two_way(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();
    // Define the ground stations.
    let ekf_num_meas = 300;
    // Set the disable time to be very low to test enable/disable sequence
    let ekf_disable_time = 3 * Unit::Minute;
    let elevation_mask = 0.0;

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Spacecraft::from(Orbit::keplerian(
        22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k,
    ));

    let mut dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        StochasticNoise::default_range_km(),
        StochasticNoise::default_doppler_km_s(),
        iau_earth,
    );
    // Set the integration time so as to generate two way measurements
    dss65_madrid.integration_time = Some(60.seconds());
    let mut dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        StochasticNoise::default_range_km(),
        StochasticNoise::default_doppler_km_s(),
        iau_earth,
    );
    dss34_canberra.integration_time = Some(60.seconds());

    // Define the tracking configurations
    let configs = BTreeMap::from([
        (dss65_madrid.name.clone(), TrkConfig::default()),
        (dss34_canberra.name.clone(), TrkConfig::default()),
    ]);

    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    // Note that we do not have Goldstone so we can test enabling and disabling the EKF.
    let mut devices = BTreeMap::new();
    devices.insert("Madrid".to_string(), dss65_madrid);
    devices.insert("Canberra".to_string(), dss34_canberra);

    // Disperse the initial state based on some orbital elements errors given from ULA Atlas 5 user guide, table 2.3.3-1 <https://www.ulalaunch.com/docs/default-source/rockets/atlasvusersguide2010a.pdf>
    // This assumes that the errors are ONE TENTH of the values given in the table. It assumes that the launch provider has provided an initial state vector, whose error is lower than the injection errors.
    // The initial covariance is computed based on the realized dispersions.
    let initial_estimate = KfEstimate::disperse_from_diag(
        initial_state,
        vec![
            StateDispersion::zero_mean(StateParameter::Inclination, 0.0025),
            StateDispersion::zero_mean(StateParameter::RAAN, 0.022),
            StateDispersion::zero_mean(StateParameter::AoP, 0.02),
        ],
        Some(0),
    )
    .unwrap();

    println!("Initial estimate:\n{}", initial_estimate);

    let initial_state_dev = initial_estimate.nominal_state;
    let (init_rss_pos_km, init_rss_vel_km_s) =
        rss_orbit_errors(&initial_state.orbit, &initial_state_dev.orbit);

    println!("Truth initial state:\n{initial_state}\n{initial_state:x}");
    println!("Filter initial state:\n{initial_state_dev}\n{initial_state_dev:x}");
    println!(
        "Initial state dev:\t{:.3} m\t{:.3} m/s\n{}",
        init_rss_pos_km * 1e3,
        init_rss_vel_km_s * 1e3,
        (initial_state.orbit - initial_state_dev.orbit).unwrap()
    );

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER, SATURN_BARYCENTER];
    let orbital_dyn = OrbitalDynamics::point_masses(bodies);
    let truth_setup = Propagator::default(SpacecraftDynamics::new(orbital_dyn));
    let (_, traj) = truth_setup
        .with(initial_state, almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(devices.clone(), traj.clone(), configs, 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // And serialize to disk
    let path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "output_data"].iter().collect();

    traj.to_parquet_simple(
        path.join("ekf_robust_two_way_traj.parquet"),
        almanac.clone(),
    )
    .unwrap();
    arc.to_parquet_simple(path.join("ekf_robust_two_way_msr.parquet"))
        .unwrap();

    println!("{arc}");

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be _nearly_ perfect because we've removed SATURN_BARYCENTER from the estimated trajectory
    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let estimator = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let setup = Propagator::default(estimator);
    let prop_est = setup.with(initial_state_dev.with_stm(), almanac.clone());

    // Define the process noise to assume an unmodeled acceleration on X, Y and Z in the ECI frame
    let sigma_q = 5e-10_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);

    let kf = KF::new(initial_estimate, process_noise);

    let trig = EkfTrigger::new(ekf_num_meas, ekf_disable_time);

    let mut odp = SpacecraftODProcess::ekf(prop_est, kf, devices, trig, None, almanac);

    // Check that exporting an empty set returns an error.
    assert!(odp
        .to_parquet(
            &arc,
            path.with_file_name("robustness_test_two_way.parquet"),
            ExportCfg::timestamped(),
        )
        .is_err());

    odp.process_arc(&arc).unwrap();

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
            &arc,
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
            "Epoch (UTC)",
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
            "Covariance X*X (Earth J2000) (km^2)",
            "Covariance X*Y (Earth J2000) (km^2)",
            "Covariance X*Z (Earth J2000) (km^2)",
            "Covariance X*Vx (Earth J2000) (km^2/s)",
            "Covariance X*Vy (Earth J2000) (km^2/s)",
            "Covariance X*Vz (Earth J2000) (km^2/s)",
            "Covariance X*Cr (Earth J2000) (km)",
            "Covariance X*Cd (Earth J2000) (km)",
            "Covariance X*Mass (Earth J2000) (km*kg)",
            "Covariance Y*Y (Earth J2000) (km^2)",
            "Covariance Y*Z (Earth J2000) (km^2)",
            "Covariance Y*Vx (Earth J2000) (km^2/s)",
            "Covariance Y*Vy (Earth J2000) (km^2/s)",
            "Covariance Y*Vz (Earth J2000) (km^2/s)",
            "Covariance Y*Cr (Earth J2000) (km)",
            "Covariance Y*Cd (Earth J2000) (km)",
            "Covariance Y*Mass (Earth J2000) (km*kg)",
            "Covariance Z*Z (Earth J2000) (km^2)",
            "Covariance Z*Vx (Earth J2000) (km^2/s)",
            "Covariance Z*Vy (Earth J2000) (km^2/s)",
            "Covariance Z*Vz (Earth J2000) (km^2/s)",
            "Covariance Z*Cr (Earth J2000) (km)",
            "Covariance Z*Cd (Earth J2000) (km)",
            "Covariance Z*Mass (Earth J2000) (km*kg)",
            "Covariance Vx*Vx (Earth J2000) (km^2/s^2)",
            "Covariance Vx*Vy (Earth J2000) (km^2/s^2)",
            "Covariance Vx*Vz (Earth J2000) (km^2/s^2)",
            "Covariance Vx*Cr (Earth J2000) (km/s)",
            "Covariance Vx*Cd (Earth J2000) (km/s)",
            "Covariance Vx*Mass (Earth J2000) (km/s*kg)",
            "Covariance Vy*Vy (Earth J2000) (km^2/s^2)",
            "Covariance Vy*Vz (Earth J2000) (km^2/s^2)",
            "Covariance Vy*Cr (Earth J2000) (km/s)",
            "Covariance Vy*Cd (Earth J2000) (km/s)",
            "Covariance Vy*Mass (Earth J2000) (km/s*kg)",
            "Covariance Vz*Vz (Earth J2000) (km^2/s^2)",
            "Covariance Vz*Cr (Earth J2000) (km/s)",
            "Covariance Vz*Cd (Earth J2000) (km/s)",
            "Covariance Vz*Mass (Earth J2000) (km/s*kg)",
            "Covariance Cr*Cr (Earth J2000) (unitless)",
            "Covariance Cr*Cd (Earth J2000) (unitless)",
            "Covariance Cr*Mass (Earth J2000) (kg^2)",
            "Covariance Cd*Cd (Earth J2000) (unitless)",
            "Covariance Cd*Mass (Earth J2000) (kg^2)",
            "Covariance Mass*Mass (Earth J2000) (kg^2)",
            "Sigma X (Earth J2000) (km)",
            "Sigma Y (Earth J2000) (km)",
            "Sigma Z (Earth J2000) (km)",
            "Sigma Vx (Earth J2000) (km/s)",
            "Sigma Vy (Earth J2000) (km/s)",
            "Sigma Vz (Earth J2000) (km/s)",
            "Sigma Cr (Earth J2000) (unitless)",
            "Sigma Cd (Earth J2000) (unitless)",
            "Sigma Mass (Earth J2000) (kg)",
            "Sigma X (RIC) (km)",
            "Sigma Y (RIC) (km)",
            "Sigma Z (RIC) (km)",
            "Sigma Vx (RIC) (km/s)",
            "Sigma Vy (RIC) (km/s)",
            "Sigma Vz (RIC) (km/s)",
        ])
        .is_ok());

    // Check that the covariance deflated
    let est = &odp.estimates[odp.estimates.len() - 1];
    let final_truth_state = traj.at(est.epoch()).unwrap();

    println!("Estimate:\n{}", est);
    println!("Truth:\n{}", final_truth_state);
    println!(
        "Delta state with truth (epoch match: {}):\n{}",
        final_truth_state.epoch() == est.epoch(),
        (final_truth_state.orbit - est.state().orbit).unwrap()
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

    assert_eq!(
        final_truth_state.epoch(),
        est.epoch(),
        "time of final EST and TRUTH epochs differ"
    );
    let delta = (est.state().orbit - final_truth_state.orbit).unwrap();
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
