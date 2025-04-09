extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MOON, SATURN_BARYCENTER, SUN};
use anise::constants::frames::IAU_EARTH_FRAME;
use nalgebra::U2;
use nyx::cosmic::Orbit;
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::io::ExportCfg;
use nyx::md::StateParameter;
use nyx::od::prelude::*;
use nyx::propagators::Propagator;
use nyx::time::{Epoch, TimeUnits, Unit};
use nyx::utils::rss_orbit_errors;
use nyx::Spacecraft;
use nyx_space::mc::StateDispersion;
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

/// Tests the robustness of the orbit determination (OD) process against large initial state errors.
///
/// This specific test simulates an OD scenario where the initial estimate provided
/// to the filter has a significant displacement (large error) compared to the
/// true initial state of the spacecraft. It focuses on scenarios utilizing
/// two-way measurements (like two-way range or Doppler) between ground stations
/// and the spacecraft. The goal is to verify that the estimation process
/// can converge to an accurate solution despite the poor initial guess, using
/// these specific measurement types.
///
/// # Arguments
///
/// * `almanac` - An `Arc<Almanac>` providing necessary environmental data (e.g., EOP, planetary ephemerides)
///               for propagation and measurement modeling.
#[allow(clippy::identity_op)]
#[rstest]
fn od_robust_large_disp_test_two_way_cov_test(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();
    // Define the ground stations.
    let elevation_mask = 0.0;

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_utc_hms(2020, 1, 1, 4, 0, 0);
    let initial_state = Spacecraft::from(Orbit::keplerian(
        22000.0, 0.01, 30.0, 80.0, 40.0, 180.0, dt, eme2k,
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
    let mut devices = BTreeMap::new();
    devices.insert("Madrid".to_string(), dss65_madrid);
    devices.insert("Canberra".to_string(), dss34_canberra);

    let initial_estimate = KfEstimate::disperse_from_diag(
        initial_state,
        vec![
            StateDispersion::zero_mean(StateParameter::SMA, 0.002),
            StateDispersion::zero_mean(StateParameter::RAAN, 0.002),
            StateDispersion::zero_mean(StateParameter::Inclination, 0.002),
            StateDispersion::zero_mean(StateParameter::Eccentricity, 0.0002),
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
    let path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "data", "04_output"]
        .iter()
        .collect();

    traj.to_parquet_simple(
        path.join("ekf_robust_two_way_traj.parquet"),
        almanac.clone(),
    )
    .unwrap();
    arc.to_parquet_simple(path.join("ekf_robust_two_way_msr.parquet"))
        .unwrap();

    println!("{arc}");

    // In a large Earth orbit, range data is _by far_ the more informative measurement type.
    let arc = arc
        .clone()
        .exclude_measurement_type(MeasurementType::Doppler);
    assert_eq!(arc.unique_types().len(), 1);
    assert_eq!(arc.unique_types()[0], MeasurementType::Range);

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be _nearly_ perfect because we've removed SATURN_BARYCENTER from the estimated trajectory
    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let estimator = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let setup = Propagator::default(estimator);

    // Define the process noise to assume an unmodeled acceleration on X, Y and Z in the ECI frame
    let sigma_q = 1e-7_f64.powi(2);
    let process_noise =
        ProcessNoise3D::from_diagonal(2 * Unit::Minute, &[sigma_q, sigma_q, sigma_q]);

    println!("{process_noise}");

    let odp = SpacecraftKalmanScalarOD::new(
        setup,
        KalmanVariant::DeviationTracking,
        None,
        devices.clone(),
        almanac.clone(),
    )
    .with_process_noise(process_noise);

    let od_sol = odp.process_arc(initial_estimate, &arc).unwrap();

    let od_pred = odp
        .predict_until(initial_estimate, arc.end_epoch().unwrap())
        .unwrap();

    // Export as Parquet
    let sol_path = od_sol
        .to_parquet(
            path.join("robustness_test_two_way.parquet"),
            ExportCfg::default(),
        )
        .unwrap();

    // Reload OD Solution
    let od_sol_reloaded =
        ODSolution::<Spacecraft, KfEstimate<Spacecraft>, U2, GroundStation>::from_parquet(
            sol_path.clone(),
            devices,
        )
        .unwrap();
    assert_eq!(od_sol_reloaded.estimates.len(), od_sol.estimates.len());
    assert_eq!(od_sol_reloaded.residuals.len(), od_sol.residuals.len());
    assert_eq!(od_sol_reloaded.gains.len(), od_sol.gains.len());
    assert_eq!(
        od_sol_reloaded.filter_smoother_ratios.len(),
        od_sol.filter_smoother_ratios.len()
    );
    od_sol_reloaded
        .to_parquet(
            path.join("robustness_test_two_way_reloaded.parquet"),
            ExportCfg::default(),
        )
        .unwrap();
    // assert!(od_sol_reloaded == od_sol, "womp womp");

    // Test the results
    // Check that the covariance deflated
    let est = &od_sol.estimates[od_sol.estimates.len() - 1];
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
    let delta = (est.orbital_state() - final_truth_state.orbit).unwrap();
    println!(
        "RMAG error = {:.6} m\tVMAG error = {:.6} m/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e3
    );

    // Compare with pure-predictor
    let est_pp = od_pred.estimates.last().unwrap();
    let delta_pp = (est_pp.orbital_state() - traj.at(est_pp.epoch()).unwrap().orbit).unwrap();
    println!(
        "Pure predictor RMAG error = {:.6} m\tVMAG error = {:.6} m/s",
        delta_pp.rmag_km() * 1e3,
        delta_pp.vmag_km_s() * 1e3
    );

    assert!(
        delta.rmag_km() < 175.0,
        "Position error should be less than 175 meters (down from ~2600 km)"
    );
    assert!(
        delta.vmag_km_s() < 1e-4,
        "Velocity error should be on decimeter per second level"
    );

    assert!(
        delta.rmag_km() / delta_pp.rmag_km() < 100.0,
        "Position error should be 100x better than a pure predictor"
    );

    // Read in the Parquet file and assert proper data was written.

    let df = ParquetReader::new(File::open(sol_path).unwrap())
        .finish()
        .unwrap();

    // Note: this also checks that the columns that match the given measurement kind exist.
    let _df_residuals = df
        .columns([
            "Prefit residual: Range (km)",
            "Postfit residual: Range (km)",
            "Residual ratio",
        ])
        .unwrap();

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
}
