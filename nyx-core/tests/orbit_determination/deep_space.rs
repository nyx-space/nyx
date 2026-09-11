extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use anise::analysis::prelude::OrbitalElement;
use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MOON, SATURN_BARYCENTER, SUN};
use nyx::Spacecraft;
use nyx::cosmic::Orbit;
use nyx::dynamics::SpacecraftDynamics;
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::io::ExportCfg;
use nyx::md::StateParameter;
use nyx::od::prelude::*;
use nyx::propagators::Propagator;
use nyx::time::{Epoch, Unit};
use nyx::utils::rss_orbit_errors;
use nyx_space::cosmic::{Mass, SRPData};
use nyx_space::mc::StateDispersion;
use std::collections::BTreeMap;
use std::env;
use std::path::PathBuf;

use anise::{constants::frames::MOON_J2000, prelude::Almanac};
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
fn od_moon_shapiro_light_time(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    // Define the ground stations.
    let elevation_mask = 5.0;

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;

    // Define state information.
    let moon_j2k = almanac.frame_info(MOON_J2000).unwrap();
    let epoch = Epoch::from_gregorian_utc_hms(2022, 7, 22, 3, 2, 1);
    let orbit = Orbit::try_keplerian_altitude(175.0, 1e-3, 51.9, 45.0, 75.0, 90.0, epoch, moon_j2k)
        .unwrap();

    let initial_state = Spacecraft::builder()
        .orbit(orbit)
        .srp(SRPData::from_area(3.21))
        .mass(Mass::from_dry_mass(159.0))
        .build();

    let mut dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        StochasticNoise::default_range_km(),
        StochasticNoise::default_doppler_km_s(),
    );
    // Set the integration time so as to generate two way measurements
    dss65_madrid.doppler_config = Some(DopplerConfig::default());
    dss65_madrid.light_time_correction = true;
    dss65_madrid.relativistic_corrections = true;
    let mut dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        StochasticNoise::default_range_km(),
        StochasticNoise::default_doppler_km_s(),
    );
    dss34_canberra.doppler_config = Some(DopplerConfig::default());
    dss34_canberra.light_time_correction = true;
    dss34_canberra.relativistic_corrections = true;

    // Define the tracking configurations
    let configs = BTreeMap::from([
        (dss65_madrid.name.clone(), TrkConfig::default()),
        (dss34_canberra.name.clone(), TrkConfig::default()),
    ]);

    let mut devices = BTreeMap::new();
    devices.insert("Madrid".to_string(), dss65_madrid);
    devices.insert("Canberra".to_string(), dss34_canberra);

    let initial_estimate = KfEstimate::from_dispersions(
        initial_state,
        vec![
            StateDispersion::zero_mean(
                StateParameter::Element(OrbitalElement::SemiMajorAxis),
                0.002,
            ),
            StateDispersion::zero_mean(StateParameter::Element(OrbitalElement::RAAN), 0.002),
            StateDispersion::zero_mean(StateParameter::Element(OrbitalElement::Inclination), 0.002),
            StateDispersion::zero_mean(
                StateParameter::Element(OrbitalElement::Eccentricity),
                0.0002,
            ),
        ],
        Some(123456),
    )
    .unwrap();

    println!("Initial estimate:\n{initial_estimate}");

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
    arc_sim.build_schedule(&almanac).unwrap();

    let arc = arc_sim.generate_measurements(&almanac).unwrap();

    // And serialize to disk
    let path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "../data", "04_output"]
        .iter()
        .collect();

    println!("{arc}");

    // Now that we have the truth data, let's start an OD and compute the estimates. We expect the
    // estimated orbit to be _nearly_ perfect because we've removed SATURN_BARYCENTER from the
    // estimated trajectory
    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let estimator = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let setup = Propagator::default(estimator);

    // Define the process noise to assume an unmodeled acceleration on X, Y and Z in the ECI frame
    let sigma_q = 1e-7_f64.powi(2);
    let process_noise =
        ProcessNoise3D::from_diagonal(&[sigma_q, sigma_q, sigma_q], 2 * Unit::Minute, None);

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

    // Export as Parquet
    od_sol
        .to_parquet(path.join("od_moon_shapiro.parquet"), ExportCfg::default())
        .unwrap();

    // Export ephemeris
    let ephem = od_sol.to_ephemeris("My Spacecraft".to_string());
    // Check that the final covariance is PSD by rotating it into another frame
    let covar_ric = ephem
        .covar_at(
            ephem.end_epoch().unwrap(),
            anise::ephemerides::ephemeris::LocalFrame::RIC,
            &almanac,
        )
        .expect("non PSD covariance?")
        .unwrap();
    println!("{covar_ric}");

    // Test the results
    // Check that the covariance deflated
    let est = &od_sol.estimates[od_sol.estimates.len() - 1];
    let final_truth_state = traj.at(est.epoch()).unwrap();

    println!("Estimate:\n{est}");
    println!("Truth:\n{final_truth_state}");
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

    assert!(
        delta.rmag_km() * 1e-3 < 75.0,
        "Position error should be less than 175 meters (down from ~2600 km)"
    );
    assert!(
        delta.vmag_km_s() < 1e-4,
        "Velocity error should be on the 10 cm/s per second level"
    );
}
