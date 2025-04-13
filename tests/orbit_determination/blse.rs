extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MOON, SUN};
use anise::constants::frames::IAU_EARTH_FRAME;
use nyx::cosmic::Orbit;
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::md::StateParameter;
use nyx::od::prelude::*;
use nyx::propagators::Propagator;
use nyx::time::{Epoch, TimeUnits, Unit};
use nyx::utils::rss_orbit_errors;
use nyx::Spacecraft;
use nyx_space::mc::StateDispersion;
use nyx_space::od::blse::BatchLeastSquares;
use std::collections::BTreeMap;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

/// Tests the robustness of the Batch least squares estimator against large initial state errors.
#[allow(clippy::identity_op)]
#[rstest]
fn blse_robust_large_disp_test(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();
    // Define the ground stations.
    let elevation_mask = 0.0;

    // BLSE usually runs over small periods of tracking data.
    let prop_time = 2 * Unit::Hour;

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_utc_hms(2020, 1, 1, 4, 0, 0);
    let initial_state = Spacecraft::from(Orbit::keplerian(
        22000.0, 0.01, 30.0, 80.0, 40.0, 45.0, dt, eme2k,
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

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
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

    // In a large Earth orbit, range data is _by far_ the more informative measurement type.
    let arc = arc
        .clone()
        .exclude_measurement_type(MeasurementType::Doppler);
    assert_eq!(arc.unique_types().len(), 1);
    assert_eq!(arc.unique_types()[0], MeasurementType::Range);

    let blse = BatchLeastSquares::builder()
        .solver(blse::BLSSolver::LevenbergMarquardt)
        .prop(truth_setup)
        .devices(devices)
        .almanac(almanac.clone())
        .build();

    let blse_solution = blse
        // .estimate(initial_estimate.nominal_state, &arc)
        .estimate(initial_state, &arc) // --> BUG: Even with the nominal state, the BLSE does not converge.
        .expect("blse should not fail");

    println!("{blse_solution}");

    assert_eq!(
        initial_state.epoch(),
        blse_solution.estimated_state.epoch(),
        "time of BLSE EST and TRUTH epochs differ"
    );
    let delta = (blse_solution.estimated_state.orbit - initial_state.orbit).unwrap();
    println!(
        "RMAG error = {:.6} m\tVMAG error = {:.6} m/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e3
    );
}
