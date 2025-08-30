extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MOON, SUN};
use anise::constants::frames::IAU_EARTH_FRAME;
use nyx::cosmic::Orbit;
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::md::StateParameter;
use nyx::od::blse::*;
use nyx::od::prelude::*;
use nyx::propagators::Propagator;
use nyx::time::Epoch;
use nyx::utils::rss_orbit_errors;
use nyx::Spacecraft;
use nyx_space::mc::StateDispersion;
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
#[case(BLSSolver::NormalEquations, 60.seconds(), 2.minutes(), false)]
#[case(BLSSolver::LevenbergMarquardt, 10.seconds(), 10.minutes(), false)]
#[case(BLSSolver::NormalEquations, 60.seconds(), 2.minutes(), true)]
#[case(BLSSolver::LevenbergMarquardt, 10.seconds(), 10.minutes(), true)]
fn blse_robust_large_disp_cov_test(
    #[case] solver: BLSSolver,
    #[case] sample: Duration,
    #[case] offset: Duration,
    #[case] disperse: bool,
    almanac: Arc<Almanac>,
) {
    let _ = pretty_env_logger::try_init();

    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();
    // Define the ground stations.
    let elevation_mask = 0.0;

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_utc_hms(2020, 1, 1, 4, 0, 0);
    let initial_state = Spacecraft::from(Orbit::keplerian(
        22000.0, 0.01, 30.0, 80.0, 40.0, 170.0, dt, eme2k,
    ));

    let prop_time = initial_state.orbit.period().unwrap();

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
            TrkConfig::from_sample_rate(sample),
        ),
        (
            dss34_canberra.name.clone(),
            TrkConfig::from_sample_rate(sample),
        ),
    ]);

    let mut devices = BTreeMap::new();
    // devices.insert("Madrid".to_string(), dss65_madrid);
    devices.insert("Canberra".to_string(), dss34_canberra);

    let initial_estimate = KfEstimate::disperse_from_diag(
        initial_state,
        vec![
            StateDispersion::zero_mean(StateParameter::SMA, 0.02),
            StateDispersion::zero_mean(StateParameter::RAAN, 0.02),
            StateDispersion::zero_mean(StateParameter::Inclination, 0.02),
            StateDispersion::zero_mean(StateParameter::Eccentricity, 0.0002),
        ],
        Some(0),
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

    let blse = BatchLeastSquares::builder()
        .solver(solver)
        .prop(truth_setup)
        .devices(devices)
        .almanac(almanac.clone())
        .build();

    let blse_solution = blse
        .estimate(
            if disperse {
                initial_estimate.nominal_state
            } else {
                initial_state
            },
            &arc.filter_by_offset(..offset),
        )
        .expect("blse should not fail");

    println!("{blse_solution}");

    assert_eq!(
        initial_state.epoch(),
        blse_solution.estimated_state.epoch(),
        "time of BLSE EST and TRUTH epochs differ"
    );
    let delta = (blse_solution.estimated_state.orbit - initial_state.orbit).unwrap();
    println!(
        "RMAG error = {:.3} m\tVMAG error = {:.3} m/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e3
    );
    let rmag_km_imp = init_rss_pos_km - delta.rmag_km();
    let vmag_km_s_imp = init_rss_vel_km_s - delta.vmag_km_s();
    println!(
        "Improvement of |R| = {:.3} m\t|V| = {:.3} m/s",
        rmag_km_imp * 1e3,
        vmag_km_s_imp * 1e3,
    );

    if disperse {
        assert!(
            rmag_km_imp > 0.0,
            "Position estimate not any better after BLSE"
        );
    } else {
        // The RMAG error should be centimeter level, roughly the noise.
        assert!(
            delta.rmag_km() < 0.1,
            "Position estimate without dispersions too large"
        );
    }

    let kf_est: KfEstimate<Spacecraft> = blse_solution.into();
    println!("{kf_est}");
}
