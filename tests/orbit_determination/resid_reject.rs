use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MOON, SATURN_BARYCENTER, SUN};
use anise::constants::frames::{EARTH_J2000, IAU_EARTH_FRAME};
use pretty_env_logger::try_init;

use rstest::*;

use nyx_space::cosmic::Orbit;
use nyx_space::dynamics::orbital::OrbitalDynamics;
use nyx_space::linalg::{Matrix2, Vector2};
use nyx_space::md::prelude::*;
use nyx_space::md::StateParameter;
use nyx_space::od::noise::GaussMarkov;
use nyx_space::od::prelude::*;
use nyx_space::propagators::{PropOpts, Propagator, RK4Fixed};
use nyx_space::time::{Epoch, TimeUnits};
use nyx_space::utils::rss_orbit_errors;
use std::collections::BTreeMap;

#[fixture]
fn epoch() -> Epoch {
    Epoch::from_gregorian_tai_at_midnight(2023, 1, 1)
}

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[fixture]
fn traj(epoch: Epoch, almanac: Arc<Almanac>) -> Traj<Spacecraft> {
    let _ = try_init().is_err();

    // Define the propagator information.
    let prop_time = 2.hours();
    let step_size = 10.0.seconds();

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let initial_state = Spacecraft::builder()
        .orbit(Orbit::keplerian(
            22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, epoch, eme2k,
        ))
        .build();

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER, SATURN_BARYCENTER];
    let orbital_dyn = OrbitalDynamics::point_masses(bodies);
    let truth_setup = Propagator::dp78(
        SpacecraftDynamics::new(orbital_dyn),
        PropOpts::with_max_step(step_size),
    );
    let (_, traj) = truth_setup
        .with(initial_state, almanac)
        .for_duration_with_traj(prop_time)
        .unwrap();

    traj
}

#[fixture]
fn devices_n_configs(
    epoch: Epoch,
    almanac: Arc<Almanac>,
) -> (Vec<GroundStation>, BTreeMap<String, TrkConfig>) {
    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();

    let elevation_mask = 0.0;

    // Define the ground stations
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
    let mut configs = BTreeMap::new();
    let cfg = TrkConfig::builder()
        .strands(vec![Strand {
            start: epoch + 60.seconds(),
            end: epoch + 2.hours(),
        }])
        .build();
    configs.insert(dss65_madrid.name.clone(), cfg.clone());
    configs.insert(dss34_canberra.name.clone(), cfg);

    (vec![dss65_madrid, dss34_canberra], configs)
}

#[fixture]
fn tracking_arc(
    traj: Traj<Spacecraft>,
    devices_n_configs: (Vec<GroundStation>, BTreeMap<String, TrkConfig>),
    almanac: Arc<Almanac>,
) -> TrackingArc<RangeDoppler> {
    let (devices, configs) = devices_n_configs;

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(devices, traj, configs, 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac).unwrap();

    println!("{arc}");

    arc
}

#[fixture]
fn initial_estimate(traj: Traj<Spacecraft>) -> KfEstimate<Spacecraft> {
    let initial_state = *(traj.first());

    let initial_estimate = KfEstimate::disperse_from_diag(
        initial_state,
        &[
            (StateParameter::SMA, 30.0),
            (StateParameter::Inclination, 0.025),
            (StateParameter::RAAN, 0.22),
        ],
        Some(10),
    );
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

    initial_estimate
}

#[rstest]
fn od_resid_reject_all_ckf_two_way(
    tracking_arc: TrackingArc<RangeDoppler>,
    initial_estimate: KfEstimate<Spacecraft>,
    devices_n_configs: (Vec<GroundStation>, BTreeMap<String, TrkConfig>),
    almanac: Arc<Almanac>,
) {
    let (devices, _configs) = devices_n_configs;

    let initial_state_dev = initial_estimate.nominal_state;

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be _nearly_ perfect because we've removed SATURN_BARYCENTER from the estimated trajectory
    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let estimator = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let setup = Propagator::new::<RK4Fixed>(estimator, PropOpts::with_fixed_step(10.seconds()));
    let prop_est = setup.with(initial_state_dev.with_stm(), almanac.clone());

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    // Define the process noise to assume an unmodeled acceleration on X, Y and Z in the ECI frame
    let sigma_q = 5e-10_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2.minutes(), &[sigma_q, sigma_q, sigma_q]);

    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    // ==> TEST <== //
    // We set up the rejection criteria to start after a single measurement and with a residual ratio of 3.0, i.e. within 3-sigmas.
    // The initial dispersion is about 20388 km, so we none of the residuals will be within those 3-sigmas.
    // Therefore, the test is to confirm that the ODP only performs time updates and zero measurement updates.

    let mut odp = ODProcess::ckf(
        prop_est,
        kf,
        Some(FltResid {
            min_accepted: 0, // Start the preprocessing filter at the first measurement because we try to filter everything out in this test
            num_sigmas: 3.0,
        }),
        almanac,
    );

    // TODO: Fix the deserialization of the measurements such that they also deserialize the integration time.
    // Without it, we're stuck having to rebuild them from scratch.
    // https://github.com/nyx-space/nyx/issues/140

    // Build the hashmap of devices from the vector using their names
    let mut devices_map = devices
        .into_iter()
        .map(|dev| (dev.name.clone(), dev))
        .collect::<BTreeMap<_, _>>();

    odp.process(
        &tracking_arc.measurements,
        &mut devices_map,
        tracking_arc.min_duration_sep().unwrap(),
    )
    .unwrap();

    for residual in odp.residuals.iter().flatten() {
        assert!(residual.rejected, "{} was not rejected!", residual.epoch);
    }

    for estimate in odp.estimates.iter() {
        assert!(
            estimate.predicted,
            "{} was not predicted!",
            estimate.epoch()
        );
    }

    // We don't check the estimation because it'll be bad since we rejected all the measurements.
}

#[rstest]
fn od_resid_reject_default_ckf_two_way(
    traj: Traj<Spacecraft>,
    tracking_arc: TrackingArc<RangeDoppler>,
    initial_estimate: KfEstimate<Spacecraft>,
    devices_n_configs: (Vec<GroundStation>, BTreeMap<String, TrkConfig>),
    almanac: Arc<Almanac>,
) {
    let (devices, _configs) = devices_n_configs;

    let initial_state_dev = initial_estimate.nominal_state;

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be _nearly_ perfect because we've removed SATURN_BARYCENTER from the estimated trajectory
    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let estimator = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let setup = Propagator::new::<RK4Fixed>(estimator, PropOpts::with_fixed_step(10.seconds()));
    let prop_est = setup.with(initial_state_dev.with_stm(), almanac.clone());

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(1e-6, 1e-3));

    // Define the process noise to assume an unmodeled acceleration on X, Y and Z in the ECI frame
    let sigma_q = 5e-10_f64.powi(2);
    let process_noise = SNC3::from_diagonal(2.minutes(), &[sigma_q, sigma_q, sigma_q]);

    let kf = KF::new(initial_estimate, process_noise, measurement_noise);

    let mut odp = ODProcess::ckf(prop_est, kf, Some(FltResid::default()), almanac);

    // TODO: Fix the deserialization of the measurements such that they also deserialize the integration time.
    // Without it, we're stuck having to rebuild them from scratch.
    // https://github.com/nyx-space/nyx/issues/140

    // Build the hashmap of devices from the vector using their names
    let mut devices_map = devices
        .into_iter()
        .map(|dev| (dev.name.clone(), dev))
        .collect::<BTreeMap<_, _>>();

    odp.process(
        &tracking_arc.measurements,
        &mut devices_map,
        tracking_arc.min_duration_sep().unwrap(),
    )
    .unwrap();

    // With the default configuration, the filter converges very fast since we have similar dynamics.

    for residual in odp.residuals.iter().flatten() {
        assert!(!residual.rejected, "{} was rejected!", residual.epoch);
    }

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

    // We start with a 20 km error and with only 120 measurements, no iteration, and no extended KF, we achieve less than 500 meters of error.
    // That's not bad!

    assert!(
        delta.rmag_km() < 0.6,
        "Position error should be less than 500 meters"
    );
    assert!(
        delta.vmag_km_s() < 1e-3,
        "Velocity error should be on meter per second level"
    );
}
