use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MOON, SATURN_BARYCENTER, SUN};
use anise::constants::frames::{EARTH_J2000, IAU_EARTH_FRAME};
use nalgebra::U2;
use nyx_space::dynamics::guidance::LocalFrame;
use pretty_env_logger::try_init;

use rand_distr::Distribution;
use rand_pcg::Pcg64Mcg;
use rstest::*;

use nyx_space::cosmic::Orbit;
use nyx_space::dynamics::orbital::OrbitalDynamics;
use nyx_space::md::prelude::*;
use nyx_space::od::prelude::*;
use nyx_space::propagators::{IntegratorMethod, IntegratorOptions, Propagator};
use nyx_space::time::{Epoch, TimeUnits};
use nyx_space::utils::rss_orbit_errors;
use std::collections::BTreeMap;
use std::path::PathBuf;

#[fixture]
fn epoch() -> Epoch {
    Epoch::from_gregorian_utc_at_midnight(2023, 1, 1)
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
    let prop_time = 8.hours();
    let step_size = 10.0.seconds();

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let initial_state = Spacecraft::builder()
        .orbit(Orbit::keplerian(
            22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, epoch, eme2k,
        ))
        .build();

    println!("{initial_state}");

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER, SATURN_BARYCENTER];
    let orbital_dyn = OrbitalDynamics::point_masses(bodies);
    let truth_setup = Propagator::dp78(
        SpacecraftDynamics::new(orbital_dyn),
        IntegratorOptions::with_max_step(step_size),
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
) -> (BTreeMap<String, GroundStation>, BTreeMap<String, TrkConfig>) {
    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();

    let elevation_mask = 0.0;

    // Define the ground stations
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
    let mut configs = BTreeMap::new();
    configs.insert(
        dss65_madrid.name.clone(),
        TrkConfig::builder()
            .strands(vec![Strand {
                start: epoch + 60.seconds(),
                end: epoch + 2.hours(),
            }])
            .build(),
    );
    configs.insert(
        dss34_canberra.name.clone(),
        TrkConfig::builder()
            .strands(vec![Strand {
                start: epoch + 4.hours(),
                end: epoch + 6.hours(),
            }])
            .build(),
    );

    let mut devices = BTreeMap::new();
    devices.insert("Madrid".to_string(), dss65_madrid);
    devices.insert("Canberra".to_string(), dss34_canberra);

    (devices, configs)
}

#[fixture]
fn tracking_arc(
    traj: Traj<Spacecraft>,
    devices_n_configs: (BTreeMap<String, GroundStation>, BTreeMap<String, TrkConfig>),
    almanac: Arc<Almanac>,
) -> TrackingDataArc {
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

    let sc = SpacecraftUncertainty::builder()
        .nominal(initial_state)
        .frame(LocalFrame::RIC)
        .x_km(1.0)
        .y_km(1.0)
        .z_km(1.0)
        .vx_km_s(0.5e-2)
        .vy_km_s(0.5e-2)
        .vz_km_s(0.5e-2)
        .build();

    let sc_estimate = sc.to_estimate().unwrap();

    // Now, let's sample from this.
    let sc_gen = sc_estimate.to_random_variable().unwrap();
    let mut rng = Pcg64Mcg::new(123); // Set the seed for reproducibility of test
    let mut initial_estimate = sc_estimate;
    initial_estimate.nominal_state = sc_gen.sample(&mut rng).state;

    println!("Initial estimate:\n{}", initial_estimate);
    println!(
        "Initial Keplerian covar:\n{:.6}",
        initial_estimate.keplerian_covar()
    );

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

#[ignore = "large state deviations to be handled in #416"]
#[rstest]
fn od_resid_reject_inflated_snc_ckf_two_way(
    traj: Traj<Spacecraft>,
    tracking_arc: TrackingDataArc,
    initial_estimate: KfEstimate<Spacecraft>,
    devices_n_configs: (BTreeMap<String, GroundStation>, BTreeMap<String, TrkConfig>),
    almanac: Arc<Almanac>,
) {
    let (devices, _configs) = devices_n_configs;

    let initial_state_dev = initial_estimate.nominal_state;

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let estimator = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let setup = Propagator::new(
        estimator,
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.seconds()),
    );
    let prop_est = setup.with(initial_state_dev.with_stm(), almanac.clone());

    // Define the process noise to assume an unmodeled acceleration on X, Y and Z in the ECI frame
    let sigma_q = 5e-9_f64.powi(2);
    let process_noise = ProcessNoise3D::from_diagonal(2.minutes(), &[sigma_q, sigma_q, sigma_q]);

    let kf =
        KF::new(initial_estimate, KalmanVariant::ReferenceUpdate).with_process_noise(process_noise);

    // ==> TEST <== //
    // We greatly inflate the SNC so that the covariance inflates tremendously. This leads to the
    // measurement noise being tremendously inflated as well, but it also means that all of
    // the measurements are accepted.
    // So we end up with an excellent estimate but an unusably high covariance.

    let mut odp = ODProcess::<_, U2, _, _, _>::new(
        prop_est,
        kf,
        devices,
        Some(ResidRejectCrit { num_sigmas: 2.0 }), // 95% to force rejections
        almanac,
    );

    let od_sol = odp.process_arc(&tracking_arc).unwrap();

    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "04_output",
        "resid_reject_inflated_snc.parquet",
    ]
    .iter()
    .collect();

    od_sol.to_parquet(path, ExportCfg::timestamped()).unwrap();

    let num_rejections = od_sol
        .residuals
        .iter()
        .flatten()
        .filter(|residual| residual.rejected)
        .count();

    // Check that the final post-fit residual isn't too bad, and definitely much better than the prefit.
    let est = &od_sol.estimates.last().unwrap();
    // BUG ? The prefit to postfit ratio seems to be a near consistent factor of two. Even with a postfit of zero, the state deviations are high: 8.4 km and 11.7 m/s.
    // This is the case both with a sequential filter and a dual filter. SNC has no effect. Neither does the seed.
    for (ith, (est, opt_resid)) in od_sol.results().enumerate() {
        if let Some(resid) = opt_resid {
            println!("RESIDUAL #{ith}");
            let truth_state = traj.at(resid.epoch).unwrap();
            let est_state = est.nominal_state() + est.state_deviation;

            let delta = (est_state.orbit - truth_state.orbit).unwrap();
            println!(
                "RMAG error = {:.6} m\tVMAG error = {:.6} m/s",
                delta.rmag_km() * 1e3,
                delta.vmag_km_s() * 1e3
            );

            println!("{truth_state:x}");
            println!("{est_state:x}");

            if ith > 50 {
                break;
            }
        }
    }
    let final_resid = &od_sol.residuals.last().unwrap().as_ref().unwrap();
    let final_truth_state = traj.at(est.epoch()).unwrap();

    println!("FINAL\n{final_resid}");

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

    assert!(
        delta.rmag_km() < 15.0,
        "Position error should be less than 15 km"
    );
    assert!(
        delta.vmag_km_s() < 2e-3,
        "Velocity error should be on meter per second level"
    );

    assert!(!final_resid.rejected, "final residual should be accepted");
    assert!(final_resid.prefit.norm() > final_resid.postfit.norm());
    assert!(final_resid.postfit.norm() < 1e-2);
    assert!(dbg!(num_rejections) < 30, "oddly high rejections");
}

#[ignore = "large state deviations to be handled in #416"]
#[rstest]
fn od_resid_reject_default_ckf_two_way_cov_test(
    traj: Traj<Spacecraft>,
    tracking_arc: TrackingDataArc,
    initial_estimate: KfEstimate<Spacecraft>,
    devices_n_configs: (BTreeMap<String, GroundStation>, BTreeMap<String, TrkConfig>),
    almanac: Arc<Almanac>,
) {
    let (devices, _configs) = devices_n_configs;

    let initial_state_dev = initial_estimate.nominal_state;

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be _nearly_ perfect because we've removed SATURN_BARYCENTER from the estimated trajectory
    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let estimator = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let setup = Propagator::new(
        estimator,
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.seconds()),
    );
    let prop_est = setup.with(initial_state_dev.with_stm(), almanac.clone());

    // Define the process noise to assume an unmodeled acceleration on X, Y and Z in the ECI frame
    let sigma_q = 5e-10_f64.powi(2);
    let process_noise = ProcessNoise3D::from_diagonal(2.minutes(), &[sigma_q, sigma_q, sigma_q]);

    let kf =
        KF::new(initial_estimate, KalmanVariant::ReferenceUpdate).with_process_noise(process_noise);

    let mut odp = ODProcess::<_, U2, _, _, _>::new(
        prop_est,
        kf,
        devices,
        Some(ResidRejectCrit::default()),
        almanac,
    );

    let od_sol = odp.process_arc(&tracking_arc).unwrap();

    // Save this result before the asserts for analysis
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "04_output",
        "resid_reject_test.parquet",
    ]
    .iter()
    .collect();
    od_sol.to_parquet(path, ExportCfg::timestamped()).unwrap();

    let num_rejections = od_sol
        .residuals
        .iter()
        .flatten()
        .filter(|residual| residual.rejected)
        .count();

    assert!(num_rejections > 220);

    // Check that the error is within the covariance.
    let est = &od_sol.estimates.last().unwrap();
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

    assert!(delta.radius_km.x < est.covar[(0, 0)].sqrt());
    assert!(delta.radius_km.y < est.covar[(1, 1)].sqrt());
    assert!(delta.radius_km.z < est.covar[(2, 2)].sqrt());

    assert!(delta.velocity_km_s.x < est.covar[(3, 3)].sqrt());
    assert!(delta.velocity_km_s.y < est.covar[(4, 4)].sqrt());
    assert!(delta.velocity_km_s.z < est.covar[(5, 5)].sqrt());
}
