extern crate nyx_space as nyx;

use anise::constants::celestial_objects::JUPITER_BARYCENTER;
use anise::constants::celestial_objects::MOON;
use anise::constants::celestial_objects::SUN;
use anise::constants::frames::IAU_EARTH_FRAME;
use nyx::od::simulator::TrackingArcSim;
use nyx::od::simulator::TrkConfig;

use self::nyx::md::prelude::*;
use self::nyx::od::prelude::*;

// Extra testing imports
use self::nyx::propagators::RK4Fixed;
use nyx::linalg::{SMatrix, SVector};
use std::collections::BTreeMap;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[allow(clippy::identity_op)]
#[rstest]
fn od_val_multi_body_ckf_perfect_stations(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();

    // Define the ground stations.
    let elevation_mask = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        StochasticNoise::MIN,
        StochasticNoise::MIN,
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        StochasticNoise::MIN,
        StochasticNoise::MIN,
        iau_earth,
    );
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        StochasticNoise::MIN,
        StochasticNoise::MIN,
        iau_earth,
    );

    // Define the tracking configurations
    let mut configs = BTreeMap::new();
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
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let orbital_dyn = OrbitalDynamics::point_masses(bodies);
    // Generate the truth data.
    let setup = Propagator::new::<RK4Fixed>(SpacecraftDynamics::new(orbital_dyn), opts);
    let mut prop = setup.with(initial_state.into(), almanac.clone());
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs, 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());
    let covar_radius_km = 1.0e-3_f64.powi(2);
    let covar_velocity_km_s = 1.0e-6_f64.powi(2);
    let init_covar = SMatrix::<f64, 9, 9>::from_diagonal(&SVector::<f64, 9>::from_iterator([
        covar_radius_km,
        covar_radius_km,
        covar_radius_km,
        covar_velocity_km_s,
        covar_velocity_km_s,
        covar_velocity_km_s,
        0.0,
        0.0,
        0.0,
    ]));

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(initial_state.into(), init_covar);

    let ckf = KF::no_snc(initial_estimate);

    let mut odp = ODProcess::<_, _, RangeDoppler, _, _, _>::ckf(prop_est, ckf, None, almanac);

    odp.process_arc::<GroundStation>(&arc).unwrap();

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
            est.state_deviation().norm() < 2e-16,
            "estimate error should be very good (perfect dynamics) ({:e})",
            est.state_deviation().norm()
        );

        last_est = Some(est);
    }

    for res in odp.residuals.iter().flatten() {
        assert!(
            res.postfit.norm() < 2e-16,
            "postfit should be zero (perfect dynamics) ({:e})",
            res
        );
    }

    let est = last_est.unwrap();
    assert!(est.state_deviation().norm() < 2e-16);
    assert!(est.covar.norm() < 1e-5);

    let delta = (est.state().orbit - final_truth.orbit).unwrap();
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 2e-16, "Position error should be zero");
    assert!(delta.vmag_km_s() < 2e-16, "Velocity error should be zero");
}

#[allow(clippy::identity_op)]
#[rstest]
fn multi_body_ckf_covar_map(almanac: Arc<Almanac>) {
    // For this test, we're only enabling one station so we can check that the covariance inflates between visibility passes.
    let _ = pretty_env_logger::try_init();

    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();
    // Define the ground stations.
    let elevation_mask = 0.0;
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        StochasticNoise::MIN,
        StochasticNoise::MIN,
        iau_earth,
    );
    // Define the tracking configurations
    let mut configs = BTreeMap::new();
    configs.insert(
        dss13_goldstone.name.clone(),
        TrkConfig::builder()
            .sampling(10.seconds())
            .scheduler(Scheduler::builder().sample_alignment(10.seconds()).build())
            .build(),
    );

    let all_stations = vec![dss13_goldstone];

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // Generate the truth data on one thread.
    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let orbital_dyn = OrbitalDynamics::point_masses(bodies);
    let setup = Propagator::new::<RK4Fixed>(SpacecraftDynamics::new(orbital_dyn), opts);
    let mut prop = setup.with(initial_state.into(), almanac.clone());

    let (_, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs, 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(Spacecraft::from(initial_state).with_stm(), almanac.clone());
    let covar_radius_km = 1.0e-3_f64.powi(2);
    let covar_velocity_km_s = 1.0e-6_f64.powi(2);
    let init_covar = SMatrix::<f64, 9, 9>::from_diagonal(&SVector::<f64, 9>::from_iterator([
        covar_radius_km,
        covar_radius_km,
        covar_radius_km,
        covar_velocity_km_s,
        covar_velocity_km_s,
        covar_velocity_km_s,
        0.0,
        0.0,
        0.0,
    ]));

    // Define the initial estimate
    let initial_estimate = KfEstimate::from_covar(initial_state.into(), init_covar);

    let ckf = KF::no_snc(initial_estimate);

    let mut odp = ODProcess::ckf(prop_est, ckf, None, almanac.clone());

    odp.process_arc::<GroundStation>(&arc).unwrap();

    let mut num_pred = 0_u32;
    for est in odp.estimates.iter() {
        if est.predicted {
            num_pred += 1;
        } else {
            // Only check that the covariance is low IF this isn't a predicted estimate
            assert!(
                est.state_deviation().norm() < 2e-16,
                "estimate error should be zero (perfect dynamics) ({:e})",
                est.state_deviation().norm()
            );
        }
        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0,
                "covar diagonal element negative @ [{}, {}]",
                i,
                i
            );
        }
    }

    // Note that we check the residuals separately from the estimates because we have many predicted estimates which do not have any associated residuals.
    for res in odp.residuals.iter().flatten() {
        assert!(
            res.postfit.norm() < 2e-16,
            "postfit should be zero (perfect dynamics) ({:e})",
            res
        );
    }

    assert!(num_pred > 0, "no predicted estimates");

    let est = odp.estimates.last().unwrap();

    println!("{:.2e}", est.state_deviation().norm());
    println!("{:.2e}", est.covar.norm());

    // Test that we can generate a navigation trajectory and search it
    let nav_traj = odp.to_traj().unwrap();
    let aop_event = Event::apoapsis();
    for found_event in nav_traj.find(&aop_event, almanac).unwrap() {
        println!("{:x}", found_event.state);
        assert!((found_event.state.orbit.ta_deg().unwrap() - 180.0).abs() < 1e-2)
    }
}
