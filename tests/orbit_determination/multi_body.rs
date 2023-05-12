extern crate nyx_space as nyx;

use nyx::od::noise::GaussMarkov;
use nyx::od::simulator::arc::TrackingArcSim;
use nyx::od::simulator::TrkConfig;

use self::nyx::md::ui::*;
use self::nyx::od::prelude::*;

// Extra testing imports
use self::nyx::linalg::{Matrix2, Matrix6, Vector2, Vector6};
use self::nyx::propagators::RK4Fixed;
use std::collections::HashMap;

#[allow(clippy::identity_op)]
#[test]
fn od_val_multi_body_ckf_perfect_stations() {
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

    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let orbital_dyn = OrbitalDynamics::point_masses(&bodies, cosm.clone());
    // Generate the truth data.
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
    let covar_radius_km = 1.0e-3_f64.powi(2);
    let covar_velocity_km_s = 1.0e-6_f64.powi(2);
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
    let measurement_noise =
        Matrix2::from_diagonal(&Vector2::new(15e-3_f64.powi(2), 1e-5_f64.powi(2)));

    let ckf = KF::no_snc(initial_estimate, measurement_noise);

    let mut odp =
        ODProcess::<_, _, RangeDoppler, _, _, _, _>::ckf(prop_est, ckf, RejectCriteria::None, cosm);

    odp.process_arc::<GroundStation>(&arc).unwrap();

    let mut wtr = csv::Writer::from_writer(io::stdout());
    let mut printed = false;
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

        if !printed {
            wtr.serialize(*est).expect("could not write to stdout");
            printed = true;
        }

        last_est = Some(est);
    }

    for res in &odp.residuals {
        assert!(
            res.postfit.norm() < 2e-16,
            "postfit should be zero (perfect dynamics) ({:e})",
            res
        );
    }

    let est = last_est.unwrap();
    assert!(est.state_deviation().norm() < 2e-16);
    assert!(est.covar.norm() < 1e-5);

    let delta = est.state() - final_truth;
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 2e-16, "Position error should be zero");
    assert!(delta.vmag_km_s() < 2e-16, "Velocity error should be zero");
}

#[allow(clippy::identity_op)]
#[test]
fn multi_body_ckf_covar_map() {
    // For this test, we're only enabling one station so we can check that the covariance inflates between visibility passes.
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::io;

    let cosm = Cosm::de438();

    let iau_earth = cosm.frame("IAU Earth");
    // Define the ground stations.
    let elevation_mask = 0.0;
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );
    // Define the tracking configurations
    let mut configs = HashMap::new();
    configs.insert(
        dss13_goldstone.name.clone(),
        TrkConfig::from_sample_rate(10.seconds()),
    );

    let all_stations = vec![dss13_goldstone];

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    // Generate the truth data on one thread.
    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let orbital_dyn = OrbitalDynamics::point_masses(&bodies, cosm.clone());
    let setup = Propagator::new::<RK4Fixed>(orbital_dyn, opts);
    let mut prop = setup.with(initial_state);

    let (_, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs, 0).unwrap();
    arc_sim.disallow_overlap(); // Prevent overlapping measurements

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let prop_est = setup.with(initial_state.with_stm());
    let covar_radius_km = 1.0e-3_f64.powi(2);
    let covar_velocity_km_s = 1.0e-6_f64.powi(2);
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
    let measurement_noise =
        Matrix2::from_diagonal(&Vector2::new(15e-3_f64.powi(2), 1e-5_f64.powi(2)));

    let ckf = KF::no_snc(initial_estimate, measurement_noise);

    let mut odp = ODProcess::ckf(prop_est, ckf, RejectCriteria::None, cosm);

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
    for res in odp.residuals.iter() {
        assert!(
            res.postfit.norm() < 2e-16,
            "postfit should be zero (perfect dynamics) ({:e})",
            res
        );
    }

    assert!(num_pred > 0, "no predicted estimates");

    let est = odp.estimates.last().unwrap();

    let mut wtr = csv::Writer::from_writer(io::stdout());
    wtr.serialize(*est).expect("could not write to stdout");

    println!("{:.2e}", est.state_deviation().norm());
    println!("{:.2e}", est.covar.norm());

    // Test that we can generate a navigation trajectory and search it
    let nav_traj = odp.to_traj().unwrap();
    let aop_event = Event::apoapsis();
    for found_event in nav_traj.find_all(&aop_event).unwrap() {
        println!("{:x}", found_event);
        assert!((found_event.ta_deg() - 180.0).abs() < 1e-2)
    }
}
