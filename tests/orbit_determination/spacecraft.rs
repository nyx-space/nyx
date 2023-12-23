extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use nyx::cosmic::{Bodies, Cosm, Orbit, Spacecraft};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::spacecraft::{SolarPressure, SpacecraftDynamics};
use nyx::linalg::{Matrix2, Matrix6, Vector2, Vector6};
use nyx::md::trajectory::ExportCfg;
use nyx::md::{Event, StateParameter};
use nyx::od::noise::GaussMarkov;
use nyx::od::prelude::*;
use nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use nyx::time::{Epoch, TimeUnits, Unit};
use std::collections::HashMap;
use std::path::PathBuf;

#[allow(clippy::identity_op)]
#[test]
fn od_val_sc_mb_srp_reals_duals_models() {
    /*
     * This tests that the state transition matrix computation is correct when multiple celestial gravities and solar radiation pressure
     * are added to the model.
     *
     * Specifically, the same dynamics are used for both the measurement generation and for the estimation.
     * However, only the estimation generation propagates the STM. When STM propagation is enabled, the code will compute
     * the dynamics using a hyperdual representation in 9 dimensions: 1 for the reals, 3 for the position partials,
     * 3 for the velocity partials, 1 for the Cr partials and 1 for the Cd partials.
     *
     * Hence, if the filter state estimation is any different from the truth data, then it means that the equations of
     * motion computed in hyperdual space differ from the ones computes in the reals.
     *
     * Thereby, this serves as a validation of the spacecraft dynamics and SRP duals implementation.
     **/
    let _ = pretty_env_logger::try_init();

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

    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let prop_time = 1 * Unit::Day;

    // Define the tracking configurations
    let mut configs = HashMap::new();
    let cfg = TrkConfig::builder()
        .strands(vec![EpochRanges {
            start: epoch,
            end: epoch + prop_time,
        }])
        .build();

    configs.insert(dss65_madrid.name.clone(), cfg.clone());
    configs.insert(dss34_canberra.name.clone(), cfg.clone());
    configs.insert(dss13_goldstone.name.clone(), cfg);

    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, epoch, eme2k);

    let dry_mass_kg = 100.0; // in kg
    let sc_area = 5.0; // m^2

    // Generate the truth data on one thread.

    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let orbital_dyn = OrbitalDynamics::point_masses(&bodies, cosm.clone());
    let sc_dynamics =
        SpacecraftDynamics::from_model(orbital_dyn, SolarPressure::default(eme2k, cosm.clone()));

    let sc_init_state = Spacecraft::from_srp_defaults(initial_state, dry_mass_kg, sc_area);

    let setup = Propagator::new::<RK4Fixed>(sc_dynamics, opts);
    let mut prop = setup.with(sc_init_state);
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Test the exporting of a spacecraft trajectory
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "output_data",
        "sc_truth_val.parquet",
    ]
    .iter()
    .collect();

    // Adding some events to the exported trajectory
    let event = Event::specific(StateParameter::Declination, 6.0, 3.0, Unit::Minute);

    let cfg = ExportCfg::from_metadata(vec![
        (
            "Dynamics".to_string(),
            "SRP, Moon, Sun, Jupiter".to_string(),
        ),
        // An `Event:` metadata will be appropriately parsed and plotted with the Nyx plotting tools.
        (
            "Event: Comms Start".to_string(),
            format!("{}", epoch + 6.minutes()),
        ),
        (
            "Event: Comms End".to_string(),
            format!("{}", epoch + 9.minutes()),
        ),
    ]);

    traj.to_parquet(path.clone(), Some(vec![&event]), cfg)
        .unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(all_stations, traj, configs, 0).unwrap();
    arc_sim.build_schedule(cosm.clone()).unwrap();

    let arc = arc_sim.generate_measurements(cosm.clone()).unwrap();

    arc.to_parquet_simple(path.with_file_name("sc_msr_arc.parquet"))
        .unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let initial_state_est = initial_state.with_stm();
    let sc_init_est = Spacecraft::from_srp_defaults(initial_state_est, dry_mass_kg, sc_area);
    // Use the same setup as earlier
    let prop_est = setup.with(sc_init_est);
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

    // Define the initial orbit estimate
    let initial_estimate = KfEstimate::from_covar(initial_state_est, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise =
        Matrix2::from_diagonal(&Vector2::new(15e-3_f64.powi(2), 1e-5_f64.powi(2)));

    let ckf = KF::no_snc(initial_estimate, measurement_noise);

    let mut odp = ODProcess::ckf(prop_est, ckf, None, cosm);

    odp.process_arc::<GroundStation>(&arc).unwrap();

    odp.to_parquet(
        path.with_file_name("spacecraft_od_results.parquet"),
        ExportCfg::timestamped(),
    )
    .unwrap();

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
            est.state_deviation().norm() < 1e-12,
            "estimate error should be zero (perfect dynamics) ({:e})",
            est.state_deviation().norm()
        );
    }

    for res in odp.residuals.iter().flatten() {
        assert!(
            res.postfit.norm() < 1e-5,
            "postfit should be zero (perfect dynamics) ({:e})",
            res
        );
    }

    let est = odp.estimates.last().unwrap();
    println!("estimate error {:.2e}", est.state_deviation().norm());
    println!("estimate covariance {:.2e}", est.covar.diagonal().norm());

    assert!(
        est.state_deviation().norm() < 1e-12,
        "estimate error should be zero (perfect dynamics) ({:e})",
        est.state_deviation().norm()
    );

    assert!(
        est.covar.diagonal().norm() < 1e-4,
        "estimate covariance norm should be small (perfect dynamics) ({:e})",
        est.covar.diagonal().norm()
    );

    let delta = est.state() - final_truth.orbit;
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 1e-9, "More than 1 micrometer error");
    assert!(delta.vmag_km_s() < 1e-9, "More than 1 micrometer/s error");
}
