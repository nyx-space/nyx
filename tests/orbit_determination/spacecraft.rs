extern crate csv;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use self::nyx::cosmic::{Bodies, Cosm, Orbit, Spacecraft};
use self::nyx::dynamics::orbital::OrbitalDynamics;
use self::nyx::dynamics::spacecraft::{SolarPressure, SpacecraftDynamics};
use self::nyx::io::formatter::NavSolutionFormatter;
use self::nyx::linalg::{Matrix2, Matrix6, Vector2, Vector6};
use self::nyx::od::ui::*;
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use self::nyx::time::{Epoch, Unit};
use std::sync::mpsc;

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
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::io;

    let cosm = Cosm::de438();

    // Define the ground stations.
    let elevation_mask = 0.0;
    let range_noise = 0.0;
    let range_rate_noise = 0.0;
    let dss65_madrid =
        GroundStation::dss65_madrid(elevation_mask, range_noise, range_rate_noise, cosm.clone());
    let dss34_canberra =
        GroundStation::dss34_canberra(elevation_mask, range_noise, range_rate_noise, cosm.clone());
    let dss13_goldstone =
        GroundStation::dss13_goldstone(elevation_mask, range_noise, range_rate_noise, cosm.clone());
    let all_stations = vec![dss65_madrid, dss34_canberra, dss13_goldstone];

    // Define the propagator information.
    let prop_time = 1 * Unit::Day;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    // Define the storages (channels for the states and a map for the measurements).
    let (truth_tx, truth_rx) = mpsc::channel();
    let mut measurements = Vec::with_capacity(10000);

    // Define state information.
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, dt, eme2k);

    let sc_dry_mass = 100.0; // in kg
    let sc_area = 5.0; // m^2

    // Generate the truth data on one thread.

    let bodies = vec![Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter];
    let orbital_dyn = OrbitalDynamics::point_masses(&bodies, cosm.clone());
    let sc_dynamics =
        SpacecraftDynamics::from_model(orbital_dyn, SolarPressure::default(eme2k, cosm.clone()));

    let sc_init_state = Spacecraft::from_srp_defaults(initial_state, sc_dry_mass, sc_area);

    let setup = Propagator::new::<RK4Fixed>(sc_dynamics, opts);
    let mut prop = setup.with(sc_init_state);
    let final_truth = prop.for_duration_with_channel(prop_time, truth_tx).unwrap();

    // Receive the states on the main thread, and populate the measurement channel.
    while let Ok(rx_sc_state) = truth_rx.try_recv() {
        for station in all_stations.iter() {
            let rx_state = rx_sc_state.orbit;
            let meas = station.measure(&rx_state).unwrap();
            if meas.visible() {
                measurements.push(meas);
                break;
            }
        }
    }

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let initial_state_est = initial_state.with_stm();
    let sc_init_est = Spacecraft::from_srp_defaults(initial_state_est, sc_dry_mass, sc_area);
    // Use the same setup as earlier
    let prop_est = setup.with(sc_init_est);
    let covar_radius = 1.0e-3_f64.powi(2);
    let covar_velocity = 1.0e-6_f64.powi(2);
    let init_covar = Matrix6::from_diagonal(&Vector6::new(
        covar_radius,
        covar_radius,
        covar_radius,
        covar_velocity,
        covar_velocity,
        covar_velocity,
    ));

    // Define the initial orbit estimate
    let initial_estimate = KfEstimate::from_covar(initial_state_est, init_covar);

    // Define the expected measurement noise (we will then expect the residuals to be within those bounds if we have correctly set up the filter)
    let measurement_noise =
        Matrix2::from_diagonal(&Vector2::new(15e-3_f64.powi(2), 1e-5_f64.powi(2)));

    let ckf = KF::no_snc(initial_estimate, measurement_noise);

    let mut odp = ODProcess::ckf(prop_est, ckf, all_stations);

    odp.process_measurements(&measurements).unwrap();

    // Initialize the formatter
    let estimate_fmtr = NavSolutionFormatter::default("sc_ckf.csv".to_owned(), cosm);

    let mut wtr = csv::Writer::from_writer(io::stdout());
    wtr.serialize(&estimate_fmtr.headers)
        .expect("could not write to stdout");
    let mut printed = false;
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

        if !printed {
            // Format the estimate
            wtr.serialize(estimate_fmtr.fmt(est))
                .expect("could not write to stdout");
            printed = true;
        }
    }

    for res in &odp.residuals {
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
        est.covar.diagonal().norm() < 1e-5,
        "estimate covariance norm should be small (perfect dynamics) ({:e})",
        est.covar.diagonal().norm()
    );

    let delta = est.state() - final_truth.orbit;
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag() * 1e3,
        delta.vmag() * 1e6
    );

    assert!(delta.rmag() < 1e-9, "More than 1 micrometer error");
    assert!(delta.vmag() < 1e-9, "More than 1 micrometer/s error");
}
