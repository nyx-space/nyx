extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MOON, SUN};
use anise::constants::frames::IAU_EARTH_FRAME;
use nyx::cosmic::{Orbit, Spacecraft};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::spacecraft::{SolarPressure, SpacecraftDynamics};
use nyx::linalg::{SMatrix, SVector};
use nyx::md::trajectory::ExportCfg;
use nyx::md::{Event, StateParameter};
use nyx::od::prelude::*;
use nyx::propagators::{IntegratorOptions, Propagator};
use nyx::time::{Epoch, TimeUnits, Unit};
use nyx_space::cosmic::{Mass, SRPData};
use nyx_space::dynamics::guidance::LocalFrame;
use nyx_space::propagators::{ErrorControl, IntegratorMethod};
use std::collections::BTreeMap;
use std::path::PathBuf;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[fixture]
fn sim_devices(almanac: Arc<Almanac>) -> BTreeMap<String, GroundStation> {
    let iau_earth = almanac.frame_info(IAU_EARTH_FRAME).unwrap();
    let elevation_mask = 0.0;
    let dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        StochasticNoise::ZERO,
        StochasticNoise::ZERO,
        iau_earth,
    );
    let dss34_canberra = GroundStation::dss34_canberra(
        elevation_mask,
        StochasticNoise::ZERO,
        StochasticNoise::ZERO,
        iau_earth,
    );
    let dss13_goldstone = GroundStation::dss13_goldstone(
        elevation_mask,
        StochasticNoise::ZERO,
        StochasticNoise::ZERO,
        iau_earth,
    );

    let mut devices = BTreeMap::new();
    devices.insert("Madrid".to_string(), dss65_madrid);
    devices.insert("Canberra".to_string(), dss34_canberra);
    devices.insert("Goldstone".to_string(), dss13_goldstone);
    devices
}

/// Devices for processing the measurement, noise may not be zero.
#[fixture]
fn proc_devices(almanac: Arc<Almanac>) -> BTreeMap<String, GroundStation> {
    let iau_earth = almanac.frame_info(IAU_EARTH_FRAME).unwrap();
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

    let mut devices = BTreeMap::new();
    devices.insert("Madrid".to_string(), dss65_madrid);
    devices.insert("Canberra".to_string(), dss34_canberra);
    devices.insert("Goldstone".to_string(), dss13_goldstone);
    devices
}

#[allow(clippy::identity_op)]
#[rstest]
fn od_val_sc_mb_srp_reals_duals_models(
    almanac: Arc<Almanac>,
    sim_devices: BTreeMap<String, GroundStation>,
    proc_devices: BTreeMap<String, GroundStation>,
) {
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

    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let prop_time = 1 * Unit::Day;

    // Define the tracking configurations
    let mut configs = BTreeMap::new();
    let cfg = TrkConfig::builder()
        .strands(vec![Strand {
            start: epoch,
            end: epoch + prop_time,
        }])
        .build();

    for name in sim_devices.keys() {
        configs.insert(name.clone(), cfg.clone());
    }

    // Define the propagator information.
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    // Define state information.
    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();
    let initial_state = Orbit::keplerian(22000.0, 0.01, 30.0, 80.0, 40.0, 0.0, epoch, eme2k);

    let dry_mass_kg = 100.0; // in kg
    let sc_area = 5.0; // m^2

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let orbital_dyn = OrbitalDynamics::point_masses(bodies);
    let sc_dynamics = SpacecraftDynamics::from_model(
        orbital_dyn,
        SolarPressure::default(eme2k, almanac.clone()).unwrap(),
    );

    let sc_init_state = Spacecraft::from_srp_defaults(initial_state, dry_mass_kg, sc_area);

    let setup = Propagator::new(sc_dynamics, IntegratorMethod::RungeKutta4, opts);
    let mut prop = setup.with(sc_init_state, almanac.clone());
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Test the exporting of a spacecraft trajectory
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "04_output",
        "sc_truth_val.parquet",
    ]
    .iter()
    .collect();

    // Adding some events to the exported trajectory
    let event = Event::specific(StateParameter::Declination, 6.0, 3.0, Unit::Minute);

    let cfg = ExportCfg::from_metadata(vec![
        (
            "Dynamics".to_string(),
            "SRP, Moon, Sun, JUPITER_BARYCENTER".to_string(),
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

    traj.to_parquet(path.clone(), Some(vec![&event]), cfg, almanac.clone())
        .unwrap();

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(sim_devices, traj, configs.clone(), 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    arc.to_parquet_simple(path.with_file_name("sc_msr_arc.parquet"))
        .unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let initial_state_est = initial_state;
    let sc_init_est =
        Spacecraft::from_srp_defaults(initial_state_est, dry_mass_kg, sc_area).with_stm();
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

    // Define the initial orbit estimate
    let initial_estimate = KfEstimate::from_covar(sc_init_est, init_covar);

    let odp = SpacecraftKalmanOD::new(
        setup,
        KalmanVariant::DeviationTracking,
        None,
        proc_devices,
        almanac,
    );

    let od_sol = odp.process_arc(initial_estimate, &arc).unwrap();

    od_sol
        .to_parquet(
            path.with_file_name("spacecraft_od_results.parquet"),
            ExportCfg::timestamped(),
        )
        .unwrap();

    for (no, est) in od_sol.estimates.iter().enumerate() {
        if no == 0 {
            // Skip the first estimate which is the initial estimate provided by user
            continue;
        }
        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0,
                "covar diagonal element negative @ [{i}, {i}]"
            );
        }
        assert!(
            est.state_deviation().norm() < 1e-12,
            "estimate error should be zero (perfect dynamics) ({:e})",
            est.state_deviation().norm()
        );
    }

    for res in od_sol.residuals.iter().flatten() {
        assert!(
            res.postfit.norm() < 1e-5,
            "postfit should be zero (perfect dynamics) ({res:e})"
        );
    }

    let est = od_sol.estimates.last().unwrap();
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

    let delta = (est.state().orbit - final_truth.orbit).unwrap();
    println!(
        "RMAG error = {:.2e} m\tVMAG error = {:.3e} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 1e-9, "More than 1 micrometer error");
    assert!(delta.vmag_km_s() < 1e-9, "More than 1 micrometer/s error");
}

#[allow(clippy::identity_op)]
#[rstest]
fn od_val_sc_srp_estimation_cov_test(
    almanac: Arc<Almanac>,
    sim_devices: BTreeMap<String, GroundStation>,
    proc_devices: BTreeMap<String, GroundStation>,
) {
    /*
     * This tests that we can correctly estimate the solar radiation pressure.
     *
     * The truth data is generated with a specific value of the Cr and then the filter is initialized with another value.
     * We expect that the estimation eventually converges onto the truth value.
     **/
    let _ = pretty_env_logger::try_init();

    let epoch = Epoch::from_gregorian_utc_at_noon(2024, 2, 29);

    // Define state information.
    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();
    // Using a GTO because Cr estimation will be more obvious.
    let initial_orbit = Orbit::keplerian(24505.9, 0.725, 7.05, 0.0, 0.0, 0.0, epoch, eme2k);

    let prop_time = initial_orbit.period().unwrap() * 5;

    let dry_mass_kg = 100.0; // in kg

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let orbital_dyn = OrbitalDynamics::point_masses(bodies);
    let sc_dynamics = SpacecraftDynamics::from_model(
        orbital_dyn,
        SolarPressure::default(eme2k, almanac.clone()).unwrap(),
    );

    let truth_cr = 1.123;

    let sc_truth = Spacecraft::builder()
        .orbit(initial_orbit)
        .srp(SRPData {
            coeff_reflectivity: truth_cr,
            area_m2: 100.0,
        })
        .mass(Mass::from_dry_mass(dry_mass_kg))
        .build();

    println!("{sc_truth}");

    let setup = Propagator::dp78(
        sc_dynamics,
        IntegratorOptions::builder()
            .max_step(30.seconds())
            .error_ctrl(ErrorControl::RSSCartesianStep)
            .build(),
    );
    let mut prop = setup.with(sc_truth, almanac.clone());
    let (_, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Test the exporting of a spacecraft trajectory
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "04_output",
        "sc_srp_truth.parquet",
    ]
    .iter()
    .collect();

    traj.to_parquet_simple(path.clone(), almanac.clone())
        .unwrap();

    // Define the tracking configurations
    let mut configs = BTreeMap::new();
    let cfg = TrkConfig::builder()
        .strands(vec![Strand {
            start: epoch,
            end: epoch + prop_time,
        }])
        .build();

    for name in sim_devices.keys() {
        configs.insert(name.clone(), cfg.clone());
    }

    let all_stations = sim_devices;

    // Simulate tracking data
    let mut arc_sim =
        TrackingArcSim::with_seed(all_stations, traj.clone(), configs.clone(), 120).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    arc.to_parquet_simple(path.with_file_name("sc_srp_msr_arc.parquet"))
        .unwrap();

    let sc_init_est = Spacecraft::builder()
        .orbit(initial_orbit)
        .srp(SRPData {
            coeff_reflectivity: 1.5,
            area_m2: 100.0,
        })
        .mass(Mass::from_dry_mass(dry_mass_kg))
        .build()
        .with_stm();

    let sc = SpacecraftUncertainty::builder()
        .nominal(sc_init_est)
        .frame(LocalFrame::RIC)
        .x_km(0.001)
        .y_km(0.001)
        .z_km(0.001)
        .vx_km_s(0.1e-6)
        .vy_km_s(0.1e-6)
        .vz_km_s(0.1e-6)
        .coeff_reflectivity(0.2)
        .build();

    // Define the initial orbit estimate
    let initial_estimate = sc.to_estimate().unwrap();

    let odp = SpacecraftKalmanOD::new(
        setup,
        KalmanVariant::ReferenceUpdate,
        None,
        proc_devices.clone(),
        almanac.clone(),
    );

    let od_sol = odp.process_arc(initial_estimate, &arc).unwrap();

    let est = od_sol.estimates.last().unwrap();

    let truth = traj.at(est.epoch()).unwrap();

    println!(
        "FINAL:\n\t{est}\nGOT:{:x}\nGOT: Cr = {:.6} +/-{:.6}\nEXP: Cr = {truth_cr:.6}\nEXP:{:x}",
        est.state().orbit,
        est.state().srp.coeff_reflectivity,
        est.covar()[(6, 6)].sqrt(),
        truth.orbit
    );

    let delta = (est.state().orbit - truth.orbit).unwrap();
    println!(
        "RMAG error = {:.6} m\tVMAG error = {:.6} mm/s",
        delta.rmag_km() * 1e3,
        delta.vmag_km_s() * 1e6
    );

    assert!(delta.rmag_km() < 2e-3, "More than 2 meter error");
    assert!(delta.vmag_km_s() < 2e-6, "More than 2 millimeter/s error");
    assert!(
        (est.state().srp.coeff_reflectivity - truth_cr).abs() < (1.5 - truth_cr),
        "Cr estimation did not improve"
    );

    for (no, est) in od_sol.estimates.iter().enumerate() {
        if no == 0 {
            // Skip the first estimate which is the initial estimate provided by user
            continue;
        }
        for i in 0..6 {
            assert!(
                est.covar[(i, i)] >= 0.0,
                "covar diagonal element negative @ [{i}, {i}]"
            );
        }
    }

    // The following are merely checks that the code works as expected, not checks on the OD.

    let pre_smooth_postfit_rms = od_sol.rms_postfit_residuals();

    od_sol
        .to_parquet("./data/04_output/od_srp_val.parquet", ExportCfg::default())
        .unwrap();

    // Reload OD Solution
    let od_sol_reloaded =
        ODSolution::from_parquet("./data/04_output/od_srp_val.parquet", proc_devices).unwrap();
    assert_eq!(od_sol_reloaded.estimates.len(), od_sol.estimates.len());
    assert_eq!(od_sol_reloaded.residuals.len(), od_sol.residuals.len());
    assert_eq!(od_sol_reloaded.gains.len(), od_sol.gains.len());
    assert_eq!(
        od_sol_reloaded.filter_smoother_ratios.len(),
        od_sol.filter_smoother_ratios.len()
    );
    assert!(od_sol_reloaded == od_sol, "womp womp");

    println!(
        "Num residuals accepted: #{}",
        od_sol.accepted_residuals().len()
    );
    println!(
        "Num residuals rejected: #{}",
        od_sol.rejected_residuals().len()
    );
    println!(
        "Percentage within +/-3: {}",
        od_sol.residual_ratio_within_threshold(3.0).unwrap()
    );
    println!("Ratios normal? {}", od_sol.is_normal(None).unwrap());

    // Regression tests, not solution test
    assert_eq!(
        od_sol.accepted_residuals().len(),
        2992,
        "all residuals should be accepted"
    );
    assert_eq!(
        od_sol.rejected_residuals().len(),
        0,
        "all residuals should be accepted"
    );
    assert!(
        (0.80..0.90).contains(&od_sol.residual_ratio_within_threshold(3.0).unwrap()),
        "expecting 80-90% of valid residual ratios"
    );
    assert!(
        !od_sol.is_normal(None).unwrap(),
        "residuals should not be normally distributed"
    );

    let od_smoothed_sol = od_sol.smooth(almanac).unwrap();

    od_smoothed_sol
        .to_parquet("./od_srp_val_smoothed.parquet", ExportCfg::default())
        .unwrap();

    println!(
        "one-pass vs smoothed postfit RMS: {pre_smooth_postfit_rms}\t{}",
        od_smoothed_sol.rms_postfit_residuals()
    );
}
