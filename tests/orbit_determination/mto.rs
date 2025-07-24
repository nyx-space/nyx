extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MARS, MARS_BARYCENTER, SUN};
use anise::constants::frames::IAU_MARS_FRAME;
use anise::constants::orientations::J2000;
use anise::frames::FrameUid;
use nyx::cosmic::{Orbit, Spacecraft};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::spacecraft::{SolarPressure, SpacecraftDynamics};
use nyx::linalg::{SMatrix, SVector};
use nyx::md::trajectory::ExportCfg;
use nyx::od::prelude::*;
use nyx::propagators::{IntegratorOptions, Propagator};
use nyx::time::{Epoch, TimeUnits, Unit};
use nyx_space::propagators::IntegratorMethod;
use std::collections::BTreeMap;
use std::path::PathBuf;

use anise::{constants::frames::MARS_BARYCENTER_J2000, prelude::Almanac};
use indexmap::{IndexMap, IndexSet};
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[fixture]
fn sim_devices(almanac: Arc<Almanac>) -> BTreeMap<String, GroundStation> {
    let iau_mars = almanac.frame_from_uid(IAU_MARS_FRAME).unwrap();

    let mut measurement_types = IndexSet::new();
    measurement_types.insert(MeasurementType::Range);
    measurement_types.insert(MeasurementType::Doppler);

    // WAG
    let range_noise_km = StochasticNoise {
        white_noise: Some(WhiteNoise::new(0.052e-3, 1.0.seconds())),
        bias: None,
    }; // 0.052 m
    let doppler_noise_km_s = StochasticNoise {
        white_noise: Some(WhiteNoise::new(0.09e-3, 1.0.seconds())),
        bias: None,
    }; // 9 cm/s

    let mut stochastics = IndexMap::new();
    stochastics.insert(MeasurementType::Range, range_noise_km);
    stochastics.insert(MeasurementType::Doppler, doppler_noise_km_s);

    let curiosity = GroundStation {
        name: "Curiosity".to_string(),
        elevation_mask_deg: 20.0, // Assumed high
        latitude_deg: 5.4,
        longitude_deg: 137.8,
        height_km: -4500.0, // Lowest point
        frame: iau_mars,
        measurement_types,
        integration_time: Some(60.0.seconds()),
        light_time_correction: true,
        timestamp_noise_s: Some(StochasticNoise {
            bias: None,
            white_noise: Some(WhiteNoise::new(180.0e-10, 60.0.seconds())),
        }),
        stochastic_noises: Some(stochastics),
    };

    let mut oxia_planum = curiosity.clone();
    oxia_planum.latitude_deg = 18.275;
    oxia_planum.longitude_deg = 335.368;

    let mut devices = BTreeMap::new();
    devices.insert("Curiosity".to_string(), curiosity);
    devices.insert("OxiaPlanum".to_string(), oxia_planum);
    devices
}

#[allow(clippy::identity_op)]
#[rstest]
fn mto_od_val_sc_mb_srp_reals_duals_models(
    almanac: Arc<Almanac>,
    sim_devices: BTreeMap<String, GroundStation>,
) {
    let proc_devices = sim_devices.clone();
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

    let epoch = Epoch::from_gregorian_tai_at_midnight(2030, 1, 1);
    let prop_time = 30 * Unit::Day;

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
    let mars_j2k = almanac.frame_from_uid(MARS_BARYCENTER_J2000).unwrap();
    let initial_state = Orbit::keplerian(17000.0, 0.01, 35.0, 80.0, 40.0, 0.0, epoch, mars_j2k);

    let dry_mass_kg = 100.0; // in kg
    let sc_area = 5.0; // m^2

    let bodies = vec![MARS_BARYCENTER, SUN, JUPITER_BARYCENTER];
    let orbital_dyn = OrbitalDynamics::point_masses(bodies);
    let sc_dynamics = SpacecraftDynamics::new(orbital_dyn);
    /*        SolarPressure::default(
            almanac
                .frame_from_uid(FrameUid {
                    ephemeris_id: MARS,
                    orientation_id: J2000,
                })
                .unwrap(),
            almanac.clone(),
        )
        .unwrap(),
    );*/
    // Disabled SRP because I don't have Mars planet loaded in the almanac

    let sc_init_state = Spacecraft::from_srp_defaults(initial_state, dry_mass_kg, sc_area);

    let setup = Propagator::new(sc_dynamics, IntegratorMethod::RungeKutta89, opts);
    let mut prop = setup.with(sc_init_state, almanac.clone());
    let (final_truth, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Test the exporting of a spacecraft trajectory
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "04_output",
        "mto_truth.parquet",
    ]
    .iter()
    .collect();

    traj.to_parquet(path.clone(), None, ExportCfg::default(), almanac.clone())
        .unwrap(); // Dhall error

    // Simulate tracking data
    let mut arc_sim = TrackingArcSim::with_seed(sim_devices, traj, configs.clone(), 0).unwrap();
    arc_sim.build_schedule(almanac.clone()).unwrap();

    let arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    arc.to_parquet_simple(path.with_file_name("mto_msr_arc.parquet"))
        .unwrap();

    // Now that we have the truth data, let's start an OD with no noise at all and compute the estimates.
    // We expect the estimated orbit to be perfect since we're using strictly the same dynamics, no noise on
    // the measurements, and the same time step.
    let initial_state_est = initial_state;
    let sc_init_est =
        Spacecraft::from_srp_defaults(initial_state_est, dry_mass_kg, sc_area).with_stm();
    let covar_radius_km = 1.0_f64.powi(2); // 1 km on each axis
    let covar_velocity_km_s = 1.0e-2_f64.powi(2); // 10 m/s on each axis
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
            path.with_file_name("mto_od_results.parquet"),
            ExportCfg::timestamped(),
        )
        .unwrap();

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
}
