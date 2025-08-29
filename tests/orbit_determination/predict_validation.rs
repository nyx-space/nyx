// Test for interlink

extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use anise::constants::celestial_objects::{EARTH, SUN};
use anise::constants::frames::{IAU_EARTH_FRAME, IAU_MOON_FRAME, MOON_J2000};
use indexmap::IndexSet;
use nyx::cosmic::Orbit;
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::md::prelude::*;
use nyx::od::prelude::*;
use nyx::propagators::Propagator;
use nyx::time::{Epoch, TimeUnits};
use nyx::Spacecraft;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;
use std::collections::BTreeMap;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

/// Ensures that a propagation in a non-J2000 frame is strictly identical to the pure predictor propagation
#[rstest]
fn val_pure_predictor(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();
    let moon_iau = almanac.frame_from_uid(IAU_MOON_FRAME).unwrap();

    let prop_time = 1.hours();

    let epoch = Epoch::from_gregorian_tai_at_midnight(2025, 1, 1);
    let nrho = Orbit::cartesian(
        166_473.631_302_239_7,
        -274_715.487_253_382_7,
        -211_233.210_176_686_7,
        0.933_451_604_520_018_4,
        0.436_775_046_841_900_9,
        -0.082_211_021_250_348_95,
        epoch,
        eme2k,
    );

    let state_luna = almanac.transform_to(nrho, MOON_J2000, None).unwrap();
    println!("Start state (dynamics: Earth, Moon, Sun gravity):\n{state_luna}");

    let bodies = vec![EARTH, SUN];
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));

    let setup = Propagator::rk89(
        dynamics,
        IntegratorOptions::builder().max_step(1.minutes()).build(),
    );

    let llo_orbit =
        Orbit::try_keplerian_altitude(110.0, 1e-4, 90.0, 0.0, 0.0, 0.0, epoch, moon_iau).unwrap();

    let llo_sc = Spacecraft::builder().orbit(llo_orbit).build();

    let mut this_prop = setup.with(llo_sc, almanac.clone());

    let (_, llo_traj) = this_prop.for_duration_with_traj(prop_time).unwrap();

    /* == Setup the OD process == */

    let mut measurement_types = IndexSet::new();
    measurement_types.insert(MeasurementType::Range);
    measurement_types.insert(MeasurementType::Doppler);

    // Devices are the transmitter, which is our NRHO vehicle.
    let mut devices = BTreeMap::new();
    devices.insert(
        "Canberra".to_string(),
        GroundStation::dss34_canberra(
            0.0,
            StochasticNoise::default_range_km(),
            StochasticNoise::default_doppler_km_s(),
            iau_earth,
        ),
    );

    let mut configs = BTreeMap::new();
    configs.insert(
        "NRHO Tx SC".to_string(),
        TrkConfig::builder()
            .strands(vec![Strand {
                start: epoch,
                end: epoch + prop_time,
            }])
            .build(),
    );

    // Run a truth OD where we estimate the LLO position
    let llo_uncertainty = SpacecraftUncertainty::builder()
        .nominal(llo_sc)
        .x_km(1.0)
        .y_km(1.0)
        .z_km(1.0)
        .vx_km_s(1e-3)
        .vy_km_s(1e-3)
        .vz_km_s(1e-3)
        .build();

    // DISABLE DISPERSIONS, cf. the bug in the estimate randomized docs.
    // Define the initial estimate, randomized, seed for reproducibility
    // let mut initial_estimate = llo_uncertainty.to_estimate_randomized(Some(0)).unwrap();
    // Inflate the covariance
    // initial_estimate.covar *= 3.0;

    let initial_estimate = llo_uncertainty.to_estimate().unwrap();

    let odp = SpacecraftKalmanOD::new(
        setup,
        KalmanVariant::DeviationTracking,
        None,
        devices,
        almanac,
    );

    let pred_truth = llo_traj.last();

    let pp = odp
        .predict_until(initial_estimate, pred_truth.epoch())
        .unwrap();

    let final_est = pp.estimates.last().unwrap().orbital_state();
    let err = final_est.ric_difference(&pred_truth.orbit).unwrap();
    assert!(err.rmag_km() < f64::EPSILON);
    assert!(err.vmag_km_s() < f64::EPSILON);
}
