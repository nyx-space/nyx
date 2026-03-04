// Test for interlink

extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use anise::constants::celestial_objects::{EARTH, SUN};
use anise::constants::frames::IAU_MOON_FRAME;
use anise::math::Vector6;
use indexmap::{IndexMap, IndexSet};
use nyx::cosmic::Orbit;
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::linalg::Const;
use nyx::md::prelude::*;
use nyx::od::interlink::InterlinkTxSpacecraft;
use nyx::od::prelude::*;
use nyx::propagators::Propagator;
use nyx::time::{Epoch, TimeUnits};
use nyx::Spacecraft;

use anise::prelude::Almanac;
use nyx_space::od::groundpnt::ground_dynamics::GroundDynamics;
use nyx_space::od::groundpnt::GroundAsset;
use rstest::*;
use std::collections::BTreeMap;
use std::path::PathBuf;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn ground_pnt_lunar_cov_test(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let moon_iau = almanac.frame_info(IAU_MOON_FRAME).unwrap();

    let epoch = Epoch::from_gregorian_utc_at_midnight(2024, 2, 29);
    /* == Propagate an LLO vehicle == */
    let llo_orbit =
        Orbit::try_keplerian_altitude(110.0, 1e-4, 35.0, 0.0, 0.0, 0.0, epoch, moon_iau).unwrap();

    let tx_llo_sc = Spacecraft::from(llo_orbit);

    println!("Start state (dynamics: Earth, Moon, Sun gravity):\n{llo_orbit}");

    let bodies = vec![EARTH, SUN];
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));

    let setup = Propagator::rk89(
        dynamics,
        IntegratorOptions::builder().max_step(0.5.minutes()).build(),
    );

    let prop_time = 10.5 * llo_orbit.period().unwrap();

    let (llo_final, mut llo_traj) = setup
        .with(tx_llo_sc, almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    llo_traj.name = Some("PNT SC".to_string());

    println!("{llo_traj}");

    /* == Setup a ground asset visible at the halfway mark of this arbitrary trajectory == */
    let (lat_deg, long_deg, _) = llo_traj
        .at(epoch + 0.55 * llo_orbit.period().unwrap())
        .unwrap()
        .orbit
        .latlongalt()
        .unwrap();

    // Set it up as immobile.
    let rover = GroundAsset::from_fixed(lat_deg, long_deg, 0.0, epoch, moon_iau);

    // Build a trajectory of the rover. It is immobile, but we demonstrate here that it can be propagated, which is required to run it through an OD.
    let rover_prop = Propagator::default(GroundDynamics {});
    let (rover_final, rover_traj) = rover_prop
        .clone()
        .with(rover, almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Check that the rover didn't move
    assert!((rover_final.latitude_deg - rover.latitude_deg).abs() < f64::EPSILON);
    assert!((rover_final.longitude_deg - rover.longitude_deg).abs() < f64::EPSILON);
    assert!((rover_final.height_km - rover.height_km).abs() < f64::EPSILON);

    /* == Check that the S/E/Z movement is properly integrated == */
    /* *** S = 1 m/s ; E = 0 m/s ; Z = 0 m/s *** */
    let rover_sb = rover.with_velocity_sez_m_s(1.0, 0.0, 0.0).unwrap();
    let rover_sb_final = Propagator::default(GroundDynamics {})
        .with(rover_sb, almanac.clone())
        .for_duration(Unit::Second * 1000)
        .unwrap();

    println!("SOUNDBOUND: {rover_sb} -> {rover_sb_final}");
    let vel_sez_m_s = rover_sb_final.velocity_sez_m_s().unwrap();

    assert!((vel_sez_m_s[0] - 1.0).abs() < 1e-12);
    assert!((vel_sez_m_s[1]).abs() < f64::EPSILON);
    assert!((vel_sez_m_s[2]).abs() < f64::EPSILON);
    // Propagated for a thousand seconds going one meter per second south bound, moving by 0.002 degrees
    assert!((rover_sb_final.latitude_deg - rover.latitude_deg).abs() < 4e-2);
    assert!((rover_sb_final.longitude_deg - rover.longitude_deg).abs() < f64::EPSILON);
    assert!((rover_sb_final.height_km - rover.height_km).abs() < f64::EPSILON);
    /* *** S = 0 m/s ; E = 1 m/s ; Z = 0 m/s *** */
    let rover_eb = rover.with_velocity_sez_m_s(0.0, 1.0, 0.0).unwrap();
    let rover_eb_final = Propagator::default(GroundDynamics {})
        .with(rover_eb, almanac.clone())
        .for_duration(Unit::Second * 1000)
        .unwrap();

    println!("EASTBOUND: {rover_eb} -> {rover_eb_final}");
    let vel_sez_m_s = rover_eb.velocity_sez_m_s().unwrap();

    assert!((vel_sez_m_s[1] - 1.0).abs() < 1e-12);
    assert!((vel_sez_m_s[0]).abs() < f64::EPSILON);
    assert!((vel_sez_m_s[2]).abs() < f64::EPSILON);

    assert!((rover_eb_final.latitude_deg - rover.latitude_deg).abs() < f64::EPSILON);
    assert!((rover_eb_final.longitude_deg - rover.longitude_deg).abs() < 4e-2);
    assert!((rover_eb_final.height_km - rover.height_km).abs() < f64::EPSILON);

    /* == Setup the interlink == */

    let mut measurement_types = IndexSet::new();
    measurement_types.insert(MeasurementType::Range);
    measurement_types.insert(MeasurementType::Doppler);

    let mut stochastics = IndexMap::new();

    let sa45_csac_allan_dev = 1e-11;

    stochastics.insert(
        MeasurementType::Range,
        StochasticNoise::from_hardware_range_km(
            sa45_csac_allan_dev,
            10.0.seconds(),
            noise::link_specific::ChipRate::StandardT4B,
            noise::link_specific::SN0::Average,
        ),
    );

    stochastics.insert(
        MeasurementType::Doppler,
        StochasticNoise::from_hardware_doppler_km_s(
            sa45_csac_allan_dev,
            10.0.seconds(),
            noise::link_specific::CarrierFreq::SBand,
            noise::link_specific::CN0::Average,
        ),
    );

    let interlink = InterlinkTxSpacecraft {
        traj: llo_traj.clone(),
        measurement_types,
        integration_time: None,
        timestamp_noise_s: None,
        ab_corr: Aberration::LT,
        stochastic_noises: Some(stochastics),
    };

    // Devices are the transmitter, which is our NRHO vehicle.
    let mut devices = BTreeMap::new();
    devices.insert("GroundAsset".to_string(), interlink);

    let mut configs = BTreeMap::new();
    configs.insert(
        "GroundAsset".to_string(),
        TrkConfig::builder()
            .strands(vec![Strand {
                start: epoch,
                end: llo_final.epoch(),
            }])
            .sampling(Unit::Second * 10)
            .build(),
    );

    let mut trk_sim =
        TrackingArcSim::with_seed(devices.clone(), rover_traj.clone(), configs, 0).unwrap();
    println!("{trk_sim}");

    let trk_data = trk_sim.generate_measurements(almanac.clone()).unwrap();
    println!("{trk_data}");

    let out = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("data/04_output/");

    trk_data
        .to_parquet_simple(out.clone().join("rover_pnt_msr.pq"))
        .unwrap();

    // Build some uncertainty for the ground asset location
    let mut rover_dispersed = rover;
    // rover_dispersed.latitude_deg += 0.5;
    // rover_dispersed.longitude_deg += 0.5;
    rover_dispersed.latitude_deg = rover_sb_final.latitude_deg;
    rover_dispersed.longitude_deg = rover_eb_final.longitude_deg;

    let asset_estimate = KfEstimate::from_diag(
        rover_dispersed,
        Vector6::from_iterator([1e-2, 1e-2, 1e-2, 0.0, 0.0, 0.0]),
    );

    let proc_devices = devices.clone();

    let odp = KalmanODProcess::<_, Const<2>, Const<3>, InterlinkTxSpacecraft>::new(
        rover_prop,
        // if ref_update {
        KalmanVariant::ReferenceUpdate,
        // } else {
        // KalmanVariant::DeviationTracking,
        // },
        None,
        // Some(ResidRejectCrit::default()),
        proc_devices,
        almanac,
    );

    let err_m = rover_dispersed.great_circle_distance_km(&rover).unwrap() * 1e3;

    println!(
        "== INIT -- Great circle dist: {err_m:.3} ==\nESTIMATE: {}\tsigmas: [{:.2e} deg, {:.2e} deg, {:.2e} m]\nTRUTH\t: {rover}",
        asset_estimate.nominal_state,
        asset_estimate.covar[(0, 0)].sqrt(),
        asset_estimate.covar[(1, 1)].sqrt(),
        asset_estimate.covar[(2, 2)].sqrt() * 1e3
    );

    let od_sol = odp.process_arc(asset_estimate, &trk_data).unwrap();

    println!("{od_sol}");

    let final_est = od_sol.estimates.last().unwrap();
    let truth = rover_traj.at(final_est.epoch()).unwrap();
    println!("Within 3 sigma: {}", final_est.within_3sigma());
    // assert!(final_est.within_3sigma(), "should be within 3 sigma");

    let final_estimate = final_est.nominal_state + final_est.state_deviation;

    let err_m = final_estimate.great_circle_distance_km(&truth).unwrap() * 1e3;

    println!(
        "== FINAL -- Great circle dist: {err_m:.3} ==\nESTIMATE: {final_estimate}\tsigmas: [{:.2e} deg, {:.2e} deg, {:.2e} m]\nTRUTH\t: {truth}",
        final_est.covar[(0, 0)].sqrt(),
        final_est.covar[(1, 1)].sqrt(),
        final_est.covar[(2, 2)].sqrt() * 1e3
    );
}
