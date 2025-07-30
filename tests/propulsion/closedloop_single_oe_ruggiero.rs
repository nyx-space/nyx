extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::nyx::cosmic::{GuidanceMode, Orbit, Spacecraft};
use self::nyx::dynamics::guidance::{Objective, Ruggiero, StateParameter, Thruster};
use self::nyx::dynamics::{OrbitalDynamics, SpacecraftDynamics};
use self::nyx::propagators::{IntegratorOptions, Propagator};
use self::nyx::time::{Epoch, Unit};

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use nyx_space::propagators::IntegratorMethod;
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn rugg_sma(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 45 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::SMA,
        42_164.0,
        1.0,
    )];

    let guid_law = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let prop_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_sma] {orbit:x}");

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let prop_usage = prop_mass - final_state.mass.prop_mass_kg;
    println!("[rugg_sma] {:x}", final_state.orbit);
    println!("[rugg_sma] prop usage: {prop_usage:.3} kg");

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((prop_usage - 21.0).abs() < 1.0);
}

#[rstest]
fn rugg_sma_regress_threshold(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.1, 0.0, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 175 * Unit::Day;

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::SMA,
        42_164.0,
        1.0,
    )];

    let prop_mass = 67.0;
    let dry_mass = 300.0;
    for (threshold, expected_prop_usage) in &[(0.9, 16.9), (0.0, 21.3)] {
        let guid_law = Ruggiero::from_ηthresholds(objectives, &[*threshold], orbit.into()).unwrap();
        let sc_state =
            Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

        let sc = SpacecraftDynamics::from_guidance_law(OrbitalDynamics::two_body(), guid_law);
        println!("[rugg_sma_regress] {orbit:x}");

        let final_state = Propagator::new(
            sc.clone(),
            IntegratorMethod::RungeKutta4,
            IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
        )
        .with(sc_state, almanac.clone())
        .for_duration(prop_time)
        .unwrap();

        let prop_usage = prop_mass - final_state.mass.prop_mass_kg;
        println!("[rugg_sma_regress] {:x}", final_state.orbit);
        println!("[rugg_sma_regress] prop usage: {prop_usage:.3} kg");

        assert!(
            sc.guidance_achieved(&final_state).unwrap(),
            "objective not achieved"
        );
        assert!((prop_usage - *expected_prop_usage).abs() < 0.5);
    }
}

#[rstest]
fn rugg_sma_decr(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(42164.0, 0.0, 0.0, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 45 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::SMA,
        24_396.0,
        1.0,
    )];

    let guid_law = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let prop_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_sma_decr] {orbit:x}");

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let prop_usage = prop_mass - final_state.mass.prop_mass_kg;
    println!("[rugg_sma_decr] {:x}", final_state.orbit);
    println!("[rugg_sma_decr] prop usage: {prop_usage:.3} kg");

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((prop_usage - 21.0).abs() < 1.0);
}

#[rstest]
fn rugg_inc(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.mean_equatorial_radius_km().unwrap() + 350.0;

    let orbit = Orbit::keplerian(sma, 0.001, 46.0, 1.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 55 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::Inclination,
        51.6,
        5e-3,
    )];

    let guid_law = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let prop_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_inc] {orbit:x}");

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let prop_usage = prop_mass - final_state.mass.prop_mass_kg;
    println!("[rugg_inc] {:x}", final_state.orbit);
    println!("[rugg_inc] prop usage: {prop_usage:.3} kg");

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((prop_usage - 25.0).abs() < 1.0);
}

#[rstest]
fn rugg_inc_threshold(almanac: Arc<Almanac>) {
    // Same inclination test as above, but with an efficiency threshold. Data comes from Figure 7 of IEPC-2011-102.

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::try_keplerian_altitude(350.0, 0.001, 46.0, 1.0, 1.0, 1.0, start_time, eme2k)
        .unwrap();

    let prop_time = 130 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::Inclination,
        51.6,
        5e-3,
    )];

    let guid_law = Ruggiero::from_ηthresholds(objectives, &[0.9], orbit.into()).unwrap();

    let prop_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_inc_threshold] {orbit:x}");

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let prop_usage = prop_mass - final_state.mass.prop_mass_kg;
    println!("[rugg_inc_threshold] {:x}", final_state.orbit);
    println!("[rugg_inc_threshold] prop usage: {prop_usage:.3} kg");

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((prop_usage - 17.0).abs() < 1.0);
}

#[rstest]
fn rugg_inc_decr(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.mean_equatorial_radius_km().unwrap() + 350.0;

    let orbit = Orbit::keplerian(sma, 0.001, 51.6, 1.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 55 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::Inclination,
        46.0,
        5e-3,
    )];

    let guid_law = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let prop_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_inc_decr] {orbit:x}");

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let prop_usage = prop_mass - final_state.mass.prop_mass_kg;
    println!("[rugg_inc_decr] {:x}", final_state.orbit);
    println!("[rugg_inc_decr] prop usage: {prop_usage:.3} kg");

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((prop_usage - 25.0).abs() < 1.0);
}

#[rstest]
fn rugg_ecc(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.mean_equatorial_radius_km().unwrap() + 9000.0;

    let orbit = Orbit::keplerian(sma, 0.01, 98.7, 0.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 30 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::Eccentricity,
        0.15,
        5e-5,
    )];

    let guid_law = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let prop_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_ecc] {orbit:x}");

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let prop_usage = prop_mass - final_state.mass.prop_mass_kg;
    println!("[rugg_ecc] {:x}", final_state.orbit);
    println!("[rugg_ecc] prop usage: {prop_usage:.3} kg");

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((prop_usage - 10.37).abs() < 1.0);
}

#[rstest]
fn rugg_ecc_regress_threshold(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.mean_equatorial_radius_km().unwrap() + 9000.0;

    let orbit = Orbit::keplerian(sma, 0.01, 98.7, 0.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 150 * Unit::Day;

    let prop_mass = 67.0;
    let dry_mass = 300.0;

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::Eccentricity,
        0.15,
        5e-5,
    )];

    for (threshold, expected_prop_usage) in &[(0.9, 8.2), (0.0, 10.37)] {
        let guid_law = Ruggiero::from_ηthresholds(objectives, &[*threshold], orbit.into()).unwrap();

        let sc_state =
            Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

        let sc = SpacecraftDynamics::from_guidance_law(OrbitalDynamics::two_body(), guid_law);
        println!("[rugg_ecc_regress] {orbit:x}");

        let final_state = Propagator::new(
            sc.clone(),
            IntegratorMethod::RungeKutta4,
            IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
        )
        .with(sc_state, almanac.clone())
        .for_duration(prop_time)
        .unwrap();

        let prop_usage = prop_mass - final_state.mass.prop_mass_kg;
        println!("[rugg_ecc_regress] {:x}", final_state.orbit);
        println!("[rugg_ecc_regress] prop usage: {prop_usage:.3} kg");

        assert!(
            sc.guidance_achieved(&final_state).unwrap(),
            "objective not achieved"
        );
        assert!((prop_usage - *expected_prop_usage).abs() < 1.0);
    }
}

#[rstest]
fn rugg_ecc_decr(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.mean_equatorial_radius_km().unwrap() + 9000.0;

    let orbit = Orbit::keplerian(sma, 0.15, 98.7, 0.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 30 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::Eccentricity,
        0.01,
        5e-5,
    )];

    let guid_law = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let prop_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_ecc_decr] {orbit:x}");

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let prop_usage = prop_mass - final_state.mass.prop_mass_kg;
    println!("[rugg_ecc_decr] {:x}", final_state.orbit);
    println!("[rugg_ecc_decr] prop usage: {prop_usage:.3} kg");

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((prop_usage - 10.37).abs() < 1.0);
}

#[rstest]
fn rugg_aop(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.mean_equatorial_radius_km().unwrap() + 900.0;

    // Note that AOP computation requires the orbit to not be equatorial or circular, hence the non-zero ECC and INC.
    let orbit = Orbit::keplerian(sma, 5e-5, 5e-3, 0.0, 178.0, 0.0, start_time, eme2k);

    // This is a very quick change because we aren't using the Ruggiero formulation for AOP change and benefit both in-plane and out of plane control.
    let prop_time = 44 * Unit::Minute + 10 * Unit::Second;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::AoP,
        183.0,
        5e-3,
    )];

    let guid_law = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let prop_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_aop] {orbit:x}");

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let prop_usage = prop_mass - final_state.mass.prop_mass_kg;
    println!("[rugg_aop] {:x}", final_state.orbit);
    println!("[rugg_aop] prop usage: {prop_usage:.3} kg");

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((prop_usage - 0.014).abs() < 1e-2);
}

#[rstest]
fn rugg_aop_decr(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.mean_equatorial_radius_km().unwrap() + 900.0;

    // Note that AOP computation requires the orbit to not be equatorial or circular, hence the non-zero ECC and INC.
    let orbit = Orbit::keplerian(sma, 5e-5, 5e-3, 0.0, 183.0, 0.0, start_time, eme2k);

    let prop_time = 44 * Unit::Minute + 10 * Unit::Second;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::AoP,
        178.0,
        5e-3,
    )];

    let guid_law = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let prop_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_aop_decr] {orbit:x}");

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let prop_usage = prop_mass - final_state.mass.prop_mass_kg;
    println!("[rugg_aop_decr] {:x}", final_state.orbit);
    println!("[rugg_aop_decr] prop usage: {prop_usage:.3} kg");

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((prop_usage - 0.014).abs() < 1e-2);
}

#[rstest]
fn rugg_raan(almanac: Arc<Almanac>) {
    use self::nyx::md::{Event, StateParameter};

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2017, 1, 1);

    let sma = eme2k.mean_equatorial_radius_km().unwrap() + 798.0;

    let orbit = Orbit::keplerian(sma, 0.00125, 98.57, 0.0, 1.0, 0.0, start_time, eme2k);

    let prop_time = 49 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(StateParameter::RAAN, 5.0, 5e-5)];

    let guid_law = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let prop_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_raan] {orbit:x}");

    let setup = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    );
    let mut prop = setup.with(sc_state, almanac.clone());
    let (final_state, traj) = prop.for_duration_with_traj(prop_time).unwrap();
    let prop_usage = prop_mass - final_state.mass.prop_mass_kg;
    println!("[rugg_raan] {:x}", final_state.orbit);
    let event = Event::new(StateParameter::RAAN, 5.0);
    println!(
        "[rugg_raan] {} => {:?}",
        event,
        traj.find(&event, None, almanac)
    );
    println!("[rugg_raan] prop usage: {prop_usage:.3} kg");

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((prop_usage - 22.189).abs() < 1.0);
}

#[rstest]
fn rugg_raan_regress_threshold(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.mean_equatorial_radius_km().unwrap() + 798.0;

    let orbit = Orbit::keplerian(sma, 0.00125, 98.57, 0.0, 1.0, 0.0, start_time, eme2k);

    let prop_time = 130 * Unit::Day;

    let prop_mass = 67.0;
    let dry_mass = 300.0;

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(StateParameter::RAAN, 5.0, 5e-5)];

    for (threshold, expected_prop_usage) in &[(0.9, 14.787), (0.0, 22.189)] {
        let guid_law = Ruggiero::from_ηthresholds(objectives, &[*threshold], orbit.into()).unwrap();

        let sc_state =
            Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

        let sc = SpacecraftDynamics::from_guidance_law(OrbitalDynamics::two_body(), guid_law);
        println!("[rugg_raan_regress] {orbit:x}");

        let final_state = Propagator::new(
            sc.clone(),
            IntegratorMethod::RungeKutta4,
            IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
        )
        .with(sc_state, almanac.clone())
        .for_duration(prop_time)
        .unwrap();

        let prop_usage = prop_mass - final_state.mass.prop_mass_kg;
        println!("[rugg_raan_regress] {:x}", final_state.orbit);
        println!("[rugg_raan_regress] prop usage: {prop_usage:.3} kg");

        assert!(
            sc.guidance_achieved(&final_state).unwrap(),
            "objective not achieved"
        );
        assert!((prop_usage - *expected_prop_usage).abs() < 1e-1);
    }
}
