extern crate nalgebra as na;
extern crate nyx_space as nyx;

use std::sync::Arc;

use self::nyx::cosmic::{GuidanceMode, Orbit, Spacecraft};
use self::nyx::dynamics::guidance::{Objective, Ruggiero, Thruster};
use self::nyx::dynamics::{OrbitalDynamics, SpacecraftDynamics};
use self::nyx::md::{Event, StateParameter};
use self::nyx::propagators::{IntegratorOptions, Propagator};
use self::nyx::time::{Epoch, Unit};

/// NOTE: Herein shows the difference between the QLaw and Ruggiero (and other control laws).
/// The Ruggiero control law takes quite some longer to converge than the QLaw.
use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use nyx_space::propagators::IntegratorMethod;
use rstest::*;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn qlaw_as_ruggiero_case_a(almanac: Arc<Almanac>) {
    // Source: AAS-2004-5089
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(398_600.433);

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(7000.0, 0.01, 0.05, 0.0, 0.0, 1.0, start_time, eme2k);

    let prop_time = 39.91 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 1.0,
        isp_s: 3100.0,
    };

    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 42_000.0, 1.0),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.01, 5e-5),
    ];

    // Events we will search later
    let events = vec![
        Event::within_tolerance(StateParameter::SMA, 42_000.0, 1.0),
        Event::within_tolerance(StateParameter::Eccentricity, 0.01, 5e-5),
    ];

    let ruggiero_ctrl = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let dry_mass = 1.0;
    let fuel_mass = 299.0;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[qlaw_as_ruggiero_case_a] {:x}", orbit);

    let setup = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    );
    let mut prop = setup.with(sc_state, almanac.clone());
    let (final_state, traj) = prop.for_duration_with_traj(prop_time).unwrap();
    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[qlaw_as_ruggiero_case_a] {:x}", final_state.orbit);
    println!("[qlaw_as_ruggiero_case_a] fuel usage: {:.3} kg", fuel_usage);
    // Find all of the events
    for e in &events {
        println!(
            "[qlaw_as_ruggiero_case_a] Found {} events of kind {}",
            traj.find(e, almanac.clone()).unwrap().len(),
            e
        );
    }

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );

    assert!((fuel_usage - 93.449).abs() < 1.0);
}

#[rstest]
fn qlaw_as_ruggiero_case_b(almanac: Arc<Almanac>) {
    // Source: AAS-2004-5089

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(24505.9, 0.725, 7.05, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 160.0 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 0.350,
        isp_s: 2000.0,
    };

    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 42_165.0, 20.0),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.001, 5e-5),
        Objective::within_tolerance(StateParameter::Inclination, 0.05, 5e-3),
    ];

    let ruggiero_ctrl = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let fuel_mass = 1999.9;
    let dry_mass = 0.1;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[qlaw_as_ruggiero_case_b] {:x}", orbit);

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[qlaw_as_ruggiero_case_b] {:x}", final_state.orbit);
    println!("[qlaw_as_ruggiero_case_b] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );

    assert!((fuel_usage - 223.515).abs() < 1.0);
}

#[rstest]
fn qlaw_as_ruggiero_case_c_cov_test(almanac: Arc<Almanac>) {
    // Source: AAS-2004-5089

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(9222.7, 0.2, 0.573, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 3.0 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 9.3,
        isp_s: 3100.0,
    };

    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 30_000.0, 1.0),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.7, 5e-5),
    ];

    let ruggiero_ctrl = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let fuel_mass = 299.9;
    let dry_mass = 0.1;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[qlaw_as_ruggiero_case_c] {:x}", orbit);

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[qlaw_as_ruggiero_case_c] {:x}", final_state.orbit);
    println!("[qlaw_as_ruggiero_case_c] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 41.742).abs() < 1.0);
}

#[rstest]
#[ignore = "https://gitlab.com/chrisrabotin/nyx/issues/103"]
fn qlaw_as_ruggiero_case_d(almanac: Arc<Almanac>) {
    // Broken: https://gitlab.com/chrisrabotin/nyx/issues/103
    // Source: AAS-2004-5089

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(24505.9, 0.725, 0.06, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 113.0 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 26_500.0, 1.0),
        Objective::within_tolerance(StateParameter::Inclination, 116.0, 5e-3),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.7, 5e-5),
        Objective::within_tolerance(StateParameter::RAAN, 360.0 - 90.0, 5e-3),
    ];

    let ruggiero_ctrl = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[qlaw_as_ruggiero_case_d] {:x}", orbit);

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[qlaw_as_ruggiero_case_d] {:x}", final_state.orbit);
    println!("[qlaw_as_ruggiero_case_d] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );

    assert!((fuel_usage - 23.0).abs() < 1.0);
}

#[rstest]
#[ignore = "https://gitlab.com/chrisrabotin/nyx/issues/103"]
fn qlaw_as_ruggiero_case_e(almanac: Arc<Almanac>) {
    // Broken: https://gitlab.com/chrisrabotin/nyx/issues/103
    // Source: AAS-2004-5089

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(24505.9, 0.725, 0.06, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 400.0 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 26_500.0, 1.0),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.7, 5e-5),
        Objective::within_tolerance(StateParameter::Inclination, 116.0, 5e-3),
        Objective::within_tolerance(StateParameter::RAAN, 270.0, 5e-3),
        Objective::within_tolerance(StateParameter::AoP, 180.0, 5e-3),
    ];

    let ruggiero_ctrl = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let fuel_mass = 1999.9;
    let dry_mass = 0.1;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[qlaw_as_ruggiero_case_e] {:x}", orbit);

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[qlaw_as_ruggiero_case_e] {:x}", final_state.orbit);
    println!("[qlaw_as_ruggiero_case_e] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );

    assert!((fuel_usage - 23.0).abs() < 1.0);
}

#[rstest]
fn qlaw_as_ruggiero_case_f(almanac: Arc<Almanac>) {
    // Source: AAS-2004-5089
    /*
        NOTE: Due to how lifetime of variables work in Rust, we need to define all of the
        components of a spacecraft before defining the spacecraft itself.
    */

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(15378.0, 0.01, 98.7, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 30.0 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    let objectives = &[Objective::new(StateParameter::Eccentricity, 0.15)];

    let ruggiero_ctrl = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[qlaw_as_ruggiero_case_f] {:x}", orbit);

    let setup = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    );
    let (final_state, traj) = setup
        .with(sc_state, almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Save as parquet
    traj.to_parquet_simple("output_data/rugg_case_f.parquet", almanac)
        .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[qlaw_as_ruggiero_case_f] {:x}", final_state.orbit);
    println!("[qlaw_as_ruggiero_case_f] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );

    assert!((fuel_usage - 10.378).abs() < 1.0);
}

#[rstest]
fn ruggiero_iepc_2011_102(almanac: Arc<Almanac>) {
    // Source: IEPC 2011 102
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.7283, 7.0, 1.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 105.0 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 42_164.0, 20.0),
        Objective::within_tolerance(StateParameter::Inclination, 0.001, 5e-3),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.011, 5e-5),
    ];

    let ruggiero_ctrl = Ruggiero::simple(objectives, orbit.into()).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[ruggiero_iepc_2011_102] {:x}", orbit);

    let final_state = Propagator::new(
        sc.clone(),
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc_state, almanac)
    .for_duration(prop_time)
    .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[ruggiero_iepc_2011_102] {:x}", final_state.orbit);
    println!("[ruggiero_iepc_2011_102] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );

    // WARNING: Paper claims this can be done with only 49kg of fuel.
    assert!((fuel_usage - 49.0).abs() < 1.0);
}
