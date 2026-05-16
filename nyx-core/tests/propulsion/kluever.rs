use anise::constants::frames::EARTH_J2000;
use anise::prelude::Almanac;
use nyx_space::cosmic::{GuidanceMode, Mass, Orbit, Spacecraft};
use nyx_space::dynamics::guidance::{Kluever, Objective, Thruster};
use nyx_space::dynamics::{OrbitalDynamics, SpacecraftDynamics};
use nyx_space::md::prelude::{OrbitalElement, StateParameter};
use nyx_space::propagators::{IntegratorOptions, Propagator};
use nyx_space::time::{Epoch, Unit};
use nyx_space::State;
use std::sync::Arc;

use rstest::*;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}
#[rstest]
fn test_qlaw_as_kluever_case_a(almanac: Arc<Almanac>) {
    let eme2k = EARTH_J2000.with_mu_km3_s2(398_600.433);

    // Same conditions as Q Law case A
    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(7000.0, 0.01, 0.05, 0.0, 0.0, 1.0, start_time, eme2k);

    let sc = Spacecraft::builder()
        .orbit(orbit)
        .mass(Mass::from_dry_and_prop_masses(1.0, 299.0))
        .thruster(Thruster {
            isp_s: 3000.0,
            thrust_N: 0.5,
        })
        .mode(GuidanceMode::Thrust)
        .build();

    // Define the objectives: Raise SMA and reduce inclination
    let objectives = &[
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::SemiMajorAxis),
            42_000.0,
            1.0,
        ),
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::Eccentricity),
            0.01,
            5e-5,
        ),
    ];
    let weights = &[1.0, 0.1];

    let kluever = Kluever::new(objectives, weights);

    println!("{kluever}");

    let orbital_dyn = OrbitalDynamics::two_body();
    let sc_dynamics = SpacecraftDynamics::from_guidance_law(orbital_dyn, kluever);

    // Propagate for a short time (e.g., 1 hour) and check that SMA increases
    let prop_time = 1.0 * Unit::Hour;
    let final_state = Propagator::default(sc_dynamics.clone())
        .with(sc, almanac.clone())
        .for_duration(prop_time)
        .unwrap();

    println!("Initial SMA: {} km", orbit.sma_km().unwrap());
    println!("Final SMA: {} km", final_state.orbit.sma_km().unwrap());
    println!("Initial Inc: {} deg", orbit.inc_deg().unwrap());
    println!("Final Inc: {} deg", final_state.orbit.inc_deg().unwrap());

    assert!(
        final_state.orbit.sma_km().unwrap() > orbit.sma_km().unwrap(),
        "SMA should have increased"
    );

    // Run full raise to compare with Ruggiero
    let prop = Propagator::rk89(
        sc_dynamics,
        IntegratorOptions::builder()
            .tolerance(1e-9)
            .min_step(1 * Unit::Second)
            .build(),
    );
    let final_state = prop
        .with(sc, almanac.clone())
        .for_duration(30 * Unit::Day)
        .unwrap();

    println!(
        "[q_law_as_kluever_case] peri: {} km\tapo: {} km",
        final_state.orbit.periapsis_km().unwrap(),
        final_state.orbit.apoapsis_km().unwrap()
    );

    let delta = final_state.epoch() - sc.epoch();

    println!(
        "[q_law_as_kluever_case_a] {delta} later SMA: {} km",
        final_state.orbit.sma_km().unwrap()
    );
    println!(
        "[q_law_as_kluever_case_a] {delta} later Inc: {} deg",
        final_state.orbit.inc_deg().unwrap()
    );
    let prop_usage = sc.mass.prop_mass_kg - final_state.mass.prop_mass_kg;
    println!("[q_law_as_kluever_case_a] prop usage: {prop_usage} kg");
}
