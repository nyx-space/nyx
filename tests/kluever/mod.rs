use std::sync::Arc;
use crate::nyx::cosmic::{GuidanceMode, Orbit, Spacecraft};
use crate::nyx::dynamics::guidance::{Kluever, Objective, Thruster};
use crate::nyx::dynamics::{OrbitalDynamics, SpacecraftDynamics};
use crate::nyx::md::prelude::{OrbitalElement, StateParameter};
use crate::nyx::propagators::{Propagator, IntegratorOptions, IntegratorMethod};
use crate::nyx::time::{Epoch, Unit};
use anise::constants::frames::EARTH_J2000;
use anise::prelude::Almanac;

#[test]
fn test_kluever_raise() {
    let almanac = Arc::new(Almanac::default());
    let eme2k = EARTH_J2000.with_mu_km3_s2(398_600.433);
    let epoch = Epoch::from_gregorian_utc_hms(2024, 2, 29, 12, 0, 0);

    // LEO start
    let orbit = Orbit::keplerian(7000.0, 0.01, 28.5, 0.0, 0.0, 0.0, epoch, eme2k);

    let sc = Spacecraft::builder()
        .orbit(orbit)
        .mass(crate::nyx::cosmic::Mass::from_dry_and_prop_masses(1000.0, 500.0))
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
            42164.0,
            10.0,
        ),
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::Inclination),
            0.0,
            0.1,
        ),
    ];
    let weights = &[1.0, 0.1];

    let kluever = Kluever::new(objectives, weights);

    let orbital_dyn = OrbitalDynamics::two_body();
    let sc_dynamics = SpacecraftDynamics::from_guidance_law(orbital_dyn, kluever);

    // Propagate for a short time (e.g., 1 hour) and check that SMA increases
    let prop_time = 1.0 * Unit::Hour;
    let final_state = Propagator::new(
        sc_dynamics,
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    )
    .with(sc, almanac.clone())
    .for_duration(prop_time)
    .unwrap();

    println!("Initial SMA: {} km", orbit.sma_km().unwrap());
    println!("Final SMA: {} km", final_state.orbit.sma_km().unwrap());
    println!("Initial Inc: {} deg", orbit.inc_deg().unwrap());
    println!("Final Inc: {} deg", final_state.orbit.inc_deg().unwrap());

    assert!(final_state.orbit.sma_km().unwrap() > orbit.sma_km().unwrap(), "SMA should have increased");
}
