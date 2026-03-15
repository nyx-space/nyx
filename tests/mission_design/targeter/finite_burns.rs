extern crate nyx_space as nyx;

use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MOON, SUN};
use hifitime::TimeUnits;
use nyx::dynamics::guidance::{LocalFrame, Maneuver, Thruster};
use nyx::dynamics::sequence::{
    AccelModels, ForceModels, GuidanceConfig, Phase, PropagatorConfig, SpacecraftSequence,
};
use nyx::dynamics::PointMasses;
use nyx::linalg::Vector3;
use nyx::md::prelude::*;
use nyx::propagators::{IntegratorMethod, IntegratorOptions};
use nyx_space::cosmic::Mass;

use crate::propagation::GMAT_EARTH_GM;
use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn thrust_dir_tgt_sma_aop_raan(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k);

    let target_delta_t = 30.seconds();

    let spacecraft = Spacecraft {
        orbit: xi_orig,
        mass: Mass::from_dry_and_prop_masses(10.0, 90.0),
        thruster: Some(Thruster {
            thrust_N: 500.0,
            isp_s: 300.0,
        }),
        mode: GuidanceMode::Thrust,
        ..Default::default()
    };

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::default_dp78(dynamics);

    // Define the objective
    let objectives = [
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::SemiMajorAxis),
            8012.176,
            0.1,
        ),
        Objective::within_tolerance(StateParameter::Element(OrbitalElement::AoP), 53.939, 1e-3),
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::RAAN),
            60.000182,
            1e-3,
        ),
    ];

    let tgt = Targeter::thrust_dir(&setup, objectives);

    println!("{tgt}");

    let achievement_epoch = orig_dt + target_delta_t;

    let solution_fd = tgt
        .try_achieve_from(spacecraft, orig_dt, achievement_epoch, almanac)
        .unwrap();

    println!("Finite differencing solution: {solution_fd}");
}

#[rstest]
fn thrust_dir_rate_tgt_sma_aop_raan(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k);

    let target_delta_t = 30.seconds();

    let spacecraft = Spacecraft {
        orbit: xi_orig,
        mass: Mass::from_dry_and_prop_masses(10.0, 90.0),
        thruster: Some(Thruster {
            thrust_N: 500.0,
            isp_s: 300.0,
        }),
        mode: GuidanceMode::Thrust,
        ..Default::default()
    };

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::default_dp78(dynamics);

    // Define the objective
    let objectives = [
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::SemiMajorAxis),
            8012.176,
            0.1,
        ),
        Objective::within_tolerance(StateParameter::Element(OrbitalElement::AoP), 53.939, 1e-2),
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::RAAN),
            60.000182,
            1e-3,
        ),
    ];

    let tgt = Targeter::thrust_dir_rate(&setup, objectives);

    println!("{tgt}");

    let achievement_epoch = orig_dt + target_delta_t;

    let solution = tgt
        .try_achieve_from(spacecraft, orig_dt, achievement_epoch, almanac.clone())
        .unwrap();

    println!("Finite differencing solution: {solution}");
    tgt.apply(&solution, almanac).unwrap();
}

#[ignore]
#[rstest]
fn thrust_profile_tgt_sma_aop_raan(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k);

    let target_delta_t = 30.seconds();

    let spacecraft = Spacecraft {
        orbit: xi_orig,
        mass: Mass::from_dry_and_prop_masses(10.0, 90.0),
        thruster: Some(Thruster {
            thrust_N: 500.0,
            isp_s: 300.0,
        }),
        mode: GuidanceMode::Thrust,
        ..Default::default()
    };

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::default_dp78(dynamics);

    // Define the objective
    let objectives = [
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::SemiMajorAxis),
            8012.176,
            0.1,
        ),
        Objective::within_tolerance(StateParameter::Element(OrbitalElement::AoP), 53.939, 1e-1),
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::RAAN),
            60.000182,
            1e-1,
        ),
    ];

    let tgt = Targeter::thrust_profile(&setup, objectives);

    println!("{tgt}");

    let achievement_epoch = orig_dt + target_delta_t;

    let solution_fd = tgt
        .try_achieve_from(spacecraft, orig_dt, achievement_epoch, almanac)
        .unwrap();

    println!("Finite differencing solution: {solution_fd}");
}

#[rstest]
fn val_tgt_finite_burn(almanac: Arc<Almanac>) {
    // In this test, we take a known finite burn solution and use the optimizer to solve for it.
    let _ = pretty_env_logger::try_init();

    let eme2k = almanac
        .frame_info(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    // Build the initial spacecraft state
    let start_time = Epoch::from_gregorian_tai_at_midnight(2002, 1, 1);
    let orbit = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_time, eme2k,
    );

    // Define the thruster
    let monoprop = Thruster {
        thrust_N: 5000.0,
        isp_s: 300.0,
    };
    let dry_mass = 1e3;
    let fuel_mass = 756.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, monoprop, GuidanceMode::Thrust);

    let prop_time = 15.0 * Unit::Second;

    // Define the dynamics
    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let orbital_dyn = OrbitalDynamics::point_masses(bodies.clone());

    // With 100% thrust: RSS errors:     pos = 3.14651e1 km      vel = 3.75245e-2 km/s

    // Define the maneuver and its schedule
    let mnvr0 = Maneuver::from_time_invariant(
        start_time + 1.seconds(),
        start_time + prop_time - 1.seconds(),
        1.0, // Full thrust
        Vector3::new(1.0, 0.0, 0.0),
        LocalFrame::Inertial,
    );

    let mut sc_seq = SpacecraftSequence::default();

    sc_seq.propagators.insert(
        "Earth".to_string(),
        PropagatorConfig {
            method: IntegratorMethod::RungeKutta89,
            options: IntegratorOptions::default(),
            accel_models: AccelModels {
                point_masses: Some(PointMasses::new(bodies)),
                gravity_field: None,
            },
            force_models: ForceModels {
                solar_pressure: None,
                drag: None,
            },
        },
    );

    sc_seq
        .thruster_sets
        .insert("Monoprop".to_string(), monoprop);

    sc_seq.seq.insert(
        start_time,
        Phase::Activity {
            name: "Initial Coast".to_string(),
            propagator: "Earth".to_string(),
            guidance: None,
            on_entry: None,
            disabled: false,
        },
    );

    sc_seq.seq.insert(
        mnvr0.start,
        Phase::Activity {
            name: "Burn".to_string(),
            propagator: "Earth".to_string(),
            guidance: Some(Box::new(GuidanceConfig::FiniteBurn {
                maneuver: mnvr0,
                thruster_model: "Monoprop".to_string(),
                disable_prop_mass: false,
            })),
            on_entry: None,
            disabled: false,
        },
    );

    sc_seq.seq.insert(
        mnvr0.end,
        Phase::Activity {
            name: "Final Coast".to_string(),
            propagator: "Earth".to_string(),
            guidance: None,
            on_entry: None,
            disabled: false,
        },
    );

    sc_seq.seq.insert(start_time + prop_time, Phase::Terminate);

    sc_seq.setup(almanac.clone()).unwrap();

    let trajectories = sc_seq.propagate(sc_state, None, almanac.clone()).unwrap();
    let sc_xf_desired = trajectories.last().unwrap().last();
    println!("started: {sc_state}\nended   :{sc_xf_desired}");

    // Build an impulsive targeter for this known solution
    let sc_no_thrust = SpacecraftDynamics::new(orbital_dyn);
    let mut prop_no_thrust = Propagator::default(sc_no_thrust);
    prop_no_thrust.set_max_step(mnvr0.duration());
    let impulsive_tgt = Targeter::delta_v(
        &prop_no_thrust,
        [
            Objective::within_tolerance(
                StateParameter::Element(OrbitalElement::X),
                sc_xf_desired.orbit.radius_km.x,
                1e-5,
            ),
            Objective::within_tolerance(
                StateParameter::Element(OrbitalElement::Y),
                sc_xf_desired.orbit.radius_km.y,
                1e-5,
            ),
            Objective::within_tolerance(
                StateParameter::Element(OrbitalElement::Z),
                sc_xf_desired.orbit.radius_km.z,
                1e-5,
            ),
        ],
    )
    .try_achieve_from(
        sc_state,
        sc_state.epoch(),
        sc_xf_desired.epoch(),
        almanac.clone(),
    )
    .unwrap();

    println!("{impulsive_tgt}");
    println!("\n\nKNOWN SOLUTION\n{mnvr0}");
}
