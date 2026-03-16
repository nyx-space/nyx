extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use anise::constants::celestial_objects::{EARTH, MOON};
use anise::constants::frames::{IAU_EARTH_FRAME, MOON_J2000};
use anise::math::Vector3;
use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use nyx::cosmic::{Orbit, Spacecraft};
use nyx::dynamics::guidance::Thruster;
use nyx::dynamics::sequence::*;
use nyx::dynamics::{Drag, SolarPressure};
use nyx::propagators::{IntegratorMethod, IntegratorOptions};
use nyx::time::{Epoch, Unit};
use nyx_space::cosmic::{Mass, SRPData};
use nyx_space::dynamics::guidance::mnvr::ImpulsiveManeuver;
use nyx_space::dynamics::guidance::{LocalFrame, Maneuver};
use nyx_space::dynamics::PointMasses;
use nyx_space::io::gravity::GravityFieldConfig;
use nyx_space::md::prelude::{Objective, OrbitalElement};
use nyx_space::md::StateParameter;
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

/// *** *** ***
/// Spacecraft Sequence is a PREDETERMINED and ABSOLUTE sequence of events.
/// Its is useful for simulating different parts of the mission given some variation in the state WITHOUT retargeting.
/// It will eventually support orbit determination thereby allowing a team to freeze the nominal design in a git-friendly
/// format and then feed it tracking data for estimation.
#[rstest]
fn spacecraft_sequence(almanac: Arc<Almanac>) {
    let _ = pel::try_init();
    // The sequence is in ABSOLUTE epochs, so let's define the starting epoch now.
    let epoch = Epoch::from_gregorian_utc_at_midnight(2010, 12, 21);

    // Build an orbit raising sequence from deployment until station.
    let mut sc_seq = SpacecraftSequence::default();
    // Build the propagators.
    sc_seq.propagators.insert(
        "Near Earth".to_string(),
        PropagatorConfig {
            method: IntegratorMethod::RungeKutta89,
            options: IntegratorOptions::default(),
            accel_models: AccelModels {
                point_masses: Some(PointMasses::new(vec![EARTH, MOON])),
                gravity_field: Some((
                    GravityFieldConfig {
                        filepath: "data/01_planetary/EGM2008_to2190_TideFree.gz".into(),
                        gunzipped: true,
                        degree: 21,
                        order: 21,
                    },
                    IAU_EARTH_FRAME.into(),
                )),
            },
            force_models: ForceModels {
                solar_pressure: None,
                drag: Some(Drag::std_atm1976(almanac.clone()).unwrap()),
            },
        },
    );

    sc_seq.propagators.insert(
        "Cislunar".to_string(),
        PropagatorConfig {
            method: IntegratorMethod::RungeKutta89,
            options: IntegratorOptions::default(),
            accel_models: AccelModels {
                point_masses: Some(PointMasses::new(vec![EARTH, MOON])),
                gravity_field: Some((
                    GravityFieldConfig {
                        filepath: "data/01_planetary/EGM2008_to2190_TideFree.gz".into(),
                        gunzipped: true,
                        degree: 8,
                        order: 8,
                    },
                    IAU_EARTH_FRAME.into(),
                )),
            },
            force_models: ForceModels {
                solar_pressure: Some(
                    SolarPressure::default_no_estimation(
                        vec![EARTH_J2000, MOON_J2000],
                        almanac.clone(),
                    )
                    .unwrap(),
                ),
                drag: None,
            },
        },
    );

    // Setup the thruster models
    sc_seq.thruster_sets.insert(
        "BiProp".to_string(),
        Thruster {
            thrust_N: 25.0,
            isp_s: 300.0,
        },
    );

    // Setup the spacecraft sequence
    sc_seq.seq.insert(
        epoch,
        Phase::Activity {
            name: "Parking orbit checkout".to_string(),
            propagator: "Near Earth".into(),
            guidance: None,
            on_entry: None,
            disabled: false, // Easy to temporarily disable a given phase
        },
    );

    sc_seq.seq.insert(
        epoch + Unit::Hour * 1.5,
        Phase::Activity {
            name: "Separation and vehicle checkout".to_string(),
            propagator: "Near Earth".to_string(),
            guidance: None,
            on_entry: Some(Box::new(DiscreteEvent::Staging {
                // Separation event is a discrete event
                impulsive_maneuver: Some(ImpulsiveManeuver {
                    dv_km_s: Vector3::new(25.0e-6, 0.0, 0.0),
                    local_frame: LocalFrame::VNC,
                }),
                decrement_properties: None,
            })),
            disabled: false,
        },
    );

    let mnvr_start = epoch + Unit::Day * 1 + Unit::Hour * 1.5;
    let mnvr_end = mnvr_start + Unit::Second * 45;

    sc_seq.seq.insert(
        mnvr_start,
        Phase::Activity {
            name: "Finite Maneuver".to_string(),
            propagator: "Cislunar".to_string(),
            guidance: Some(Box::new(GuidanceConfig::FiniteBurn {
                maneuver: Maneuver::from_time_invariant(
                    mnvr_start,
                    mnvr_end,
                    1.0,
                    Vector3::new(1.0, 0.0, 0.0),
                    LocalFrame::VNC,
                ),
                thruster_model: "BiProp".to_string(),
                disable_prop_mass: false,
            })),
            on_entry: None,
            disabled: false,
        },
    );

    // All sequences MUST end with a Terminate phase.
    sc_seq.seq.insert(epoch + Unit::Day * 30, Phase::Terminate);

    // Test Dhall serialization
    println!(
        "{}",
        serde_dhall::serialize(&sc_seq.propagators["Near Earth"])
            .static_type_annotation()
            .to_string()
            .unwrap()
    );
    // Initialize the propagators.
    sc_seq.setup(almanac.clone()).unwrap();

    // Set up the initial state

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let orbit =
        Orbit::try_keplerian_altitude(300.0, 2e-4, 28.5, 10.0, 0.0, 0.0, epoch, eme2k).unwrap();

    let sc = Spacecraft::builder()
        .srp(SRPData {
            area_m2: 16.0,
            coeff_reflectivity: 1.2,
        })
        .mass(Mass::from_dry_and_prop_masses(300.0, 250.0))
        .orbit(orbit)
        .build();

    println!("{sc:x}");

    // Propagate, returning a vector of trajectories, one per phase.
    let trajectories = sc_seq.propagate(sc, None, almanac.clone()).unwrap();
    // Minus one because there is no trajectory for the terminate phase.
    assert_eq!(trajectories.len(), sc_seq.seq.len() - 1);
}

#[rstest]
#[case(GuidanceConfig::Ruggiero {
    thruster_model: "HET".to_string(),
    disable_prop_mass: false,
    objectives: vec![
        (
            Objective::new(
                StateParameter::Element(OrbitalElement::SemiMajorAxis),
                7_300.0,
            ),
            0.0,
        ),
        (
            Objective::new(StateParameter::Element(OrbitalElement::Eccentricity), 1e-4),
            0.0,
        ),
    ],
    max_eclipse_prct: Some(0.5),
},
    GuidanceConfig::Ruggiero {
    thruster_model: "HET".to_string(),
    disable_prop_mass: false,
    objectives: vec![
        (
            Objective::new(
                StateParameter::Element(OrbitalElement::SemiMajorAxis),
                8_000.0,
            ),
            0.0,
        ),
        (
            Objective::new(StateParameter::Element(OrbitalElement::Eccentricity), 1e-4),
            0.0,
        ),
        (
            Objective::new(StateParameter::Element(OrbitalElement::Inclination), 35.0),
            0.0,
        ),
    ],
    max_eclipse_prct: Some(0.5),}, "Ruggiero")]
#[case(GuidanceConfig::Kluever {
    thruster_model: "HET".to_string(),
    disable_prop_mass: false,
    objectives: vec![
        (
            Objective::new(
                StateParameter::Element(OrbitalElement::SemiMajorAxis),
                7_300.0,
            ),
            1.0,
        ),
        (
            Objective::new(StateParameter::Element(OrbitalElement::Eccentricity), 1e-4),
            1.0,
        ),
    ],
    max_eclipse_prct: Some(0.5),
},
    GuidanceConfig::Kluever {
    thruster_model: "HET".to_string(),
    disable_prop_mass: false,
    objectives: vec![
        (
            Objective::new(
                StateParameter::Element(OrbitalElement::SemiMajorAxis),
                8_000.0,
            ),
            1.0,
        ),
        (
            Objective::new(StateParameter::Element(OrbitalElement::Eccentricity), 1e-4),
            1.0,
        ),
        (
            Objective::new(StateParameter::Element(OrbitalElement::Inclination), 35.0),
            1.0,
        ),
    ],
    max_eclipse_prct: Some(0.5),}, "Kluever")]
fn spacecraft_low_thrust_orbit_raise(
    #[case] first_raise: GuidanceConfig,
    #[case] second_raise: GuidanceConfig,
    #[case] alias: &'static str,
    almanac: Arc<Almanac>,
) {
    let _ = pel::try_init();
    // The sequence is in ABSOLUTE epochs, so let's define the starting epoch now.
    let epoch = Epoch::from_gregorian_utc_at_midnight(2010, 12, 21);

    // Build an orbit raising sequence from deployment until station.
    let mut sc_seq = SpacecraftSequence::default();
    // Build the propagators.
    sc_seq.propagators.insert(
        "Earth MPOP".to_string(),
        PropagatorConfig {
            method: IntegratorMethod::DormandPrince78,
            options: IntegratorOptions::builder()
                .min_step(Unit::Second * 1)
                .build(),
            accel_models: AccelModels {
                point_masses: Some(PointMasses::new(vec![EARTH, MOON])),
                gravity_field: Some((
                    GravityFieldConfig {
                        filepath: "data/01_planetary/EGM2008_to2190_TideFree.gz".into(),
                        gunzipped: true,
                        degree: 8,
                        order: 8,
                    },
                    IAU_EARTH_FRAME.into(),
                )),
            },
            force_models: ForceModels {
                solar_pressure: None,
                drag: Some(Drag::std_atm1976(almanac.clone()).unwrap()),
            },
        },
    );

    // Setup the thruster models
    sc_seq.thruster_sets.insert(
        "HET".to_string(),
        Thruster {
            thrust_N: 0.5,
            isp_s: 3000.0,
        },
    );

    // Setup the spacecraft sequence
    sc_seq.seq.insert(
        epoch,
        Phase::Activity {
            name: "Parking orbit checkout".to_string(),
            propagator: "Earth MPOP".into(),
            guidance: None,
            on_entry: None,
            disabled: false,
        },
    );

    sc_seq.seq.insert(
        epoch + Unit::Hour * 1.5,
        Phase::Activity {
            name: "Separation and vehicle checkout".to_string(),
            propagator: "Earth MPOP".to_string(),
            guidance: None,
            on_entry: Some(Box::new(DiscreteEvent::Staging {
                impulsive_maneuver: Some(ImpulsiveManeuver {
                    dv_km_s: Vector3::new(25.0e-6, 0.0, 0.0),
                    local_frame: LocalFrame::VNC,
                }),
                decrement_properties: None,
            })),
            disabled: false,
        },
    );

    sc_seq.seq.insert(
        epoch + Unit::Day * 1 + Unit::Hour * 1.5,
        Phase::Activity {
            name: format!("Raise to higher checkout for {alias}"),
            propagator: "Earth MPOP".to_string(),
            guidance: Some(Box::new(first_raise)),
            on_entry: None,
            disabled: false,
        },
    );

    sc_seq.seq.insert(
        epoch + Unit::Day * 30,
        Phase::Activity {
            name: format!("Raise to station for {alias}"),
            propagator: "Earth MPOP".to_string(),
            guidance: Some(Box::new(second_raise)),
            on_entry: None,
            disabled: false,
        },
    );

    // All sequences MUST end with a Terminate phase.
    sc_seq.seq.insert(epoch + Unit::Day * 90, Phase::Terminate);

    // Initialize the propagators.
    sc_seq.setup(almanac.clone()).unwrap();

    // Set up the initial state

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let orbit =
        Orbit::try_keplerian_altitude(300.0, 2e-4, 28.5, 10.0, 0.0, 0.0, epoch, eme2k).unwrap();

    let sc = Spacecraft::builder()
        .srp(SRPData {
            area_m2: 16.0,
            coeff_reflectivity: 1.2,
        })
        .mass(Mass::from_dry_and_prop_masses(300.0, 250.0))
        .orbit(orbit)
        .build();

    println!("{sc:x}");

    // Propagate, returning a vector of trajectories, one per phase.
    let trajectories = sc_seq.propagate(sc, None, almanac.clone()).unwrap();
    // Minus one because there is no trajectory for the terminate phase.
    assert_eq!(trajectories.len(), sc_seq.seq.len() - 1);
}
