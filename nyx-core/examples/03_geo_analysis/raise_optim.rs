#![doc = include_str!("./README.md")]
extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use anise::{
    almanac::{Almanac, metaload::MetaFile},
    constants::{
        celestial_objects::{MOON, SUN},
        frames::{EARTH_J2000, IAU_EARTH_FRAME},
    },
};
use hifitime::{Epoch, TimeUnits, Unit};
use log::info;
use nyx::{
    Spacecraft,
    cosmic::{GuidanceMode, Mass, MetaAlmanac, Orbit, SRPData},
    dynamics::{
        GravityField, OrbitalDynamics, SolarPressure, SpacecraftDynamics,
        guidance::{Ruggiero, Thruster},
    },
    io::gravity::GravityFieldData,
    md::prelude::{Objective, OrbitalElement, StateParameter},
    propagators::{ErrorControl, IntegratorOptions, Propagator},
};
use radiate::problem::EngineProblem;
use radiate::*;
use std::{error::Error, sync::Arc};

// Shared state struct for the fitness evaluation to avoid reading files thousands of times
struct SharedState {
    almanac: Arc<Almanac>,
    harmonics: Arc<GravityField>,
    srp_dyn: Arc<SolarPressure>,
}

impl SharedState {
    fn new() -> Result<Self, Box<dyn Error>> {
        let almanac = Arc::new(MetaAlmanac::latest().map_err(Box::new)?);

        let mut jgm3_meta = MetaFile {
            uri: "http://public-data.nyxspace.com/nyx/models/JGM3.cof.gz".to_string(),
            crc32: Some(0xF446F027),
        };
        jgm3_meta.process(true)?;

        let harmonics = GravityField::new(GravityFieldData::from_cof(
            &jgm3_meta.uri,
            4,
            4,
            true,
            almanac.frame_info(IAU_EARTH_FRAME)?,
        )?);
        let srp_dyn = SolarPressure::default_flux(EARTH_J2000, &almanac)?;

        Ok(Self {
            almanac,
            harmonics,
            srp_dyn,
        })
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    pel::init();

    /* let (prop_usage_kg, penalty) = evaluate_weights(
        &[0.22033301, 0.8512096, 0.49421895],
        60.0,
        Arc::new(SharedState::new()?),
    )
    .unwrap();

    println!("Best weight prop usage = {prop_usage_kg:.3} kg \t penalty = {penalty:.3}"); */

    // Set up shared state (read large files only once!)
    let shared_state = Arc::new(SharedState::new()?);

    // Set up the genetic algorithm optimization
    let codec = FloatCodec::vector(3, 0.1_f32..1.0_f32); // 3 weights for SMA, Ecc, Inc
    let problem = EngineProblem {
        objective: radiate::Objective::Multi(vec![Optimize::Minimize, Optimize::Minimize]), // NSGA2 Multi Objective
        codec: Arc::new(codec),
        fitness_fn: Some(Arc::new(move |weights: Vec<f32>| {
            // Full 60 days propagation for evaluating the actual performance, but running fast due to shared state
            let (prop_usage, penalty) =
                evaluate_weights(&weights, 60.0, shared_state.clone()).unwrap_or((1e6, 1e6));
            Score::from(vec![prop_usage as f32, penalty as f32])
        })),
        raw_fitness_fn: None,
    };

    let mut engine = GeneticEngine::<FloatChromosome<f32>, Vec<f32>>::builder()
        .population_size(20)
        .parallel()
        .multi_objective(vec![Optimize::Minimize, Optimize::Minimize])
        .problem(problem)
        .survivor_selector(NSGA2Selector::new())
        .build();

    // Wrap the engine with the UI
    let final_generation = engine.run(|generation: &Generation<FloatChromosome<f32>, Vec<f32>>| {
        let scores = generation.score().as_slice();
        println!(
            "[ {:?} ]: Best Score: Prop usage {:.3} kg, Penalty {:.3}",
            generation.index(),
            scores[0],
            scores[1]
        );
        generation.index() >= 5
    });

    let best_weights = final_generation
        .value()
        .iter()
        .map(|w| format!("W: = {w}"))
        .collect::<Vec<String>>()
        .join(", ");
    let best_score = final_generation
        .score()
        .iter()
        .enumerate()
        .map(|(i, w)| format!("S[{i}]: = {w}"))
        .collect::<Vec<String>>()
        .join(", ");
    println!("Optimization finished. Best weights: [{best_weights}] -> Best score: [{best_score}]");

    // Evaluate these weights.
    let best_weights: Vec<f32> = final_generation.value().to_vec();

    let (prop_usage_kg, penalty) =
        evaluate_weights(&best_weights, 60.0, Arc::new(SharedState::new()?)).unwrap();

    println!("Best weight prop usage = {prop_usage_kg:.3} kg \t penalty = {penalty:.3}");

    Ok(())
}

fn evaluate_weights(
    weights: &[f32],
    prop_time_days: f64,
    state: Arc<SharedState>,
) -> Result<(f64, f64), Box<dyn Error>> {
    let ηthresholds: Vec<f64> = weights.iter().map(|w| *w as f64).collect();

    let eme2k = state.almanac.frame_info(EARTH_J2000).unwrap();
    let epoch = Epoch::from_gregorian_utc_hms(2024, 2, 29, 12, 13, 14);

    let orbit = Orbit::keplerian(24505.9, 0.725, 7.05, 0.0, 0.0, 0.0, epoch, eme2k);

    let sc = Spacecraft::builder()
        .orbit(orbit)
        .mass(Mass::from_dry_and_prop_masses(1000.0, 1000.0))
        .srp(SRPData::from_area(3.0 * 6.0))
        .thruster(Thruster {
            isp_s: 4435.0,
            thrust_N: 0.472,
        })
        .mode(GuidanceMode::Thrust)
        .build();

    let prop_time = prop_time_days * Unit::Day;

    let objectives = &[
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::SemiMajorAxis),
            30_000.0,
            20.0,
        ),
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::Eccentricity),
            0.001,
            5e-5,
        ),
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::Inclination),
            0.05,
            1e-2,
        ),
    ];

    let ctrl = Ruggiero::from_ηthresholds(objectives, &ηthresholds, sc)?;

    let mut orbital_dyn = OrbitalDynamics::point_masses(vec![MOON, SUN]);
    orbital_dyn.accel_models.push(state.harmonics.clone());

    let sc_dynamics = SpacecraftDynamics::from_model(orbital_dyn, state.srp_dyn.clone())
        .with_guidance_law(ctrl.clone());

    let (final_state, _traj) = Propagator::rk89(
        sc_dynamics.clone(),
        IntegratorOptions::builder()
            .min_step(10.0_f64.seconds())
            .tolerance(1e-8)
            .error_ctrl(ErrorControl::RSSCartesianStep)
            .build(),
    )
    .with(sc, state.almanac.clone())
    .for_duration_with_traj(prop_time)?;

    let prop_usage = sc.mass.prop_mass_kg - final_state.mass.prop_mass_kg;

    let mut penalty = 0.0;
    for obj in objectives {
        let (achieved, error) = obj.assess(&final_state)?;
        if !achieved {
            penalty += error.abs();
        }
        info!("{obj} error: {error:.3}, achieved? {achieved}");
    }

    info!("{ηthresholds:?} -> {prop_usage:.3} kg\tpenalty = {penalty:.3}");

    Ok((prop_usage, penalty * 1000.0))
}
