#![doc = include_str!("./README.md")]
extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use anise::{
    almanac::{metaload::MetaFile, Almanac},
    constants::{
        celestial_objects::{MOON, SUN},
        frames::{EARTH_J2000, IAU_EARTH_FRAME},
    },
};
use hifitime::{Epoch, TimeUnits, Unit};
use log::info;
use nyx::{
    cosmic::{GuidanceMode, Mass, MetaAlmanac, Orbit, SRPData},
    dynamics::{
        guidance::{Ruggiero, Thruster},
        GravityField, OrbitalDynamics, SolarPressure, SpacecraftDynamics,
    },
    io::gravity::GravityFieldData,
    md::prelude::{Objective, OrbitalElement, StateParameter},
    propagators::{ErrorControl, IntegratorOptions, Propagator},
    Spacecraft,
};
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

        let harmonics = GravityField::from_stor(
            almanac.frame_info(IAU_EARTH_FRAME)?,
            GravityFieldData::from_cof(&jgm3_meta.uri, 8, 8, true)?,
        );
        let srp_dyn = SolarPressure::default_flux(EARTH_J2000, almanac.clone())?;

        Ok(Self {
            almanac,
            harmonics,
            srp_dyn,
        })
    }
}

fn main() -> Result<(), Box<dyn Error>> {
    pel::init();

    // Set up shared state (read large files only once!)
    let shared_state = Arc::new(SharedState::new()?);

    // Set up the genetic algorithm optimization
    let codec = FloatCodec::vector(3, 0.1_f32..1.0_f32); // 3 weights for SMA, Ecc, Inc
    let problem = EngineProblem {
        objective: radiate::Objective::Multi(vec![Optimize::Minimize, Optimize::Minimize]), // NSGA2 Multi Objective
        codec: Arc::new(codec.clone()),
        fitness_fn: Some(Arc::new(move |weights: Vec<f32>| {
            // Full 60 days propagation for evaluating the actual performance, but running fast due to shared state
            let (prop_usage, penalty) =
                evaluate_weights(&weights, 60.0, shared_state.clone()).unwrap_or((1e6, 1e6));
            Score::from(vec![prop_usage as f32, penalty as f32])
        })),
        raw_fitness_fn: None,
    };

    let mut engine = GeneticEngine::<FloatChromosome<f32>, Vec<f32>>::builder()
        .codec(codec)
        .executor(Executor::FixedSizedWorkerPool(8))
        .problem(problem)
        .population_size(20)
        .survivor_selector(NSGA2Selector::new())
        .build();

    let result = engine.run(|generation: &Generation<FloatChromosome<f32>, Vec<f32>>| {
        let scores = generation.score().as_slice();
        println!(
            "[ {:?} ]: Best Score: Prop usage {:.3} kg, Penalty {:.3}",
            generation.index(),
            scores[0],
            scores[1]
        );
        generation.index() >= 10 // Max generations
    });

    let best_as_string = result
        .value()
        .iter()
        .map(|w| w.to_string())
        .collect::<Vec<String>>()
        .join(", ");
    println!("Optimization finished. Best weights: [{}]", best_as_string);

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

    // let kluever_ctrl = Kluever::from_max_eclipse(objectives, &weights_f64, 0.2);
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
    }

    info!("{ηthresholds:?} -> {prop_usage:.3} kg\tpenalty = {penalty:.3}");

    Ok((prop_usage, penalty * 1000.0))
}
