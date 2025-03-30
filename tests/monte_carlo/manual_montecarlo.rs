extern crate nyx_space as nyx;
extern crate rand;
extern crate rand_distr;
extern crate rayon;

use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MOON, SUN};
use anise::constants::frames::IAU_EARTH_FRAME;
use nyx::cosmic::Orbit;
use nyx::dynamics::{Harmonics, SpacecraftDynamics};
use nyx::dynamics::{OrbitalDynamics, PointMasses};
use nyx::io::gravity::*;
use nyx::time::{Epoch, Unit};
use nyx::State;
use nyx::{propagators::*, Spacecraft};
use rand::thread_rng;
use rand_distr::{Distribution, Normal};
use rayon::prelude::*;
use std::sync::Arc;
use std::time::Instant as StdInstant;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[allow(clippy::identity_op)]
#[rstest]
fn multi_thread_monte_carlo_demo(almanac: Arc<Almanac>) {
    /*
    In this demo, we'll be running a 100 runs with the same dynamics and end state, but with a slightly variation in eccentricity.
    */
    extern crate pretty_env_logger;
    let _ = pretty_env_logger::try_init();

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();

    let earth_sph_harm =
        HarmonicsMem::from_cof("data/01_planetary/JGM3.cof.gz", 70, 70, true).unwrap();
    let harmonics = Harmonics::from_stor(iau_earth, earth_sph_harm);

    let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
    let state = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);

    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::new(vec![
        PointMasses::new(vec![SUN, MOON, JUPITER_BARYCENTER]),
        harmonics,
    ]));

    // We need to wrap the propagator setup in an Arc to enable multithreading.
    let setup = Arc::new(Propagator::default_dp78(orbital_dyn));

    // Around 1 km of error
    let sma_dist = Normal::new(0.0, 1.0).unwrap();

    // Generate all 100 initial states
    let init_states: Vec<Spacecraft> = sma_dist
        .sample_iter(&mut thread_rng())
        .take(100)
        .map(|delta_sma| {
            Spacecraft::from(
                state
                    .with_sma_km(state.sma_km().unwrap() + delta_sma)
                    .unwrap(),
            )
        })
        .collect();

    let prop_time = 1 * Unit::Day;
    let start = StdInstant::now();
    let end_epoch = dt + prop_time;
    init_states
        .par_iter()
        .for_each_with((setup, almanac), |(setup, almanac), state| {
            let final_state = setup
                .with(*state, almanac.clone())
                .for_duration(prop_time)
                .unwrap();
            assert_eq!(end_epoch, final_state.epoch());
        });

    let clock_time = StdInstant::now() - start;
    println!(
        "Propagated {} states in {} seconds",
        init_states.len(),
        clock_time.as_secs_f64()
    );
}
