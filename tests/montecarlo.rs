extern crate nyx_space as nyx;
extern crate rand;
extern crate rand_distr;

use nyx::celestia::{Bodies, Cosm, Orbit};
use nyx::dynamics::Harmonics;
use nyx::dynamics::{Dynamics, OrbitalDynamics};
use nyx::io::gravity::*;
use nyx::propagators::*;
use nyx::time::{Epoch, TimeUnit};
use rand::thread_rng;
use rand_distr::{Distribution, Normal};

// #[allow(clippy::identity_op)]
// #[test]
// fn multi_thread_monte_carlo_demo() {
//     /*
//     In this demo, we'll be running a 100 runs with the same dynamics and end state, but with a slightly variation in eccentricity.
//     */
//     extern crate pretty_env_logger;
//     if pretty_env_logger::try_init().is_err() {
//         println!("could not init env_logger");
//     }
//     use nyx::dynamics::Harmonics;
//     use nyx::io::gravity::*;

//     let cosm = Cosm::de438();
//     let eme2k = cosm.frame("EME2000");
//     let iau_earth = cosm.frame("IAU Earth");

//     let earth_sph_harm = HarmonicsMem::from_cof("data/JGM3.cof.gz", 70, 70, true).unwrap();
//     let harmonics = Harmonics::from_stor(iau_earth, earth_sph_harm, &cosm);

//     let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
//     let state = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);

//     // TODO: Change point masses to be a MOVED vector.
//     let mut orbital_dyn = OrbitalDynamics::point_masses(
//         state.frame,
//         &vec![Bodies::Sun, Bodies::Luna, Bodies::JupiterBarycenter],
//         &cosm,
//     );
//     orbital_dyn.add_model(harmonics);

//     let setup = Propagator::default(&orbital_dyn);

//     let mut threads = vec![];
//     let mut final_states: Vec<Orbit> = vec![];

//     // Around 1 km of error
//     let sma_dist = Normal::new(0.0, 1.0).unwrap();

//     // Generate all 100 initial states
//     let init_states: Vec<Orbit> = sma_dist
//         .sample_iter(&mut thread_rng())
//         .take(100)
//         .map(|delta_sma| state.with_sma(state.sma() + delta_sma))
//         .collect();

//     for state in init_states {
//         threads.push(std::thread::spawn(move || {
//             let setp = setup.clone();
//             setup.with(state).for_duration(1 * TimeUnit::Day).unwrap()
//         }));
//     }

//     for t in threads {
//         final_states.push(t.join().unwrap());
//     }

//     println!("Collected {} states", final_states.len());

//     // for _ in (0..100) {
//     //     let mut this_state = state;
//     //     this_state.sma += sma_dist.sample_iter(rng: R)
//     // }
//     // .with(state)
//     // .for_duration(1 * TimeUnit::Day)
//     // .unwrap();
// }
