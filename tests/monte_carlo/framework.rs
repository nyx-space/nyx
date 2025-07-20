extern crate nyx_space as nyx;

use nyx::mc::*;
use nyx::md::prelude::*;

use anise::{
    constants::{
        celestial_objects::{JUPITER_BARYCENTER, MOON, SUN},
        frames::EARTH_J2000,
    },
    prelude::Almanac,
};
use rstest::*;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn test_monte_carlo_epoch_cov_test(almanac: Arc<Almanac>) {
    extern crate pretty_env_logger;
    let _ = pretty_env_logger::try_init();

    // Build the initial state

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
    let state = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);

    // Build the state generator using a Gaussian distribution (you may use any distribution from rand_distr)
    // 5% error on SMA and 5% on Eccentricity
    let nominal_state = Spacecraft::from(state);

    let random_state = MvnSpacecraft::new(
        nominal_state,
        vec![
            StateDispersion::zero_mean(StateParameter::SMA, 0.05),
            StateDispersion::zero_mean(StateParameter::Eccentricity, 0.05),
        ],
    )
    .unwrap();

    // Set up the dynamics
    let orbital_dyn = SpacecraftDynamics::new(OrbitalDynamics::new(vec![PointMasses::new(vec![
        SUN,
        MOON,
        JUPITER_BARYCENTER,
    ])]));

    let prop = Propagator::default_dp78(orbital_dyn);

    // Setup the Monte Carlo

    let my_mc = MonteCarlo {
        random_state,
        nominal_state,
        seed: Some(0),
        scenario: "test_monte_carlo_epoch".to_string(),
    };

    let rslts = my_mc.run_until_epoch(prop, almanac.clone(), dt + 1.0_f64 * Unit::Day, 10);

    let average_sma_dispersion = rslts
        .dispersion_values_of(StateParameter::SMA)
        .unwrap()
        .amean()
        .unwrap();

    let average_sma = rslts
        .every_value_of(StateParameter::SMA, 5.minutes(), None)
        .amean()
        .unwrap();

    let average_initial_sma = rslts
        .first_values_of(StateParameter::SMA, None)
        .amean()
        .unwrap();

    let average_final_sma = rslts
        .last_values_of(StateParameter::SMA, None)
        .amean()
        .unwrap();

    println!("Average SMA dispersion = {average_sma_dispersion} km");
    println!("Average initial SMA = {average_initial_sma} km");
    println!("Average final SMA = {average_final_sma} km");
    println!("Average SMA = {average_sma} km");
}
