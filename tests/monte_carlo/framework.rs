extern crate nyx_space as nyx;

use nyx::mc::*;
use nyx::md::ui::*;

#[test]
fn test_monte_carlo_epoch() {
    extern crate pretty_env_logger;
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();

    // Build the innitial state

    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
    let state = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);

    // Build the state generator using a Gaussian distribution (you may use any distribution from rand_distr)
    // 5% error on SMA and 5% on Eccentricity
    let generator = GaussianGenerator::from_std_dev_prcts(
        state,
        &[
            (StateParameter::SMA, 0.05),
            (StateParameter::Eccentricity, 0.05),
        ],
    )
    .unwrap();

    // Set up the dynamics
    let orbital_dyn = OrbitalDynamics::new(vec![PointMasses::new(
        &[Bodies::Sun, Bodies::Luna, Bodies::JupiterBarycenter],
        cosm,
    )]);

    let prop = Propagator::default_dp78(orbital_dyn);

    // Setup the Monte Carlo

    let my_mc = MonteCarlo {
        generator,
        seed: 0,
        scenario: "test_monte_carlo_epoch".to_string(),
    };

    let rslts = my_mc.run_until_epoch(prop, dt + 1.0_f64 * Unit::Day, 100);

    let average_sma_dispersion = rslts
        .dispersion_values_of(&StateParameter::SMA)
        .unwrap()
        .amean()
        .unwrap();

    let average_sma = rslts
        .every_value_of(&StateParameter::SMA, 5.minutes(), None)
        .amean()
        .unwrap();

    let average_initial_sma = rslts
        .first_values_of(&StateParameter::SMA, None)
        .amean()
        .unwrap();

    let average_final_sma = rslts
        .last_values_of(&StateParameter::SMA, None)
        .amean()
        .unwrap();

    println!("Average SMA dispersion = {} km", average_sma_dispersion);
    println!("Average initial SMA = {} km", average_initial_sma);
    println!("Average final SMA = {} km", average_final_sma);
    println!("Average SMA = {} km", average_sma);
}
