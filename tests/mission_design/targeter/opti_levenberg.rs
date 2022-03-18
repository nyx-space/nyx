extern crate nyx_space as nyx;

use nyx::md::optimizer::*;
use nyx::md::ui::*;

// Semi major axis

#[test]
fn tgt_levenberg_sma_from_apo() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 180.0, orig_dt, eme2k);

    let target_delta_t: Duration = xi_orig.period() / 2.0;

    println!("Period: {} s", xi_orig.period().in_seconds() / 2.0);

    let spacecraft = Spacecraft::from_srp_defaults(xi_orig, 100.0, 0.0);

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::default_dp78(dynamics);

    // Try to increase SMA
    let xf_desired_sma = 8_100.0;
    let xf_desired_ecc = 0.40;
    let xf_desired_aop = 60.0;

    // Define the objective
    let objectives = [
        Objective::new(StateParameter::SMA, xf_desired_sma),
        Objective::new(StateParameter::AoP, xf_desired_aop),
        Objective::new(StateParameter::Eccentricity, xf_desired_ecc),
    ];

    let tgt = Optimizer::delta_v(&setup, objectives);

    println!("{}", tgt);

    tgt.minimize(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    // let gmat_sol = 0.05312024615278713;
    // // GMAT validation
    // assert!(
    //     dbg!(solution_fd.correction.norm() - gmat_sol).abs() < 1e-6,
    //     "Finite differencing result different from GMAT (greater than 1 mm/s)."
    // );

    // // Check that the solutions nearly match
    // println!(
    //     "GMAT validation - tgt_sma_from_apo: Δv = {:.3} m/s\terr = {:.6} m/s (better = {})",
    //     solution_fd.correction.norm() * 1e3,
    //     (solution_fd.correction.norm() - gmat_sol).abs() * 1e3,
    //     solution_fd.correction.norm() < gmat_sol
    // );
}
