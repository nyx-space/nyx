extern crate nyx_space as nyx;

use nyx::md::targeter::*;
use nyx::md::ui::*;

#[test]
fn tgt_sma_from_apo() {
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
    let setup = Propagator::default(dynamics);

    // Try to increase SMA
    let xf_desired = Orbit::keplerian(
        8_100.0,
        0.2,
        30.0,
        60.0,
        60.0,
        180.0,
        orig_dt + target_delta_t,
        eme2k,
    );

    // Define the objective
    let objectives = vec![Objective::new(StateParameter::SMA, xf_desired.sma())];

    let tgt = Targeter::delta_v(Arc::new(&setup), objectives);

    println!("{}", tgt);

    // let solution = tgt
    //     .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
    //     .unwrap();

    let solution_fd = tgt
        .try_achieve_from_with_guess(
            spacecraft,
            &[0.0, 0.0, 0.0],
            orig_dt,
            orig_dt + target_delta_t,
        )
        .unwrap();

    println!("Finite differencing solution: {}", solution_fd);

    let solution = tgt
        .try_achieve_from_with_guess_dual(
            spacecraft,
            &[0.0, 0.0, 0.0],
            orig_dt,
            orig_dt + target_delta_t,
        )
        .unwrap();

    println!("{}", solution);
}
