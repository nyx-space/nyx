extern crate nyx_space as nyx;

use nyx::dynamics::guidance::Thruster;
use nyx::md::optimizer::*;
use nyx::md::ui::*;

#[test]
fn fb_tgt_sma_ecc() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k);

    let target_delta_t: Duration = xi_orig.period() / 2.0;

    let spacecraft = Spacecraft {
        orbit: xi_orig,
        dry_mass_kg: 10.0,
        fuel_mass_kg: 90.0,
        thruster: Some(Thruster {
            thrust: 500.0,
            isp: 300.0,
        }),
        mode: GuidanceMode::Thrust,
        ..Default::default()
    };

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::default(dynamics);

    // Define the objective
    let objectives = [
        Objective::within_tolerance(StateParameter::Eccentricity, 0.4, 1e-5),
        Objective::within_tolerance(StateParameter::SMA, 8100.0, 0.1),
    ];

    // The variables in this targeter
    let variables = [
        Variable::from(Vary::MnvrAlpha).with_initial_guess(-0.3021017411736592_f64.to_radians()),
        // Variable::from(Vary::MnvrAlphaDot).with_initial_guess(45.0),
        Variable::from(Vary::MnvrAlphaDDot)
            .with_initial_guess(-2.1098425649685995_f64.to_radians()),
        Variable::from(Vary::MnvrBeta).with_initial_guess(0.3530352682197084_f64.to_radians()),
        // Variable::from(Vary::MnvrBetaDot).with_initial_guess(45.0),
        Variable::from(Vary::MnvrBetaDDot)
            .with_initial_guess(4.152947118658474e-7_f64.to_radians()),
        // Variable::from(Vary::Duration).with_initial_guess(5.0),
    ];

    let tgt = Optimizer::new(&setup, variables, objectives);

    println!("{}", tgt);

    let achievement_epoch = orig_dt + target_delta_t;

    let solution_fd = tgt
        .try_achieve_from(spacecraft, orig_dt, achievement_epoch)
        .unwrap();

    println!("Finite differencing solution: {}", solution_fd);

    let gmat_sol = 3.1160765514523914;
    println!(
        "GMAT validation - tgt_sma_from_peri: Δv = {:.3} m/s\terr = {:.6} m/s",
        solution_fd.correction.norm() * 1e3,
        (solution_fd.correction.norm() - gmat_sol).abs() * 1e3
    );
    // GMAT validation
    assert!(
        (solution_fd.correction.norm() - gmat_sol).abs() < 1e-6
            || solution_fd.correction.norm() < gmat_sol,
        "Finite differencing result different from GMAT and greater!"
    );
}
