extern crate nyx_space as nyx;

use nyx::md::targeter::*;
use nyx::md::prelude::*;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn tgt_vnc_c3_decl(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k);

    let target_delta_t: Duration = xi_orig.period().unwrap() / 2.0;

    println!("Period: {} s", xi_orig.period().unwrap().to_seconds() / 2.0);

    let spacecraft = Spacecraft::from_srp_defaults(xi_orig, 100.0, 0.0);

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::default(dynamics);

    // Define the objective
    let objectives = [
        Objective::within_tolerance(StateParameter::Declination, 5.0, 0.1),
        Objective::within_tolerance(StateParameter::C3, -5.0, 0.5),
    ];

    let tgt = Targeter::vnc(&setup, objectives);

    println!("{}", tgt);

    let solution_fd = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t, almanac)
        .unwrap();

    println!("Finite differencing solution: {}", solution_fd);

    let gmat_sol = 2.385704523944014;
    println!(
        "GMAT validation - tgt_sma_from_peri: Δv = {:.3} m/s\terr = {:.6} m/s",
        solution_fd.correction.norm() * 1e3,
        (solution_fd.correction.norm() - gmat_sol).abs() * 1e3
    );
    // GMAT validation
    assert!(
        (solution_fd.correction.norm() - gmat_sol).abs() < 6e-3,
        "Finite differencing result different from GMAT (greater than 6 m/s)."
    );
}

#[rstest]
fn tgt_vnc_sma_ecc(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k);

    let target_delta_t: Duration = xi_orig.period().unwrap() / 2.0;

    println!("Period: {} s", xi_orig.period().unwrap().to_seconds() / 2.0);

    let spacecraft = Spacecraft::from_srp_defaults(xi_orig, 100.0, 0.0);

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::default(dynamics);

    // Define the objective
    let objectives = [
        Objective::within_tolerance(StateParameter::Eccentricity, 0.4, 1e-5),
        Objective::within_tolerance(StateParameter::SMA, 8100.0, 0.1),
    ];

    let tgt = Targeter::vnc_with_components(
        &setup,
        [
            Variable {
                component: Vary::VelocityX,
                max_step: 0.5,
                ..Default::default()
            },
            Variable {
                component: Vary::VelocityY,
                max_step: 0.5,
                ..Default::default()
            },
            Variable {
                component: Vary::VelocityZ,
                max_step: 0.5,
                ..Default::default()
            },
        ],
        objectives,
    );

    println!("{}", tgt);

    let solution_fd = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t, almanac)
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
