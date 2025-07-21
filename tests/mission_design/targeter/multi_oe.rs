extern crate nyx_space as nyx;

use nyx::dynamics::guidance::Thruster;
use nyx::md::prelude::*;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use nyx_space::cosmic::Mass;
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn tgt_c3_decl(almanac: Arc<Almanac>) {
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

    let tgt = Targeter::delta_v(&setup, objectives);

    println!("{tgt}");

    let solution_fd = tgt
        .try_achieve_from(
            spacecraft,
            orig_dt,
            orig_dt + target_delta_t,
            almanac.clone(),
        )
        .unwrap();

    println!("Finite differencing solution: {solution_fd}");

    tgt.apply(&solution_fd, almanac).unwrap();

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
fn conv_tgt_sma_ecc(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k);

    let target_delta_t: Duration = xi_orig.period().unwrap() / 2.0;

    println!("Period: {} s", xi_orig.period().unwrap().to_seconds() / 2.0);

    let spacecraft = Spacecraft {
        orbit: xi_orig,
        mass: Mass {
            dry_mass_kg: 10.0,
            prop_mass_kg: 90.0,
            extra_mass_kg: 0.0,
        },
        thruster: Some(Thruster {
            thrust_N: 500.0,
            isp_s: 300.0,
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

    let tgt = Targeter::new(
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

    println!("{tgt}");

    let achievement_epoch = orig_dt + target_delta_t;

    let solution_fd = tgt
        .try_achieve_from(spacecraft, orig_dt, achievement_epoch, almanac)
        .unwrap();

    println!("Finite differencing solution: {solution_fd}");

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

    /* *** */
    /* Convert to a finite burn and make sure this converges */
    /* *** */

    // Optimizer::convert_impulsive_mnvr(solution_fd.corrected_state, solution_fd.correction, &setup)
    //     .unwrap();

    // let mut setup = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    // setup.set_tolerance(1e-3);
    // let mnvr = Mnvr::impulsive_to_finite(
    //     orig_dt,
    //     solution_fd.correction,
    //     spacecraft,
    //     &setup,
    //     Frame::Inertial,
    // )
    // .unwrap();
    // println!("CONVERGED ON {}", mnvr);

    // // Propagate for the duration of the burn
    // setup.dynamics = setup.dynamics.with_ctrl(Arc::new(mnvr));
    // // Propagate until maneuver start
    // let pre_mnvr = setup
    //     .with(spacecraft)
    //     .until_epoch(mnvr.start - mnvr.duration())
    //     .unwrap();
    // // Use a small step through the maneuver
    // setup.set_max_step(mnvr.end - mnvr.start);
    // let xf = setup
    //     .with(pre_mnvr.with_guidance_mode(GuidanceMode::Custom(0)))
    //     .until_epoch(achievement_epoch)
    //     .unwrap();
    // println!("{:x}", xf);
}

#[rstest]
fn tgt_hd_sma_ecc_cov_test(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k);

    let target_delta_t: Duration = xi_orig.period().unwrap() / 20.0;

    let spacecraft = Spacecraft::from_srp_defaults(xi_orig, 100.0, 0.0).with_stm();

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::default_dp78(dynamics);

    // Define the objective
    let objectives = [
        Objective::within_tolerance(StateParameter::Eccentricity, 0.4, 1e-5),
        Objective::within_tolerance(StateParameter::SMA, 8100.0, 0.1),
    ];

    let tgt = Targeter::new(
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

    println!("{tgt}");

    let solution_fd = tgt
        .try_achieve_dual(spacecraft, orig_dt, orig_dt + target_delta_t, almanac)
        .unwrap();

    println!("Finite differencing solution: {solution_fd}");

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
