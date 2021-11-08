extern crate nyx_space as nyx;

use nyx::md::targeter::*;
use nyx::md::ui::*;

// Semi major axis

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
    let xf_desired_sma = 8_100.0;

    // Define the objective
    let objectives = vec![Objective::new(StateParameter::SMA, xf_desired_sma)];

    let tgt = Targeter::delta_v(&setup, objectives);

    println!("{}", tgt);

    let solution_fd = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("Finite differencing solution: {}", solution_fd);

    let gmat_sol = 0.05312024615278713;
    // GMAT validation
    assert!(
        dbg!(solution_fd.correction.norm() - gmat_sol).abs() < 1e-6,
        "Finite differencing result different from GMAT (greater than 1 mm/s)."
    );

    // Check that the solutions nearly match
    println!(
        "GMAT validation - tgt_sma_from_apo: Δv = {:.3} m/s\terr = {:.6} m/s",
        solution_fd.correction.norm() * 1e3,
        (solution_fd.correction.norm() - gmat_sol).abs() * 1e3
    );
}

#[test]
fn tgt_sma_from_peri_fd() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k);

    let target_delta_t: Duration = xi_orig.period() / 20.0;

    println!("Period: {} s", xi_orig.period().in_seconds() / 2.0);

    let spacecraft = Spacecraft::from_srp_defaults(xi_orig, 100.0, 0.0);

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter],
        cosm,
    ));
    let setup = Propagator::default(dynamics);

    // Try to increase SMA
    let xf_desired_sma = 8_100.0;

    // Define the objective
    let objectives = vec![Objective::new(StateParameter::SMA, xf_desired_sma)];

    let tgt = Targeter::delta_v(&setup, objectives);

    println!("{}", tgt);

    let solution_fd = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("Finite differencing solution: {}", solution_fd);

    let gmat_sol = 0.03550369448069638;
    println!(
        "GMAT validation - tgt_sma_from_peri: Δv = {:.3} m/s\terr = {:.6} m/s",
        solution_fd.correction.norm() * 1e3,
        (solution_fd.correction.norm() - gmat_sol).abs() * 1e3
    );
    // GMAT validation
    assert!(
        (solution_fd.correction.norm() - gmat_sol).abs() < 1e-6,
        "Finite differencing result different from GMAT (greater than 1 mm/s)."
    );
}

#[test]
fn tgt_hd_sma_from_peri() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k).with_stm();

    let target_delta_t: Duration = xi_orig.period() / 40.0;

    println!("Period: {} s", xi_orig.period().in_seconds() / 2.0);

    let spacecraft = Spacecraft::new(xi_orig, 100.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter],
        cosm,
    ));
    let setup = Propagator::default(dynamics);

    // Try to increase SMA
    let xf_desired_sma = 8_100.0;

    // Define the objective
    let objectives = vec![Objective::new(StateParameter::SMA, xf_desired_sma)];

    let mut tgt = Targeter::delta_v(&setup, objectives);
    tgt.iterations = 5;

    println!("{}", tgt);

    let solution_fd = tgt
        .try_achieve_dual(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("Finite differencing solution: {}", solution_fd);

    let gmat_sol = 0.03550369448069638;
    println!(
        "GMAT validation - tgt_sma_from_peri: Δv = {:.3} m/s\terr = {:.6} m/s",
        solution_fd.correction.norm() * 1e3,
        (solution_fd.correction.norm() - gmat_sol).abs() * 1e3
    );
    // GMAT validation
    assert!(
        (solution_fd.correction.norm() - gmat_sol).abs() < 1e-6,
        "Finite differencing result different from GMAT (greater than 1 mm/s)."
    );
}

#[test]
fn orbit_stm_chk() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k);

    // let target_delta_t: Duration = xi_orig.period() / 2.0;
    let target_delta_t = 100.0 * TimeUnit::Second;

    println!("Period: {} s", xi_orig.period().in_seconds() / 2.0);

    // let spacecraft = Spacecraft::from_srp_defaults(xi_orig, 100.0, 0.0);

    let dynamics = OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter],
        cosm,
    );
    let setup = Propagator::default(dynamics);
    let mut prop_instance = setup.with(xi_orig.with_stm());

    let achievement_epoch = orig_dt + target_delta_t;

    loop {
        let prev_state = prop_instance.state.to_cartesian_vec();
        prop_instance.single_step().unwrap();
        if prop_instance.state.epoch() > achievement_epoch {
            // Go backward if we've done too far
            prop_instance.until_epoch(achievement_epoch).unwrap();
        }
        let stm_k_kp1 = prop_instance.state.stm().unwrap();
        println!(
            "{}=>err = {}",
            stm_k_kp1,
            stm_k_kp1 * prev_state - prop_instance.state.to_cartesian_vec()
        );
        // traj_stm *= stm_k_kp1;
        // And reset the STM
        prop_instance.state.reset_stm();
        if prop_instance.state.epoch() == achievement_epoch {
            break;
        }
    }
}

// Eccentricity
#[test]
fn tgt_ecc_from_apo() {
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

    let xf_desired_ecc = 0.4;

    let tgt = Targeter::new(
        &setup,
        vec![
            Variable {
                component: Vary::VelocityX,
                max_step: 5.0,
                ..Default::default()
            },
            Variable {
                component: Vary::VelocityY,
                max_step: 5.0,
                ..Default::default()
            },
            Variable {
                component: Vary::VelocityZ,
                max_step: 5.0,
                ..Default::default()
            },
        ],
        vec![Objective::new(StateParameter::Eccentricity, xf_desired_ecc)],
    );

    println!("{}", tgt);

    let solution_fd = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("Finite differencing solution: {}", solution_fd);

    let gmat_sol = 0.7721483022815125;
    println!(
        "GMAT validation - tgt_ecc_from_apo: Δv = {:.3} m/s\terr = {:.6} m/s",
        solution_fd.correction.norm() * 1e3,
        (solution_fd.correction.norm() - gmat_sol).abs() * 1e3
    );
    // GMAT validation
    assert!(
        (solution_fd.correction.norm() - gmat_sol).abs() < 1e-3,
        "Finite differencing result different from GMAT (greater than 1 m/s)."
    );
}

#[test]
fn tgt_ecc_from_peri() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k);

    let target_delta_t: Duration = xi_orig.period() / 2.0;

    println!("Period: {} s", xi_orig.period().in_seconds() / 2.0);

    let spacecraft = Spacecraft::from_srp_defaults(xi_orig, 100.0, 0.0);

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter],
        cosm,
    ));
    let setup = Propagator::default(dynamics);

    let xf_desired_ecc = 0.4;

    let tgt = Targeter::new(
        &setup,
        vec![
            Variable {
                component: Vary::VelocityX,
                max_step: 5.0,
                ..Default::default()
            },
            Variable {
                component: Vary::VelocityY,
                max_step: 5.0,
                ..Default::default()
            },
            Variable {
                component: Vary::VelocityZ,
                max_step: 5.0,
                ..Default::default()
            },
        ],
        vec![Objective::new(StateParameter::Eccentricity, xf_desired_ecc)],
    );

    println!("{}", tgt);

    let solution_fd = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("Finite differencing solution: {}", solution_fd);

    let gmat_sol = 0.6926746704643234;
    println!(
        "GMAT validation - tgt_ecc_from_peri: Δv = {:.3} m/s\terr = {:.6} m/s",
        solution_fd.correction.norm() * 1e3,
        (solution_fd.correction.norm() - gmat_sol).abs() * 1e3
    );
    // GMAT validation
    assert!(
        (solution_fd.correction.norm() - gmat_sol).abs() < 1e-3,
        "Finite differencing result different from GMAT (greater than 1 m/s)."
    );
}

// RAAN
#[test]
fn tgt_raan_from_apo() {
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

    let xf_desired_raan = 65.0;

    // Define the objective
    let objectives = vec![Objective::new(StateParameter::RAAN, xf_desired_raan)];

    let tgt = Targeter::delta_v(&setup, objectives);

    println!("{}", tgt);

    let solution_fd = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("Finite differencing solution: {}", solution_fd);

    let gmat_sol = 0.30344716711198855;
    println!(
        "GMAT validation - tgt_raan_from_apo: Δv = {:.3} m/s\terr = {:.6} m/s",
        solution_fd.correction.norm() * 1e3,
        (solution_fd.correction.norm() - gmat_sol).abs() * 1e3
    );
    // GMAT validation
    assert!(
        (solution_fd.correction.norm() - gmat_sol).abs() < 1e-3,
        "Finite differencing result different from GMAT (greater than 1 m/s)."
    );
}

#[test]
fn tgt_raan_from_peri() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k);

    let target_delta_t: Duration = xi_orig.period() / 2.0;

    println!("Period: {} s", xi_orig.period().in_seconds() / 2.0);

    let spacecraft = Spacecraft::from_srp_defaults(xi_orig, 100.0, 0.0);

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::default(dynamics);

    let xf_desired_raan = 65.0;

    // Define the objective
    let objectives = vec![Objective::new(StateParameter::RAAN, xf_desired_raan)];

    let tgt = Targeter::new(
        &setup,
        vec![
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
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("Finite differencing solution: {}", solution_fd);

    let gmat_sol = 0.45110541873478793;
    println!(
        "GMAT validation - tgt_raan_from_peri: Δv = {:.3} m/s\terr = {:.6} m/s",
        solution_fd.correction.norm() * 1e3,
        (solution_fd.correction.norm() - gmat_sol).abs() * 1e3
    );
    // GMAT validation
    assert!(
        (solution_fd.correction.norm() - gmat_sol).abs() < 6e-3,
        "Finite differencing result different from GMAT (greater than 6 m/s)."
    );
}

// AoP
#[test]
fn tgt_aop_from_apo() {
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

    let xf_desired_aop = 65.0;

    // Define the objective
    let objectives = vec![Objective::new(StateParameter::AoP, xf_desired_aop)];

    let tgt = Targeter::delta_v(&setup, objectives);

    println!("{}", tgt);

    let solution_fd = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("Finite differencing solution: {}", solution_fd);

    let gmat_sol = 0.11772316331182386;
    println!(
        "GMAT validation - tgt_aop_from_apo: Δv = {:.3} m/s\terr = {:.6} m/s",
        solution_fd.correction.norm() * 1e3,
        (solution_fd.correction.norm() - gmat_sol).abs() * 1e3
    );
    // GMAT validation
    assert!(
        (solution_fd.correction.norm() - gmat_sol).abs() < 1e-3,
        "Finite differencing result different from GMAT (greater than 1 m/s)."
    );
}

#[test]
fn tgt_aop_from_peri() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 0.0, orig_dt, eme2k);

    let target_delta_t: Duration = xi_orig.period() / 2.0;

    println!("Period: {} s", xi_orig.period().in_seconds() / 2.0);

    let spacecraft = Spacecraft::from_srp_defaults(xi_orig, 100.0, 0.0);

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::default(dynamics);

    let xf_desired_aop = 65.0;

    // Define the objective
    let objectives = vec![Objective::new(StateParameter::AoP, xf_desired_aop)];

    let tgt = Targeter::delta_v(&setup, objectives);

    println!("{}", tgt);

    let solution_fd = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("Finite differencing solution: {}", solution_fd);

    let gmat_sol = 0.12197875695918228;
    println!(
        "GMAT validation - tgt_aop_from_peri: Δv = {:.3} m/s\terr = {:.6} m/s",
        solution_fd.correction.norm() * 1e3,
        (solution_fd.correction.norm() - gmat_sol).abs() * 1e3
    );
    // GMAT validation
    assert!(
        (solution_fd.correction.norm() - gmat_sol).abs() < 6e-3,
        "Finite differencing result different from GMAT (greater than 6 m/s)."
    );
}
