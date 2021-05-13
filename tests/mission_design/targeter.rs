extern crate nyx_space as nyx;

use nyx::md::targeter::*;
use nyx::md::ui::*;

#[test]
fn tgt_basic_position() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 180.0, orig_dt, eme2k);

    let target_delta_t: Duration = xi_orig.period() / 2.0;

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

    // Define the objectives
    let objectives = vec![
        Objective::within_tolerance(StateParameter::X, xf_desired.x, 0.1),
        Objective::within_tolerance(StateParameter::Y, xf_desired.y, 0.1),
        Objective::within_tolerance(StateParameter::Z, xf_desired.z, 0.1),
    ];

    // let tgt = Targeter::delta_v(Arc::new(&setup), objectives);
    let tgt = Targeter::new(
        Arc::new(&setup),
        vec![Vary::VelocityV, Vary::VelocityC],
        objectives,
    );

    println!("{}", tgt);

    let solution = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("{}", solution);
}

#[test]
fn tgt_basic_sma() {
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

    let solution = tgt
        .try_achieve_from_with_guess_fd(
            spacecraft,
            &[0.0, 0.0, 0.0],
            orig_dt,
            orig_dt + target_delta_t,
        )
        .unwrap();

    println!("{}", solution);
}

#[test]
fn tgt_position_sma() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 180.0, orig_dt, eme2k);

    let target_delta_t: Duration = xi_orig.period() / 2.0;

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

    // Define the objectives
    let objectives = vec![Objective::within_tolerance(
        StateParameter::SMA,
        xf_desired.sma(),
        0.1,
    )];

    let tgt = Targeter::delta_v(Arc::new(&setup), objectives);

    println!("{}", tgt);

    let solution = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("{}", solution);

    tgt.apply(solution).unwrap();

    // As expected, the further out we are, the better the less delta-V is needed to match a B-Plane
    // assert!((solution.correction.norm() - 7563.095e-3).abs() < 1e-6);
}

#[test]
fn tgt_c3_ra_decl_velocity() {
    // TODO: Reubild this in GMAT see if it works there
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");
    let luna = cosm.frame("luna");

    let orig_dt = Epoch::from_gregorian_utc_at_midnight(2020, 1, 1);

    let xi_orig = Orbit::keplerian(8_000.0, 0.2, 30.0, 60.0, 60.0, 180.0, orig_dt, eme2k);
    let xi_moon = cosm.frame_chg(&xi_orig, luna);

    let spacecraft = Spacecraft::from_srp_defaults(xi_moon, 100.0, 0.0);

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::default(dynamics);

    // Define the objective
    let objectives = vec![
        // Objective::within_tolerance(StateParameter::C3, -2.0, 0.5),
        Objective::within_tolerance(StateParameter::RightAscension, 1.0, 0.1),
        Objective::within_tolerance(StateParameter::Declination, 2.0, 0.1),
    ];

    let tgt = Targeter::delta_r(Arc::new(&setup), objectives);
    println!("{}", tgt);

    let solution = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + 4 * TimeUnit::Day)
        .unwrap();

    println!("{}", solution);
}

#[test]
fn tgt_b_plane_sanity() {
    // Rebuild the "in-place" targeting from the B-Plane test of `try_achieve`

    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    // This is a reproduction of the B-plane computation from the `Ex_LunarTransfer.script` file from GMAT
    let cosm = Cosm::de438_gmat();
    // Define the epoch
    let epoch = Epoch::from_gregorian_utc_at_midnight(2016, 1, 1);

    // Hyperbolic orbit
    let orbit = Orbit::cartesian(
        546507.344255845,
        -527978.380486028,
        531109.066836708,
        -4.9220589268733,
        5.36316523097915,
        -5.22166308425181,
        epoch,
        cosm.frame("EME2000"),
    );

    let prop = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter],
        cosm,
    )));

    let spacecraft = Spacecraft::from_srp_defaults(orbit, 100.0, 0.0);

    let b_plane_tgt = BPlaneTarget::from_bt_br(13135.7982982557, 5022.26511510685);

    let tgt = Targeter::delta_v(Arc::new(&prop), b_plane_tgt.to_objectives());

    let sol = tgt.try_achieve_from(spacecraft, epoch, epoch).unwrap();

    println!("{}", sol);

    // Note that we allow for slightly larger error than the other in-place correction
    assert!((sol.correction[0] - -0.25386251697606466).abs() < 1e-6);
    assert!((sol.correction[1] - -0.18774460089778605).abs() < 1e-6);
    assert!((sol.correction[2] - 0.046145009839345504).abs() < 1e-6);

    tgt.apply(sol).unwrap();
}

#[test]
fn tgt_b_plane_legit() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    // This is a reproduction of the B-plane computation from the `Ex_LunarTransfer.script` file from GMAT
    let cosm = Cosm::de438_gmat();
    // Grab the frame
    let eme2k = cosm.frame("EME2000");
    let luna = cosm.frame("Luna");
    // Define the epoch
    let epoch = Epoch::from_gregorian_utc(2014, 7, 22, 11, 29, 10, 811_000);
    // Define the initial orbit in EME2k but convert it to Moon J2000
    let orbit = Orbit::cartesian(
        -137380.1984338506,
        75679.87867537055,
        21487.63875187856,
        -0.2324532014235503,
        -0.4462753967758019,
        0.08561205662877103,
        epoch,
        eme2k,
    );

    let prop = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter],
        cosm.clone(),
    )));

    // let loi_epoch = Epoch::from_str("2014-07-28 22:08:02.448000000 TAI").unwrap();
    // let loi_epoch = epoch + 556697 * TimeUnit::Second;

    let spacecraft = Spacecraft::from_srp_defaults(orbit, 1000.0, 1.0);
    println!("{}", spacecraft);

    // Propagate to periapsis
    let periapse_spacecraft = prop
        .with(spacecraft)
        .until_event(1 * orbit.period(), &Event::periapsis(), 1)
        .unwrap()
        .0;

    // Convert to the Moon J2000 frame
    // let orbit_moon = cosm.frame_chg(&periapse_spacecraft.orbit, luna);
    // let periapse_spacecraft_moon = spacecraft.with_orbit(orbit_moon);

    // println!(
    //     "{}\n{}",
    //     periapse_spacecraft.orbit, periapse_spacecraft_moon.orbit
    // );

    // let b_plane_tgt = BPlaneTarget::from_bt_br(391732.3347895856, 104579.9942274809);
    let b_plane_tgt = BPlaneTarget::from_bt_br(15_000.4, 4_000.6);

    let tgt = Targeter::delta_v_in_frame(
        Arc::new(&prop),
        b_plane_tgt.to_objectives_with_tolerance(3.0),
        luna,
        cosm.clone(),
    );

    let tcm_epoch = periapse_spacecraft.epoch();
    let loi_epoch = tcm_epoch + 556697 * TimeUnit::Second;

    let sol = tgt
        .try_achieve_from(periapse_spacecraft, tcm_epoch, loi_epoch)
        .unwrap();

    println!("{}", sol);

    assert!((sol.correction.norm() - 43.197e-3).abs() < 1e-6);

    tgt.apply(sol).unwrap();
}

#[test]
fn tgt_b_plane_with_propagation() {
    // Rebuild the "in-place" targeting from the B-Plane test of `try_achieve`
    // But perform a backward propagation to make sure that applying the linearization works.

    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    // This is a reproduction of the B-plane computation from the `Ex_LunarTransfer.script` file from GMAT
    let cosm = Cosm::de438_gmat();
    // Define the epoch
    let epoch = Epoch::from_gregorian_utc_at_midnight(2016, 1, 1);

    // Hyperbolic orbit
    let orbit = Orbit::cartesian(
        546507.344255845,
        -527978.380486028,
        531109.066836708,
        -4.9220589268733,
        5.36316523097915,
        -5.22166308425181,
        epoch,
        cosm.frame("EME2000"),
    );

    let prop = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter],
        cosm,
    )));

    let spacecraft = Spacecraft::from_srp_defaults(orbit, 100.0, 0.0);

    let prior_sc = prop
        .with(spacecraft)
        .for_duration(-12 * TimeUnit::Hour)
        .unwrap();

    let b_plane_tgt = BPlaneTarget::from_bt_br(13135.7982982557, 5022.26511510685);

    let tgt = Targeter::delta_v(Arc::new(&prop), b_plane_tgt.to_objectives());

    let sol = tgt
        .try_achieve_from(prior_sc, prior_sc.epoch(), epoch)
        .unwrap();

    println!("{}", sol);

    // As expected, the further out we are, the better the less delta-V is needed to match a B-Plane
    assert!((sol.correction.norm() - 225.379e-3).abs() < 1e-6);

    tgt.apply(sol).unwrap();
}

#[test]
fn tgt_b_plane_remove_this() {
    // Rebuild the "in-place" targeting from the B-Plane test of `try_achieve`
    // But perform a backward propagation to make sure that applying the linearization works.

    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    // This is a reproduction of the B-plane computation from the `Ex_LunarTransfer.script` file from GMAT
    let cosm = Cosm::de438_gmat();
    // Define the epoch
    let epoch = Epoch::from_gregorian_utc_at_midnight(2022, 11, 29);

    // Hyperbolic orbit
    let orbit = Orbit::cartesian(
        4395.725834,
        -8831.846344,
        -5422.661606,
        7.919679,
        -1.783247,
        -1.689868,
        epoch,
        cosm.frame("EME2000"),
    );

    let luna = cosm.frame("Luna");
    let orbit_moon = cosm.frame_chg(&orbit, luna);

    let prop = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter],
        cosm,
    )));

    let spacecraft = Spacecraft::from_srp_defaults(orbit_moon, 100.0, 0.0);

    let b_plane_tgt = BPlaneTarget::from_bt_br(306.550415207, -5639.9403447);

    let tgt = Targeter::delta_v(Arc::new(&prop), b_plane_tgt.to_objectives());

    let tcm1_epoch = epoch + 128349.0 * TimeUnit::Second;
    let loi_epoch = tcm1_epoch + 322559.0 * TimeUnit::Second;

    let sol = tgt
        .try_achieve_from(spacecraft, tcm1_epoch, loi_epoch)
        .unwrap();

    println!("{}", sol);

    // As expected, the further out we are, the better the less delta-V is needed to match a B-Plane
    // assert!((sol.correction.norm() - 225.387e-3).abs() < 1e-3);

    tgt.apply(sol).unwrap();
}
