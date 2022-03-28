extern crate nyx_space as nyx;

use nyx::md::optimizer::*;
use nyx::md::ui::*;

#[test]
fn tgt_b_plane_earth_gravity_assist_no_propagation() {
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

    let prop = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter],
        cosm,
    )));

    let spacecraft = Spacecraft::from_srp_defaults(orbit, 100.0, 0.0);

    let b_plane_tgt = BPlaneTarget::from_bt_br(13135.7982982557, 5022.26511510685);

    let tgt = Optimizer::delta_v(&prop, b_plane_tgt.to_objectives());

    let sol = tgt.try_achieve_from(spacecraft, epoch, epoch).unwrap();

    println!("{}", sol);

    // This is the exact GMAT data from EarthGA.script
    let gmat_sol = 0.31909814507892165;
    println!(
        "GMAT validation - tgt_b_plane_earth_gravity_assist: Δv = {:.3} m/s\terr = {:.6} m/s",
        sol.correction.norm() * 1e3,
        (sol.correction.norm() - gmat_sol).abs() * 1e3
    );
    // GMAT validation
    assert!(
        (sol.correction.norm() - gmat_sol).abs() < 1e-3,
        "Finite differencing result different from GMAT by over 1 m/s"
    );

    tgt.apply(&sol).unwrap();
}

#[allow(clippy::identity_op)]
#[test]
#[ignore = "https://gitlab.com/nyx-space/nyx/-/issues/212"]
fn tgt_b_plane_lunar_transfer() {
    // WARNING: This test is ignored until https://gitlab.com/nyx-space/nyx/-/issues/212
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

    let prop = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter],
        cosm.clone(),
    )));

    let spacecraft = Spacecraft::from_srp_defaults(orbit, 1000.0, 0.0);
    println!("{}", spacecraft);

    // Propagate to periapsis
    let periapse_spacecraft = prop
        .with(spacecraft)
        .until_nth_event(1 * orbit.period(), &Event::periapsis(), 1)
        .unwrap()
        .0;

    let b_plane_tgt = BPlaneTarget::from_bt_br(15_000.4, 4_000.6);

    // GMAT truth with central differencing: 1.15740867962, -0.576350387399, 0.632247251449
    // GMAT truth with forward differencing: 1.33490412071, -0.5447988683, 1.77697094604 (approach currently in Nyx)

    let tgt = Optimizer::in_frame(
        &prop,
        [
            Variable {
                component: Vary::VelocityX,
                min_value: -3.0,
                max_value: 3.0,
                perturbation: 0.0001,
                max_step: 0.5,
                ..Default::default()
            },
            Variable {
                component: Vary::VelocityY,
                min_value: -3.0,
                max_value: 3.0,
                perturbation: 0.0001,
                max_step: 0.5,
                ..Default::default()
            },
            Variable {
                component: Vary::VelocityZ,
                min_value: -3.0,
                max_value: 3.0,
                perturbation: 0.0001,
                max_step: 0.5,
                ..Default::default()
            },
        ],
        b_plane_tgt.to_objectives_with_tolerance(3.0),
        luna,
        cosm.clone(),
    );

    let tcm_epoch = periapse_spacecraft.epoch();
    let loi_epoch = tcm_epoch + 556697 * Unit::Second;

    let sol = tgt
        .try_achieve_from(periapse_spacecraft, tcm_epoch, loi_epoch)
        .unwrap();

    println!("{}", sol);
    let gmat_sol = 2.2883182823767747;
    // GMAT validation
    assert!(
        (sol.correction.norm() - gmat_sol).abs() < 1e-6,
        "Finite differencing result different from GMAT (greater than 1 mm/s)."
    );

    // Check that the solutions nearly match
    println!(
        "GMAT validation - tgt_b_plane_lunar_transfer: Δv = {:.3} m/s\terr = {:.6} m/s",
        sol.correction.norm() * 1e3,
        (sol.correction.norm() - gmat_sol).abs() * 1e3
    );

    // Check that the solution works with the same dynamics.
    tgt.apply(&sol).unwrap();
}

#[test]
fn tgt_b_plane_earth_gravity_assist_with_propagation() {
    // Rebuild the "tgt_b_plane_earth_gravity_assist" scenario but with a propagation and applying the dv earlier

    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

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

    let prop = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter],
        cosm,
    )));

    let spacecraft = Spacecraft::from_srp_defaults(orbit, 100.0, 0.0);

    let prior_sc = prop
        .with(spacecraft)
        .for_duration(-12 * Unit::Hour)
        .unwrap();

    let b_plane_tgt = BPlaneTarget::from_bt_br(13135.7982982557, 5022.26511510685);

    let tgt = Optimizer::delta_v(&prop, b_plane_tgt.to_objectives());

    let sol = tgt
        .try_achieve_from(prior_sc, prior_sc.epoch(), epoch)
        .unwrap();

    println!("{}", sol);

    // As expected, the further out we are, the better the less delta-V is needed to match a B-Plane
    assert!(sol.correction.norm() <= 225.309e-3);

    tgt.apply(&sol).unwrap();
}
