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
        Objective {
            parameter: StateParameter::X,
            desired_value: xf_desired.x,
            tolerance: 0.1,
        },
        Objective {
            parameter: StateParameter::Y,
            desired_value: xf_desired.y,
            tolerance: 0.1,
        },
        Objective {
            parameter: StateParameter::Z,
            desired_value: xf_desired.z,
            tolerance: 0.1,
        },
    ];

    let tgt = Targeter::delta_v(Arc::new(&setup), objectives);

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

    let solution = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
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
    let objectives = vec![Objective {
        parameter: StateParameter::SMA,
        desired_value: xf_desired.sma(),
        tolerance: 0.1,
    }];

    let tgt = Targeter {
        prop: Arc::new(&setup),
        objectives,
        corrector: Corrector::Position,
        iterations: 50,
    };

    println!("{}", tgt);

    let solution = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + target_delta_t)
        .unwrap();

    println!("{}", solution);
}

#[test]
fn tgt_c3_ra_decl_velocity() {
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
        Objective::new(StateParameter::C3, -2.0),
        Objective::new(StateParameter::RightAscension, 0.0),
        Objective::new(StateParameter::Declination, 0.0),
    ];

    let tgt = Targeter::delta_v(Arc::new(&setup), objectives);

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
    // Grab the frame
    let luna = cosm.frame("Luna");
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
    // Propagate until periapse
    let prop = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter],
        cosm.clone(),
    )));

    let orbit_moon = cosm.frame_chg(&orbit, luna);
    let spacecraft = Spacecraft::from_srp_defaults(orbit_moon, 100.0, 0.0);

    let b_plane_tgt = BPlaneTarget::from_b_plane(5022.26511510685, 13135.7982982557);

    let tgt = Targeter::delta_v(Arc::new(&prop), b_plane_tgt.to_objectives());

    let sol = tgt.try_achieve_from(spacecraft, epoch, epoch).unwrap();

    println!("{}", sol);

    tgt.apply(sol).unwrap();
}

#[test]
fn tgt_b_plane_legit() {
    use std::str::FromStr;
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
    // Propagate until periapse
    let prop = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun, Bodies::JupiterBarycenter],
        cosm.clone(),
    )));

    let loi_epoch = Epoch::from_str("2014-07-22 22:48:32.066000000 TAI").unwrap();

    let orbit_moon = cosm.frame_chg(&orbit, luna);
    let spacecraft = Spacecraft::from_srp_defaults(orbit_moon, 100.0, 0.0);

    let b_plane_tgt = BPlaneTarget::from_b_plane(104579.9942274809, 391732.3347895856);

    let tgt = Targeter::delta_v(Arc::new(&prop), b_plane_tgt.to_objectives());

    let sol = tgt.try_achieve_from(spacecraft, epoch, loi_epoch).unwrap();

    println!("{}", sol);

    tgt.apply(sol).unwrap();
}
