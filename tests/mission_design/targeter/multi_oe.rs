extern crate nyx_space as nyx;

use nyx::md::targeter::*;
use nyx::md::ui::*;

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
        Objective::within_tolerance(StateParameter::C3, -2.0, 0.5),
        Objective::within_tolerance(StateParameter::RightAscension, 1.0, 0.1),
        Objective::within_tolerance(StateParameter::Declination, 2.0, 0.1),
    ];

    let tgt = Targeter::delta_v(Arc::new(&setup), objectives);
    println!("{}", tgt);

    let solution = tgt
        .try_achieve_from(spacecraft, orig_dt, orig_dt + 4 * TimeUnit::Day)
        .unwrap();

    println!("{}", solution);
}
