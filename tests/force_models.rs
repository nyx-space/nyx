extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use hifitime::{Epoch, SECONDS_PER_DAY};
use na::Vector6;
use nyx::celestia::{bodies, Cosm, State};
use nyx::dynamics::celestial::CelestialDynamics;
use nyx::dynamics::drag::ExpEarthDrag;
use nyx::dynamics::solarpressure::SolarPressure;
use nyx::dynamics::spacecraft::Spacecraft;
use nyx::dynamics::thrustctrl::NoThrustControl;
use nyx::dynamics::Dynamics;
use nyx::propagators::{PropOpts, Propagator};
use nyx::utils::rss_state_errors;

#[test]
fn srp_earth() {
    let mut cosm = Cosm::from_xb("./de438s");
    cosm.mut_gm_for_frame_id(bodies::EARTH, 398_600.441_5);
    let earth = cosm.frame_by_id(bodies::EARTH);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = State::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, earth);

    let prop_time = 24.0 * SECONDS_PER_DAY;

    // Define the dynamics
    let dynamics = CelestialDynamics::two_body(orbit);

    let shadow_bodies = vec![earth.exb_id()];

    let srp = SolarPressure::default(1.0, shadow_bodies, &cosm);

    let dry_mass = 300.0;

    let mut sc = Spacecraft::<NoThrustControl>::with_srp(dynamics, srp, dry_mass);
    println!("{:o}", orbit);

    let mut prop = Propagator::default(&mut sc, &PropOpts::default());
    prop.until_time_elapsed(prop_time);

    let final_state = prop.dynamics.state();
    println!("{}", final_state);

    // GMAT result
    let rslt = Vector6::new(
        -10_256.012_795_218_48,
        -22_135.875_068_323_23,
        0.000_486_861_760_106_139_9,
        3.667_559_854_760_713,
        -1.699_197_984_179_455,
        -5.729_183_400_502_337e-8,
    );

    let (err_r, err_v) = rss_state_errors(&final_state.orbit.to_cartesian_vec(), &rslt);
    println!("{:e}\t{:e}", err_r, err_v);
    // Note that we have quite large SRP differences with GMAT compared to the other errors.
    // Cf. VALIDATION.MD for details.
    assert!(err_r < 1.9, "position error too large for SRP");
    assert!(err_v < 3e-4, "velocity error too large for SRP");
}

#[test]
fn drag_earth() {
    let mut cosm = Cosm::from_xb("./de438s");
    cosm.mut_gm_for_frame_id(bodies::EARTH, 398_600.441_5);
    let earth = cosm.frame_by_id(bodies::EARTH);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = State::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, earth);

    let prop_time = 24.0 * SECONDS_PER_DAY;

    // Define the dynamics
    let dynamics = CelestialDynamics::two_body(orbit);

    let shadow_bodies = vec![earth.exb_id()];

    let srp = SolarPressure::default(1.0, shadow_bodies, &cosm);
    let drag = ExpEarthDrag {
        sc_area: 1.0,
        cd: 2.2,
        cosm: &cosm,
    };

    let dry_mass = 300.0;

    let mut sc = Spacecraft::<NoThrustControl>::with_srp(dynamics, srp, dry_mass);
    // Add the drag model to the spacecraft
    sc.exp_drag = Some(drag);
    println!("{:o}", orbit);

    let mut prop = Propagator::default(&mut sc, &PropOpts::default());
    prop.until_time_elapsed(prop_time);

    let final_state = prop.dynamics.state();
    println!("{}", final_state);
    println!("{}", final_state.orbit);
}
