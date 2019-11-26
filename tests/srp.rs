extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use hifitime::{Epoch, SECONDS_PER_DAY};
use na::Vector6;
use nyx::celestia::{bodies, Cosm, Geoid, State};
use nyx::dynamics::celestial::CelestialDynamics;
use nyx::dynamics::solarpressure::SolarPressure;
use nyx::dynamics::spacecraft::Spacecraft;
use nyx::dynamics::thrustctrl::NoThrustControl;
use nyx::dynamics::Dynamics;
use nyx::propagators::{PropOpts, Propagator, RK89};
use nyx::utils::rss_state_errors;

#[test]
fn srp_earth() {
    let mut cosm = Cosm::from_xb("./de438s");
    cosm.mut_gm_for_geoid_id(bodies::EARTH, 398_600.441_5);
    let earth = cosm.geoid_from_id(bodies::EARTH);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = State::<Geoid>::from_keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, earth);

    let prop_time = 24.0 * SECONDS_PER_DAY;

    // Define the dynamics
    let mut dynamics = CelestialDynamics::two_body(orbit);

    let shadow_bodies = vec![earth];

    let mut srp = SolarPressure::default(1.0, shadow_bodies, &cosm);

    let dry_mass = 300.0;

    let mut sc = Spacecraft::<NoThrustControl>::with_srp(dynamics, srp, dry_mass);
    println!("{:o}", orbit);

    let mut prop = Propagator::new::<RK89>(&mut sc, &PropOpts::default());
    prop.until_time_elapsed(prop_time);

    let final_state = prop.dynamics.state();
    println!("{}", final_state);

    // GMAT result
    let rslt = Vector6::new(
        -10256.01279521848,
        -22135.87506832323,
        0.0004868617601061399,
        3.667559854760713,
        -1.699197984179455,
        -5.729183400502337e-08,
    );

    let (err_r, err_v) = rss_state_errors(&final_state.orbit.to_cartesian_vec(), &rslt);
    println!("{:e}\t{:e}", err_r, err_v);
    // Note that we have quite large SRP differences with GMAT compared to the other errors.
    // Cf. VALIDATION.MD for details.
    assert!(err_r < 5e-4, "position error too large for SRP");
    assert!(err_v < 1e-7, "velocity error too large for SRP");
}
