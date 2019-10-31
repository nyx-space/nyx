extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::hifitime::Epoch;
use self::na::Vector3;
use self::nyx::celestia::{bodies, Cosm, Geoid, State};
use self::nyx::dynamics::celestial::CelestialDynamics;
use self::nyx::dynamics::propulsion::{Propulsion, Thruster};
use self::nyx::dynamics::spacecraft::Spacecraft;
use self::nyx::dynamics::thrustctrl::{FiniteBurns, Mnvr};
use self::nyx::dynamics::Dynamics;
use self::nyx::propagators::error_ctrl::RSSStatePV;
use self::nyx::propagators::{PropOpts, Propagator, RK89};
use self::nyx::utils::rss_state_errors;

#[test]
fn transfer_schedule_no_depl() {
    /*
        NOTE: Due to how lifetime of variables work in Rust, we need to define all of the
        components of a spacecraft before defining the spacecraft itself.
    */

    let mut cosm = Cosm::from_xb("./de438s");
    // Modify GMs to match GMAT's
    cosm.mut_gm_for_geoid_id(bodies::EARTH, 398_600.441_5);
    cosm.mut_gm_for_geoid_id(bodies::EARTH_MOON, 4_902.800_582_147_8);
    cosm.mut_gm_for_geoid_id(bodies::JUPITER_BARYCENTER, 126_712_767.857_80);
    cosm.mut_gm_for_geoid_id(bodies::SUN, 132_712_440_017.99);
    let earth = cosm.geoid_from_id(bodies::EARTH);

    let start_time = Epoch::from_gregorian_tai_at_midnight(2002, 1, 1);

    let orbit = State::<Geoid>::from_cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_time, earth,
    );

    let prop_time = 3000.0;

    let mut end_time = start_time;
    end_time.mut_add_secs(prop_time);

    let rslt = State::<Geoid>::from_cartesian(
        4_172.396_780_515_64f64,
        436.944_560_056_202_8,
        -6_518.328_156_815_674,
        -3.979_610_765_995_537,
        5.540_316_900_333_103,
        -2.207_082_771_390_863,
        end_time,
        earth,
    );

    // Define the dynamics
    let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
    let mut dynamics = CelestialDynamics::new(orbit, bodies, &cosm);

    // Define the thruster
    let biprop = vec![Thruster {
        thrust: 10.0,
        isp: 300.0,
    }];

    // Define the maneuver and its schedule
    let mnvr0 = Mnvr {
        start: Epoch::from_gregorian_tai_at_midnight(2002, 1, 1),
        end: end_time,
        thrust_lvl: 1.0, // Full thrust
        vector: Vector3::new(1.0, 0.0, 0.0),
    };

    let mut schedule = FiniteBurns::from_mnvrs(vec![mnvr0], earth);
    let dry_mass = 1e3;
    let fuel_mass = 756.0;

    let mut prop_subsys =
        Propulsion::new(&mut schedule, dry_mass, fuel_mass, orbit.dt, biprop, false);

    let mut sc = Spacecraft::with_prop(&mut dynamics, &mut prop_subsys);

    let mut prop = Propagator::new::<RK89>(&mut sc, &PropOpts::with_fixed_step(10.0));
    prop.until_time_elapsed(prop_time);

    // Compute the errors
    let (err_r, err_v) = rss_state_errors(
        &prop.dynamics.celestial.state_vector(),
        &rslt.to_cartesian_vec(),
    );
    println!("Absolute errors");
    let delta = prop.dynamics.celestial.state_vector() - rslt.to_cartesian_vec();
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} km\tvel = {:.5e} km/s",
        err_r, err_v,
    );

    assert!(
        err_r < 2e-10,
        format!("finite burn position wrong: {:.5e}", err_r)
    );
    assert!(
        err_v < 1e-13,
        format!("finite burn velocity wrong: {:.5e}", err_v)
    );

    // Ensure that there was no change in fuel mass since tank depletion was off
    assert!(
        (prop.dynamics.prop.as_ref().unwrap().fuel_mass - fuel_mass).abs() < std::f64::EPSILON,
        "incorrect fuel mass"
    );
}

#[test]
fn transfer_schedule_depl() {
    /*
        NOTE: Due to how lifetime of variables work in Rust, we need to define all of the
        components of a spacecraft before defining the spacecraft itself.
    */

    let mut cosm = Cosm::from_xb("./de438s");
    // Modify GMs to match GMAT's
    cosm.mut_gm_for_geoid_id(bodies::EARTH, 398_600.441_5);
    cosm.mut_gm_for_geoid_id(bodies::EARTH_MOON, 4_902.800_582_147_8);
    cosm.mut_gm_for_geoid_id(bodies::JUPITER_BARYCENTER, 126_712_767.857_80);
    cosm.mut_gm_for_geoid_id(bodies::SUN, 132_712_440_017.99);
    let earth = cosm.geoid_from_id(bodies::EARTH);

    let start_time = Epoch::from_gregorian_tai_at_midnight(2002, 1, 1);

    let orbit = State::<Geoid>::from_cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_time, earth,
    );

    let prop_time = 3000.0;

    let mut end_time = start_time;
    end_time.mut_add_secs(prop_time);

    let rslt = State::<Geoid>::from_cartesian(
        4_172.433_936_615_18,
        436.936_159_720_413,
        -6_518.368_821_953_345,
        -3.979_569_721_967_499,
        5.540_321_146_839_762,
        -2.207_146_819_283_441,
        end_time,
        earth,
    );

    // Define the dynamics
    let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
    let mut dynamics = CelestialDynamics::new(orbit, bodies, &cosm);

    // Define the thruster
    let biprop = vec![Thruster {
        thrust: 10.0,
        isp: 300.0,
    }];

    // Define the maneuver and its schedule
    let mnvr0 = Mnvr {
        start: Epoch::from_gregorian_tai_at_midnight(2002, 1, 1),
        end: end_time,
        thrust_lvl: 1.0, // Full thrust
        vector: Vector3::new(1.0, 0.0, 0.0),
    };

    let mut schedule = FiniteBurns::from_mnvrs(vec![mnvr0], earth);
    let dry_mass = 1e3;
    let fuel_mass = 756.0;

    let mut prop_subsys =
        Propulsion::new(&mut schedule, dry_mass, fuel_mass, orbit.dt, biprop, true);

    let mut sc = Spacecraft::with_prop(&mut dynamics, &mut prop_subsys);

    let mut prop = Propagator::new::<RK89>(&mut sc, &PropOpts::with_fixed_step(10.0));
    prop.until_time_elapsed(prop_time);

    // Compute the errors
    let (err_r, err_v) = rss_state_errors(
        &prop.dynamics.celestial.state_vector(),
        &rslt.to_cartesian_vec(),
    );
    println!("Absolute errors");
    let delta = prop.dynamics.celestial.state_vector() - rslt.to_cartesian_vec();
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} km\tvel = {:.5e} km/s",
        err_r, err_v,
    );

    assert!(
        err_r < 2e-10,
        format!("finite burn position wrong: {:.5e}", err_r)
    );
    assert!(
        err_v < 1e-13,
        format!("finite burn velocity wrong: {:.5e}", err_v)
    );
    assert!(
        ((prop.dynamics.prop.as_ref().unwrap().fuel_mass - 745.802_837_870_161).abs()) < 2e-10,
        "incorrect fuel mass"
    );
}
