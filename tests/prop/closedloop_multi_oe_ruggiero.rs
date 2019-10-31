extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::hifitime::{Epoch, SECONDS_PER_DAY};
use self::nyx::celestia::{bodies, Cosm, Geoid, State};
use self::nyx::dynamics::celestial::CelestialDynamics;
use self::nyx::dynamics::propulsion::{Propulsion, Thruster};
use self::nyx::dynamics::spacecraft::Spacecraft;
use self::nyx::dynamics::thrustctrl::{Achieve, Ruggiero};
use self::nyx::dynamics::Dynamics;
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};

/// NOTE: Herein shows the difference between the QLaw and Ruggiero (and other control laws).
/// The Ruggiero control law takes quite some longer to converge than the QLaw.

#[test]
fn qlaw_as_ruggiero_case_a() {
    // Source: AAS-2004-5089

    let mut cosm = Cosm::from_xb("./de438s");
    cosm.mut_gm_for_geoid_id(bodies::EARTH, 398_600.433);
    let earth = cosm.geoid_from_id(bodies::EARTH);

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = earth.equatorial_radius + 1000.0;

    let orbit = State::<Geoid>::from_keplerian(sma, 0.01, 0.05, 0.0, 0.0, 1.0, start_time, earth);

    let prop_time = 39.91 * SECONDS_PER_DAY;

    // Define the dynamics
    let mut dynamics = CelestialDynamics::two_body(orbit);

    // Define the thruster
    let lowt = vec![Thruster {
        thrust: 1.0,
        isp: 3100.0,
    }];

    // Define the objectives
    let objectives = vec![
        Achieve::Sma {
            target: 42164.0,
            tol: 1.0,
        },
        Achieve::Ecc {
            target: 0.01,
            tol: 5e-5,
        },
    ];

    let mut ruggiero = Ruggiero::new(objectives, orbit);

    let dry_mass = 1.0;
    let fuel_mass = 299.0;

    let mut prop_subsys = Propulsion::new(&mut ruggiero, fuel_mass, orbit.dt, lowt, true);

    let mut sc = Spacecraft::with_prop(&mut dynamics, &mut prop_subsys, dry_mass);
    println!("{:o}", orbit);

    let mut prop = Propagator::new::<RK4Fixed>(&mut sc, &PropOpts::with_fixed_step(10.0));
    prop.until_time_elapsed(prop_time);
    let final_state = prop.dynamics.celestial.state();
    let fuel_usage = fuel_mass - sc.prop.unwrap().fuel_mass;
    println!("{:o}", final_state);
    println!("fuel usage: {:.3} kg", fuel_usage);

    assert!(ruggiero.achieved(&final_state), "objective not achieved");

    assert!((fuel_usage - 113.0).abs() < 1.0);
}

#[test]
fn qlaw_as_ruggiero_case_b() {
    // Source: AAS-2004-5089
    let cosm = Cosm::from_xb("./de438s");
    let earth = cosm.geoid_from_id(bodies::EARTH);

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit =
        State::<Geoid>::from_keplerian(24505.9, 0.725, 7.05, 0.0, 0.0, 0.0, start_time, earth);

    let prop_time = 160.0 * SECONDS_PER_DAY;

    // Define the dynamics
    let mut dynamics = CelestialDynamics::two_body(orbit);

    // Define the thruster
    let lowt = vec![Thruster {
        thrust: 0.350,
        isp: 2000.0,
    }];

    // Define the objectives
    let objectives = vec![
        Achieve::Sma {
            target: 42165.0,
            tol: 20.0,
        },
        Achieve::Ecc {
            target: 0.001,
            tol: 5e-5,
        },
        Achieve::Inc {
            target: 0.05,
            tol: 5e-3,
        },
    ];

    let mut ruggiero = Ruggiero::new(objectives, orbit);

    let fuel_mass = 1999.9;
    let dry_mass = 0.1;

    let mut prop_subsys = Propulsion::new(&mut ruggiero, fuel_mass, orbit.dt, lowt, true);

    let mut sc = Spacecraft::with_prop(&mut dynamics, &mut prop_subsys, dry_mass);
    println!("{:o}", orbit);

    let mut prop = Propagator::new::<RK4Fixed>(&mut sc, &PropOpts::with_fixed_step(10.0));
    prop.until_time_elapsed(prop_time);

    let final_state = prop.dynamics.celestial.state();
    let fuel_usage = fuel_mass - sc.prop.unwrap().fuel_mass;
    println!("{:o}", final_state);
    println!("fuel usage: {:.3} kg", fuel_usage);

    assert!(ruggiero.achieved(&final_state), "objective not achieved");

    assert!((fuel_usage - 247.0).abs() < 1.0);
}

#[test]
fn qlaw_as_ruggiero_case_c() {
    // Source: AAS-2004-5089
    let cosm = Cosm::from_xb("./de438s");
    let earth = cosm.geoid_from_id(bodies::EARTH);

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit =
        State::<Geoid>::from_keplerian(9222.7, 0.2, 0.573, 0.0, 0.0, 0.0, start_time, earth);

    let prop_time = 3.0 * SECONDS_PER_DAY;

    // Define the dynamics
    let mut dynamics = CelestialDynamics::two_body(orbit);

    // Define the thruster
    let lowt = vec![Thruster {
        thrust: 9.3,
        isp: 3100.0,
    }];

    // Define the objectives
    let objectives = vec![
        Achieve::Sma {
            target: 30000.0,
            tol: 1.0,
        },
        Achieve::Ecc {
            target: 0.7,
            tol: 5e-5,
        },
    ];

    let mut ruggiero = Ruggiero::new(objectives, orbit);

    let fuel_mass = 299.9;
    let dry_mass = 0.1;

    let mut prop_subsys = Propulsion::new(&mut ruggiero, fuel_mass, orbit.dt, lowt, true);

    let mut sc = Spacecraft::with_prop(&mut dynamics, &mut prop_subsys, dry_mass);
    println!("{:o}", orbit);

    let mut prop = Propagator::new::<RK4Fixed>(&mut sc, &PropOpts::with_fixed_step(10.0));
    prop.until_time_elapsed(prop_time);

    let final_state = prop.dynamics.celestial.state();
    let fuel_usage = fuel_mass - sc.prop.unwrap().fuel_mass;
    println!("{:o}", final_state);
    println!("fuel usage: {:.3} kg", fuel_usage);

    assert!(ruggiero.achieved(&final_state), "objective not achieved");

    assert!((fuel_usage - 64.0).abs() < 1.0);
}

#[test]
#[ignore]
fn qlaw_as_ruggiero_case_d() {
    // BUG due to broken RAAN correction: https://gitlab.com/chrisrabotin/nyx/issues/83
    // Source: AAS-2004-5089
    let cosm = Cosm::from_xb("./de438s");
    let earth = cosm.geoid_from_id(bodies::EARTH);

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit =
        State::<Geoid>::from_keplerian(24505.9, 0.725, 0.06, 0.0, 0.0, 0.0, start_time, earth);

    let prop_time = 113.0 * SECONDS_PER_DAY;

    // Define the dynamics
    let mut dynamics = CelestialDynamics::two_body(orbit);

    // Define the thruster
    let lowt = vec![Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    }];

    // Define the objectives
    let objectives = vec![
        Achieve::Sma {
            target: 26500.0,
            tol: 1.0,
        },
        Achieve::Inc {
            target: 116.0,
            tol: 5e-3,
        },
        Achieve::Ecc {
            target: 0.7,
            tol: 5e-5,
        },
        Achieve::Raan {
            target: 360.0 - 90.0,
            tol: 5e-3,
        },
    ];

    let mut ruggiero = Ruggiero::new(objectives, orbit);

    let fuel_mass = 67.0;
    let dry_mass = 300.0;

    let mut prop_subsys = Propulsion::new(&mut ruggiero, fuel_mass, orbit.dt, lowt, true);

    let mut sc = Spacecraft::with_prop(&mut dynamics, &mut prop_subsys, dry_mass);
    println!("{:o}", orbit);

    let mut prop = Propagator::new::<RK4Fixed>(&mut sc, &PropOpts::with_fixed_step(10.0));
    prop.until_time_elapsed(prop_time);

    let final_state = prop.dynamics.celestial.state();
    let fuel_usage = fuel_mass - sc.prop.unwrap().fuel_mass;
    println!("{:o}", final_state);
    println!("fuel usage: {:.3} kg", fuel_usage);

    assert!(ruggiero.achieved(&final_state), "objective not achieved");

    assert!((fuel_usage - 23.0).abs() < 1.0);
}

#[test]
#[ignore]
fn qlaw_as_ruggiero_case_e() {
    // BUG due to broken RAAN correction: https://gitlab.com/chrisrabotin/nyx/issues/83
    // Source: AAS-2004-5089
    let cosm = Cosm::from_xb("./de438s");
    let earth = cosm.geoid_from_id(bodies::EARTH);

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit =
        State::<Geoid>::from_keplerian(24505.9, 0.725, 0.06, 0.0, 0.0, 0.0, start_time, earth);

    let prop_time = 400.0 * SECONDS_PER_DAY;

    // Define the dynamics
    let mut dynamics = CelestialDynamics::two_body(orbit);

    // Define the thruster
    let lowt = vec![Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    }];

    // Define the objectives
    let objectives = vec![
        Achieve::Sma {
            target: 26500.0,
            tol: 1.0,
        },
        Achieve::Ecc {
            target: 0.7,
            tol: 5e-5,
        },
        Achieve::Inc {
            target: 116.0,
            tol: 5e-3,
        },
        Achieve::Raan {
            target: 270.0,
            tol: 5e-3,
        },
        Achieve::Aop {
            target: 180.0,
            tol: 5e-3,
        },
    ];

    let mut ruggiero = Ruggiero::new(objectives, orbit);

    let fuel_mass = 1999.9;
    let dry_mass = 0.1;

    let mut prop_subsys = Propulsion::new(&mut ruggiero, fuel_mass, orbit.dt, lowt, true);

    let mut sc = Spacecraft::with_prop(&mut dynamics, &mut prop_subsys, dry_mass);
    println!("{:o}", orbit);

    let mut prop = Propagator::new::<RK4Fixed>(&mut sc, &PropOpts::with_fixed_step(10.0));
    prop.until_time_elapsed(prop_time);

    let final_state = prop.dynamics.celestial.state();
    let fuel_usage = fuel_mass - sc.prop.unwrap().fuel_mass;
    println!("{:o}", final_state);
    println!("fuel usage: {:.3} kg", fuel_usage);

    assert!(ruggiero.achieved(&final_state), "objective not achieved");

    assert!((fuel_usage - 23.0).abs() < 1.0);
}

#[test]
fn qlaw_as_ruggiero_case_f() {
    // Source: AAS-2004-5089
    /*
        NOTE: Due to how lifetime of variables work in Rust, we need to define all of the
        components of a spacecraft before defining the spacecraft itself.
    */

    let cosm = Cosm::from_xb("./de438s");
    let earth = cosm.geoid_from_id(bodies::EARTH);

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit =
        State::<Geoid>::from_keplerian(15378.0, 0.01, 98.7, 0.0, 0.0, 0.0, start_time, earth);

    let prop_time = 30.0 * SECONDS_PER_DAY;

    // Define the dynamics
    let mut dynamics = CelestialDynamics::two_body(orbit);

    // Define the thruster
    let lowt = vec![Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    }];

    // Define the objectives
    let objectives = vec![Achieve::Ecc {
        target: 0.15,
        tol: 1e-5,
    }];

    let mut ruggiero = Ruggiero::new(objectives, orbit);

    let fuel_mass = 67.0;
    let dry_mass = 300.0;

    let mut prop_subsys = Propulsion::new(&mut ruggiero, fuel_mass, orbit.dt, lowt, true);

    let mut sc = Spacecraft::with_prop(&mut dynamics, &mut prop_subsys, dry_mass);
    println!("{:o}", orbit);

    let mut prop = Propagator::new::<RK4Fixed>(&mut sc, &PropOpts::with_fixed_step(10.0));
    prop.until_time_elapsed(prop_time);

    let final_state = prop.dynamics.celestial.state();
    let fuel_usage = fuel_mass - sc.prop.unwrap().fuel_mass;
    println!("{:o}", final_state);
    println!("fuel usage: {:.3} kg", fuel_usage);

    assert!(ruggiero.achieved(&final_state), "objective not achieved");

    assert!((fuel_usage - 14.0).abs() < 1.0);
}

#[test]
#[ignore]
fn ruggiero_iepc_2011_102() {
    // Source: IEPC 2011 102
    // WARNING: The paper claims that only 103.9 days are needed. I could not reproduce this, their here nor in another tool.
    let cosm = Cosm::from_xb("./de438s");
    let earth = cosm.geoid_from_id(bodies::EARTH);

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit =
        State::<Geoid>::from_keplerian(24396.0, 0.7283, 7.0, 1.0, 1.0, 1.0, start_time, earth);

    // let prop_time = 104.0 * SECONDS_PER_DAY;
    let prop_time = 140.0 * SECONDS_PER_DAY;
    // let prop_time = 3404.0 * 3600.0 + 21. * 60.;

    // Define the dynamics
    let mut dynamics = CelestialDynamics::two_body(orbit);

    // Define the thruster
    let lowt = vec![Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    }];

    // Define the objectives
    let objectives = vec![
        Achieve::Sma {
            target: 42164.0,
            tol: 20.0,
        },
        Achieve::Inc {
            target: 0.001,
            tol: 5e-3,
        },
        Achieve::Ecc {
            target: 0.001,
            tol: 5e-5,
        },
    ];

    let mut ruggiero = Ruggiero::new(objectives, orbit);

    let fuel_mass = 67.0;
    let dry_mass = 300.0;

    let mut prop_subsys = Propulsion::new(&mut ruggiero, fuel_mass, orbit.dt, lowt, true);

    let mut sc = Spacecraft::with_prop(&mut dynamics, &mut prop_subsys, dry_mass);
    println!("{:o}", orbit);

    let mut prop = Propagator::new::<RK4Fixed>(&mut sc, &PropOpts::with_fixed_step(10.0));
    prop.until_time_elapsed(prop_time);

    let final_state = prop.dynamics.celestial.state();
    let fuel_usage = fuel_mass - sc.prop.unwrap().fuel_mass;
    println!("{:o}", final_state);
    println!("fuel usage: {:.3} kg", fuel_usage);

    assert!(ruggiero.achieved(&final_state), "objective not achieved");

    // WARNING: Paper claims this can be done with only 49kg of fuel.
    assert!((fuel_usage - 49.0).abs() < 1.0);
}
