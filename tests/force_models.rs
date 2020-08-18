extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use hifitime::{Epoch, SECONDS_PER_DAY};
use na::Vector6;
use nyx::celestia::{Cosm, State};
use nyx::dynamics::drag::Drag;
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::solarpressure::SolarPressure;
use nyx::dynamics::spacecraft::Spacecraft;
use nyx::propagators::{PropOpts, Propagator};
use nyx::utils::rss_errors;

#[test]
fn srp_earth() {
    let mut cosm = Cosm::de438();
    cosm.mut_gm_for_frame("EME2000", 398_600.441_5);
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = State::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24.0 * SECONDS_PER_DAY;

    // Define the dynamics
    let dynamics = OrbitalDynamics::two_body(orbit);

    let shadow_bodies = vec![eme2k];

    let srp = SolarPressure::default(1.0, shadow_bodies, &cosm);

    let dry_mass = 300.0;

    let mut sc = Spacecraft::new(dynamics, dry_mass);
    // Add the SRP model to the spacecraft
    sc.add_model(Box::new(srp));
    println!("{:o}", orbit);

    let mut prop = Propagator::default(&mut sc, &PropOpts::default());
    prop.until_time_elapsed(prop_time);

    let final_state = prop.state();
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

    let (err_r, err_v) = rss_errors(&final_state.orbit.to_cartesian_vec(), &rslt);
    println!("{:e}\t{:e}", err_r, err_v);
    // Note that we have quite large SRP differences with GMAT compared to the other errors.
    // Cf. VALIDATION.MD for details.
    assert!(err_r < 1.9, "position error too large for SRP");
    assert!(err_v < 3e-4, "velocity error too large for SRP");
}

#[test]
fn exp_drag_earth() {
    let mut cosm = Cosm::de438();
    cosm.mut_gm_for_frame("EME2000", 398_600.441_5);
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = State::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24.0 * SECONDS_PER_DAY;

    // Define the dynamics
    let dynamics = OrbitalDynamics::two_body(orbit);

    let shadow_bodies = vec![eme2k];

    let srp = SolarPressure::default(1.0, shadow_bodies, &cosm);
    let drag = Drag::earth_exp(1.0, 2.0, &cosm);

    let dry_mass = 300.0;

    let mut sc = Spacecraft::new(dynamics, dry_mass);
    // Add the SRP model to the spacecraft
    sc.add_model(Box::new(srp));
    // Add the drag model to the spacecraft
    sc.add_model(Box::new(drag));
    println!("{:o}", orbit);

    let mut prop = Propagator::default(&mut sc, &PropOpts::default());
    prop.until_time_elapsed(prop_time);

    let final_state = prop.state();
    println!("{}", final_state);
    println!("{}", final_state.orbit);
}

#[test]
fn std_atm_drag_earth() {
    let mut cosm = Cosm::de438();
    cosm.mut_gm_for_frame("EME2000", 398_600.441_5);
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = State::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24.0 * SECONDS_PER_DAY;

    // Define the dynamics
    let dynamics = OrbitalDynamics::two_body(orbit);

    let shadow_bodies = vec![eme2k];

    let srp = SolarPressure::default(1.0, shadow_bodies, &cosm);
    let drag = Drag::std_atm1976(1.0, 2.0, &cosm);

    let dry_mass = 300.0;

    let mut sc = Spacecraft::new(dynamics, dry_mass);
    // Add the SRP model to the spacecraft
    sc.add_model(Box::new(srp));
    // Add the drag model to the spacecraft
    sc.add_model(Box::new(drag));
    println!("{:o}", orbit);

    let mut prop = Propagator::default(&mut sc, &PropOpts::default());
    prop.until_time_elapsed(prop_time);

    let final_state = prop.state();
    println!("{}", final_state);
    println!("{}", final_state.orbit);

    /*
    Test: compared with exponential drag model, and the final states are similar:

    exp_drag_earth
    [Earth J2000] 2000-01-25T00:00:00 TAI   sma = 24394.167595 km   ecc = 0.000019  inc = 0.000001 deg      raan = 299.937993 deg   aop = 264.180317 deg    ta = 42.152440 deg      300 kg
    [Earth J2000] 2000-01-25T00:00:00 TAI   position = [-9816.442834, -22331.499423, -0.000477] km  velocity = [3.700562, -1.626744, 0.000000] km/s

    std_atm_drag_earth
    [Earth J2000] 2000-01-25T00:00:00 TAI   sma = 24396.000574 km   ecc = 0.000020  inc = 0.000001 deg      raan = 299.400845 deg   aop = 264.478503 deg    ta = 41.265314 deg      300 kg
    [Earth J2000] 2000-01-25T00:00:00 TAI   position = [-10254.183112, -22135.911958, -0.000484] km velocity = [3.667742, -1.699095, 0.000000] km/s

    */
}

#[test]
fn std_atm_drag_earth_low() {
    let mut cosm = Cosm::de438();
    cosm.mut_gm_for_frame("EME2000", 398_600.441_5);
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = State::keplerian(
        eme2k.equatorial_radius() + 300.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        dt,
        eme2k,
    );

    let prop_time = 24.0 * SECONDS_PER_DAY;

    // Define the dynamics
    let dynamics = OrbitalDynamics::two_body(orbit);

    let shadow_bodies = vec![eme2k];

    let srp = SolarPressure::default(1.0, shadow_bodies, &cosm);
    let drag = Drag::std_atm1976(1.0, 2.0, &cosm);

    let dry_mass = 300.0;

    let mut sc = Spacecraft::new(dynamics, dry_mass);
    // Add the SRP model to the spacecraft
    sc.add_model(Box::new(srp));
    // Add the drag model to the spacecraft
    sc.add_model(Box::new(drag));
    println!("{:o}", orbit);

    let mut prop = Propagator::default(&mut sc, &PropOpts::default());
    prop.until_time_elapsed(prop_time);

    let final_state = prop.state();
    println!("{}", final_state);
    println!("{}", final_state.orbit);

    /*
    Test: compared with exponential drag model, and the final states are similar:

    exp_drag_earth
    [Earth J2000] 2000-01-25T00:00:00 TAI   sma = 24394.167595 km   ecc = 0.000019  inc = 0.000001 deg      raan = 299.937993 deg   aop = 264.180317 deg    ta = 42.152440 deg      300 kg
    [Earth J2000] 2000-01-25T00:00:00 TAI   position = [-9816.442834, -22331.499423, -0.000477] km  velocity = [3.700562, -1.626744, 0.000000] km/s

    std_atm_drag_earth
    [Earth J2000] 2000-01-25T00:00:00 TAI   sma = 24396.000574 km   ecc = 0.000020  inc = 0.000001 deg      raan = 299.400845 deg   aop = 264.478503 deg    ta = 41.265314 deg      300 kg
    [Earth J2000] 2000-01-25T00:00:00 TAI   position = [-10254.183112, -22135.911958, -0.000484] km velocity = [3.667742, -1.699095, 0.000000] km/s

    */
}
