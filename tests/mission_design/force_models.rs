extern crate nyx_space as nyx;

use nyx::celestia::{Cosm, Orbit, SpacecraftState};
use nyx::dimensions::Vector6;
use nyx::dynamics::{Drag, OrbitalDynamics, SolarPressure, Spacecraft};
use nyx::propagators::Propagator;
use nyx::time::{Epoch, TimeUnit};
use nyx::utils::rss_errors;

#[test]
fn srp_earth() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24 * TimeUnit::Day;

    // Define the dynamics

    let shadow_bodies = vec![eme2k];

    let srp = SolarPressure::default(1.0, shadow_bodies, cosm);

    let dry_mass = 300.0;

    // Create a spacecraft with SRP model
    let sc_dyn = Spacecraft::with_model(OrbitalDynamics::two_body(), srp);
    println!("{:o}", orbit);

    let sc = SpacecraftState::from_srp_defaults(orbit, dry_mass, 1.0);

    let setup = Propagator::default(sc_dyn);
    let mut prop = setup.with(sc);
    let final_state = prop.for_duration(prop_time).unwrap();

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
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24 * TimeUnit::Day;

    // Define the dynamics

    let shadow_bodies = vec![eme2k];

    let srp = SolarPressure::default(1.0, shadow_bodies, cosm.clone());
    let drag = Drag::earth_exp(1.0, 2.0, cosm);

    let dry_mass = 300.0;

    // Build a spacecraft with SRP and Drag enabled.
    let sc_dyn = Spacecraft::with_models(OrbitalDynamics::two_body(), vec![srp, drag]);
    println!("{:o}", orbit);

    let sc = SpacecraftState::from_srp_defaults(orbit, dry_mass, 1.0).with_drag(1.0, 2.0);

    let setup = Propagator::default(sc_dyn);
    let mut prop = setup.with(sc);
    prop.for_duration(prop_time).unwrap();

    let final_state = prop.state;
    println!("{}", final_state);
    println!("{}", final_state.orbit);
}

#[test]
fn std_atm_drag_earth() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24 * TimeUnit::Day;

    // Define the dynamics

    let shadow_bodies = vec![eme2k];

    let srp = SolarPressure::default(1.0, shadow_bodies, cosm.clone());
    let drag = Drag::std_atm1976(1.0, 2.0, cosm);

    let dry_mass = 300.0;

    // Build a spacecraft with SRP and Drag enabled.
    let sc_dyn = Spacecraft::with_models(OrbitalDynamics::two_body(), vec![srp, drag]);
    println!("{:o}", orbit);

    let sc = SpacecraftState::from_srp_defaults(orbit, dry_mass, 1.0).with_drag(1.0, 2.0);

    let setup = Propagator::default(sc_dyn);
    let mut prop = setup.with(sc);
    prop.for_duration(prop_time).unwrap();

    let final_state = prop.state;
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
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(
        eme2k.equatorial_radius() + 300.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        dt,
        eme2k,
    );

    let prop_time = 24 * TimeUnit::Day;

    // Define the dynamics

    let shadow_bodies = vec![eme2k];

    let srp = SolarPressure::default(1.0, shadow_bodies, cosm.clone());
    let drag = Drag::std_atm1976(1.0, 2.0, cosm);

    let dry_mass = 300.0;

    // Build a spacecraft with SRP and Drag enabled.
    let sc_dyn = Spacecraft::with_models(OrbitalDynamics::two_body(), vec![srp, drag]);
    println!("{:o}", orbit);

    let sc = SpacecraftState::from_srp_defaults(orbit, dry_mass, 1.0).with_drag(1.0, 2.0);

    let setup = Propagator::default(sc_dyn);
    let mut prop = setup.with(sc);
    prop.for_duration(prop_time).unwrap();

    let final_state = prop.state;
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
