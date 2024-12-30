extern crate nyx_space as nyx;
use self::nyx::cosmic::{GuidanceMode, Orbit, Spacecraft};
use self::nyx::dynamics::guidance::{FiniteBurns, Maneuver, Thruster};
use self::nyx::dynamics::{OrbitalDynamics, SpacecraftDynamics};
use self::nyx::linalg::Vector3;
use self::nyx::propagators::{IntegratorOptions, Propagator};
use self::nyx::time::{Epoch, Unit};
use self::nyx::utils::rss_orbit_vec_errors;
use crate::propagation::GMAT_EARTH_GM;
use nyx::dynamics::guidance::LocalFrame;
use std::sync::Arc;

use anise::constants::celestial_objects::{JUPITER_BARYCENTER, MOON, SUN};
use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn val_transfer_schedule_no_depl(almanac: Arc<Almanac>) {
    /*
        NOTE: Due to how lifetime of variables work in Rust, we need to define all of the
        components of a spacecraft before defining the spacecraft itself.
    */

    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    // Build the initial spacecraft state
    let start_time = Epoch::from_gregorian_tai_at_midnight(2002, 1, 1);
    let orbit = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_time, eme2k,
    );

    // Define the thruster
    let monoprop = Thruster {
        thrust_N: 10.0,
        isp_s: 300.0,
    };
    let dry_mass = 1e3;
    let prop_mass = 756.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, prop_mass, monoprop, GuidanceMode::Coast);

    let prop_time = 50.0 * Unit::Minute;

    let end_time = start_time + prop_time;

    // Define the dynamics
    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let orbital_dyn = OrbitalDynamics::point_masses(bodies);

    // Define the maneuver and its schedule
    let mnvr0 = Maneuver::from_time_invariant(
        Epoch::from_gregorian_tai_at_midnight(2002, 1, 1),
        end_time,
        1.0, // Full thrust
        Vector3::new(1.0, 0.0, 0.0),
        LocalFrame::VNC,
    );

    let schedule = FiniteBurns::from_mnvrs(vec![mnvr0]);

    // And create the spacecraft with that controller
    // Disable prop mass decrement
    let sc = SpacecraftDynamics::from_guidance_law_no_decr(orbital_dyn, schedule);
    // Setup a propagator, and propagate for that duration
    // NOTE: We specify the use an RK89 to match the GMAT setup.
    let final_state = Propagator::rk89(sc, IntegratorOptions::with_fixed_step(10.0 * Unit::Second))
        .with(sc_state, almanac)
        .for_duration(prop_time)
        .unwrap();

    // Compute the errors
    let rslt = Orbit::cartesian(
        4_172.396_780_515_64f64,
        436.944_560_056_202_8,
        -6_518.328_156_815_674,
        -3.979_610_765_995_537,
        5.540_316_900_333_103,
        -2.207_082_771_390_863,
        end_time,
        eme2k,
    );

    let (err_r, err_v) = rss_orbit_vec_errors(
        &final_state.orbit.to_cartesian_pos_vel(),
        &rslt.to_cartesian_pos_vel(),
    );
    println!("Absolute errors");
    let delta = final_state.orbit.to_cartesian_pos_vel() - rslt.to_cartesian_pos_vel();
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} km\tvel = {:.5e} km/s",
        err_r, err_v,
    );

    assert!(err_r < 5e-10, "finite burn position wrong: {:.5e}", err_r);
    assert!(err_v < 6e-13, "finite burn velocity wrong: {:.5e}", err_v);

    // Ensure that there was no change in prop mass since tank depletion was off
    assert!(
        (final_state.mass.prop_mass_kg - prop_mass).abs() < f64::EPSILON,
        "incorrect prop mass"
    );
}

#[rstest]
fn val_transfer_schedule_depl_cov_test(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    // Build the initial spacecraft state
    let start_time = Epoch::from_gregorian_tai_at_midnight(2002, 1, 1);
    let orbit = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_time, eme2k,
    );

    // Define the thruster
    let monoprop = Thruster {
        thrust_N: 10.0,
        isp_s: 300.0,
    };
    let dry_mass = 1e3;
    let prop_mass = 756.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, prop_mass, monoprop, GuidanceMode::Coast);

    let prop_time = 50.0 * Unit::Minute;

    let end_time = start_time + prop_time;

    // Define the dynamics
    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let orbital_dyn = OrbitalDynamics::point_masses(bodies);

    // With 100% thrust: RSS errors:     pos = 3.14651e1 km      vel = 3.75245e-2 km/s

    // Define the maneuver and its schedule
    let mnvr0 = Maneuver::from_time_invariant(
        Epoch::from_gregorian_tai_at_midnight(2002, 1, 1),
        end_time,
        1.0, // Full thrust
        Vector3::new(1.0, 0.0, 0.0),
        LocalFrame::VNC,
    );

    let schedule = FiniteBurns::from_mnvrs(vec![mnvr0]);

    // And create the spacecraft with that controller
    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, schedule);
    // Setup a propagator, and propagate for that duration
    // NOTE: We specify the use an RK89 to match the GMAT setup.
    let setup = Propagator::rk89(sc, IntegratorOptions::with_fixed_step(10.0 * Unit::Second));
    let final_state = setup
        .with(sc_state, almanac.clone())
        .for_duration(prop_time)
        .unwrap();

    // Compute the errors
    let rslt = Orbit::cartesian(
        4_172.433_936_615_18,
        436.936_159_720_413,
        -6_518.368_821_953_345,
        -3.979_569_721_967_499,
        5.540_321_146_839_762,
        -2.207_146_819_283_441,
        end_time,
        eme2k,
    );

    let rslt_prop_mass = 745.802_837_870_161;

    let (err_r, err_v) = rss_orbit_vec_errors(
        &final_state.orbit.to_cartesian_pos_vel(),
        &rslt.to_cartesian_pos_vel(),
    );
    println!("Absolute errors");
    let delta = final_state.orbit.to_cartesian_pos_vel() - rslt.to_cartesian_pos_vel();
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} km\tvel = {:.5e} km/s",
        err_r, err_v,
    );

    assert!(err_r < 5e-10, "finite burn position wrong: {:.5e}", err_r);
    assert!(err_v < 5e-13, "finite burn velocity wrong: {:.5e}", err_v);

    let delta_prop_mass = (final_state.mass.prop_mass_kg - rslt_prop_mass).abs();
    println!("Absolute prop mass error: {:.0e} kg", delta_prop_mass);
    assert!(delta_prop_mass < 2e-10, "incorrect prop mass");

    // Now, test that backward propagation of maneuvers also works.
    let backward_state = setup
        .with(final_state, almanac)
        .for_duration(-prop_time)
        .unwrap();
    println!("Reached: {}\nWanted:  {}", backward_state, sc_state);

    let (err_r, err_v) = rss_orbit_vec_errors(
        &backward_state.orbit.to_cartesian_pos_vel(),
        &sc_state.orbit.to_cartesian_pos_vel(),
    );
    println!("Backprop Absolute errors");
    let delta = backward_state.orbit.to_cartesian_pos_vel() - sc_state.orbit.to_cartesian_pos_vel();
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} km\tvel = {:.5e} km/s",
        err_r, err_v,
    );

    assert!(
        err_r < 1.0,
        "finite burn backprop position wrong: {:.5e}",
        err_r
    );
    assert!(
        err_v < 1e-3,
        "finite burn backprop velocity wrong: {:.5e}",
        err_v
    );
}

#[rstest]
fn val_transfer_single_maneuver_depl_cov_test(almanac: Arc<Almanac>) {
    /* This is the same test as val_transfer_schedule_depl but uses the maneuver directly as the guidance law. It should work in the same way. */

    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    // Build the initial spacecraft state
    let start_time = Epoch::from_gregorian_tai_at_midnight(2002, 1, 1);
    let orbit = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_time, eme2k,
    );

    // Define the thruster
    let monoprop = Thruster {
        thrust_N: 10.0,
        isp_s: 300.0,
    };
    let dry_mass_kg = 1e3;
    let prop_mass_kg = 756.0;
    let sc_state = Spacecraft::from_thruster(
        orbit,
        dry_mass_kg,
        prop_mass_kg,
        monoprop,
        GuidanceMode::Coast,
    );

    let prop_time = 50.0 * Unit::Minute;

    let end_time = start_time + prop_time;

    // Define the dynamics
    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let orbital_dyn = OrbitalDynamics::point_masses(bodies);

    // With 100% thrust: RSS errors:     pos = 3.14651e1 km      vel = 3.75245e-2 km/s

    // Define the maneuver and its schedule
    let mnvr0 = Maneuver::from_time_invariant(
        Epoch::from_gregorian_tai_at_midnight(2002, 1, 1),
        end_time,
        1.0, // Full thrust
        Vector3::new(1.0, 0.0, 0.0),
        LocalFrame::VNC,
    );

    // And create the spacecraft with that controller
    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, Arc::new(mnvr0));
    // Setup a propagator, and propagate for that duration
    // NOTE: We specify the use an RK89 to match the GMAT setup.
    let setup = Propagator::rk89(sc, IntegratorOptions::with_fixed_step(10.0 * Unit::Second));
    let final_state = setup
        .with(sc_state, almanac.clone())
        .for_duration(prop_time)
        .unwrap();

    // Compute the errors
    let rslt = Orbit::cartesian(
        4_172.433_936_615_18,
        436.936_159_720_413,
        -6_518.368_821_953_345,
        -3.979_569_721_967_499,
        5.540_321_146_839_762,
        -2.207_146_819_283_441,
        end_time,
        eme2k,
    );

    let rslt_prop_mass = 745.802_837_870_161;

    let (err_r, err_v) = rss_orbit_vec_errors(
        &final_state.orbit.to_cartesian_pos_vel(),
        &rslt.to_cartesian_pos_vel(),
    );
    println!("Absolute errors");
    let delta = final_state.orbit.to_cartesian_pos_vel() - rslt.to_cartesian_pos_vel();
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} km\tvel = {:.5e} km/s",
        err_r, err_v,
    );

    assert!(err_r < 5e-10, "finite burn position wrong: {:.5e}", err_r);
    assert!(err_v < 5e-13, "finite burn velocity wrong: {:.5e}", err_v);

    let delta_prop_mass = (final_state.mass.prop_mass_kg - rslt_prop_mass).abs();
    println!("Absolute prop mass error: {:.0e} kg", delta_prop_mass);
    assert!(delta_prop_mass < 2e-10, "incorrect prop mass");

    // Now, test that backward propagation of maneuvers also works.
    let backward_state = setup
        .with(final_state, almanac)
        .for_duration(-prop_time)
        .unwrap();
    println!("Reached: {}\nWanted:  {}", backward_state, sc_state);

    let (err_r, err_v) = rss_orbit_vec_errors(
        &backward_state.orbit.to_cartesian_pos_vel(),
        &sc_state.orbit.to_cartesian_pos_vel(),
    );
    println!("Backprop Absolute errors");
    let delta = backward_state.orbit.to_cartesian_pos_vel() - sc_state.orbit.to_cartesian_pos_vel();
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} km\tvel = {:.5e} km/s",
        err_r, err_v,
    );

    assert!(
        err_r < 1.0,
        "finite burn backprop position wrong: {:.5e}",
        err_r
    );
    assert!(
        err_v < 1e-3,
        "finite burn backprop velocity wrong: {:.5e}",
        err_v
    );
}
