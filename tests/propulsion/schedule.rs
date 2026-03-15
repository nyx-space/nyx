extern crate nyx_space as nyx;
use self::nyx::cosmic::{GuidanceMode, Orbit, Spacecraft};
use self::nyx::dynamics::guidance::{Maneuver, Thruster};
use self::nyx::dynamics::sequence::{
    AccelModels, ForceModels, GuidanceConfig, Phase, PropagatorConfig, SpacecraftSequence,
};
use self::nyx::dynamics::{OrbitalDynamics, PointMasses, SpacecraftDynamics};
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
    let eme2k = almanac
        .frame_info(EARTH_J2000)
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

    // Define the maneuver and its schedule
    let mnvr0 = Maneuver::from_time_invariant(
        Epoch::from_gregorian_tai_at_midnight(2002, 1, 1),
        end_time,
        1.0, // Full thrust
        Vector3::new(1.0, 0.0, 0.0),
        LocalFrame::VNC,
    );

    let mut sc_seq = SpacecraftSequence::default();

    sc_seq.propagators.insert(
        "Earth".to_string(),
        PropagatorConfig {
            method: nyx::propagators::IntegratorMethod::RungeKutta89,
            options: IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
            accel_models: AccelModels {
                point_masses: Some(PointMasses::new(vec![MOON, SUN, JUPITER_BARYCENTER])),
                gravity_field: None,
            },
            force_models: ForceModels {
                solar_pressure: None,
                drag: None,
            },
        },
    );

    sc_seq
        .thruster_sets
        .insert("Monoprop".to_string(), monoprop);

    sc_seq.seq.insert(
        start_time,
        Phase::Activity {
            name: "Burn".to_string(),
            propagator: "Earth".to_string(),
            guidance: Some(Box::new(GuidanceConfig::FiniteBurn {
                maneuver: mnvr0,
                thruster_model: "Monoprop".to_string(),
                disable_prop_mass: true,
            })),
            on_entry: None,
            disabled: false,
        },
    );

    sc_seq.seq.insert(end_time, Phase::Terminate);

    sc_seq.setup(almanac.clone()).unwrap();

    let trajectories = sc_seq.propagate(sc_state, None, almanac).unwrap();
    let final_state = trajectories.last().unwrap().last();

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

    println!("RSS errors:\tpos = {err_r:.5e} km\tvel = {err_v:.5e} km/s",);

    assert!(err_r < 5e-8, "finite burn position wrong: {err_r:.5e}");
    assert!(err_v < 2e-8, "finite burn velocity wrong: {err_v:.5e}");

    // Ensure that there was no change in prop mass since tank depletion was off
    assert!(
        (final_state.mass.prop_mass_kg - prop_mass).abs() < f64::EPSILON,
        "incorrect prop mass"
    );
}

#[rstest]
fn val_transfer_schedule_depl(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_info(EARTH_J2000)
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
    let orbital_dyn = OrbitalDynamics::point_masses(bodies.clone());

    // With 100% thrust: RSS errors:     pos = 3.14651e1 km      vel = 3.75245e-2 km/s

    // Define the maneuver and its schedule
    let mnvr0 = Maneuver::from_time_invariant(
        Epoch::from_gregorian_tai_at_midnight(2002, 1, 1),
        end_time,
        1.0, // Full thrust
        Vector3::new(1.0, 0.0, 0.0),
        LocalFrame::VNC,
    );

    let mut sc_seq = SpacecraftSequence::default();

    sc_seq.propagators.insert(
        "Earth".to_string(),
        PropagatorConfig {
            method: nyx::propagators::IntegratorMethod::RungeKutta89,
            options: IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
            accel_models: AccelModels {
                point_masses: Some(PointMasses::new(bodies)),
                gravity_field: None,
            },
            force_models: ForceModels {
                solar_pressure: None,
                drag: None,
            },
        },
    );

    sc_seq
        .thruster_sets
        .insert("Monoprop".to_string(), monoprop);

    sc_seq.seq.insert(
        start_time,
        Phase::Activity {
            name: "Burn".to_string(),
            propagator: "Earth".to_string(),
            guidance: Some(Box::new(GuidanceConfig::FiniteBurn {
                maneuver: mnvr0,
                thruster_model: "Monoprop".to_string(),
                disable_prop_mass: false,
            })),
            on_entry: None,
            disabled: false,
        },
    );

    sc_seq.seq.insert(end_time, Phase::Terminate);

    sc_seq.setup(almanac.clone()).unwrap();

    let trajectories = sc_seq.propagate(sc_state, None, almanac.clone()).unwrap();
    let final_state = trajectories.last().unwrap().last();

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

    println!("RSS errors:\tpos = {err_r:.5e} km\tvel = {err_v:.5e} km/s",);

    assert!(err_r < 5e-8, "finite burn position wrong: {err_r:.5e}");
    assert!(err_v < 2e-8, "finite burn velocity wrong: {err_v:.5e}");

    let delta_prop_mass = (final_state.mass.prop_mass_kg - rslt_prop_mass).abs();
    println!("Absolute prop mass error: {delta_prop_mass:.0e} kg");
    assert!(delta_prop_mass < 1e-5, "incorrect prop mass");

    // Now, test that backward propagation of maneuvers also works.
    // SpacecraftSequence does not directly support backwards propagation of sequences right now,
    // so we set up a reverse propagator manually with the maneuver.
    let backward_sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, Arc::new(mnvr0));
    let setup = Propagator::rk89(
        backward_sc,
        IntegratorOptions::with_fixed_step(10.0 * Unit::Second),
    );

    let backward_state = setup
        .with(*final_state, almanac)
        .for_duration(-prop_time)
        .unwrap();
    println!("Reached: {backward_state}\nWanted:  {sc_state}");

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

    println!("RSS errors:\tpos = {err_r:.5e} km\tvel = {err_v:.5e} km/s",);

    assert!(
        err_r < 1.0,
        "finite burn backprop position wrong: {err_r:.5e}"
    );
    assert!(
        err_v < 1e-3,
        "finite burn backprop velocity wrong: {err_v:.5e}"
    );
}

#[rstest]
fn val_transfer_single_maneuver_depl(almanac: Arc<Almanac>) {
    /* This is the same test as val_transfer_schedule_depl but uses the maneuver directly as the guidance law. It should work in the same way. */

    let eme2k = almanac
        .frame_info(EARTH_J2000)
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

    println!("RSS errors:\tpos = {err_r:.5e} km\tvel = {err_v:.5e} km/s",);

    assert!(err_r < 5e-10, "finite burn position wrong: {err_r:.5e}");
    assert!(err_v < 5e-13, "finite burn velocity wrong: {err_v:.5e}");

    let delta_prop_mass = (final_state.mass.prop_mass_kg - rslt_prop_mass).abs();
    println!("Absolute prop mass error: {delta_prop_mass:.0e} kg");
    assert!(delta_prop_mass < 2e-10, "incorrect prop mass");

    // Now, test that backward propagation of maneuvers also works.
    let backward_state = setup
        .with(final_state, almanac)
        .for_duration(-prop_time)
        .unwrap();
    println!("Reached: {backward_state}\nWanted:  {sc_state}");

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

    println!("RSS errors:\tpos = {err_r:.5e} km\tvel = {err_v:.5e} km/s",);

    assert!(
        err_r < 1.0,
        "finite burn backprop position wrong: {err_r:.5e}"
    );
    assert!(
        err_v < 1e-3,
        "finite burn backprop velocity wrong: {err_v:.5e}"
    );
}

#[rstest]
fn finite_burns_respects_gaps_between_maneuvers(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_info(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let start_time = Epoch::from_gregorian_tai_at_midnight(2002, 1, 1);
    let orbit = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_time, eme2k,
    );

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

    let mnvr0 = Maneuver::from_time_invariant(
        start_time,
        start_time + 60.0 * Unit::Second,
        1.0,
        Vector3::new(1.0, 0.0, 0.0),
        LocalFrame::VNC,
    );

    // notice the 120 second gap here

    let mnvr1 = Maneuver::from_time_invariant(
        start_time + 180.0 * Unit::Second,
        start_time + 240.0 * Unit::Second,
        1.0,
        Vector3::new(1.0, 0.0, 0.0),
        LocalFrame::VNC,
    );

    let mut sc_seq = SpacecraftSequence::default();

    sc_seq.propagators.insert(
        "Earth".to_string(),
        PropagatorConfig {
            method: nyx::propagators::IntegratorMethod::RungeKutta89,
            options: IntegratorOptions::with_fixed_step(1.0 * Unit::Second),
            accel_models: AccelModels {
                point_masses: Some(PointMasses::new(vec![MOON, SUN, JUPITER_BARYCENTER])),
                gravity_field: None,
            },
            force_models: ForceModels {
                solar_pressure: None,
                drag: None,
            },
        },
    );

    sc_seq
        .thruster_sets
        .insert("Monoprop".to_string(), monoprop);

    sc_seq.seq.insert(
        start_time,
        Phase::Activity {
            name: "Initial Coast".to_string(),
            propagator: "Earth".to_string(),
            guidance: None,
            on_entry: None,
            disabled: false,
        },
    );

    sc_seq.seq.insert(
        mnvr0.start,
        Phase::Activity {
            name: "Burn 1".to_string(),
            propagator: "Earth".to_string(),
            guidance: Some(Box::new(GuidanceConfig::FiniteBurn {
                maneuver: mnvr0,
                thruster_model: "Monoprop".to_string(),
                disable_prop_mass: false,
            })),
            on_entry: None,
            disabled: false,
        },
    );

    sc_seq.seq.insert(
        mnvr0.end,
        Phase::Activity {
            name: "Intermediate Coast".to_string(),
            propagator: "Earth".to_string(),
            guidance: None,
            on_entry: None,
            disabled: false,
        },
    );

    sc_seq.seq.insert(
        mnvr1.start,
        Phase::Activity {
            name: "Burn 2".to_string(),
            propagator: "Earth".to_string(),
            guidance: Some(Box::new(GuidanceConfig::FiniteBurn {
                maneuver: mnvr1,
                thruster_model: "Monoprop".to_string(),
                disable_prop_mass: false,
            })),
            on_entry: None,
            disabled: false,
        },
    );

    sc_seq.seq.insert(
        mnvr1.end,
        Phase::Activity {
            name: "Final Coast".to_string(),
            propagator: "Earth".to_string(),
            guidance: None,
            on_entry: None,
            disabled: false,
        },
    );

    sc_seq
        .seq
        .insert(start_time + 300.0 * Unit::Second, Phase::Terminate);

    sc_seq.setup(almanac.clone()).unwrap();

    let trajectories = sc_seq.propagate(sc_state, None, almanac).unwrap();

    let get_mass_at = |epoch| {
        trajectories
            .iter()
            .find_map(|t| t.at(epoch).ok())
            .unwrap()
            .mass
            .prop_mass_kg
    };

    let m_during_first = get_mass_at(start_time + 30.0 * Unit::Second);
    let m_after_first = get_mass_at(start_time + 70.0 * Unit::Second);
    let m_in_gap = get_mass_at(start_time + 150.0 * Unit::Second);
    let m_during_second = get_mass_at(start_time + 220.0 * Unit::Second);
    let m_after_second = get_mass_at(start_time + 250.0 * Unit::Second);
    let m_after_second_2 = get_mass_at(start_time + 280.0 * Unit::Second);

    assert!(
        m_after_first < m_during_first,
        "first burn did not consume propellant"
    );
    assert!(
        (m_in_gap - m_after_first).abs() < 1e-9,
        "propellant changed in burn gap"
    );
    assert!(
        m_during_second < m_in_gap,
        "second burn did not consume propellant"
    );
    assert!(
        m_after_second < m_during_second,
        "second burn did not consume propellant"
    );
    assert!(
        (m_after_second_2 - m_after_second).abs() < 1e-9,
        "propellant changed after burn"
    );
}
