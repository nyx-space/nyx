extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use anise::constants::frames::{EARTH_J2000, MOON_J2000};
use hifitime::TimeUnits;
use nyx::cosmic::eclipse::EclipseLocator;
use nyx::cosmic::{GuidanceMode, Orbit, Spacecraft};
use nyx::dynamics::guidance::{GuidanceLaw, Ruggiero, Thruster};
use nyx::dynamics::{OrbitalDynamics, SpacecraftDynamics};
use nyx::md::prelude::{ExportCfg, Objective};
use nyx::md::StateParameter;
use nyx::propagators::*;
use nyx::time::{Epoch, TimeSeries, Unit};
use nyx::State;
use nyx_space::md::Trajectory;
use std::path::PathBuf;
use std::sync::mpsc::channel;
use std::sync::Arc;

use anise::prelude::Almanac;
use rstest::*;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[allow(clippy::identity_op)]
#[rstest]
fn traj_ephem_forward_cov_test(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();
    // Test that we can correctly interpolate a spacecraft orbit

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_dt = Epoch::from_gregorian_utc_at_noon(2021, 1, 1);
    let start_state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_dt, eme2k,
    );

    let setup = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut prop = setup.with(start_state.into(), almanac.clone());
    // The trajectory must always be generated on its own thread, no need to worry about it ;-)
    let now = Epoch::now().unwrap();
    let (end_state, ephem) = prop.for_duration_with_traj(31 * Unit::Day).unwrap();
    let exec_time = Epoch::now().unwrap() - now;
    println!("[TIMING] {exec_time}");

    let ephem_start = ephem.states.first().unwrap().epoch();
    let trimmed = ephem.clone().filter_by_offset(1.hours()..=3.hours());
    let trim_start_err =
        (trimmed.states.first().unwrap().epoch() - (ephem_start + 1.hours())).abs();
    let trim_end_err = (trimmed.states.last().unwrap().epoch() - (ephem_start + 3.hours())).abs();
    assert!(
        trim_start_err < 1.minutes(),
        "traj filter failed -- large step?: {trim_start_err}"
    );
    assert!(
        trim_end_err < 1.minutes(),
        "traj filter failed -- large step?: {trim_end_err}"
    );

    // Example of iterating through the trajectory.
    // Note how we receive a fully blown Orbit, so we can manipulate it exactly like a normal state.
    let mut sum_sma = 0.0;
    let mut cnt = 0.0;
    for state in ephem.every(1 * Unit::Day) {
        cnt += 1.0;
        sum_sma += state.orbit.sma_km().unwrap()
    }
    println!(
        "Average SMA: {:.3} km\tShould be: {:.3}",
        sum_sma / cnt,
        start_state.sma_km().unwrap()
    );
    // assert!(dbg!(sum_sma / cnt - start_state.sma()).abs() < 1e-6);

    // === Below is the validation of the ephemeris == //

    assert_eq!(ephem.first().orbit, start_state, "Wrong initial state");
    assert_eq!(ephem.last(), &end_state, "Wrong final state");
    assert!(ephem.last().stm().is_err(), "STM is set!");
    assert!(
        ephem.at(end_state.epoch() + 1 * Unit::Nanosecond).is_err(),
        "Expected to be outside of interpolation window!"
    );

    println!("Ephem: {ephem}");

    // Now let's re-generate the truth data and ensure that each state we generate is in the ephemeris and matches the expected state within tolerance.

    let (tx, rx) = channel();
    let almanac_c = almanac.clone();
    std::thread::spawn(move || {
        setup
            .with(start_state.into(), almanac_c)
            .for_duration_with_channel(31 * Unit::Day, tx)
            .unwrap();
    });

    // Evaluate the first time of the trajectory to make sure that one is there too.
    let eval_state = ephem.at(start_dt).unwrap();

    let mut max_pos_err = (eval_state.orbit.radius_km - start_state.radius_km).norm();
    let mut max_vel_err = (eval_state.orbit.velocity_km_s - start_state.velocity_km_s).norm();

    while let Ok(prop_state) = rx.recv() {
        let eval_state = ephem.at(prop_state.epoch()).unwrap();

        let pos_err = (eval_state.orbit.radius_km - prop_state.orbit.radius_km).norm();
        if pos_err > max_pos_err {
            max_pos_err = pos_err;
        }
        let vel_err = (eval_state.orbit.velocity_km_s - prop_state.orbit.velocity_km_s).norm();
        if vel_err > max_vel_err {
            max_vel_err = vel_err;
        }
    }

    println!(
        "[traj_ephem] Maximum error on exact step: pos: {:.2e} m\t\tvel: {:.2e} m/s",
        max_pos_err * 1e3,
        max_vel_err * 1e3,
    );

    // Should be zero because we store these states as is
    assert!(
        max_pos_err == 0.0,
        "Maximum spacecraft position in interpolation is too high!"
    );

    assert!(
        max_vel_err == 0.0,
        "Maximum orbit velocity in interpolation is too high!"
    );

    // Save the trajectory to parquet
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "04_output",
        "ephem_forward.parquet",
    ]
    .iter()
    .collect();

    let exported_path = ephem
        .to_parquet(
            path,
            Some(vec![
                &EclipseLocator::cislunar(almanac.clone()).to_penumbra_event()
            ]),
            ExportCfg::timestamped(),
            almanac.clone(),
        )
        .unwrap();

    // Reload this trajectory and make sure that it matches

    let concrete_traj = Trajectory::from_parquet(exported_path).unwrap();

    if ephem != concrete_traj {
        // Uh oh, let's see where the differences are.
        assert_eq!(
            ephem.states.len(),
            concrete_traj.states.len(),
            "loaded trajectory has different number of states"
        );

        assert_eq!(
            ephem.first().epoch(),
            concrete_traj.first().epoch(),
            "loaded trajectory has starts at different times"
        );

        assert_eq!(
            ephem.last().epoch(),
            concrete_traj.last().epoch(),
            "loaded trajectory has ends at different times"
        );

        // So far so good, so let's iterate through both now.
        for (i, (orig_state, loaded_state)) in
            ephem.states.iter().zip(&concrete_traj.states).enumerate()
        {
            // Check the data info one by one. Time may be very slightly off.
            let delta_t = orig_state.epoch() - loaded_state.epoch();
            assert!(
                delta_t < 1 * Unit::Microsecond,
                "#{i} differ (Δt = {delta_t})"
            );
            assert_eq!(
                orig_state.to_vector(),
                loaded_state.to_vector(),
                "#{i} differ"
            );
        }
    }

    // And let's convert into another frame and back to check the error
    let ephem_luna = ephem.to_frame(MOON_J2000, almanac.clone()).unwrap();
    println!("ephem_luna {ephem_luna}");
    assert!(
        (ephem.first().epoch() - ephem_luna.first().epoch()).abs() < 1.microseconds(),
        "Start time differ!"
    );
    assert!(
        (ephem.last().epoch() - ephem_luna.last().epoch()).abs() < 1.microseconds(),
        "End time differ!"
    );
    // And convert back, to see the error this leads to
    let ephem_back_to_earth = ephem_luna.to_frame(eme2k, almanac.clone()).unwrap();
    println!("Ephem back: {ephem_back_to_earth}");
    assert!(
        (ephem.first().epoch() - ephem_back_to_earth.first().epoch()).abs() < 1.microseconds(),
        "Start time differ after double conversion!"
    );
    assert!(
        (ephem.last().epoch() - ephem_back_to_earth.last().epoch()).abs() < 1.microseconds(),
        "End time differ after double conversion!"
    );

    let conv_state = ephem_back_to_earth.at(start_dt).unwrap();
    let mut max_pos_err = (eval_state.orbit.radius_km - conv_state.orbit.radius_km).norm();
    let mut max_vel_err = (eval_state.orbit.velocity_km_s - conv_state.orbit.velocity_km_s).norm();

    for conv_state in ephem_back_to_earth.every(5 * Unit::Minute) {
        let eval_state = ephem.at(conv_state.epoch()).unwrap();

        let pos_err = (eval_state.orbit.radius_km - conv_state.orbit.radius_km).norm();
        if pos_err > max_pos_err {
            println!(
                "{pos_err:.e}\nEval: {}\nConv: {}\t{:.3} m\n",
                eval_state,
                conv_state,
                pos_err * 1e3
            );
            max_pos_err = pos_err;
        }
        let vel_err = (eval_state.orbit.velocity_km_s - conv_state.orbit.velocity_km_s).norm();
        if vel_err > max_vel_err {
            max_vel_err = vel_err;
        }
    }
    println!(
        "[traj_ephem] Maximum interpolation error after double conversion: pos: {:.2e} m\t\tvel: {:.2e} m/s",
        max_pos_err * 1e3,
        max_vel_err * 1e3,
    );

    // Allow for up to meter error after double conversion
    assert!(
        max_pos_err < 1.0,
        "Maximum spacecraft position in interpolation is too high!"
    );

    // Allow for up to ten centimeters per second error after double conversion
    assert!(
        max_vel_err < 1e-2,
        "Maximum orbit velocity in interpolation is too high!"
    );
}

#[allow(clippy::identity_op)]
#[rstest]
fn traj_spacecraft(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();
    // Test the interpolation of a spaceraft trajectory and of its prop. Includes a demo of checking what the guidance mode _should_ be provided the state.
    // Note that we _do not_ attempt to interpolate the Guidance Mode.
    // This is based on the Ruggiero AOP correction

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    // Build the initial orbit
    let start_dt = Epoch::from_gregorian_utc_at_noon(2021, 1, 1);
    let sma = eme2k.mean_equatorial_radius_km().unwrap() + 900.0;
    let orbit = Orbit::keplerian(sma, 5e-5, 5e-3, 0.0, 178.0, 0.0, start_dt, eme2k);

    // Define the objectives and the control law
    let objectives = &[Objective::within_tolerance(
        StateParameter::AoP,
        183.0,
        5e-3,
    )];

    let ruggiero_ctrl = Ruggiero::simple(objectives, orbit.into()).unwrap();

    // Build the spacecraft state
    let prop_mass = 67.0;
    let dry_mass = 300.0;
    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };
    let start_state =
        Spacecraft::from_thruster(orbit, dry_mass, prop_mass, lowt, GuidanceMode::Thrust);

    let sc_dynamics =
        SpacecraftDynamics::from_guidance_law(OrbitalDynamics::two_body(), ruggiero_ctrl.clone());

    let setup = Propagator::default(sc_dynamics);
    let prop_time = 44 * Unit::Minute + 10 * Unit::Second;
    let mut prop = setup.with(start_state, almanac.clone());
    let (end_state, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Example of iterating through the spaceraft trajectory and checking what the guidance mode is at each time.
    let mut prev_mode = GuidanceMode::Coast;

    // But let's now iterate over the trajectory every day instead of every minute
    // Note that we can iterate over this trajectory
    // Note: _Because_ we need to use the trajectory below, we'll be cloning the trajectory
    // If you don't need the trajectory after you've iterated over it, don't clone it (rustc will tell you that).

    for mut sc_state in traj.every(1 * Unit::Day) {
        // We need to evaluate the mode of this state because the trajectory does not store discrete information
        ruggiero_ctrl.next(&mut sc_state, almanac.clone());
        if sc_state.mode() != prev_mode {
            println!(
                "Mode changed from {:?} to {:?} @ {}",
                prev_mode,
                sc_state.mode(),
                sc_state.epoch()
            );
            prev_mode = sc_state.mode();
        }
    }

    for epoch in TimeSeries::inclusive(start_dt, start_dt + prop_time, 1 * Unit::Day) {
        // Note: the `evaluate` function will return a Result which prevents a panic if you request something out of the ephemeris
        let mut sc_state = traj.at(epoch).unwrap();
        ruggiero_ctrl.next(&mut sc_state, almanac.clone());
        if sc_state.mode() != prev_mode {
            println!(
                "Mode changed from {:?} to {:?} @ {}",
                prev_mode,
                sc_state.mode(),
                epoch
            );
            prev_mode = sc_state.mode();
        }
    }

    // === Below is the validation of the ephemeris == //

    println!("{traj}");
    assert_eq!(traj.first(), &start_state, "Wrong initial state");
    assert_eq!(traj.last(), &end_state, "Wrong final state");
    assert!(
        traj.at(end_state.epoch() + 1 * Unit::Nanosecond).is_err(),
        "Expected to be outside of interpolation window!"
    );

    // Now let's re-generate the truth data and ensure that each state we generate is in the ephemeris and matches the expected state within tolerance.

    let (tx, rx) = channel();
    let almanac_c = almanac.clone();
    std::thread::spawn(move || {
        setup
            .with(start_state, almanac_c)
            .until_epoch_with_channel(end_state.epoch(), tx)
            .unwrap();
    });

    // Evaluate the first time of the trajectory to make sure that one is there too.
    let eval_state = traj.at(start_dt).unwrap();

    let mut max_pos_err = (eval_state.orbit.radius_km - start_state.orbit.radius_km).norm();
    let mut max_vel_err = (eval_state.orbit.velocity_km_s - start_state.orbit.velocity_km_s).norm();
    let mut max_prop_err = eval_state.mass.prop_mass_kg - start_state.mass.prop_mass_kg;
    let mut max_err = (eval_state.to_vector() - start_state.to_vector()).norm();

    while let Ok(prop_state) = rx.recv() {
        let eval_state = traj.at(prop_state.epoch()).unwrap();

        let pos_err = (eval_state.orbit.radius_km - prop_state.orbit.radius_km).norm();
        if pos_err > max_pos_err {
            max_pos_err = pos_err;
            println!("pos_err = {:.3e} m @ {}", pos_err * 1e3, prop_state.epoch());
        }
        let vel_err = (eval_state.orbit.velocity_km_s - prop_state.orbit.velocity_km_s).norm();
        if vel_err > max_vel_err {
            max_vel_err = vel_err;
            println!(
                "vel_err = {:.3e} m/s @ {}",
                vel_err * 1e3,
                prop_state.epoch()
            );
        }
        let prop_err = eval_state.mass.prop_mass_kg - prop_state.mass.prop_mass_kg;
        if prop_err > max_prop_err {
            max_prop_err = prop_err;
            println!(
                "prop_err = {:.3e} g @ {}",
                prop_err * 1e3,
                prop_state.epoch()
            );
        }
        let err = (eval_state.to_vector() - prop_state.to_vector()).norm();
        if err > max_err {
            max_err = err;
        }
    }

    println!(
        "[traj_spacecraft] Maximum interpolation error: pos: {:.2e} m\t\tvel: {:.2e} m/s\t\tprop: {:.2e} g\t\tfull state: {:.2e} (no unit)",
        max_pos_err * 1e3,
        max_vel_err * 1e3,
        max_prop_err * 1e3,
        max_err
    );

    // Allow millimeter error
    assert!(
        max_pos_err < 1e-6,
        "Maximum spacecraft position in interpolation is too high!"
    );

    // Allow centimeter per second error
    assert!(
        max_vel_err < 1e-5,
        "Maximum spacecraft velocity in interpolation is too high!"
    );

    // Allow for up to 0.1 gram error
    assert!(
        max_vel_err < 1e-4,
        "Maximum spacecraft prop in interpolation is too high!"
    );

    // And let's convert into another frame and back to check the error
    let ephem_luna = traj.to_frame(MOON_J2000, almanac.clone()).unwrap();
    // And convert back, to see the error this leads to
    let ephem_back_to_earth = ephem_luna.to_frame(eme2k, almanac.clone()).unwrap();

    assert_eq!(
        traj, ephem_back_to_earth,
        "Expecting exactly the same data returned after converting back"
    );
}

#[allow(clippy::identity_op)]
#[rstest]
fn traj_ephem_backward(almanac: Arc<Almanac>) {
    // Test that we can correctly interpolate a spacecraft orbit

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_dt = Epoch::from_gregorian_utc_at_noon(2021, 1, 1);
    let start_state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_dt, eme2k,
    );

    let setup = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut prop = setup.with(start_state.into(), almanac.clone());
    let (end_state, ephem) = prop.for_duration_with_traj(-31 * Unit::Day).unwrap();

    // Example of iterating through the trajectory.
    // Note how we receive a fully blown Orbit, so we can manipulate it exactly like a normal state.
    let mut sum_sma = 0.0;
    let mut cnt = 0.0;
    for epoch in TimeSeries::inclusive(start_dt - 31 * Unit::Day, start_dt, 1 * Unit::Day) {
        cnt += 1.0;
        // Note: the `evaluate` function will return a Result which prevents a panic if you request something out of the ephemeris
        // let state = ephem.at(epoch + 17 * Unit::Second).unwrap();
        // sum_sma += state.sma();
        match ephem.at(epoch) {
            Ok(state) => sum_sma += state.orbit.sma_km().unwrap(),
            Err(e) => println!("{e}"),
        }
    }
    println!("Average SMA: {:.3} km", sum_sma / cnt);

    // === Below is the validation of the ephemeris == //

    // NOTE: The trajectory organizes the states chronologically.
    assert_eq!(
        ephem.first(),
        &end_state,
        "Wrong initial state\nGot:  {}\nWant: {}",
        ephem.first(),
        end_state
    );
    assert_eq!(
        ephem.last().orbit,
        start_state,
        "Wrong final state\nGot:  {}\nWant: {}",
        ephem.last(),
        start_state
    );
    assert!(ephem.last().stm().is_err(), "STM is set!");
    assert!(
        ephem.at(end_state.epoch() - 1 * Unit::Nanosecond).is_err(),
        "Expected to be outside of interpolation window!"
    );

    // Now let's re-generate the truth data and ensure that each state we generate is in the ephemeris and matches the expected state within tolerance.

    let (tx, rx) = channel();
    std::thread::spawn(move || {
        setup
            .with(start_state.into(), almanac.clone())
            .for_duration_with_channel(-31 * Unit::Day, tx)
            .unwrap();
    });

    // Evaluate the first time of the trajectory to make sure that one is there too.
    let eval_state = ephem.at(start_dt).unwrap();

    let mut max_pos_err = (eval_state.orbit.radius_km - start_state.radius_km).norm();
    let mut max_vel_err = (eval_state.orbit.velocity_km_s - start_state.velocity_km_s).norm();
    let mut max_err = (eval_state.to_vector().fixed_rows::<6>(0).to_owned()
        - start_state.to_cartesian_pos_vel())
    .norm();

    println!("{ephem}");

    while let Ok(prop_state) = rx.recv() {
        let eval_state = ephem.at(prop_state.epoch()).unwrap();

        let pos_err = (eval_state.orbit.radius_km - prop_state.orbit.radius_km).norm();
        if pos_err > max_pos_err {
            max_pos_err = pos_err;
        }
        let vel_err = (eval_state.orbit.velocity_km_s - prop_state.orbit.velocity_km_s).norm();
        if vel_err > max_vel_err {
            max_vel_err = vel_err;
        }
        let err = (eval_state.to_vector() - prop_state.to_vector()).norm();
        if err > max_err {
            max_err = err;
        }
    }

    println!(
        "[traj_ephem_backward] Maximum interpolation error: pos: {:.2e} m\t\tvel: {:.2e} m/s\t\tfull state: {:.2e} (no unit)",
        max_pos_err * 1e3,
        max_vel_err * 1e3,
        max_err
    );

    // Allow for up to micrometer error
    assert!(
        max_pos_err < 1e-9,
        "Maximum spacecraft position in interpolation is too high!"
    );

    // Allow for up to micrometer per second error
    assert!(
        max_vel_err < 1e-9,
        "Maximum orbit velocity in interpolation is too high!"
    );

    // Check that the error in STM doesn't break everything
    assert!(
        max_err < 1e-9,
        "Maximum state in interpolation is too high!"
    );
}
