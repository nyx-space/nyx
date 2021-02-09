extern crate nyx_space as nyx;

use nyx::celestia::{Cosm, GuidanceMode, Orbit, SpacecraftState};
use nyx::dynamics::thrustctrl::{Achieve, Ruggiero, ThrustControl, Thruster};
use nyx::dynamics::{OrbitalDynamics, Spacecraft};
use nyx::md::{Ephemeris, ScTraj};
use nyx::propagators::*;
use nyx::time::{Epoch, TimeSeries, TimeUnit};
use nyx::{State, TimeTagged};
use std::sync::mpsc::channel;

#[allow(clippy::identity_op)]
#[test]
fn traj_ephem() {
    // Test that we can correctly interpolate a spacecraft orbit
    let (tx, rx) = channel();
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_dt = Epoch::from_gregorian_utc_at_noon(2021, 1, 1);
    let start_state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_dt, eme2k,
    );

    // The trajectory must always be generated on its own thread.
    let ephem_thread = std::thread::spawn(move || Ephemeris::new(start_state, rx));

    let setup = Propagator::default(OrbitalDynamics::two_body());
    let mut prop = setup.with(start_state.with_stm()).with_tx(tx);
    let end_state = prop.for_duration(31 * TimeUnit::Day).unwrap();

    // Retrieve the ephemeris by unwrapping it twice:
    // + The first is to unwrap the thread (i.e. assume the thread has not failed)
    // + The second assumes that the generation of the ephemeris didn't fail.
    let ephem = ephem_thread.join().unwrap().unwrap();

    // Example of iterating through the trajectory.
    // Note how we receive a fully blown Orbit, so we can manipulate it exactly like a normal state.
    let mut sum_sma = 0.0;
    let mut cnt = 0.0;
    for epoch in TimeSeries::inclusive(start_dt, start_dt + 31 * TimeUnit::Day, 1 * TimeUnit::Day) {
        cnt += 1.0;
        // Note: the `evaluate` function will return a Result which prevents a panic if you request something out of the ephemeris
        // let state = ephem.evaluate(epoch + 17 * TimeUnit::Second).unwrap();
        // sum_sma += state.sma();
        match ephem.evaluate(epoch) {
            Ok(state) => sum_sma += state.sma(),
            Err(e) => println!("{}", e),
        }
    }
    println!("Average SMA: {:.3} km", sum_sma / cnt);

    // === Below is the validation of the ephemeris == //

    assert_eq!(
        ephem.segments.len(),
        1010,
        "Wrong number of expected segments"
    );

    assert_eq!(ephem.first(), start_state, "Wrong initial state");
    assert_eq!(ephem.last(), end_state, "Wrong final state");
    assert!(ephem.last().stm().norm() > 0.0, "STM is not set!");
    assert!(
        ephem
            .evaluate(end_state.dt + 1 * TimeUnit::Nanosecond)
            .is_err(),
        "Expected to be outside of interpolation window!"
    );

    // Now let's re-generate the truth data and ensure that each state we generate is in the ephemeris and matches the expected state within tolerance.

    let (tx, rx) = channel();
    std::thread::spawn(move || {
        let mut prop = setup.with(start_state).with_tx(tx);
        prop.for_duration(31 * TimeUnit::Day).unwrap();
    });

    // Evaluate the first time of the trajectory to make sure that one is there too.
    let eval_state = ephem.evaluate(start_dt).unwrap();

    let mut max_pos_err = (eval_state.radius() - start_state.radius()).norm();
    let mut max_vel_err = (eval_state.velocity() - start_state.velocity()).norm();
    let mut max_err = (eval_state.as_vector().unwrap() - start_state.as_vector().unwrap()).norm();

    while let Ok(prop_state) = rx.recv() {
        let eval_state = ephem.evaluate(prop_state.dt).unwrap();

        let pos_err = (eval_state.radius() - prop_state.radius()).norm();
        if pos_err > max_pos_err {
            max_pos_err = pos_err;
        }
        let vel_err = (eval_state.velocity() - prop_state.velocity()).norm();
        if vel_err > max_vel_err {
            max_vel_err = vel_err;
        }
        let err = (eval_state.as_vector().unwrap() - prop_state.as_vector().unwrap()).norm();
        if err > max_err {
            max_err = err;
        }
    }

    println!(
        "[traj_ephem] Maximum interpolation error: pos: {:.2e} m\t\tvel: {:.2e} m/s\t\tfull state: {:.2e} (no unit)",
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

#[allow(clippy::identity_op)]
#[test]
fn traj_spacecraft() {
    // Test the interpolation of a spaceraft trajectory and of its fuel. Includes a demo of checking what the guidance mode _should_ be provided the state.
    // Note that we _do not_ attempt to interpolate the Guidance Mode.
    // This is based on the Ruggiero AOP correction
    let (tx, rx) = channel();
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    // Build the initial orbit
    let start_dt = Epoch::from_gregorian_utc_at_noon(2021, 1, 1);
    let sma = eme2k.equatorial_radius() + 900.0;
    let orbit = Orbit::keplerian(sma, 5e-5, 5e-3, 0.0, 178.0, 0.0, start_dt, eme2k);

    // Define the objectives and the control law
    let objectives = vec![Achieve::Aop {
        target: 183.0,
        tol: 5e-3,
    }];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit);

    // Build the spacecraft state
    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    // Define the thruster
    let lowt = Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    };
    let start_state =
        SpacecraftState::with_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc_dynamics = Spacecraft::with_ctrl(OrbitalDynamics::two_body(), ruggiero_ctrl.clone());

    // The trajectory must always be generated on its own thread.
    let ephem_thread = std::thread::spawn(move || ScTraj::new(start_state, rx));

    let setup = Propagator::default(sc_dynamics);
    let prop_time = 44 * TimeUnit::Minute + 10 * TimeUnit::Second;
    let mut prop = setup.with(start_state).with_tx(tx);
    let end_state = prop.for_duration(prop_time).unwrap();

    // Retrieve the ephemeris by unwrapping it twice:
    // + The first is to unwrap the thread (i.e. assume the thread has not failed)
    // + The second assumes that the generation of the ephemeris didn't fail.
    let ephem = ephem_thread.join().unwrap().unwrap();

    // Example of iterating through the spaceraft trajectory and checking what the guidance mode is at each time.
    let mut prev_mode = GuidanceMode::Coast;
    for epoch in TimeSeries::inclusive(start_dt, start_dt + prop_time, 1 * TimeUnit::Day) {
        // Note: the `evaluate` function will return a Result which prevents a panic if you request something out of the ephemeris
        let sc_state = ephem.evaluate(epoch).unwrap();
        let mode_then = ruggiero_ctrl.next(&sc_state);
        if mode_then != prev_mode {
            println!(
                "Mode changed from {:?} to {:?} @ {}",
                prev_mode, mode_then, epoch
            );
            prev_mode = mode_then;
        }
    }

    // === Below is the validation of the ephemeris == //

    assert_eq!(ephem.segments.len(), 3, "Wrong number of expected segments");

    assert_eq!(ephem.first(), start_state, "Wrong initial state");
    assert_eq!(ephem.last(), end_state, "Wrong final state");
    assert!(
        ephem
            .evaluate(end_state.epoch() + 1 * TimeUnit::Nanosecond)
            .is_err(),
        "Expected to be outside of interpolation window!"
    );

    // Now let's re-generate the truth data and ensure that each state we generate is in the ephemeris and matches the expected state within tolerance.

    let (tx, rx) = channel();
    std::thread::spawn(move || {
        let mut prop = setup.with(start_state).with_tx(tx);
        prop.for_duration(prop_time).unwrap();
    });

    // Evaluate the first time of the trajectory to make sure that one is there too.
    let eval_state = ephem.evaluate(start_dt).unwrap();

    let mut max_pos_err = (eval_state.orbit.radius() - start_state.orbit.radius()).norm();
    let mut max_vel_err = (eval_state.orbit.velocity() - start_state.orbit.velocity()).norm();
    let mut max_fuel_err = eval_state.fuel_mass_kg - start_state.fuel_mass_kg;
    let mut max_err = (eval_state.as_vector().unwrap() - start_state.as_vector().unwrap()).norm();

    while let Ok(prop_state) = rx.recv() {
        let eval_state = ephem.evaluate(prop_state.epoch()).unwrap();

        let pos_err = (eval_state.orbit.radius() - prop_state.orbit.radius()).norm();
        if pos_err > max_pos_err {
            max_pos_err = pos_err;
        }
        let vel_err = (eval_state.orbit.velocity() - prop_state.orbit.velocity()).norm();
        if vel_err > max_vel_err {
            max_vel_err = vel_err;
        }
        let fuel_err = eval_state.fuel_mass_kg - prop_state.fuel_mass_kg;
        if fuel_err > max_fuel_err {
            max_fuel_err = fuel_err;
        }
        let err = (eval_state.as_vector().unwrap() - prop_state.as_vector().unwrap()).norm();
        if err > max_err {
            max_err = err;
        }
    }

    println!(
        "[traj_spacecraft] Maximum interpolation error: pos: {:.2e} m\t\tvel: {:.2e} m/s\t\tfuel: {:.2e} g\t\tfull state: {:.2e} (no unit)",
        max_pos_err * 1e3,
        max_vel_err * 1e3,
        max_fuel_err * 1e3,
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
        "Maximum spacecraft velocity in interpolation is too high!"
    );

    // Allow for up to microgram error
    assert!(
        max_vel_err < 1e-9,
        "Maximum spacecraft fuel in interpolation is too high!"
    );
}
