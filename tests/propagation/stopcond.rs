extern crate nalgebra as na;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use hifitime::J2000_OFFSET;
use na::Vector3;
use nyx::cosmic::{Bodies, Cosm, Orbit};
use nyx::dynamics::guidance::{FiniteBurns, Mnvr, Thruster};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::md::{Event, EventEvaluator, StateParameter};
use nyx::propagators::error_ctrl::RSSCartesianStep;
use nyx::propagators::{PropOpts, Propagator};
use nyx::time::{Epoch, TimeUnits, Unit};
use nyx::{Spacecraft, State};

#[test]
fn stop_cond_3rd_apo() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.01, start_dt, eme2k,
    );

    let period = state.period();

    // Track how many times we've passed by that TA again
    let apo_event = Event::apoapsis(); // Special event shortcut!

    let setup = Propagator::default(OrbitalDynamics::two_body());
    let mut prop = setup.with(state);
    // Propagate for at five orbital periods so we know we've passed the third one
    // NOTE: We start counting at ZERO, so finding the 3rd means grabbing the second found.
    let (third_apo, traj) = prop.until_nth_event(5 * period, &apo_event, 2).unwrap();

    let events = traj.find_all(&apo_event).unwrap();
    let mut prev_event_match = events[0].epoch();
    for event_match in events.iter().skip(1) {
        let delta_period = event_match.epoch() - prev_event_match - period;
        assert!(delta_period.abs() < 50.microseconds(), "in two body dyn, event finding should be extremely precise, instead time error of {delta_period}");
        prev_event_match = event_match.epoch();
    }

    let min_epoch = start_dt + 2.0 * period;
    let max_epoch = start_dt + 3.0 * period;

    println!("{}\t{}\t\t{:x}", min_epoch, max_epoch, third_apo);
    // Confirm that this is the third apoapse event which is found
    // We use a weird check because it actually converged on a time that's 0.00042 nanoseconds _after_ the max time
    assert!(
        (third_apo.epoch - min_epoch) >= 1.nanoseconds(),
        "Found apoapse is {} before min epoch",
        third_apo.epoch - min_epoch
    );
    assert!(
        (third_apo.epoch - max_epoch) <= 1.nanoseconds(),
        "Found apoapse is {} after max epoch",
        third_apo.epoch - max_epoch
    );

    assert!(
        (180.0 - third_apo.ta_deg()).abs() < 1e-3,
        "converged, yet convergence criteria not met"
    );
}

#[test]
fn stop_cond_3rd_peri() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.01, start_dt, eme2k,
    );

    let period = state.period();

    // Track how many times we've passed by that TA again
    let peri_event = Event::periapsis(); // Special event shortcut!

    let setup = Propagator::default(OrbitalDynamics::two_body());
    let mut prop = setup.with(state);
    // Propagate for at four orbital periods so we know we've passed the third one
    // NOTE: We're fetching the 3rd item because the initial state is actually at periapse,
    // which the event finder will find.
    let (third_peri, traj) = prop.until_nth_event(5 * period, &peri_event, 2).unwrap();

    let events = traj.find_all(&peri_event).unwrap();
    let mut prev_event_match = events[0].epoch();
    for event_match in events.iter().skip(1) {
        let delta_period = event_match.epoch() - prev_event_match - period;
        assert!(delta_period.abs() < 50.microseconds(), "in two body dyn, event finding should be extremely precise, instead time error of {delta_period}");
        prev_event_match = event_match.epoch();
    }

    let min_epoch = start_dt + 2.0 * period;
    let max_epoch = start_dt + 3.0 * period;

    println!("{}\t{}\t\t{:x}", min_epoch, max_epoch, third_peri);
    // Confirm that this is the third apoapse event which is found
    // We use a weird check because it actually converged on a time that's 0.00042 nanoseconds _after_ the max time
    assert!(
        (third_peri.epoch - min_epoch) >= 1.nanoseconds(),
        "Found apoapse is {} before min epoch",
        third_peri.epoch - min_epoch
    );
    assert!(
        (third_peri.epoch - max_epoch) <= 1.nanoseconds(),
        "Found apoapse is {} after max epoch",
        third_peri.epoch - max_epoch
    );

    assert!(
        third_peri.ta_deg().abs() < 1e-1 || (360.0 - third_peri.ta_deg().abs() < 1e-1),
        "converged, yet convergence criteria not met"
    );
}

#[test]
fn stop_cond_nrho_apo() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }
    use std::time::Instant;
    // The following test technically works, but the transformation of thousands of states
    // into another frame is quite slow...
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");
    let luna = cosm.frame("Luna");

    let dt = Epoch::from_gregorian_tai(2021, 5, 29, 19, 51, 16, 852_000);
    let state = Orbit::cartesian(
        166_473.631_302_239_7,
        -274_715.487_253_382_7,
        -211_233.210_176_686_7,
        0.933_451_604_520_018_4,
        0.436_775_046_841_900_9,
        -0.082_211_021_250_348_95,
        dt,
        eme2k,
    );

    let state_luna = cosm.frame_chg(&state, luna);
    println!(
        "Start state (dynamics: Earth, Moon, Sun gravity):\n{}",
        state_luna
    );

    let bodies = vec![Bodies::Earth, Bodies::Sun];
    let dynamics = OrbitalDynamics::point_masses(&bodies, cosm.clone());

    let setup = Propagator::rk89(
        dynamics,
        PropOpts::with_adaptive_step_s(1.0, 60.0, 1e-6, RSSCartesianStep {}),
    );

    // NOTE: Here, we will propagate for the maximum duration in the original frame
    // Then convert that trajectory into the other frame, and perform the search there.
    // We can only do that for spacecraft and orbit trajectories since those have a frame.
    let prop_time = 0.5 * state_luna.period();
    let start = Instant::now();
    let (orbit, traj) = setup.with(state).for_duration_with_traj(prop_time).unwrap();

    let end_prop = Instant::now();
    println!(
        "Propagated for {} in {} ms:\n{:x}\n{}\n",
        prop_time,
        (end_prop - start).as_millis(),
        orbit,
        traj
    );

    // Create the event
    let near_apo_event = Event::new(StateParameter::TrueAnomaly, 172.0);

    // Convert this trajectory into the Luna frame
    let traj_luna = traj.to_frame(luna, cosm).unwrap();
    let end_conv = Instant::now();
    println!(
        "Converted EME2000 trajectory into Moon J2000 in {} ms\nFrom: {}\nTo  : {}",
        (end_conv - end_prop).as_millis(),
        traj,
        traj_luna
    );
    assert!(
        (traj.first().epoch() - traj_luna.first().epoch()).abs() < 1.milliseconds(),
        "First epoch of converted trajectories do not match"
    );
    assert!(
        (traj.last().epoch() - traj_luna.last().epoch()).abs() < 1.milliseconds(),
        "Last epoch of converted trajectories do not match"
    );

    // Now, find all of the requested events
    let events = traj_luna.find_all(&near_apo_event).unwrap();
    println!(
        "Found all {} events in {} ms",
        near_apo_event,
        (Instant::now() - end_conv).as_millis()
    );
    for event_state in &events {
        let delta_t = event_state.epoch() - dt;
        println!("{} after start:\n{:x}", delta_t, event_state);
        assert!((event_state.ta_deg() - 172.0).abs() < near_apo_event.value_precision);
    }
}

#[test]
fn line_of_nodes() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_dt, eme2k,
    );

    let period = state.period();

    let lon_event = Event::new(StateParameter::GeodeticLongitude, 0.0);

    let setup = Propagator::default(OrbitalDynamics::two_body());
    let mut prop = setup.with(state);
    let (lon_state, _) = prop.until_event(3 * period, &lon_event).unwrap();
    println!(
        "{:x} => longitude = {} degrees",
        lon_state,
        lon_state.geodetic_longitude_deg()
    );

    assert!(
        lon_state.geodetic_longitude_deg().abs() < lon_event.value_precision,
        "converged, yet convergence criteria not met"
    );
}

#[test]
fn latitude() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_dt, eme2k,
    );

    let period = state.period();

    let lat_event = Event::new(StateParameter::GeodeticLatitude, 2.0);

    let setup = Propagator::default_dp78(OrbitalDynamics::two_body());
    let mut prop = setup.with(state);
    let (lon_state, _) = prop.until_event(3 * period, &lat_event).unwrap();
    println!(
        "{:x} => latitude = {} degrees",
        lon_state,
        lon_state.geodetic_latitude_deg()
    );

    assert!(
        (2.0 - lon_state.geodetic_latitude_deg()).abs() < lat_event.value_precision,
        "converged, yet convergence criteria not met"
    );
}

#[test]
fn event_and_combination() {
    /// Event combinations cannot be implemented with a brent solver (the approache used by Nyx to find events).
    /// Instead, two events must be sought for, one after another.
    use nyx::cosmic::Frame;
    use nyx::dynamics::GuidanceMode;

    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    // Setup a scenario
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");
    let moonj2k = cosm.frame("Moon J2000");

    let epoch = Epoch::now().unwrap();
    // We're at periapse of a GTO
    let orbit = Orbit::keplerian_altitude(42_165.0, 0.7, 30.0, 45.0, 45.0, 0.01, epoch, eme2k);

    let sc = Spacecraft::from_thruster(
        orbit,
        100.0,
        50.0,
        Thruster {
            isp_s: 300.0,
            thrust_N: 50.0,
        },
        GuidanceMode::Thrust,
    );

    println!(
        "{sc}\tinitial c3 = {}",
        sc.value(StateParameter::C3).unwrap()
    );

    // Thrust in the +X direction continuously
    let burn = FiniteBurns::from_mnvrs(vec![Mnvr::from_time_invariant(
        epoch + 1.minutes(),
        epoch + 15.minutes(),
        1.0,
        Vector3::x(),
        Frame::VNC,
    )]);

    let orbital_dyn = OrbitalDynamics::two_body();
    let dynamics = SpacecraftDynamics::from_guidance_law(orbital_dyn, burn);

    let setup = Propagator::default(dynamics);
    let mut prop = setup.with(sc);

    // First, propagate until apoapsis
    let (sc_apo, traj) = prop
        .until_event(orbit.period() * 4.0, &Event::apoapsis())
        .unwrap();

    // Check that the fuel always decreases or stays constant
    let mut cur_fuel = traj.states[0].fuel_mass_kg;
    for state in traj.states.iter().skip(1) {
        assert!(
            state.fuel_mass_kg - cur_fuel <= 1e-6, // Check that fuel never increases, at least a mg level
            "{cur_fuel} > {}",
            state.fuel_mass_kg
        );
        cur_fuel = state.fuel_mass_kg;
    }

    // Convert the trajectory to the Moon frame
    let traj_moon = traj.to_frame(moonj2k, cosm).unwrap();

    let sc_moon_apo = traj_moon.at(sc_apo.epoch()).unwrap();

    println!(
        "Earth Apoapse\n{:x}\tc3 = {} km^2/s^2\n{:x}\tdecl = {} deg",
        sc_apo,
        sc_apo.value(StateParameter::C3).unwrap(),
        sc_moon_apo,
        sc_moon_apo.value(StateParameter::Declination).unwrap()
    );

    println!(
        "End of prop\n{:x}\tc3 = {} km^2/s^2\n{:x}\tdecl = {} deg",
        traj.last(),
        traj.last().value(StateParameter::C3).unwrap(),
        traj_moon.last(),
        traj_moon.last().value(StateParameter::Declination).unwrap()
    );

    // Now let's find when the declination with the Moon is zero.
    // Within one minute and with a precision of 3.0 degrees.
    // NOTE: We're unwrapping here, so if the event isn't found, this will cause the test to fail.
    let event = Event::specific(StateParameter::Declination, 6.0, 3.0, Unit::Minute);
    let mut decl_deg = 0.0;
    if let Ok(matching_states) = traj_moon.find_all(&event) {
        for sc_decl_zero in matching_states {
            decl_deg = sc_decl_zero.value(StateParameter::Declination).unwrap();
            println!(
                "{event}: {} => decl = {} deg",
                event.eval_string(&sc_decl_zero),
                decl_deg,
            );
            assert!((decl_deg - 6.0).abs() < 3.0);
        }

        // We should be able to find a similar event with a tighter bound too.
        if let Ok(tighter_states) = traj_moon.find_all(&Event::specific(
            StateParameter::Declination,
            decl_deg,
            1.0,
            Unit::Minute,
        )) {
            for sc_decl_zero in tighter_states {
                let found_decl_deg = sc_decl_zero.value(StateParameter::Declination).unwrap();
                println!("{sc_decl_zero:x} => decl = {} deg", found_decl_deg);
                assert!((decl_deg - found_decl_deg).abs() < 1.0);
            }
        }
    }
}
