extern crate nalgebra as na;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use std::sync::Arc;

use anise::constants::celestial_objects::{EARTH, SUN};
use anise::constants::frames::{EARTH_J2000, MOON_J2000};
use anise::prelude::Almanac;
use hifitime::JD_J2000;
use nalgebra::Vector3;
use nyx::cosmic::Orbit;
use nyx::dynamics::guidance::{FiniteBurns, LocalFrame, Mnvr, Thruster};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::md::{Event, StateParameter};
use nyx::propagators::{IntegratorOptions, Propagator};
use nyx::time::{Epoch, TimeUnits, Unit};
use nyx::{Spacecraft, State};

use nyx_space::propagators::ErrorControl;
use rstest::*;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn stop_cond_3rd_apo(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_dt = Epoch::from_mjd_tai(JD_J2000);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.01, start_dt, eme2k,
    );

    let period = state.period().unwrap();

    // Track how many times we've passed by that TA again
    let apo_event = Event::apoapsis(); // Special event shortcut!

    let setup = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut prop = setup.with(state.into(), almanac.clone());
    // Propagate for at five orbital periods so we know we've passed the third one
    // NOTE: We start counting at ZERO, so finding the 3rd means grabbing the second found.
    let (third_apo, traj) = prop.until_nth_event(5 * period, &apo_event, 2).unwrap();

    let events = traj.find(&apo_event, almanac).unwrap();
    let mut prev_event_match = events[0].state.epoch();
    for event_match in events.iter().skip(1) {
        let delta_period = event_match.state.epoch() - prev_event_match - period;
        assert!(delta_period.abs() < 10.milliseconds(), "in two body dyn, event finding should be extremely precise, instead time error of {delta_period}");
        prev_event_match = event_match.state.epoch();
    }

    let min_epoch = start_dt + 2.0 * period;
    let max_epoch = start_dt + 3.0 * period;

    println!("{}\t{}\t\t{:x}", min_epoch, max_epoch, third_apo);
    // Confirm that this is the third apoapse event which is found
    // We use a weird check because it actually converged on a time that's 0.00042 nanoseconds _after_ the max time
    assert!(
        (third_apo.orbit.epoch - min_epoch) >= 1.nanoseconds(),
        "Found apoapse is {} before min epoch",
        third_apo.orbit.epoch - min_epoch
    );
    assert!(
        (third_apo.orbit.epoch - max_epoch) <= 1.nanoseconds(),
        "Found apoapse is {} after max epoch",
        third_apo.orbit.epoch - max_epoch
    );

    assert!(
        (180.0 - third_apo.orbit.ta_deg().unwrap()).abs() < 1e-3,
        "converged, yet convergence criteria not met"
    );
}

#[rstest]
fn stop_cond_3rd_peri(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_dt = Epoch::from_mjd_tai(JD_J2000);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.01, start_dt, eme2k,
    );

    let period = state.period().unwrap();

    // Track how many times we've passed by that TA again
    let peri_event = Event::periapsis(); // Special event shortcut!

    let setup = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut prop = setup.with(state.into(), almanac.clone());
    // Propagate for at four orbital periods so we know we've passed the third one
    // NOTE: We're fetching the 3rd item because the initial state is actually at periapse,
    // which the event finder will find.
    let (third_peri, traj) = prop.until_nth_event(5 * period, &peri_event, 2).unwrap();

    let events = traj.find(&peri_event, almanac).unwrap();
    let mut prev_event_match = events[0].state.epoch();
    for event_match in events.iter().skip(1) {
        let delta_period = event_match.state.epoch() - prev_event_match - period;
        assert!(delta_period.abs() < 10.milliseconds(), "in two body dyn, event finding should be extremely precise, instead time error of {delta_period}");
        prev_event_match = event_match.state.epoch();
    }

    let min_epoch = start_dt + 2.0 * period;
    let max_epoch = start_dt + 3.0 * period;

    println!("{}\t{}\t\t{:x}", min_epoch, max_epoch, third_peri);
    // Confirm that this is the third apoapse event which is found
    // We use a weird check because it actually converged on a time that's 0.00042 nanoseconds _after_ the max time
    assert!(
        (third_peri.orbit.epoch - min_epoch) >= 1.nanoseconds(),
        "Found apoapse is {} before min epoch",
        third_peri.orbit.epoch - min_epoch
    );
    assert!(
        (third_peri.orbit.epoch - max_epoch) <= 1.nanoseconds(),
        "Found apoapse is {} after max epoch",
        third_peri.orbit.epoch - max_epoch
    );

    assert!(
        third_peri.orbit.ta_deg().unwrap().abs() < 1e-1
            || (360.0 - third_peri.orbit.ta_deg().unwrap().abs() < 1e-1),
        "converged, yet convergence criteria not met"
    );
}

#[rstest]
fn stop_cond_nrho_apo(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();
    use std::time::Instant;
    // The following test technically works, but the transformation of thousands of states
    // into another frame is quite slow...

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

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

    let state_luna = almanac.transform_to(state, MOON_J2000, None).unwrap();
    println!(
        "Start state (dynamics: Earth, Moon, Sun gravity):\n{}",
        state_luna
    );

    let bodies = vec![EARTH, SUN];
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));

    let setup = Propagator::rk89(
        dynamics,
        IntegratorOptions::with_adaptive_step_s(1.0, 60.0, 1e-6, ErrorControl::RSSCartesianStep),
    );

    // NOTE: Here, we will propagate for the maximum duration in the original frame
    // Then convert that trajectory into the other frame, and perform the search there.
    // We can only do that for spacecraft and orbit trajectories since those have a frame.
    let prop_time = 0.5 * state_luna.period().unwrap();
    let start = Instant::now();
    let (orbit, traj) = setup
        .with(Spacecraft::builder().orbit(state).build(), almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

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
    let traj_luna = traj.to_frame(MOON_J2000, almanac.clone()).unwrap();
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
    let events = traj_luna.find(&near_apo_event, almanac).unwrap();
    println!(
        "Found all {} events in {} ms",
        near_apo_event,
        (Instant::now() - end_conv).as_millis()
    );
    for event in &events {
        let event_state = event.state;
        let delta_t = event_state.epoch() - dt;
        println!("{delta_t} after start:\n{event_state:x}");
        assert!(
            (event_state.orbit.ta_deg().unwrap() - 172.0).abs() < near_apo_event.value_precision
        );
    }
}

#[rstest]
fn line_of_nodes(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_dt = Epoch::from_mjd_tai(JD_J2000);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_dt, eme2k,
    );

    let period = state.period().unwrap();

    let lon_event = Event::new(StateParameter::Longitude, 0.0);

    let setup = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut prop = setup.with(state.into(), almanac);
    let (lon_state, _) = prop.until_event(3 * period, &lon_event).unwrap();
    println!(
        "{:x} => longitude = {} degrees",
        lon_state,
        lon_state.orbit.longitude_deg()
    );

    assert!(
        lon_state.orbit.longitude_deg().abs() < lon_event.value_precision,
        "converged, yet convergence criteria not met"
    );
}

#[rstest]
fn latitude(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_dt = Epoch::from_mjd_tai(JD_J2000);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_dt, eme2k,
    );

    let period = state.period().unwrap();

    let lat_event = Event::new(StateParameter::Latitude, 2.0);

    let setup = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut prop = setup.with(state.into(), almanac);
    let (lon_state, _) = prop.until_event(3 * period, &lat_event).unwrap();
    println!(
        "{:x} => latitude = {} degrees",
        lon_state,
        lon_state.orbit.latitude_deg().unwrap()
    );

    assert!(
        (2.0 - lon_state.orbit.latitude_deg().unwrap()).abs() < lat_event.value_precision,
        "converged, yet convergence criteria not met"
    );
}

#[rstest]
fn event_and_combination(almanac: Arc<Almanac>) {
    /// Event combinations cannot be implemented with a brent solver (the approache used by Nyx to find events).
    /// Instead, two events must be sought for, one after another.
    use nyx::dynamics::GuidanceMode;

    let _ = pretty_env_logger::try_init();

    // Setup a scenario

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let epoch = Epoch::now().unwrap();
    // We're at periapse of a GTO
    let orbit =
        Orbit::try_keplerian_altitude(42_165.0, 0.7, 30.0, 45.0, 45.0, 0.01, epoch, eme2k).unwrap();

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
        LocalFrame::VNC,
    )]);

    let dynamics = SpacecraftDynamics::from_guidance_law(OrbitalDynamics::two_body(), burn);

    let setup = Propagator::default(dynamics);
    let mut prop = setup.with(sc, almanac.clone());

    // First, propagate until apoapsis
    let (sc_apo, traj) = prop
        .until_event(orbit.period().unwrap() * 4.0, &Event::apoapsis())
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
    let traj_moon = traj.to_frame(MOON_J2000, almanac.clone()).unwrap();

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
    if let Ok(matching_states) = traj_moon.find(&event, almanac.clone()) {
        for sc_decl_zero in matching_states {
            decl_deg = sc_decl_zero
                .state
                .value(StateParameter::Declination)
                .unwrap();
            println!("{sc_decl_zero} => decl = {} deg", decl_deg,);
            assert!((decl_deg - 6.0).abs() < 3.0);
        }

        // We should be able to find a similar event with a tighter bound too.
        if let Ok(tighter_states) = traj_moon.find(
            &Event::specific(StateParameter::Declination, decl_deg, 1.0, Unit::Minute),
            almanac,
        ) {
            for sc_decl_zero in tighter_states {
                let found_decl_deg = sc_decl_zero
                    .state
                    .value(StateParameter::Declination)
                    .unwrap();
                println!("{sc_decl_zero} => decl = {} deg", found_decl_deg);
                assert!((decl_deg - found_decl_deg).abs() < 1.0);
            }
        }
    }
}
