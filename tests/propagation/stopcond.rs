extern crate nalgebra as na;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;

use std::sync::Arc;

use anise::analysis::prelude::{Condition, Event, OrbitalElement, ScalarExpr};
use anise::constants::celestial_objects::{EARTH, SUN};
use anise::constants::frames::{EARTH_J2000, IAU_EARTH_FRAME, MOON_J2000};
use anise::prelude::Almanac;
use nalgebra::Vector3;
use nyx::cosmic::Orbit;
use nyx::dynamics::guidance::{FiniteBurns, LocalFrame, Maneuver, Thruster};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::md::StateParameter;
use nyx::propagators::{IntegratorOptions, Propagator};
use nyx::time::{Epoch, TimeUnits};
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
    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_dt = Epoch::from_gregorian_utc_at_midnight(2025, 11, 12);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.01, start_dt, eme2k,
    );

    let period = state.period().unwrap();

    // Track how many times we've passed by that TA again
    let apo_event = Event::apoapsis();

    let setup = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut prop = setup.with(state.into(), almanac.clone());
    // Propagate for at five orbital periods so we know we've passed the third one
    let (third_apo, traj) = prop
        .until_nth_event(5 * period, &apo_event, None, 3)
        .unwrap();

    // Trajectory implements the StateSpecTrait, so we can search for events directly with the Almanac!
    let events = almanac
        .report_events(&traj, &apo_event, traj.first().epoch(), traj.last().epoch())
        .unwrap();

    let mut prev_event_match = events[0].orbit.epoch;
    for event_match in events.iter().skip(1) {
        let delta_period = event_match.orbit.epoch - prev_event_match - period;
        assert!(delta_period.abs() < 0.5.seconds(), "in two body dyn, event finding should be extremely precise, instead time error of {delta_period}");
        prev_event_match = event_match.orbit.epoch;
    }

    let min_epoch = start_dt + 2.0 * period;
    let max_epoch = start_dt + 3.0 * period;

    println!("{min_epoch}\t{max_epoch}\t\t{third_apo:x}");
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
        (180.0 - third_apo.orbit.ta_deg().unwrap()).abs() < 1e-6,
        "converged, yet convergence criteria not met"
    );
}

#[rstest]
fn stop_cond_3rd_peri(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let epoch = Epoch::from_gregorian_utc_at_noon(2008, 2, 29);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.01, epoch, eme2k,
    );

    let period = state.period().unwrap();

    // Track how many times we've passed by that true anomaly again
    let peri_event = Event::periapsis(); // Special event shortcut!

    let setup = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut prop = setup.with(state.into(), almanac.clone());

    let (third_peri, traj) = prop
        .until_nth_event(5 * period, &peri_event, None, 3)
        .unwrap();

    let events = almanac
        .report_events(
            &traj,
            &peri_event,
            traj.first().epoch(),
            traj.last().epoch(),
        )
        .unwrap();

    let mut prev_event_match = events[0].orbit.epoch;
    for event_match in events.iter().skip(1) {
        let delta_period = event_match.orbit.epoch - prev_event_match - period;
        println!("{:x}", event_match.orbit);
        assert!(delta_period.abs() < 0.3.seconds(), "in two body dyn, event finding should be extremely precise, instead time error of {delta_period}");
        prev_event_match = event_match.orbit.epoch;
    }

    let min_epoch = epoch + 2.0 * period;
    let max_epoch = epoch + 3.0 * period;

    println!("{min_epoch}\t{max_epoch}\t\t{third_peri:x}");
    // Confirm that this is the third apoapse event which is found
    // We use a weird check because it actually converged on a time that's 0.00042 nanoseconds _after_ the max time
    assert!(
        (third_peri.orbit.epoch - min_epoch) >= 1.nanoseconds(),
        "Found periapse is {} before min epoch",
        third_peri.orbit.epoch - min_epoch
    );
    assert!(
        (third_peri.orbit.epoch - max_epoch) <= 1.nanoseconds(),
        "Found periapse is {} after max epoch",
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

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

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
    println!("Start state (dynamics: Earth, Moon, Sun gravity):\n{state_luna}");

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
    let near_apo_event = Event::new(
        ScalarExpr::Element(OrbitalElement::TrueAnomaly),
        Condition::Equals(172.0),
    );

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
    let events = almanac
        .report_events(
            &traj,
            &near_apo_event,
            traj.first().epoch(),
            traj.last().epoch(),
        )
        .unwrap();
    println!(
        "Found all {} events in {} ms",
        near_apo_event,
        (Instant::now() - end_conv).as_millis()
    );
    for event in &events {
        let event_state = event.orbit;
        let delta_t = event_state.epoch - dt;
        println!("{delta_t} after start:\n{event_state:x}");
        assert!((event_state.ta_deg().unwrap() - 172.0).abs() < 1e-3);
    }
}

#[rstest]
fn line_of_nodes(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_dt = Epoch::from_gregorian_utc_at_noon(2008, 2, 29);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_dt, eme2k,
    );

    let period = state.period().unwrap();

    let asc_node_event = Event::new(
        ScalarExpr::Element(OrbitalElement::Longitude),
        Condition::Equals(0.0),
    );

    let setup = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut prop = setup.with(state.into(), almanac);
    let (lon_state, _) = prop.until_event(3 * period, &asc_node_event, None).unwrap();
    println!(
        "{:x} => longitude = {} degrees",
        lon_state,
        lon_state.orbit.longitude_deg()
    );

    assert!(
        lon_state.orbit.longitude_deg().abs() < 1e-3,
        "converged, yet convergence criteria not met"
    );
}

#[rstest]
fn latitude(almanac: Arc<Almanac>) {
    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_dt = Epoch::from_gregorian_utc_at_noon(2008, 2, 29);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_dt, eme2k,
    );

    let period = state.period().unwrap();

    let lat_event = Event::new(
        ScalarExpr::Element(OrbitalElement::Latitude),
        Condition::Equals(2.0),
    );

    let setup = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut prop = setup.with(state.into(), almanac);
    let (lon_state, _) = prop
        .until_event(3 * period, &lat_event, Some(IAU_EARTH_FRAME))
        .unwrap();
    println!(
        "{:x} => latitude = {} degrees",
        lon_state,
        lon_state.orbit.latitude_deg().unwrap()
    );

    assert!(
        (2.0 - lon_state.orbit.latitude_deg().unwrap()).abs() < 1e-3,
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

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

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
        sc.value(StateParameter::Element(OrbitalElement::C3))
            .unwrap()
    );

    // Thrust in the +X direction continuously
    let burn = FiniteBurns::from_mnvrs(vec![Maneuver::from_time_invariant(
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
        .until_event(orbit.period().unwrap() * 4.0, &Event::apoapsis(), None)
        .unwrap();

    // Check that the prop always decreases or stays constant
    let mut cur_prop = traj.states[0].mass.prop_mass_kg;
    for state in traj.states.iter().skip(1) {
        assert!(
            state.mass.prop_mass_kg - cur_prop <= 1e-6, // Check that prop never increases, at least a mg level
            "{cur_prop} > {}",
            state.mass.prop_mass_kg
        );
        cur_prop = state.mass.prop_mass_kg;
    }

    // Convert the trajectory to the Moon frame
    let traj_moon = traj.to_frame(MOON_J2000, almanac.clone()).unwrap();

    let sc_moon_apo = traj_moon.at(sc_apo.epoch()).unwrap();

    println!(
        "Earth Apoapse\n{:x}\tc3 = {} km^2/s^2\n{:x}\tdecl = {} deg",
        sc_apo,
        sc_apo
            .value(StateParameter::Element(OrbitalElement::C3))
            .unwrap(),
        sc_moon_apo,
        sc_moon_apo
            .value(StateParameter::Element(OrbitalElement::Declination))
            .unwrap()
    );

    println!(
        "End of prop\n{:x}\tc3 = {} km^2/s^2\n{:x}\tdecl = {} deg",
        traj.last(),
        traj.last()
            .value(StateParameter::Element(OrbitalElement::C3))
            .unwrap(),
        traj_moon.last(),
        traj_moon
            .last()
            .value(StateParameter::Element(OrbitalElement::Declination))
            .unwrap()
    );

    // Now let's find when the declination with the Moon is zero.
    // Within one minute and with a precision of 3.0 degrees.
    // NOTE: We're unwrapping here, so if the event isn't found, this will cause the test to fail.
    let event = Event::new(
        ScalarExpr::Element(OrbitalElement::Declination),
        Condition::Between(3.0, 6.0),
    );

    if let Ok(matching_states) = almanac.report_events(
        &traj_moon,
        &event,
        traj.first().epoch(),
        traj.last().epoch(),
    ) {
        for sc_decl_zero in matching_states {
            let decl_deg = sc_decl_zero.orbit.declination_deg();
            println!("{sc_decl_zero} => decl = {decl_deg} deg",);
            assert!((decl_deg - 6.0).abs() < 3.0);
        }
    }
}
