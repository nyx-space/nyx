extern crate nyx_space as nyx;
use std::{fmt::Write, sync::Arc};

use anise::{
    constants::frames::{EARTH_J2000, IAU_EARTH_FRAME, SUN_J2000},
    prelude::Almanac,
};
use rstest::*;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn event_tracker_true_anomaly(almanac: Arc<Almanac>) {
    use nyx::cosmic::eclipse::{EclipseLocator, EclipseState};
    use nyx::md::prelude::*;
    use nyx::od::GroundStation;

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let dt = Epoch::from_gregorian_tai_at_noon(2020, 1, 1);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );

    let prop_time = state.period().unwrap();

    // Track how many times we've passed by that TA again
    let peri_event = Event::periapsis(); // Special event shortcut!
    let apo_event = Event::apoapsis(); // Special event shortcut!
    let ta_event0 = Event::new(StateParameter::TrueAnomaly, 35.1);
    let ta_event1 = Event::new(StateParameter::TrueAnomaly, 235.1);

    let events = vec![peri_event, apo_event, ta_event0, ta_event1];

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::rk89(dynamics, PropOpts::with_tolerance(1e-9));
    let mut prop = setup.with(state.into(), almanac.clone());
    let (_, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Find all of the events
    for e in &events {
        let found_events = traj.find(e, almanac.clone()).unwrap();
        let pretty = found_events
            .iter()
            .fold(String::new(), |mut output, orbit_event| {
                let _ = writeln!(
                    output,
                    "{:x}\tevent value: {}",
                    orbit_event.state, orbit_event.value
                );
                output
            });
        println!("[ta_tracker] {} =>\n{}", e, pretty);
    }

    // Find all eclipses!
    let e_loc = EclipseLocator {
        light_source: SUN_J2000,
        shadow_bodies: vec![eme2k],
    };

    // Adding this print to confirm that the penumbra calculation continuously increases and then decreases.
    let mut e_state = EclipseState::Umbra;
    // Also see what is the max elevation of this spacecraft over the Grand Canyon
    let gc = GroundStation::from_point(
        "Grand Canyon".to_string(),
        36.0544,
        112.1402,
        0.0,
        IAU_EARTH_FRAME,
    );
    let mut min_el = std::f64::INFINITY;
    let mut max_el = std::f64::NEG_INFINITY;
    let mut min_dt = dt;
    let mut max_dt = dt;
    for state in traj.every(10 * Unit::Second) {
        let new_e_state = e_loc.compute(state.orbit, almanac.clone()).unwrap();
        if e_state != new_e_state {
            println!("{:x}\t{}", state, new_e_state);
            e_state = new_e_state;
        }

        // Compute the elevation
        let aer = gc.azimuth_elevation_of(state.orbit, &almanac).unwrap();
        let elevation = aer.elevation_deg;
        if elevation > max_el {
            max_el = elevation;
            max_dt = state.epoch();
        }

        if elevation < min_el {
            min_el = elevation;
            min_dt = state.epoch();
        }
    }

    println!("Min elevation {} degrees @ {}", min_el, min_dt);
    println!("Max elevation {} degrees @ {}", max_el, max_dt);

    let umbra_event_loc = e_loc.to_umbra_event();
    let umbra_events = traj.find(&umbra_event_loc, almanac.clone()).unwrap();

    let pretty = umbra_events
        .iter()
        .fold(String::new(), |mut output, orbit_event| {
            let orbit = orbit_event.state.orbit;
            let _ = writeln!(
                output,
                "{:x}\tevent value: {}\t(-10s: {}\t+10s: {})",
                orbit,
                &e_loc.compute(orbit, almanac.clone()).unwrap(),
                &e_loc
                    .compute(
                        traj.at(orbit.epoch - 10 * Unit::Second).unwrap().orbit,
                        almanac.clone()
                    )
                    .unwrap(),
                &e_loc
                    .compute(
                        traj.at(orbit.epoch + 10 * Unit::Second).unwrap().orbit,
                        almanac.clone()
                    )
                    .unwrap()
            );
            output
        });
    println!("[eclipses] {} =>\n{}", umbra_event_loc, pretty);

    let penumbra_event_loc = e_loc.to_penumbra_event();
    let penumbra_events = traj.find(&penumbra_event_loc, almanac.clone()).unwrap();

    let pretty = penumbra_events
        .iter()
        .fold(String::new(), |mut output, orbit_event| {
            let orbit = orbit_event.state.orbit;
            let _ = writeln!(
                output,
                "{:x}\tevent value: {}\t(-10s: {}\t+10s: {})",
                orbit,
                &e_loc.compute(orbit, almanac.clone()).unwrap(),
                &e_loc
                    .compute(
                        traj.at(orbit.epoch - 10 * Unit::Second).unwrap().orbit,
                        almanac.clone()
                    )
                    .unwrap(),
                &e_loc
                    .compute(
                        traj.at(orbit.epoch + 10 * Unit::Second).unwrap().orbit,
                        almanac.clone()
                    )
                    .unwrap()
            );
            output
        });
    println!("[eclipses] {} =>\n{}", penumbra_event_loc, pretty);
}
