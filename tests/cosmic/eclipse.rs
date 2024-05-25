extern crate nyx_space as nyx;

use anise::constants::celestial_objects::{JUPITER, SUN};
use anise::constants::frames::SUN_J2000;
use nyx::cosmic::eclipse::{EclipseLocator, EclipseState};
use nyx::cosmic::Orbit;
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::propagators::{PropOpts, Propagator};
use nyx::time::{Epoch, Unit};
use std::sync::{mpsc, Arc};
use std::thread;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn leo_sun_earth_eclipses(almanac: Arc<Almanac>) {
    let prop_time = 2.0 * Unit::Day;

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let leo = Orbit::keplerian(6778.0, 0.1, 60.0, 0.0, 0.0, 0.0, start_time, eme2k);

    let (truth_tx, truth_rx) = mpsc::channel();

    let bodies = vec![SUN, JUPITER];

    let almanac_c = almanac.clone();
    thread::spawn(move || {
        let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
        let setup = Propagator::rk89(dynamics, PropOpts::with_fixed_step_s(60.0));

        setup
            .with(leo.into(), almanac_c)
            .for_duration_with_channel(prop_time, truth_tx)
            .unwrap();
    });

    // Initialize the EclipseLocator
    let e_loc = EclipseLocator {
        light_source: SUN_J2000,
        shadow_bodies: vec![eme2k],
    };

    // Receive the states on the main thread.
    let mut prev_eclipse_state = EclipseState::Umbra;
    let mut cnt_changes = 0;
    while let Ok(rx_state) = truth_rx.recv() {
        let new_eclipse_state = e_loc.compute(rx_state.orbit, almanac.clone()).unwrap();
        if new_eclipse_state != prev_eclipse_state {
            println!("{:.6} now in {:?}", rx_state.orbit.epoch, new_eclipse_state);
            prev_eclipse_state = new_eclipse_state;
            cnt_changes += 1;
        }
    }

    assert_eq!(cnt_changes, 68, "wrong number of eclipse state changes");
}

#[rstest]
fn geo_sun_earth_eclipses(almanac: Arc<Almanac>) {
    let prop_time = 2 * Unit::Day;

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    // GEO are in shadow or near shadow during the equinoxes.
    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 3, 19);

    let geo = Orbit::keplerian(42000.0, 0.1, 0.1, 0.0, 0.0, 0.0, start_time, eme2k);

    let (truth_tx, truth_rx) = mpsc::channel();

    let bodies = vec![SUN, JUPITER];

    let almanac_c = almanac.clone();

    thread::spawn(move || {
        let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
        let setup = Propagator::rk89(dynamics, PropOpts::with_fixed_step_s(60.0));

        setup
            .with(geo.into(), almanac_c)
            .for_duration_with_channel(prop_time, truth_tx)
            .unwrap();
    });

    // Initialize the EclipseLocator
    let e_loc = EclipseLocator {
        light_source: SUN_J2000,
        shadow_bodies: vec![eme2k],
    };

    // Receive the states on the main thread.
    let mut prev_eclipse_state = EclipseState::Umbra;
    let mut cnt_changes = 0;
    while let Ok(rx_state) = truth_rx.recv() {
        let new_eclipse_state = e_loc.compute(rx_state.orbit, almanac.clone()).unwrap();
        if new_eclipse_state != prev_eclipse_state {
            println!("{:.6} now in {:?}", rx_state.orbit.epoch, new_eclipse_state);
            prev_eclipse_state = new_eclipse_state;
            cnt_changes += 1;
        }
    }

    assert_eq!(cnt_changes, 15, "wrong number of eclipse state changes");
}
