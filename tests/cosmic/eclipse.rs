extern crate nyx_space as nyx;

use anise::astro::Occultation;
use anise::constants::celestial_objects::{JUPITER_BARYCENTER, SUN};
use anise::constants::frames::SUN_J2000;
use nyx::cosmic::eclipse::EclipseLocator;
use nyx::cosmic::Orbit;
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::propagators::{IntegratorOptions, Propagator};
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

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let leo = Orbit::keplerian(6778.0, 0.1, 60.0, 0.0, 0.0, 0.0, start_time, eme2k);

    let (truth_tx, truth_rx) = mpsc::channel();

    let bodies = vec![SUN, JUPITER_BARYCENTER];

    let almanac_c = almanac.clone();
    thread::spawn(move || {
        let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
        let setup = Propagator::rk89(dynamics, IntegratorOptions::with_fixed_step_s(60.0));

        setup
            .with(leo.into(), almanac_c)
            .for_duration_with_channel(prop_time, truth_tx)
            .unwrap();
    });

    // Initialize the EclipseLocator
    let e_loc = EclipseLocator {
        light_source: almanac.frame_info(SUN_J2000).unwrap(),
        shadow_bodies: vec![eme2k],
    };

    // Receive the states on the main thread.
    let mut prev_eclipse_state: Option<Occultation> = None;
    let mut cnt_changes = 0;
    while let Ok(rx_state) = truth_rx.recv() {
        let new_eclipse_state = e_loc.compute(rx_state.orbit, almanac.clone()).unwrap();
        if let Some(prev_state) = prev_eclipse_state {
            if new_eclipse_state.percentage != prev_state.percentage {
                println!("{:.6} now in {}", rx_state.orbit.epoch, new_eclipse_state);
                prev_eclipse_state = Some(new_eclipse_state);
                cnt_changes += 1;
            }
        } else {
            prev_eclipse_state = Some(new_eclipse_state);
        }
    }

    assert_eq!(cnt_changes, 68, "wrong number of eclipse state changes");
}

#[rstest]
fn geo_sun_earth_eclipses(almanac: Arc<Almanac>) {
    let prop_time = 2 * Unit::Day;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    // GEO are in shadow or near shadow during the equinoxes.
    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 3, 19);

    let geo = Orbit::keplerian(42000.0, 0.1, 0.1, 0.0, 0.0, 0.0, start_time, eme2k);

    let (truth_tx, truth_rx) = mpsc::channel();

    let bodies = vec![SUN, JUPITER_BARYCENTER];

    let almanac_c = almanac.clone();

    thread::spawn(move || {
        let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
        let setup = Propagator::rk89(dynamics, IntegratorOptions::with_fixed_step_s(60.0));

        setup
            .with(geo.into(), almanac_c)
            .for_duration_with_channel(prop_time, truth_tx)
            .unwrap();
    });

    // Initialize the EclipseLocator
    let e_loc = EclipseLocator {
        light_source: almanac.frame_info(SUN_J2000).unwrap(),
        shadow_bodies: vec![eme2k],
    };

    // Receive the states on the main thread.
    let mut prev_eclipse_state: Option<Occultation> = None;
    let mut cnt_changes = 0;
    while let Ok(rx_state) = truth_rx.recv() {
        let new_eclipse_state = e_loc.compute(rx_state.orbit, almanac.clone()).unwrap();
        if let Some(prev_state) = prev_eclipse_state {
            if new_eclipse_state.percentage != prev_state.percentage {
                println!("{:.6} now in {}", rx_state.orbit.epoch, new_eclipse_state);
                prev_eclipse_state = Some(new_eclipse_state);
                cnt_changes += 1;
            }
        } else {
            prev_eclipse_state = Some(new_eclipse_state);
        }
    }

    assert_eq!(cnt_changes, 14, "wrong number of eclipse state changes");
}
