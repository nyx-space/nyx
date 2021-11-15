extern crate nyx_space as nyx;

use nyx::cosmic::eclipse::{EclipseLocator, EclipseState};
use nyx::cosmic::{Bodies, Cosm, Orbit};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::propagators::{PropOpts, Propagator};
use nyx::time::{Epoch, TimeUnit};
use std::sync::mpsc;
use std::thread;

#[test]
fn leo_sun_earth_eclipses() {
    let prop_time = 2.0 * TimeUnit::Day;

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let leo = Orbit::keplerian(6778.0, 0.1, 60.0, 0.0, 0.0, 0.0, start_time, eme2k);

    let (truth_tx, truth_rx) = mpsc::channel();

    let bodies = vec![Bodies::Sun, Bodies::JupiterBarycenter];

    let cosmc = cosm.clone();
    thread::spawn(move || {
        let dynamics = OrbitalDynamics::point_masses(&bodies, cosmc);
        let setup = Propagator::rk89(dynamics, PropOpts::with_fixed_step_s(60.0));

        setup
            .with(leo)
            .for_duration_with_channel(prop_time, truth_tx)
            .unwrap();
    });

    // Initialize the EclipseLocator
    let e_loc = EclipseLocator {
        light_source: cosm.frame("Sun J2000"),
        shadow_bodies: vec![eme2k],
        cosm,
    };

    // Receive the states on the main thread.
    let mut prev_eclipse_state = EclipseState::Umbra;
    let mut cnt_changes = 0;
    while let Ok(rx_state) = truth_rx.recv() {
        let new_eclipse_state = e_loc.compute(&rx_state);
        if new_eclipse_state != prev_eclipse_state {
            println!(
                "{:.6} now in {:?}",
                rx_state.dt.as_jde_tai_days(),
                new_eclipse_state
            );
            prev_eclipse_state = new_eclipse_state;
            cnt_changes += 1;
        }
    }

    assert_eq!(cnt_changes, 68, "wrong number of eclipse state changes");
}

#[test]
fn geo_sun_earth_eclipses() {
    let prop_time = 2 * TimeUnit::Day;

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    // GEO are in shadow or near shadow during the equinoxes.
    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 3, 19);

    let geo = Orbit::keplerian(42000.0, 0.1, 0.1, 0.0, 0.0, 0.0, start_time, eme2k);

    let (truth_tx, truth_rx) = mpsc::channel();

    let bodies = vec![Bodies::Sun, Bodies::JupiterBarycenter];

    thread::spawn(move || {
        let cosm = Cosm::de438();
        let dynamics = OrbitalDynamics::point_masses(&bodies, cosm);
        let setup = Propagator::rk89(dynamics, PropOpts::with_fixed_step_s(60.0));

        setup
            .with(geo)
            .for_duration_with_channel(prop_time, truth_tx)
            .unwrap();
    });

    // Initialize the EclipseLocator
    let e_loc = EclipseLocator {
        light_source: cosm.frame("Sun J2000"),
        shadow_bodies: vec![eme2k],
        cosm,
    };

    // Receive the states on the main thread.
    let mut prev_eclipse_state = EclipseState::Umbra;
    let mut cnt_changes = 0;
    while let Ok(rx_state) = truth_rx.recv() {
        let new_eclipse_state = e_loc.compute(&rx_state);
        if new_eclipse_state != prev_eclipse_state {
            println!(
                "{:.6} now in {:?}",
                rx_state.dt.as_jde_tai_days(),
                new_eclipse_state
            );
            prev_eclipse_state = new_eclipse_state;
            cnt_changes += 1;
        }
    }

    assert_eq!(cnt_changes, 15, "wrong number of eclipse state changes");
}
