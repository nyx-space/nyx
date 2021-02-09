extern crate nyx_space as nyx;

use nyx::celestia::{Cosm, Orbit};
use nyx::dynamics::OrbitalDynamics;
use nyx::md::{Ephemeris, ScTraj};
use nyx::propagators::*;
use nyx::time::{Epoch, TimeSeries, TimeUnit};
use nyx::State;
use std::sync::mpsc::channel;

#[test]
fn traj_ephem() {
    let (tx, rx) = channel();
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_gregorian_utc_at_noon(2021, 1, 1);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );

    let ephem_thread = std::thread::spawn(move || Ephemeris::new(state, rx));

    println!("Init: {}", state);

    let dynamics = OrbitalDynamics::two_body();
    let setup = Propagator::rk89(dynamics, PropOpts::<RSSStepPV>::default());
    let mut prop = setup.with(state).with_tx(tx);
    prop.for_duration(31 * TimeUnit::Day).unwrap();

    // Retrieve the ephemeris by unwrapping it twice:
    // + The first is to unwrap the thread (i.e. assume the thread has not failed)
    // + The second assumes that the generation of the ephemeris didn't fail.
    let ephem = ephem_thread.join().unwrap().unwrap();

    assert_eq!(
        ephem.segments.len(),
        1009,
        "Wrong number of expected segments"
    );

    assert_eq!(ephem.start_state, state, "Wrong initial state");

    // Now let's re-generate the truth data and ensure that each state we generate is in the ephemeris and matches the expected state within tolerance.

    let (tx, rx) = channel();
    std::thread::spawn(move || {
        let mut prop = setup.with(state).with_tx(tx);
        prop.for_duration(31 * TimeUnit::Day).unwrap();
    });

    let mut max_err = std::f64::NEG_INFINITY;
    let mut max_pos_err = std::f64::NEG_INFINITY;
    let mut max_vel_err = std::f64::NEG_INFINITY;

    while let Ok(prop_state) = rx.recv() {
        match ephem.evaluate(prop_state.dt) {
            Ok(eval_state) => {
                let pos_err = (eval_state.radius() - prop_state.radius()).norm();
                if pos_err > max_pos_err {
                    max_pos_err = pos_err;
                }
                let vel_err = (eval_state.velocity() - prop_state.velocity()).norm();
                if vel_err > max_vel_err {
                    max_vel_err = vel_err;
                }
                let err =
                    (eval_state.as_vector().unwrap() - prop_state.as_vector().unwrap()).norm();
                if err > max_err {
                    max_err = err;
                }
            }
            Err(e) => println!("Err: {}", e),
        }
    }

    println!(
        "[traj_ephem] Maximum interpolation error: pos: {:.2e} m\t\tvel: {:.2e} m/s\t\tfull state: {:.2e} (no unit)",
        max_pos_err * 1e3,
        max_vel_err * 1e3,
        max_err
    );

    assert!(
        max_err < 1e-10,
        "Maximum orbit in interpolation is too high!"
    );

    // Print a bunch of states throughout the orbit
    // for epoch in TimeSeries::inclusive(dt, dt + 31 * TimeUnit::Day, 1 * TimeUnit::Day) {
    //     match ephem.evaluate(epoch) {
    //         Ok(state) => println!("{}", state),
    //         Err(e) => println!("{} -- Oh no {}", epoch, e),
    //     };
    // }
}
