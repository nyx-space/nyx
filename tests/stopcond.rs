extern crate hifitime;
extern crate nalgebra as na;

extern crate nyx_space as nyx;

use hifitime::{Epoch, J2000_OFFSET};
use nyx::celestia::{bodies, Cosm, Geoid, State};
use nyx::dynamics::celestial::CelestialDynamics;
use nyx::propagators::error_ctrl::RSSStepPV;
use nyx::propagators::events::{EventKind, OrbitalEvent, StopCondition};
use nyx::propagators::{PropOpts, Propagator, RK89};

#[test]
fn stop_cond_3rd_apo() {
    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(bodies::EARTH);

    let dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let state = State::<Geoid>::from_cartesian(
        -2436.45,
        -2436.45,
        6891.037,
        5.088_611,
        -5.088_611,
        0.0,
        dt,
        earth_geoid,
    );

    let period = state.period();

    // Track how many times we've passed by that TA again
    let apo_event = OrbitalEvent::new(EventKind::Apoapse);
    let condition = StopCondition::after_hits(apo_event, 3, 4.0 * period, 1e-6);

    let mut dynamics = CelestialDynamics::two_body(state);

    let mut prop = Propagator::new::<RK89>(
        &mut dynamics,
        &PropOpts::with_adaptive_step(1.0, 60.0, 1e-9, RSSStepPV {}),
    );

    let rslt = prop.until_event(condition);

    // Check how many times we have found that event
    println!("{}", prop.event_trackers);
    let orbit = rslt.expect("condition should have been found");
    println!("{:o}", orbit);
    // Confirm that this is the third apoapse event which is found
    assert!(
        orbit.dt - dt < 3.0 * period && orbit.dt - dt >= 2.0 * period,
        "converged on the wrong apoapse"
    );
    assert!(
        (180.0 - orbit.ta()) < 1e-6,
        "converged, yet convergence critera not met"
    );
}

#[test]
fn nrho_apo() {
    let cosm = Cosm::from_xb("./de438s");
    let earth = cosm.geoid_from_id(bodies::EARTH);
    let luna = cosm.geoid_from_id(bodies::EARTH_MOON);

    let dt = Epoch::from_gregorian_tai(2021, 5, 29, 19, 51, 16, 852_000);
    let state = State::<Geoid>::from_cartesian(
        166_473.631_302_239_7,
        -274_715.487_253_382_7,
        -211_233.210_176_686_7,
        0.933_451_604_520_018_4,
        0.436_775_046_841_900_9,
        -0.082_211_021_250_348_95,
        dt,
        earth,
    );

    let state_luna = cosm.frame_chg(&state, luna);

    // Track how many times we've passed by that TA again
    // let apo_event = OrbitalEvent::in_frame(EventKind::Apoapse, luna, &cosm);
    let apo_event = OrbitalEvent::new(EventKind::Apoapse);
    let condition = StopCondition::new(apo_event, 2.0 * 86_400.0, 1e-1);

    let mut dynamics = CelestialDynamics::new(state_luna, vec![bodies::EARTH], &cosm);

    let mut prop = Propagator::new::<RK89>(
        &mut dynamics,
        &PropOpts::with_adaptive_step(1.0, 60.0, 1e-9, RSSStepPV {}),
    );

    let rslt = prop.until_event(condition);

    // Check how many times we have found that event
    println!("{}", prop.event_trackers);
    let orbit = rslt.expect("condition should have been found");
    println!("Luna: {:o}", orbit);
    let rslt_eme = cosm.frame_chg(&orbit, earth);
    println!("EME2k: {:o}", rslt_eme);
    println!("EME2k: {}", rslt_eme);
    let delta_t = orbit.dt - dt;
    println!(
        "Found {} seconds / {} days after",
        delta_t,
        delta_t / 86_400.0
    );
}
