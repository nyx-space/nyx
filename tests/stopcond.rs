extern crate nalgebra as na;

extern crate nyx_space as nyx;

use nyx::celestia::{bodies, Cosm, Orbit};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::propagators::error_ctrl::RSSStepPV;
use nyx::propagators::events::{EventKind, OrbitalEvent, StopCondition};
use nyx::propagators::{PropOpts, Propagator};
use nyx::time::{Epoch, TimeUnit, J2000_OFFSET};

#[ignore]
#[test]
fn stop_cond_3rd_apo() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );

    let period = state.period();

    // Track how many times we've passed by that TA again
    let apo_event = OrbitalEvent::new(EventKind::Apoapse);
    let condition = StopCondition::after_hits(apo_event, 3, 4.0 * period, 1e-10);

    let dynamics = OrbitalDynamics::two_body();

    let setup = Propagator::rk89(
        &dynamics,
        PropOpts::with_adaptive_step_s(1.0, 60.0, 1e-9, RSSStepPV {}),
    );

    let mut prop = setup.with(state);

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
        (180.0 - orbit.ta()).abs() < 1e-3,
        "converged, yet convergence critera not met"
    );
}

#[ignore]
#[test]
fn stop_cond_3rd_peri() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );

    let period = state.period();

    // Track how many times we've passed by that TA again
    let apo_event = OrbitalEvent::new(EventKind::Periapse);
    let condition = StopCondition::after_hits(apo_event, 3, 4.0 * period, 1e-10);

    let dynamics = OrbitalDynamics::two_body();

    let setup = Propagator::rk89(
        &dynamics,
        PropOpts::with_adaptive_step_s(1.0, 60.0, 1e-9, RSSStepPV {}),
    );

    let mut prop = setup.with(state);

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
        orbit.ta().abs() < 1e-1 || (360.0 - orbit.ta().abs() < 1e-1),
        "converged, yet convergence critera not met"
    );
}

#[ignore]
#[test]
fn nrho_apo() {
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

    // Track how many times we've passed by that TA again
    // let apo_event = OrbitalEvent::in_frame(EventKind::Apoapse, luna, &cosm);
    let apo_event = OrbitalEvent::new(EventKind::Apoapse);
    let condition = StopCondition::new(apo_event, 2 * TimeUnit::Day, 1e-1);

    let bodies = vec![bodies::EARTH, bodies::SUN];
    let dynamics = OrbitalDynamics::point_masses(luna, &bodies, &cosm);

    let setup = Propagator::rk89(
        &dynamics,
        PropOpts::with_adaptive_step_s(1.0, 60.0, 1e-9, RSSStepPV {}),
    );

    let mut prop = setup.with(state_luna);

    let rslt = prop.until_event(condition);

    // Check how many times we have found that event
    let orbit = rslt.expect("condition should have been found");
    println!("Luna: {:o}", orbit);
    assert!((orbit.ta() - 180.0).abs() < 0.1);
    let rslt_eme = cosm.frame_chg(&orbit, eme2k);
    println!("EME2k: {}", rslt_eme);
    let delta_t = orbit.dt - dt;
    println!("Found {}days after", delta_t);
}
