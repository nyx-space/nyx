extern crate hifitime;
extern crate nalgebra as na;

extern crate nyx_space as nyx;

#[test]
fn stop_cond_3rd_apo() {
    use hifitime::{Epoch, J2000_OFFSET};
    use nyx::celestia::{bodies, Cosm, Geoid, State};
    use nyx::dynamics::celestial::CelestialDynamics;
    use nyx::propagators::error_ctrl::RSSStepPV;
    use nyx::propagators::events::{EventKind, OrbitalEvent, StopCondition};
    use nyx::propagators::*;

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
