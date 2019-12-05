extern crate hifitime;
extern crate nalgebra as na;

extern crate nyx_space as nyx;

#[test]
fn event_tracker_true_anomaly() {
    use hifitime::{Epoch, J2000_OFFSET};
    use nyx::celestia::{bodies, Cosm, Geoid, State};
    use nyx::dynamics::celestial::CelestialDynamics;
    use nyx::propagators::error_ctrl::RSSStepPV;
    use nyx::propagators::events::{EventTrackers, StateEvent, StateEventKind};
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

    let prop_time = 5.0 * state.period();

    // Track how many times we've passed by that TA again
    let ta_event = StateEvent {
        kind: StateEventKind::TA(0.0),
    };

    let tracker = EventTrackers::from_event(Box::new(ta_event));

    let mut dynamics = CelestialDynamics::two_body(state);

    let mut prop = Propagator::new::<RK89>(
        &mut dynamics,
        &PropOpts::with_adaptive_step(1.0, 60.0, 1e-9, RSSStepPV {}),
    );
    prop.event_trackers = tracker;
    prop.until_time_elapsed(prop_time);

    // Check how many times we have found that event
    println!("{:?}", prop.event_trackers);
}
