extern crate hifitime;
extern crate nalgebra as na;

extern crate nyx_space as nyx;

#[test]
fn event_tracker_true_anomaly() {
    use hifitime::{Epoch, J2000_OFFSET};
    use nyx::celestia::{bodies, Cosm, OrbitState};
    use nyx::dynamics::celestial::CelestialDynamics;
    use nyx::propagators::error_ctrl::RSSStepPV;
    use nyx::propagators::events::{EventKind, EventTrackers, OrbitalEvent};
    use nyx::propagators::*;

    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.frame_by_id(bodies::EARTH);

    let dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let state = OrbitState::cartesian(
        -2436.45,
        -2436.45,
        6891.037,
        5.088_611,
        -5.088_611,
        0.0,
        dt,
        earth_geoid,
    );

    let prop_time = state.period();

    // Track how many times we've passed by that TA again
    let peri_event = OrbitalEvent::new(EventKind::Periapse);
    let apo_event = OrbitalEvent::new(EventKind::Apoapse);
    let ta_event0 = OrbitalEvent::new(EventKind::TA(35.1));
    let ta_event1 = OrbitalEvent::new(EventKind::TA(235.1));

    let tracker = EventTrackers::from_events(vec![peri_event, apo_event, ta_event0, ta_event1]);

    let mut dynamics = CelestialDynamics::two_body(state);

    let mut prop = Propagator::default(
        &mut dynamics,
        &PropOpts::with_adaptive_step(1.0, 60.0, 1e-9, RSSStepPV {}),
    );
    prop.event_trackers = tracker;
    prop.until_time_elapsed(prop_time);

    // Check how many times we have found that event
    println!("{}", prop.event_trackers);
}
