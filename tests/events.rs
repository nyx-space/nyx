extern crate nyx_space as nyx;

#[test]
fn event_tracker_true_anomaly() {
    use nyx::celestia::{Cosm, Orbit};
    use nyx::dynamics::orbital::OrbitalDynamics;
    use nyx::md::events::{Event, StateParameter};
    use nyx::propagators::error_ctrl::RSSStepPV;
    use nyx::propagators::*;
    use nyx::time::{Epoch, J2000_OFFSET};

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );

    let prop_time = state.period();

    // Track how many times we've passed by that TA again
    let peri_event = Event::new(StateParameter::Periapsis, 0.0);
    let apo_event = Event::new(StateParameter::Apoapsis, 0.0);
    let ta_event0 = Event::new(StateParameter::TrueAnomaly, 35.1);
    let ta_event1 = Event::new(StateParameter::TrueAnomaly, 235.1);

    let events = vec![peri_event, apo_event, ta_event0, ta_event1];

    let dynamics = OrbitalDynamics::two_body();
    let setup = Propagator::rk89(
        dynamics,
        PropOpts::with_adaptive_step_s(1.0, 60.0, 1e-9, RSSStepPV {}),
    );
    let mut prop = setup.with(state);
    let (_, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Find all of the events
    for e in &events {
        println!("[ta_tracker] {} => {:?}", e, traj.find_all(e));
    }
}
