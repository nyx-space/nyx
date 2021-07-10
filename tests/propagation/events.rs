extern crate nyx_space as nyx;

#[test]
fn event_tracker_true_anomaly() {
    use nyx::cosmic::eclipse::{EclipseLocator, EclipseState};
    use nyx::md::ui::*;
    use nyx::md::EventEvaluator; // Only needed because we're manually calling e.eval
    use nyx::od::measurement::GroundStation;

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_gregorian_tai_at_noon(2020, 1, 1);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );

    let prop_time = state.period();

    // Track how many times we've passed by that TA again
    let peri_event = Event::periapsis(); // Special event shortcut!
    let apo_event = Event::apoapsis(); // Special event shortcut!
    let ta_event0 = Event::new(StateParameter::TrueAnomaly, 35.1);
    let ta_event1 = Event::new(StateParameter::TrueAnomaly, 235.1);

    let events = vec![peri_event, apo_event, ta_event0, ta_event1];

    let dynamics = OrbitalDynamics::two_body();
    let setup = Propagator::rk89(dynamics, PropOpts::with_tolerance(1e-9));
    let mut prop = setup.with(state);
    let (_, traj) = prop.for_duration_with_traj(prop_time).unwrap();

    // Find all of the events
    for e in &events {
        let found_events = traj.find_all(e).unwrap();
        let pretty = found_events
            .iter()
            .map(|orbit| format!("{:x}\tevent value: {}\n", orbit, e.eval(orbit)))
            .collect::<String>();
        println!("[ta_tracker] {} =>\n{}", e, pretty);
    }

    // Find all eclipses!
    let e_loc = EclipseLocator {
        light_source: cosm.frame("Sun J2000"),
        shadow_bodies: vec![cosm.frame("EME2000")],
        cosm: cosm.clone(),
    };

    // Adding this print to confirm that the penumbra calculation continuously increases and then decreases.
    let mut e_state = EclipseState::Umbra;
    // Also see what is the max elevation of this spacecraft over the Grand Canyon
    let gc = GroundStation::from_point(
        "Grand Canyon".to_string(),
        36.0544,
        112.1402,
        0.0,
        cosm.frame("IAU Earth"),
        cosm.clone(),
    );
    let mut min_el = std::f64::INFINITY;
    let mut max_el = std::f64::NEG_INFINITY;
    let mut min_dt = dt;
    let mut max_dt = dt;
    for state in traj.every(10 * TimeUnit::Second) {
        let new_e_state = e_loc.compute(&state);
        if e_state != new_e_state {
            println!("{:x}\t{}", state, new_e_state);
            e_state = new_e_state;
        }

        // Compute the elevation
        let (elevation, _, _) = gc.elevation_of(&state);
        if elevation > max_el {
            max_el = elevation;
            max_dt = state.epoch();
        }

        if elevation < min_el {
            min_el = elevation;
            min_dt = state.epoch();
        }
    }

    println!("Min elevation {} degrees @ {}", min_el, min_dt);
    println!("Max elevation {} degrees @ {}", max_el, max_dt);

    let umbra_event_loc = e_loc.to_umbra_event();
    let umbra_events = traj.find_all(&umbra_event_loc).unwrap();

    let pretty = umbra_events
        .iter()
        .map(|orbit| {
            format!(
                "{:x}\tevent value: {}\t(-10s: {}\t+10s: {})\n",
                orbit,
                &e_loc.compute(orbit),
                &e_loc.compute(&traj.at(orbit.epoch() - 10 * TimeUnit::Second).unwrap()),
                &e_loc.compute(&traj.at(orbit.epoch() + 10 * TimeUnit::Second).unwrap())
            )
        })
        .collect::<String>();
    println!("[eclipses] {} =>\n{}", umbra_event_loc, pretty);

    let penumbra_event_loc = e_loc.to_penumbra_event();
    let penumbra_events = traj.find_all(&penumbra_event_loc).unwrap();

    let pretty = penumbra_events
        .iter()
        .map(|orbit| {
            format!(
                "{:x}\tevent value: {}\t(-10s: {}\t+10s: {})\n",
                orbit,
                &e_loc.compute(orbit),
                &e_loc.compute(&traj.at(orbit.epoch() - 10 * TimeUnit::Second).unwrap()),
                &e_loc.compute(&traj.at(orbit.epoch() + 10 * TimeUnit::Second).unwrap())
            )
        })
        .collect::<String>();
    println!("[eclipses] {} =>\n{}", penumbra_event_loc, pretty);
}
