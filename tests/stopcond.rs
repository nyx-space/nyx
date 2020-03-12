extern crate hifitime;
extern crate nalgebra as na;

extern crate nyx_space as nyx;

use hifitime::{Epoch, J2000_OFFSET};
use nyx::celestia::{bodies, Cosm, OrbitState};
use nyx::dynamics::celestial::CelestialDynamics;
use nyx::propagators::error_ctrl::RSSStepPV;
use nyx::propagators::events::{EventKind, OrbitalEvent, StopCondition};
use nyx::propagators::{PropOpts, Propagator};

#[test]
fn stop_cond_3rd_apo() {
    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(bodies::EARTH);

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

    let period = state.period();

    // Track how many times we've passed by that TA again
    let apo_event = OrbitalEvent::new(EventKind::Apoapse);
    let condition = StopCondition::after_hits(apo_event, 3, 4.0 * period, 1e-6);

    let mut dynamics = CelestialDynamics::two_body(state);

    let mut prop = Propagator::default(
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
    let state = OrbitState::cartesian(
        166_473.631_302_239_7,
        -274_715.487_253_382_7,
        -211_233.210_176_686_7,
        0.933_451_604_520_018_4,
        0.436_775_046_841_900_9,
        -0.082_211_021_250_348_95,
        dt,
        earth,
    );

    // Note that this expected state was generated using SRP and a lunar gravity field
    // Hence, we allow for a greater error since these are not modeled here.
    let expect = OrbitState::cartesian(
        266_375.578_868_798,
        -213_365.467_957_944,
        -203_571.279_542_228,
        0.741_790_420_281_572,
        0.588_200_782_187_968,
        0.202_695_184_897_909,
        dt + 118_753.007_910_251,
        earth,
    );

    let state_luna = cosm.frame_chg(&state, luna);

    // Track how many times we've passed by that TA again
    // let apo_event = OrbitalEvent::in_frame(EventKind::Apoapse, luna, &cosm);
    let apo_event = OrbitalEvent::new(EventKind::Apoapse);
    let condition = StopCondition::new(apo_event, 2.0 * 86_400.0, 1e-1);

    let mut dynamics = CelestialDynamics::new(state_luna, vec![bodies::EARTH, bodies::SUN], &cosm);

    let mut prop = Propagator::default(
        &mut dynamics,
        &PropOpts::with_adaptive_step(1.0, 60.0, 1e-9, RSSStepPV {}),
    );

    let rslt = prop.until_event(condition);

    // Check how many times we have found that event
    let orbit = rslt.expect("condition should have been found");
    println!("Luna: {:o}", orbit);
    let rslt_eme = cosm.frame_chg(&orbit, earth);
    println!("EME2k: {}", rslt_eme);
    assert!((rslt_eme - expect).rmag() < 1.0);
    assert!((rslt_eme - expect).vmag() < 1e-5);
    let delta_t = orbit.dt - dt;
    println!(
        "Found {} seconds / {} days after",
        delta_t,
        delta_t / 86_400.0
    );
}
