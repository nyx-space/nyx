extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::hifitime::Epoch;
use self::na::Vector3;
use self::nyx::celestia::{bodies, Cosm, OrbitState};
use self::nyx::dynamics::celestial::CelestialDynamics;
use self::nyx::dynamics::deltavctrl::{InstantBurns, Mnvr};
use self::nyx::dynamics::missionarc::MissionArc;
use self::nyx::dynamics::Dynamics;
use self::nyx::propagators::{PropOpts, Propagator};

#[test]
fn arc_example() {
    // This is an example of two delta Vs happening one after another. Other implementations of
    // a DeltaVctrl would be much more useful than this schedule of maneuvers.
    let cosm = Cosm::from_xb("./de438s");
    let earth = cosm.geoid_from_id(bodies::EARTH);

    let start_time = Epoch::from_gregorian_tai_at_midnight(2002, 1, 1);

    let orbit = OrbitState::keplerian(6678.0, 0.0, 0.1, 60.0, 30.0, 0.0, start_time, earth);

    let prop_time = 86_400.0;

    let mut mnvr1_dt = start_time;
    mnvr1_dt.mut_add_secs(2.0 * 3600.0);

    let mut end_time = start_time;
    end_time.mut_add_secs(prop_time);

    // Define the dynamics
    let celestial = CelestialDynamics::two_body(orbit);

    let mnvr0 = Mnvr::instantaneous(start_time, Vector3::new(2.42, 0.0, 0.0));
    let mnvr1 = Mnvr::instantaneous(mnvr1_dt, Vector3::new(-1.46, 0.0, 0.0));

    let schedule = InstantBurns::from_mnvrs(vec![mnvr0, mnvr1]);

    let mut arc = MissionArc::new(celestial, schedule);
    let mut prop = Propagator::default(&mut arc, &PropOpts::with_fixed_step(10.0));
    prop.until_time_elapsed(prop_time);

    println!("final state: {:o}", prop.dynamics.celestial.state());
}
