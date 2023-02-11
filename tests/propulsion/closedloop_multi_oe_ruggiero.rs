extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::nyx::cosmic::{Cosm, GuidanceMode, Orbit, Spacecraft};
use self::nyx::dynamics::guidance::{Objective, Ruggiero, Thruster};
use self::nyx::dynamics::{OrbitalDynamics, SpacecraftDynamics};
use self::nyx::md::{Event, StateParameter};
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use self::nyx::time::{Epoch, Unit};

/// NOTE: Herein shows the difference between the QLaw and Ruggiero (and other control laws).
/// The Ruggiero control law takes quite some longer to converge than the QLaw.

#[test]
fn qlaw_as_ruggiero_case_a() {
    // Source: AAS-2004-5089

    let mut cosm = Cosm::de438_raw();
    cosm.frame_mut_gm("EME2000", 398_600.433);
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(7000.0, 0.01, 0.05, 0.0, 0.0, 1.0, start_time, eme2k);

    let prop_time = 39.91 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 1.0,
        isp_s: 3100.0,
    };

    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 42_000.0, 1.0),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.01, 5e-5),
    ];

    // Events we will search later
    let events = vec![
        Event::within_tolerance(StateParameter::SMA, 42_000.0, 1.0),
        Event::within_tolerance(StateParameter::Eccentricity, 0.01, 5e-5),
    ];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit).unwrap();

    let dry_mass = 1.0;
    let fuel_mass = 299.0;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[qlaw_as_ruggiero_case_a] {:x}", orbit);

    let setup =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second));
    let mut prop = setup.with(sc_state);
    let (final_state, traj) = prop.for_duration_with_traj(prop_time).unwrap();
    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[qlaw_as_ruggiero_case_a] {:x}", final_state.orbit);
    println!("[qlaw_as_ruggiero_case_a] fuel usage: {:.3} kg", fuel_usage);
    // Find all of the events
    for e in &events {
        println!(
            "[qlaw_as_ruggiero_case_a] Found {} events of kind {}",
            traj.find_all(e).unwrap().len(),
            e
        );
    }

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );

    assert!((fuel_usage - 93.449).abs() < 1.0);
}

#[test]
fn qlaw_as_ruggiero_case_b() {
    // Source: AAS-2004-5089
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(24505.9, 0.725, 7.05, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 160.0 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 0.350,
        isp_s: 2000.0,
    };

    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 42_165.0, 20.0),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.001, 5e-5),
        Objective::within_tolerance(StateParameter::Inclination, 0.05, 5e-3),
    ];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 1999.9;
    let dry_mass = 0.1;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[qlaw_as_ruggiero_case_b] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[qlaw_as_ruggiero_case_b] {:x}", final_state.orbit);
    println!("[qlaw_as_ruggiero_case_b] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );

    assert!((fuel_usage - 223.515).abs() < 1.0);
}

#[test]
fn qlaw_as_ruggiero_case_c() {
    // Source: AAS-2004-5089
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(9222.7, 0.2, 0.573, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 3.0 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 9.3,
        isp_s: 3100.0,
    };

    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 30_000.0, 1.0),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.7, 5e-5),
    ];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 299.9;
    let dry_mass = 0.1;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[qlaw_as_ruggiero_case_c] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[qlaw_as_ruggiero_case_c] {:x}", final_state.orbit);
    println!("[qlaw_as_ruggiero_case_c] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 41.742).abs() < 1.0);
}

#[test]
#[ignore = "https://gitlab.com/chrisrabotin/nyx/issues/103"]
fn qlaw_as_ruggiero_case_d() {
    // Broken: https://gitlab.com/chrisrabotin/nyx/issues/103
    // Source: AAS-2004-5089
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(24505.9, 0.725, 0.06, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 113.0 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 26_500.0, 1.0),
        Objective::within_tolerance(StateParameter::Inclination, 116.0, 5e-3),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.7, 5e-5),
        Objective::within_tolerance(StateParameter::RAAN, 360.0 - 90.0, 5e-3),
    ];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[qlaw_as_ruggiero_case_d] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[qlaw_as_ruggiero_case_d] {:x}", final_state.orbit);
    println!("[qlaw_as_ruggiero_case_d] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );

    assert!((fuel_usage - 23.0).abs() < 1.0);
}

#[test]
#[ignore = "https://gitlab.com/chrisrabotin/nyx/issues/103"]
fn qlaw_as_ruggiero_case_e() {
    // Broken: https://gitlab.com/chrisrabotin/nyx/issues/103
    // Source: AAS-2004-5089
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(24505.9, 0.725, 0.06, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 400.0 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 26_500.0, 1.0),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.7, 5e-5),
        Objective::within_tolerance(StateParameter::Inclination, 116.0, 5e-3),
        Objective::within_tolerance(StateParameter::RAAN, 270.0, 5e-3),
        Objective::within_tolerance(StateParameter::AoP, 180.0, 5e-3),
    ];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 1999.9;
    let dry_mass = 0.1;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[qlaw_as_ruggiero_case_e] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[qlaw_as_ruggiero_case_e] {:x}", final_state.orbit);
    println!("[qlaw_as_ruggiero_case_e] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );

    assert!((fuel_usage - 23.0).abs() < 1.0);
}

#[test]
fn qlaw_as_ruggiero_case_f() {
    // Source: AAS-2004-5089
    /*
        NOTE: Due to how lifetime of variables work in Rust, we need to define all of the
        components of a spacecraft before defining the spacecraft itself.
    */

    // We'll export this trajectory as a POC. Adding the needed crates here
    extern crate csv;
    use std::sync::mpsc;
    use std::sync::mpsc::{Receiver, Sender};
    use std::thread;

    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(15378.0, 0.01, 98.7, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 30.0 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    let objectives = &[Objective::new(StateParameter::Eccentricity, 0.15)];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[qlaw_as_ruggiero_case_f] {:x}", orbit);

    let (tx, rx): (Sender<Spacecraft>, Receiver<Spacecraft>) = mpsc::channel();

    // Set up the writing channel
    thread::spawn(move || {
        let mut wtr = csv::Writer::from_path("./rugg_case_f.csv").expect("could not create file");
        while let Ok(rx_state) = rx.recv() {
            // Serialize only the orbital state
            wtr.serialize(rx_state.orbit)
                .expect("could not serialize state");
        }
    });

    let setup =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second));
    let final_state = setup
        .with(sc_state)
        .for_duration_with_channel(prop_time, tx)
        .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[qlaw_as_ruggiero_case_f] {:x}", final_state.orbit);
    println!("[qlaw_as_ruggiero_case_f] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );

    assert!((fuel_usage - 10.378).abs() < 1.0);
}

#[test]
fn ruggiero_iepc_2011_102() {
    // Source: IEPC 2011 102
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.7283, 7.0, 1.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 105.0 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 42_164.0, 20.0),
        Objective::within_tolerance(StateParameter::Inclination, 0.001, 5e-3),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.011, 5e-5),
    ];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;

    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, ruggiero_ctrl);
    println!("[ruggiero_iepc_2011_102] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[ruggiero_iepc_2011_102] {:x}", final_state.orbit);
    println!("[ruggiero_iepc_2011_102] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );

    // WARNING: Paper claims this can be done with only 49kg of fuel.
    assert!((fuel_usage - 49.0).abs() < 1.0);
}
