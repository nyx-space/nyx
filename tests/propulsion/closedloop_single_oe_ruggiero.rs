extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::nyx::cosmic::{Cosm, GuidanceMode, Orbit, Spacecraft};
use self::nyx::dynamics::guidance::{Achieve, Ruggiero, Thruster};
use self::nyx::dynamics::{OrbitalDynamics, SpacecraftDynamics};
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use self::nyx::time::{Epoch, TimeUnit};

#[test]
fn rugg_sma() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 45 * TimeUnit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    };

    // Define the objectives
    let objectives = vec![Achieve::Sma {
        target: 42164.0,
        tol: 1.0,
    }];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit);

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_ctrl(orbital_dyn, ruggiero_ctrl);
    println!("[rugg_sma] {:x}", orbit);

    let final_state = Propagator::new::<RK4Fixed>(
        sc.clone(),
        PropOpts::with_fixed_step(10.0 * TimeUnit::Second),
    )
    .with(sc_state)
    .for_duration(prop_time)
    .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_sma] {:x}", final_state.orbit);
    println!("[rugg_sma] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.ctrl_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 21.0).abs() < 1.0);
}

#[test]
fn rugg_sma_decr() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(42164.0, 0.0, 0.0, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 45 * TimeUnit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    };

    // Define the objectives
    let objectives = vec![Achieve::Sma {
        target: 24396.0,
        tol: 1.0,
    }];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit);

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_ctrl(orbital_dyn, ruggiero_ctrl);
    println!("[rugg_sma_decr] {:x}", orbit);

    let final_state = Propagator::new::<RK4Fixed>(
        sc.clone(),
        PropOpts::with_fixed_step(10.0 * TimeUnit::Second),
    )
    .with(sc_state)
    .for_duration(prop_time)
    .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_sma_decr] {:x}", final_state.orbit);
    println!("[rugg_sma_decr] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.ctrl_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 21.0).abs() < 1.0);
}

#[test]
fn rugg_inc() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.equatorial_radius() + 350.0;

    let orbit = Orbit::keplerian(sma, 0.001, 46.0, 1.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 55 * TimeUnit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    };

    // Define the objectives
    let objectives = vec![Achieve::Inc {
        target: 51.6,
        tol: 5e-3,
    }];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit);

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_ctrl(orbital_dyn, ruggiero_ctrl);
    println!("[rugg_inc] {:x}", orbit);

    let final_state = Propagator::new::<RK4Fixed>(
        sc.clone(),
        PropOpts::with_fixed_step(10.0 * TimeUnit::Second),
    )
    .with(sc_state)
    .for_duration(prop_time)
    .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_inc] {:x}", final_state.orbit);
    println!("[rugg_inc] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.ctrl_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 25.0).abs() < 1.0);
}

#[test]
fn rugg_inc_decr() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.equatorial_radius() + 350.0;

    let orbit = Orbit::keplerian(sma, 0.001, 51.6, 1.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 55 * TimeUnit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    };

    // Define the objectives
    let objectives = vec![Achieve::Inc {
        target: 46.0,
        tol: 5e-3,
    }];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit);

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_ctrl(orbital_dyn, ruggiero_ctrl);
    println!("[rugg_inc_decr] {:x}", orbit);

    let final_state = Propagator::new::<RK4Fixed>(
        sc.clone(),
        PropOpts::with_fixed_step(10.0 * TimeUnit::Second),
    )
    .with(sc_state)
    .for_duration(prop_time)
    .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_inc_decr] {:x}", final_state.orbit);
    println!("[rugg_inc_decr] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.ctrl_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 25.0).abs() < 1.0);
}

#[test]
fn rugg_ecc() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.equatorial_radius() + 9000.0;

    let orbit = Orbit::keplerian(sma, 0.01, 98.7, 0.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 30 * TimeUnit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    };

    // Define the objectives
    let objectives = vec![Achieve::Ecc {
        target: 0.15,
        tol: 5e-5,
    }];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit);

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_ctrl(orbital_dyn, ruggiero_ctrl);
    println!("[rugg_ecc] {:x}", orbit);

    let final_state = Propagator::new::<RK4Fixed>(
        sc.clone(),
        PropOpts::with_fixed_step(10.0 * TimeUnit::Second),
    )
    .with(sc_state)
    .for_duration(prop_time)
    .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_ecc] {:x}", final_state.orbit);
    println!("[rugg_ecc] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.ctrl_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 10.37).abs() < 1.0);
}

#[test]
fn rugg_ecc_decr() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.equatorial_radius() + 9000.0;

    let orbit = Orbit::keplerian(sma, 0.15, 98.7, 0.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 30 * TimeUnit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    };

    // Define the objectives
    let objectives = vec![Achieve::Ecc {
        target: 0.01,
        tol: 5e-5,
    }];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit);

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_ctrl(orbital_dyn, ruggiero_ctrl);
    println!("[rugg_ecc_decr] {:x}", orbit);

    let final_state = Propagator::new::<RK4Fixed>(
        sc.clone(),
        PropOpts::with_fixed_step(10.0 * TimeUnit::Second),
    )
    .with(sc_state)
    .for_duration(prop_time)
    .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_ecc_decr] {:x}", final_state.orbit);
    println!("[rugg_ecc_decr] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.ctrl_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 10.37).abs() < 1.0);
}

#[test]
fn rugg_aop() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.equatorial_radius() + 900.0;

    // Note that AOP computation requires the orbit to not be equatorial or circular, hence the non-zero ECC and INC.
    let orbit = Orbit::keplerian(sma, 5e-5, 5e-3, 0.0, 178.0, 0.0, start_time, eme2k);

    // This is a very quick change because we aren't using the Ruggiero formulation for AOP change and benefit both in-plane and out of plane control.
    let prop_time = 44 * TimeUnit::Minute + 10 * TimeUnit::Second;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    };

    // Define the objectives
    let objectives = vec![Achieve::Aop {
        target: 183.0,
        tol: 5e-3,
    }];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit);

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_ctrl(orbital_dyn, ruggiero_ctrl);
    println!("[rugg_aop] {:x}", orbit);

    let final_state = Propagator::new::<RK4Fixed>(
        sc.clone(),
        PropOpts::with_fixed_step(10.0 * TimeUnit::Second),
    )
    .with(sc_state)
    .for_duration(prop_time)
    .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_aop] {:x}", final_state.orbit);
    println!("[rugg_aop] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.ctrl_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 0.014).abs() < 1e-2);
}

#[test]
fn rugg_aop_decr() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.equatorial_radius() + 900.0;

    // Note that AOP computation requires the orbit to not be equatorial or circular, hence the non-zero ECC and INC.
    let orbit = Orbit::keplerian(sma, 5e-5, 5e-3, 0.0, 183.0, 0.0, start_time, eme2k);

    let prop_time = 44 * TimeUnit::Minute + 10 * TimeUnit::Second;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    };

    // Define the objectives
    let objectives = vec![Achieve::Aop {
        target: 178.0,
        tol: 5e-3,
    }];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit);

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_ctrl(orbital_dyn, ruggiero_ctrl);
    println!("[rugg_aop_decr] {:x}", orbit);

    let final_state = Propagator::new::<RK4Fixed>(
        sc.clone(),
        PropOpts::with_fixed_step(10.0 * TimeUnit::Second),
    )
    .with(sc_state)
    .for_duration(prop_time)
    .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_aop_decr] {:x}", final_state.orbit);
    println!("[rugg_aop_decr] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.ctrl_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 0.014).abs() < 1e-2);
}

#[test]
fn rugg_raan() {
    use self::nyx::md::{Event, StateParameter};
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2017, 1, 1);

    let sma = eme2k.equatorial_radius() + 798.0;

    let orbit = Orbit::keplerian(sma, 0.00125, 98.57, 0.0, 1.0, 0.0, start_time, eme2k);

    let prop_time = 49 * TimeUnit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust: 89e-3,
        isp: 1650.0,
    };

    // Define the objectives
    let objectives = vec![Achieve::Raan {
        target: 5.0,
        tol: 5e-3,
    }];

    let ruggiero_ctrl = Ruggiero::new(objectives, orbit);

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_ctrl(orbital_dyn, ruggiero_ctrl);
    println!("[rugg_raan] {:x}", orbit);

    let setup = Propagator::new::<RK4Fixed>(
        sc.clone(),
        PropOpts::with_fixed_step(10.0 * TimeUnit::Second),
    );
    let mut prop = setup.with(sc_state);
    let (final_state, traj) = prop.for_duration_with_traj(prop_time).unwrap();
    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_raan] {:x}", final_state.orbit);
    let event = Event::new(StateParameter::RAAN, 5.0);
    println!("[rugg_raan] {} => {:?}", event, traj.find_all(&event));
    println!("[rugg_raan] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.ctrl_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 22.189).abs() < 1.0);
}
