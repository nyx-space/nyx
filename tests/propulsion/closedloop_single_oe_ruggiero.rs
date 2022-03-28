extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::nyx::cosmic::{Cosm, GuidanceMode, Orbit, Spacecraft};
use self::nyx::dynamics::guidance::{Objective, Ruggiero, StateParameter, Thruster};
use self::nyx::dynamics::{OrbitalDynamics, SpacecraftDynamics};
use self::nyx::propagators::{PropOpts, Propagator, RK4Fixed};
use self::nyx::time::{Epoch, Unit};

#[test]
fn rugg_sma() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 45 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::SMA,
        42_164.0,
        1.0,
    )];

    let guid_law = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_sma] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_sma] {:x}", final_state.orbit);
    println!("[rugg_sma] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 21.0).abs() < 1.0);
}

#[test]
fn rugg_sma_regress_threshold() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.1, 0.0, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 175 * Unit::Day;

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::SMA,
        42_164.0,
        1.0,
    )];

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    for (threshold, expected_fuel_usage) in &[(0.9, 16.9), (0.0, 21.3)] {
        let guid_law = Ruggiero::with_ηthresholds(objectives, &[*threshold], orbit).unwrap();
        let sc_state =
            Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

        let sc = SpacecraftDynamics::from_guidance_law(OrbitalDynamics::two_body(), guid_law);
        println!("[rugg_sma_regress] {:x}", orbit);

        let final_state =
            Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
                .with(sc_state)
                .for_duration(prop_time)
                .unwrap();

        let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
        println!("[rugg_sma_regress] {:x}", final_state.orbit);
        println!("[rugg_sma_regress] fuel usage: {:.3} kg", fuel_usage);

        assert!(
            sc.guidance_achieved(&final_state).unwrap(),
            "objective not achieved"
        );
        assert!((fuel_usage - *expected_fuel_usage).abs() < 0.5);
    }
}

#[test]
fn rugg_sma_decr() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian(42164.0, 0.0, 0.0, 0.0, 0.0, 0.0, start_time, eme2k);

    let prop_time = 45 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::SMA,
        24_396.0,
        1.0,
    )];

    let guid_law = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_sma_decr] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_sma_decr] {:x}", final_state.orbit);
    println!("[rugg_sma_decr] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
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

    let prop_time = 55 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::Inclination,
        51.6,
        5e-3,
    )];

    let guid_law = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_inc] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_inc] {:x}", final_state.orbit);
    println!("[rugg_inc] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 25.0).abs() < 1.0);
}

#[test]
fn rugg_inc_threshold() {
    // Same inclination test as above, but with an efficiency threshold. Data comes from Figure 7 of IEPC-2011-102.
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let orbit = Orbit::keplerian_altitude(350.0, 0.001, 46.0, 1.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 130 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::Inclination,
        51.6,
        5e-3,
    )];

    let guid_law = Ruggiero::with_ηthresholds(objectives, &[0.9], orbit).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_inc_threshold] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_inc_threshold] {:x}", final_state.orbit);
    println!("[rugg_inc_threshold] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 17.0).abs() < 1.0);
}

#[test]
fn rugg_inc_decr() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.equatorial_radius() + 350.0;

    let orbit = Orbit::keplerian(sma, 0.001, 51.6, 1.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 55 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::Inclination,
        46.0,
        5e-3,
    )];

    let guid_law = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_inc_decr] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_inc_decr] {:x}", final_state.orbit);
    println!("[rugg_inc_decr] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
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

    let prop_time = 30 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::Eccentricity,
        0.15,
        5e-5,
    )];

    let guid_law = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_ecc] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_ecc] {:x}", final_state.orbit);
    println!("[rugg_ecc] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 10.37).abs() < 1.0);
}

#[test]
fn rugg_ecc_regress_threshold() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.equatorial_radius() + 9000.0;

    let orbit = Orbit::keplerian(sma, 0.01, 98.7, 0.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 150 * Unit::Day;

    let fuel_mass = 67.0;
    let dry_mass = 300.0;

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::Eccentricity,
        0.15,
        5e-5,
    )];

    for (threshold, expected_fuel_usage) in &[(0.9, 8.2), (0.0, 10.37)] {
        let guid_law = Ruggiero::with_ηthresholds(objectives, &[*threshold], orbit).unwrap();

        let sc_state =
            Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

        let sc = SpacecraftDynamics::from_guidance_law(OrbitalDynamics::two_body(), guid_law);
        println!("[rugg_ecc_regress] {:x}", orbit);

        let final_state =
            Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
                .with(sc_state)
                .for_duration(prop_time)
                .unwrap();

        let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
        println!("[rugg_ecc_regress] {:x}", final_state.orbit);
        println!("[rugg_ecc_regress] fuel usage: {:.3} kg", fuel_usage);

        assert!(
            sc.guidance_achieved(&final_state).unwrap(),
            "objective not achieved"
        );
        assert!((fuel_usage - *expected_fuel_usage).abs() < 1.0);
    }
}

#[test]
fn rugg_ecc_decr() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.equatorial_radius() + 9000.0;

    let orbit = Orbit::keplerian(sma, 0.15, 98.7, 0.0, 1.0, 1.0, start_time, eme2k);

    let prop_time = 30 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::Eccentricity,
        0.01,
        5e-5,
    )];

    let guid_law = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_ecc_decr] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_ecc_decr] {:x}", final_state.orbit);
    println!("[rugg_ecc_decr] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
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
    let prop_time = 44 * Unit::Minute + 10 * Unit::Second;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::AoP,
        183.0,
        5e-3,
    )];

    let guid_law = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_aop] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_aop] {:x}", final_state.orbit);
    println!("[rugg_aop] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
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

    let prop_time = 44 * Unit::Minute + 10 * Unit::Second;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(
        StateParameter::AoP,
        178.0,
        5e-3,
    )];

    let guid_law = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_aop_decr] {:x}", orbit);

    let final_state =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
            .with(sc_state)
            .for_duration(prop_time)
            .unwrap();

    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_aop_decr] {:x}", final_state.orbit);
    println!("[rugg_aop_decr] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
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

    let prop_time = 49 * Unit::Day;

    // Define the dynamics
    let orbital_dyn = OrbitalDynamics::two_body();

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(StateParameter::RAAN, 5.0, 5e-5)];

    let guid_law = Ruggiero::new(objectives, orbit).unwrap();

    let fuel_mass = 67.0;
    let dry_mass = 300.0;
    let sc_state =
        Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

    let sc = SpacecraftDynamics::from_guidance_law(orbital_dyn, guid_law);
    println!("[rugg_raan] {:x}", orbit);

    let setup =
        Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second));
    let mut prop = setup.with(sc_state);
    let (final_state, traj) = prop.for_duration_with_traj(prop_time).unwrap();
    let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
    println!("[rugg_raan] {:x}", final_state.orbit);
    let event = Event::new(StateParameter::RAAN, 5.0);
    println!("[rugg_raan] {} => {:?}", event, traj.find_all(&event));
    println!("[rugg_raan] fuel usage: {:.3} kg", fuel_usage);

    assert!(
        sc.guidance_achieved(&final_state).unwrap(),
        "objective not achieved"
    );
    assert!((fuel_usage - 22.189).abs() < 1.0);
}

#[test]
fn rugg_raan_regress_threshold() {
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let sma = eme2k.equatorial_radius() + 798.0;

    let orbit = Orbit::keplerian(sma, 0.00125, 98.57, 0.0, 1.0, 0.0, start_time, eme2k);

    let prop_time = 130 * Unit::Day;

    let fuel_mass = 67.0;
    let dry_mass = 300.0;

    // Define the thruster
    let lowt = Thruster {
        thrust_N: 89e-3,
        isp_s: 1650.0,
    };

    // Define the objectives
    let objectives = &[Objective::within_tolerance(StateParameter::RAAN, 5.0, 5e-5)];

    for (threshold, expected_fuel_usage) in &[(0.9, 14.787), (0.0, 22.189)] {
        let guid_law = Ruggiero::with_ηthresholds(objectives, &[*threshold], orbit).unwrap();

        let sc_state =
            Spacecraft::from_thruster(orbit, dry_mass, fuel_mass, lowt, GuidanceMode::Thrust);

        let sc = SpacecraftDynamics::from_guidance_law(OrbitalDynamics::two_body(), guid_law);
        println!("[rugg_raan_regress] {:x}", orbit);

        let final_state =
            Propagator::new::<RK4Fixed>(sc.clone(), PropOpts::with_fixed_step(10.0 * Unit::Second))
                .with(sc_state)
                .for_duration(prop_time)
                .unwrap();

        let fuel_usage = fuel_mass - final_state.fuel_mass_kg;
        println!("[rugg_raan_regress] {:x}", final_state.orbit);
        println!("[rugg_raan_regress] fuel usage: {:.3} kg", fuel_usage);

        assert!(
            sc.guidance_achieved(&final_state).unwrap(),
            "objective not achieved"
        );
        assert!((fuel_usage - *expected_fuel_usage).abs() < 1e-1);
    }
}
