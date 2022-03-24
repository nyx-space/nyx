extern crate nyx_space as nyx;

use nyx::dynamics::guidance::Thruster;
use nyx::md::opti::multishoot::*;
use nyx::md::ui::*;

#[test]
fn alt_orbit_raising() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let iau_earth = cosm.frame("IAU_Earth");

    /* Define the parking orbit */
    let epoch = Epoch::from_gregorian_utc_at_noon(2022, 03, 04);
    let start = Orbit::keplerian_altitude(300.0, 0.01, 30.0, 90.0, 90.0, 60.0, epoch, eme2k);

    /* Build the spacecraft -- really only the mass is needed here */
    let sc = Spacecraft {
        orbit: start,
        dry_mass_kg: 100.0,
        fuel_mass_kg: 500.0,
        thruster: Some(Thruster {
            thrust_N: 150.0,
            isp_s: 300.0,
        }),
        ext: GuidanceMode::Thrust,

        ..Default::default()
    };

    /* Define the target orbit */
    let target = Orbit::keplerian_altitude(
        1500.0,
        0.01,
        30.0,
        90.0,
        90.0,
        60.0,
        epoch + 2 * start.period(),
        eme2k,
    );

    /* Define the multiple shooting parameters */
    let node_count = 30; // We're targeting 30 minutes in the future, so using 30 nodes

    let prop = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut opti = MultipleShooting::linear_altitude_heuristic(
        sc,
        target,
        node_count,
        iau_earth,
        &prop,
        cosm.clone(),
    )
    .unwrap();

    // Check that all nodes are above the surface
    println!("Initial nodes\nNode no,X (km),Y (km),Z (km),Epoch:GregorianUtc");
    for (i, node) in opti.targets.iter().enumerate() {
        println!(
            "{}, {}, {}, {}, '{}'",
            i,
            node.x,
            node.y,
            node.z,
            node.epoch.as_gregorian_utc_str()
        );
    }

    let multishoot_sol = opti.solve(CostFunction::MinimumFuel).unwrap();

    println!("Final nodes\nNode no,X (km),Y (km),Z (km),Epoch:GregorianUtc");
    for (i, node) in opti.targets.iter().enumerate() {
        println!(
            "{}, {}, {}, {}, '{}'",
            i,
            node.x,
            node.y,
            node.z,
            node.epoch.as_gregorian_utc_str()
        );
    }

    let all_trajectories = multishoot_sol.build_trajectories(&prop).unwrap();

    let mut full_traj = all_trajectories[0].clone();

    for (i, traj) in all_trajectories.iter().enumerate() {
        traj.to_csv_with_step(
            &format!("multishoot_to_node_{}.csv", i),
            2 * Unit::Second,
            cosm.clone(),
        )
        .unwrap();
        if i > 0 {
            full_traj += traj;
        }
    }

    assert_eq!(
        full_traj.first().epoch(),
        all_trajectories[0].first().epoch(),
        "Initial epochs differ"
    );
    assert_eq!(
        full_traj.last().epoch(),
        all_trajectories.last().unwrap().last().epoch(),
        "Final epochs differ: {} != {}",
        full_traj.last().epoch(),
        all_trajectories.last().unwrap().last().epoch(),
    );

    let solution = &multishoot_sol.solutions[node_count - 1];
    let sc_sol = solution.achieved_state;

    println!("{}", multishoot_sol);

    // Compute the total delta-v of the solution
    let mut dv_ms = 0.0;
    for sol in &multishoot_sol.solutions {
        dv_ms += sol.correction.norm() * 1e3;
    }
    println!(
        "Multiple shooting solution requires a total of {:.3} m/s",
        dv_ms
    );

    // Propagate the initial orbit too
    prop.with(sc)
        .for_duration_with_traj(start.period())
        .unwrap()
        .1
        .to_csv_with_step(
            &"multishoot_start.csv".to_string(),
            2 * Unit::Second,
            cosm.clone(),
        )
        .unwrap();

    // Propagate the initial orbit too
    prop.with(sc.with_orbit(target))
        .for_duration_with_traj(target.period())
        .unwrap()
        .1
        .to_csv_with_step(
            &"multishoot_target.csv".to_string(),
            2 * Unit::Second,
            cosm.clone(),
        )
        .unwrap();

    // Just propagate this spacecraft for one orbit for plotting
    let (_, end_traj) = prop
        .with(sc_sol)
        .for_duration_with_traj(2 * Unit::Hour)
        .unwrap();

    end_traj
        .to_csv_with_step(
            &"multishoot_to_end.csv".to_string(),
            2 * Unit::Second,
            cosm.clone(),
        )
        .unwrap();

    // Check that error is 50km or less. That isn't great, but I blame that on the scenario and the final node being optimized.
    let achieved_geoheight = cosm
        .frame_chg(
            &multishoot_sol
                .solutions
                .last()
                .unwrap()
                .achieved_state
                .orbit,
            iau_earth,
        )
        .geodetic_height();
    let target_geoheight = cosm.frame_chg(&target, iau_earth).geodetic_height();
    assert!(
        (achieved_geoheight - target_geoheight).abs() < 1e-3,
        "Geodetic height achieved greater than 1 m above goal"
    );
}

#[ignore = "Does not have a valid test case yet"]
#[test]
fn vmag_orbit_raising() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let iau_earth = cosm.frame("IAU_Earth");

    /* Define the parking orbit */
    let epoch = Epoch::from_gregorian_utc_at_noon(2022, 03, 04);
    let start = Orbit::keplerian_altitude(300.0, 0.01, 30.0, 90.0, 90.0, 60.0, epoch, eme2k);

    /* Build the spacecraft -- really only the mass is needed here */
    let sc = Spacecraft {
        orbit: start,
        dry_mass_kg: 100.0,
        fuel_mass_kg: 500.0,
        thruster: Some(Thruster {
            thrust_N: 150.0,
            isp_s: 300.0,
        }),
        ext: GuidanceMode::Thrust,

        ..Default::default()
    };

    /* Define the target orbit */
    let target = Orbit::keplerian_altitude(
        1500.0,
        0.01,
        30.0,
        90.0,
        90.0,
        60.0,
        epoch + 0.05 * start.period(),
        eme2k,
    );

    /* Define the multiple shooting parameters */
    let node_count = 300;

    let prop = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut opti = MultipleShooting::equidistant_nodes(sc, target, node_count, &prop).unwrap();

    // Check that all nodes are above the surface
    println!("Initial nodes\nNode no,X (km),Y (km),Z (km),Epoch:GregorianUtc");
    for (i, node) in opti.targets.iter().enumerate() {
        println!(
            "{}, {}, {}, {}, {}, '{}'",
            i,
            node.x,
            node.y,
            node.z,
            node.vmag,
            node.epoch.as_gregorian_utc_str()
        );
    }

    let multishoot_sol = opti.solve(CostFunction::MinimumFuel).unwrap();

    println!("Final nodes\nNode no,X (km),Y (km),Z (km),Epoch:GregorianUtc");
    for (i, node) in opti.targets.iter().enumerate() {
        println!(
            "{}, {}, {}, {}, '{}'",
            i,
            node.x,
            node.y,
            node.z,
            node.epoch.as_gregorian_utc_str()
        );
    }

    let all_trajectories = multishoot_sol.build_trajectories(&prop).unwrap();

    let mut full_traj = all_trajectories[0].clone();

    for (i, traj) in all_trajectories.iter().enumerate() {
        traj.to_csv_with_step(
            &format!("multishoot_to_node_{}.csv", i),
            2 * Unit::Second,
            cosm.clone(),
        )
        .unwrap();
        if i > 0 {
            full_traj += traj;
        }
    }

    assert_eq!(
        full_traj.first().epoch(),
        all_trajectories[0].first().epoch(),
        "Initial epochs differ"
    );
    assert_eq!(
        full_traj.last().epoch(),
        all_trajectories.last().unwrap().last().epoch(),
        "Final epochs differ: {} != {}",
        full_traj.last().epoch(),
        all_trajectories.last().unwrap().last().epoch(),
    );

    assert_eq!(
        full_traj.first().epoch(),
        all_trajectories[0].first().epoch(),
        "Initial epochs differ"
    );
    assert_eq!(
        full_traj.last().epoch(),
        all_trajectories.last().unwrap().last().epoch(),
        "Final epochs differ: {} != {}",
        full_traj.last().epoch(),
        all_trajectories.last().unwrap().last().epoch(),
    );

    let solution = &multishoot_sol.solutions[node_count - 1];
    let sc_sol = solution.achieved_state;

    println!("{}", multishoot_sol);

    // Compute the total delta-v of the solution
    let mut dv_ms = 0.0;
    for sol in &multishoot_sol.solutions {
        dv_ms += sol.correction.norm() * 1e3;
    }
    println!(
        "Multiple shooting solution requires a total of {:.3} m/s",
        dv_ms
    );

    assert!((dv_ms - 735.9).abs() < 0.1, "Wrong total DV");

    // Propagate the initial orbit too
    prop.with(sc)
        .for_duration_with_traj(start.period())
        .unwrap()
        .1
        .to_csv_with_step(
            &"multishoot_start.csv".to_string(),
            2 * Unit::Second,
            cosm.clone(),
        )
        .unwrap();

    // Propagate the initial orbit too
    prop.with(sc.with_orbit(target))
        .for_duration_with_traj(target.period())
        .unwrap()
        .1
        .to_csv_with_step(
            &"multishoot_target.csv".to_string(),
            2 * Unit::Second,
            cosm.clone(),
        )
        .unwrap();

    // Just propagate this spacecraft for one orbit for plotting
    let (_, end_traj) = prop
        .with(sc_sol)
        .for_duration_with_traj(2 * Unit::Hour)
        .unwrap();

    end_traj
        .to_csv_with_step(
            &"multishoot_to_end.csv".to_string(),
            2 * Unit::Second,
            cosm.clone(),
        )
        .unwrap();

    // Check that error is 50km or less. That isn't great, but I blame that on the scenario and the final node being optimized.
    let achieved_geoheight = cosm
        .frame_chg(
            &multishoot_sol
                .solutions
                .last()
                .unwrap()
                .achieved_state
                .orbit,
            iau_earth,
        )
        .geodetic_height();
    let target_geoheight = cosm.frame_chg(&target, iau_earth).geodetic_height();
    assert!(
        (achieved_geoheight - target_geoheight).abs() < 1e-3,
        "Geodetic height achieved greater than 1 m above goal"
    );

    for (i, traj) in multishoot_sol
        .build_trajectories(&prop)
        .unwrap()
        .iter()
        .enumerate()
    {
        traj.to_csv_with_step(
            &format!("multishoot_to_node_{}.csv", i),
            2 * Unit::Second,
            cosm.clone(),
        )
        .unwrap();
    }
}
