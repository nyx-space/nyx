extern crate nyx_space as nyx;

use std::path::PathBuf;

use nyx::dynamics::guidance::Thruster;
use nyx::md::opti::multishoot::*;
use nyx::md::prelude::*;

use anise::{
    constants::{
        frames::{EARTH_J2000, IAU_EARTH_FRAME},
        usual_planetary_constants::MEAN_EARTH_ANGULAR_VELOCITY_DEG_S,
    },
    prelude::Almanac,
};
use nyx_space::cosmic::Mass;
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn alt_orbit_raising_cov_test(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();

    /* Define the parking orbit */
    let epoch = Epoch::from_gregorian_utc_at_noon(2022, 3, 4);
    let start =
        Orbit::try_keplerian_altitude(300.0, 0.01, 30.0, 90.0, 90.0, 60.0, epoch, eme2k).unwrap();

    /* Build the spacecraft -- really only the mass is needed here */
    let sc = Spacecraft {
        orbit: start,
        mass: Mass::from_dry_and_prop_masses(100.0, 500.0),
        thruster: Some(Thruster {
            thrust_N: 150.0,
            isp_s: 300.0,
        }),
        mode: GuidanceMode::Thrust,

        ..Default::default()
    };

    /* Define the target orbit */
    let target = Orbit::try_keplerian_altitude(
        1500.0,
        0.01,
        30.0,
        90.0,
        90.0,
        60.0,
        epoch + 2 * start.period().unwrap(),
        eme2k,
    )
    .unwrap();

    /* Define the multiple shooting parameters */
    let node_count = 30; // We're targeting 30 minutes in the future, so using 30 nodes

    let prop = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut opti = MultipleShooting::linear_altitude_heuristic(
        sc,
        target,
        node_count,
        MEAN_EARTH_ANGULAR_VELOCITY_DEG_S,
        iau_earth,
        &prop,
        almanac.clone(),
    )
    .unwrap();

    // Check that all nodes are above the surface
    println!("Initial nodes\nNode no,X (km),Y (km),Z (km),Epoch:GregorianUtc");
    for (i, node) in opti.targets.iter().enumerate() {
        println!(
            "{}, {}, {}, {}, '{}'",
            i, node.x, node.y, node.z, node.epoch
        );
    }

    let multishoot_sol = opti
        .solve(CostFunction::MinimumFuel, almanac.clone())
        .unwrap();

    println!("Final nodes\nNode no,X (km),Y (km),Z (km),Epoch:GregorianUtc");
    for (i, node) in opti.targets.iter().enumerate() {
        println!(
            "{}, {}, {}, {}, '{}'",
            i, node.x, node.y, node.z, node.epoch
        );
    }

    let all_trajectories = multishoot_sol
        .build_trajectories(&prop, almanac.clone())
        .unwrap();

    let mut full_traj = all_trajectories[0].clone();

    let output_path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "output_data"].iter().collect();

    for (i, traj) in all_trajectories.iter().enumerate() {
        traj.to_parquet_with_step(
            output_path.join(format!("multishoot_to_node_{i}.parquet")),
            2 * Unit::Second,
            almanac.clone(),
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
    prop.with(sc, almanac.clone())
        .for_duration_with_traj(start.period().unwrap())
        .unwrap()
        .1
        .to_parquet_with_step(
            output_path.join("multishoot_start.parquet"),
            2 * Unit::Second,
            almanac.clone(),
        )
        .unwrap();

    // Propagate the initial orbit too
    prop.with(sc.with_orbit(target), almanac.clone())
        .for_duration_with_traj(target.period().unwrap())
        .unwrap()
        .1
        .to_parquet_with_step(
            output_path.join("multishoot_target.parquet"),
            2 * Unit::Second,
            almanac.clone(),
        )
        .unwrap();

    // Just propagate this spacecraft for one orbit for plotting
    let (_, end_traj) = prop
        .with(sc_sol, almanac.clone())
        .for_duration_with_traj(2 * Unit::Hour)
        .unwrap();

    end_traj
        .to_parquet_with_step(
            output_path.join("multishoot_to_end.parquet"),
            2 * Unit::Second,
            almanac.clone(),
        )
        .unwrap();

    // Check that error is 50km or less. That isn't great, but I blame that on the scenario and the final node being optimized.
    let achieved_geoheight = almanac
        .transform_to(
            multishoot_sol
                .solutions
                .last()
                .unwrap()
                .achieved_state
                .orbit,
            iau_earth,
            None,
        )
        .unwrap()
        .height_km()
        .unwrap();
    let target_geoheight = almanac
        .transform_to(target, iau_earth, None)
        .unwrap()
        .height_km()
        .unwrap();
    assert!(
        (achieved_geoheight - target_geoheight).abs() < 1e-3,
        "Geodetic height achieved greater than 1 m above goal"
    );
}

#[ignore = "Does not have a valid test case yet"]
#[rstest]
fn vmag_orbit_raising(almanac: Arc<Almanac>) {
    let _ = pretty_env_logger::try_init();

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();

    /* Define the parking orbit */
    let epoch = Epoch::from_gregorian_utc_at_noon(2022, 3, 4);
    let start =
        Orbit::try_keplerian_altitude(300.0, 0.01, 30.0, 90.0, 90.0, 60.0, epoch, eme2k).unwrap();

    /* Build the spacecraft -- really only the mass is needed here */
    let sc = Spacecraft {
        orbit: start,
        mass: Mass::from_dry_and_prop_masses(100.0, 500.0),
        thruster: Some(Thruster {
            thrust_N: 150.0,
            isp_s: 300.0,
        }),
        mode: GuidanceMode::Thrust,

        ..Default::default()
    };

    /* Define the target orbit */
    let target = Orbit::try_keplerian_altitude(
        1500.0,
        0.01,
        30.0,
        90.0,
        90.0,
        60.0,
        epoch + 0.05 * start.period().unwrap(),
        eme2k,
    )
    .unwrap();

    /* Define the multiple shooting parameters */
    let node_count = 300;

    let prop = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::two_body()));
    let mut opti = MultipleShooting::equidistant_nodes(sc, target, node_count, &prop).unwrap();

    // Check that all nodes are above the surface
    println!("Initial nodes\nNode no,X (km),Y (km),Z (km),Epoch:GregorianUtc");
    for (i, node) in opti.targets.iter().enumerate() {
        println!(
            "{}, {}, {}, {}, {}, '{}'",
            i, node.x, node.y, node.z, node.vmag, node.epoch
        );
    }

    let multishoot_sol = opti
        .solve(CostFunction::MinimumFuel, almanac.clone())
        .unwrap();

    println!("Final nodes\nNode no,X (km),Y (km),Z (km),Epoch:GregorianUtc");
    for (i, node) in opti.targets.iter().enumerate() {
        println!(
            "{}, {}, {}, {}, '{}'",
            i, node.x, node.y, node.z, node.epoch
        );
    }

    let all_trajectories = multishoot_sol
        .build_trajectories(&prop, almanac.clone())
        .unwrap();

    let mut full_traj = all_trajectories[0].clone();

    let output_path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "output_data"].iter().collect();

    for (i, traj) in all_trajectories.iter().enumerate() {
        traj.to_parquet_with_step(
            output_path.join(format!("multishoot_to_node_{i}.parquet")),
            2 * Unit::Second,
            almanac.clone(),
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
    prop.with(sc, almanac.clone())
        .for_duration_with_traj(start.period().unwrap())
        .unwrap()
        .1
        .to_parquet_with_step(
            output_path.join("multishoot_start.parquet"),
            2 * Unit::Second,
            almanac.clone(),
        )
        .unwrap();

    // Propagate the initial orbit too
    prop.with(sc.with_orbit(target), almanac.clone())
        .for_duration_with_traj(target.period().unwrap())
        .unwrap()
        .1
        .to_parquet_with_step(
            output_path.join("multishoot_to_target.parquet"),
            2 * Unit::Second,
            almanac.clone(),
        )
        .unwrap();

    // Just propagate this spacecraft for one orbit for plotting
    let (_, end_traj) = prop
        .with(sc_sol, almanac.clone())
        .for_duration_with_traj(2 * Unit::Hour)
        .unwrap();

    end_traj
        .to_parquet_with_step(
            output_path.join("multishoot_end.parquet"),
            2 * Unit::Second,
            almanac.clone(),
        )
        .unwrap();

    // Check that error is 50km or less. That isn't great, but I blame that on the scenario and the final node being optimized.
    let achieved_geoheight = almanac
        .transform_to(
            multishoot_sol
                .solutions
                .last()
                .unwrap()
                .achieved_state
                .orbit,
            iau_earth,
            None,
        )
        .unwrap()
        .height_km()
        .unwrap();
    let target_geoheight = almanac
        .transform_to(target, iau_earth, None)
        .unwrap()
        .height_km()
        .unwrap();
    assert!(
        (achieved_geoheight - target_geoheight).abs() < 1e-3,
        "Geodetic height achieved greater than 1 m above goal"
    );

    for (i, traj) in multishoot_sol
        .build_trajectories(&prop, almanac.clone())
        .unwrap()
        .iter()
        .enumerate()
    {
        traj.to_parquet_with_step(
            format!("multishoot_to_node_{}.parquet", i),
            2 * Unit::Second,
            almanac.clone(),
        )
        .unwrap();
    }
}
