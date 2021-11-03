extern crate nyx_space as nyx;

use nyx::dynamics::guidance::Thruster;
use nyx::md::ui::*;
use nyx::opti::multishoot::*;

// TODO -- Add an Earth orbit to Lunar orbit multiple shooting algorithm and see if that works

#[test]
fn landing_demo() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    const SITE_LAT_DEG: f64 = -86.798;
    const SITE_LONG_DEG: f64 = -21.150;
    const SITE_HEIGHT_KM: f64 = 0.4;
    const ALTITUDE_BUFFER_KM: f64 = 0.5;
    let cosm = Cosm::de438_gmat();
    let moonj2k = cosm.frame("Luna");
    let moon_iau = cosm.frame("IAU Moon");

    /* *** */
    /* Define the landing epoch, landing site, and initial orbit. */
    /* *** */
    // Landing epoch
    let e_landing = Epoch::from_gregorian_utc(2023, 11, 25, 14, 19, 26, 789000000);
    // let e_landing = Epoch::from_str("2023-11-25T14:19:26.789000 UTC").unwrap();
    // let vertical_landing_duration = 2 * TimeUnit::Minute;
    // let e_pre_landing = e_landing - vertical_landing_duration;
    // Landing site with buffer
    let ls_buf = Orbit::from_geodesic(
        SITE_LAT_DEG,
        SITE_LONG_DEG,
        SITE_HEIGHT_KM + ALTITUDE_BUFFER_KM,
        e_landing,
        moon_iau,
    );
    println!(
        "Target landing site (with buffer):\n{}\nlong = {:.2} deg\tlat = {:.2} deg\talt = {:.3}km\t",
        ls_buf,
        ls_buf.geodetic_longitude(),
        ls_buf.geodetic_latitude(),
        ls_buf.geodetic_height()
    );

    let ls = Orbit::from_geodesic(
        SITE_LAT_DEG,
        SITE_LONG_DEG,
        SITE_HEIGHT_KM,
        e_landing,
        moon_iau,
    );

    // DRM state in Moon J2000
    let start_orbit = Orbit::cartesian(
        -77.19951157,
        637.0714101,
        -1639.05629855,
        0.72861834,
        1.47489041,
        0.48984635,
        e_landing,
        moonj2k,
    );

    let start_orbit_iau = cosm.frame_chg(&start_orbit, moon_iau);
    println!(
        "DRM end\n{}\nlong = {:.2} deg\tlat = {:.2} deg\talt = {:.3}km\t",
        start_orbit_iau,
        start_orbit_iau.geodetic_longitude(),
        start_orbit_iau.geodetic_latitude(),
        start_orbit_iau.geodetic_height()
    );

    let thruster = Thruster {
        isp: 300.0,
        thrust: 6.0 * 667.233,
    };
    // Masses from δCDR
    let dry_mass_kg = 1100.0;
    let fuel_mass_kg = 1479.4 - 44.353;
    let xl1 = Spacecraft::from_thruster(
        start_orbit,
        dry_mass_kg,
        fuel_mass_kg,
        thruster,
        GuidanceMode::Coast,
    );

    /* *** */
    /* Run the differential corrector for the initial guess of the velocity vector. */
    /* *** */
    let prop = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));

    let pdi_start = prop.with(xl1).for_duration(-8 * TimeUnit::Minute).unwrap();

    let pdi_start_fixed = cosm.frame_chg(&pdi_start.orbit, cosm.frame("IAU Moon"));

    let num_nodes = 8 * 3;
    let ls_luna = cosm.frame_chg(&ls, moonj2k);

    // Convert the landing site into the same frame as the spacecraft and use that as targeting values
    let ls_luna_buf = cosm.frame_chg(&ls_buf, moonj2k);

    println!("PDI start: {}", pdi_start);
    println!(
        "PDI Start: |r| = {:.4} km\t|v| = {:.4} km/s\taltitude = {:.4} km\tlanding altitude = {:.4} km\nPDI lat = {:.2} deg\tlong = {:.2} deg",
        pdi_start.orbit.rmag(),
        pdi_start.orbit.vmag(),
        pdi_start_fixed.geodetic_height(),
        ls_buf.geodetic_height(),
        pdi_start_fixed.geodetic_latitude(),
        pdi_start_fixed.geodetic_longitude()
    );

    println!("LANDING SITE: {}", ls_luna);
    println!(
        "LANDING SITE slant angle: φ = {} deg",
        pdi_start
            .orbit
            .r_hat()
            .dot(&ls_luna.r_hat())
            .acos()
            .to_degrees()
    );

    // And run the multiple shooting algorithm
    let mut opti = MultipleShooting::linear_heuristic(
        pdi_start,
        ls_luna_buf,
        num_nodes,
        cosm.frame("IAU Moon"),
        &prop,
        cosm.clone(),
    )
    .unwrap();

    // Check that all nodes are above the surface
    println!("Initial nodes\nNode no,X (km),Y (km),Z (km),Epoch:GregorianUtc");
    for (i, node) in opti.nodes.iter().enumerate() {
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
    for (i, node) in opti.nodes.iter().enumerate() {
        println!(
            "{}, {}, {}, {}, '{}'",
            i,
            node.x,
            node.y,
            node.z,
            node.epoch.as_gregorian_utc_str()
        );
    }

    for (i, traj) in multishoot_sol
        .build_trajectories(&prop)
        .unwrap()
        .iter()
        .enumerate()
    {
        traj.to_csv_with_step(
            &format!("multishoot_to_node_{}.csv", i),
            2 * TimeUnit::Second,
            cosm.clone(),
        )
        .unwrap();
    }
    // And propagate the final state too
    let (_, traj) = prop
        .with(multishoot_sol.solutions.last().unwrap().achieved)
        .for_duration_with_traj(2 * TimeUnit::Hour)
        .unwrap();

    traj.to_csv_with_step(
        &"multishoot_to_end.csv".to_string(),
        2 * TimeUnit::Second,
        cosm.clone(),
    )
    .unwrap();

    let solution = &multishoot_sol.solutions[num_nodes - 1];
    let sc_near_ls = solution.achieved;
    // Compute the total delta-v of the solution
    let mut dv_ms = 0.0;
    for sol in &multishoot_sol.solutions {
        dv_ms += sol.correction.norm() * 1e3;
    }
    println!(
        "Multiple shooting solution requires a total of {} m/s",
        dv_ms
    );

    println!("{}", multishoot_sol);

    println!(
        "SC above landing site: {}\taltitude = {:.3} km\t\tlanding site altitude: = {:.3} km",
        sc_near_ls,
        sc_near_ls.orbit.geodetic_height(),
        ls_luna.geodetic_height()
    );
    // Check r_miss and v_miss
    let r_miss = (sc_near_ls.orbit.radius() - ls_luna_buf.radius()).norm();
    let v_miss = (sc_near_ls.orbit.velocity() - ls_luna_buf.velocity()).norm();
    println!(
        "\nABOVE\n{}\n\tr_miss = {:.1} m\tv_miss = {:.1} m/s\n\theight = {:.3} km (want: {:.3} km)\tlat = {:.3} deg (want: {:.3} deg)\tlong = {:.3} deg (want: {:.3} deg)",
        sc_near_ls,
        r_miss * 1e3,
        v_miss * 1e3,
        sc_near_ls.orbit.geodetic_height(),
        ls_luna.geodetic_height(),
        sc_near_ls.orbit.geodetic_latitude(),
        ls_luna.geodetic_latitude(),
        sc_near_ls.orbit.geodetic_longitude(),
        ls_luna.geodetic_longitude(),
    );

    // Just propagate this spacecraft until lunar impact
    let interception = prop
        .with(sc_near_ls)
        .until_event(
            2 * TimeUnit::Hour,
            &Event::within_tolerance(StateParameter::Rmag, ls_luna.rmag(), 1e-3),
            0,
        )
        .unwrap()
        .0;

    // Now, try to land with zero velocity and check again the accuracy
    // let objectives = vec![Objective::new(StateParameter::Vmag, 1.0e-3)];
    // let tgt = Targeter::delta_v(&prop, objectives);
    // let sc_landing_tgt = tgt
    //     .try_achieve_from(sc_near_ls, sc_near_ls.epoch(), e_landing)
    //     .unwrap();
    // Check r_miss and v_miss
    let r_miss = (interception.orbit.radius() - ls_luna.radius()).norm();
    let v_miss = (interception.orbit.velocity() - ls_luna.velocity()).norm();
    println!(
        "\nFINALLY\n{}\n\tr_miss = {:.1} m\tv_miss = {:.1} m/s",
        interception,
        r_miss * 1e3,
        v_miss * 1e3
    );
}
