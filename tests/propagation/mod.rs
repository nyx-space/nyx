mod events;
mod propagators;
mod stm;
mod stopcond;
mod trajectory;

#[test]
fn orbit_finder() {
    const SITE_LAT_DEG: f64 = -86.794;
    const SITE_LONG_DEG: f64 = -21.186;
    const SITE_HEIGHT_KM: f64 = 0.4;

    extern crate nyx_space as nyx;
    use nyx::cosmic::eclipse::{EclipseLocator, EclipseState};
    use nyx::cosmic::{Bodies, Cosm, LightTimeCalc, Orbit};
    use nyx::dimensions::Vector3;
    use nyx::md::ui::*;
    use nyx::od::ui::GroundStation;
    use nyx::time::Epoch;

    let cosm = Cosm::de438_gmat();

    // Approx landing epoch
    // let landing_epoch = Epoch::from_gregorian_tai_at_noon(2023, 11, 25);
    // From DRM
    // 2023 Nov 25 14:20:35.972000 TDB
    let landing_epoch = Epoch::from_gregorian_tai_hms(2023, 11, 25, 14, 20, 36);

    let ls = GroundStation::from_point(
        "Hayworth".to_string(),
        SITE_LAT_DEG,
        SITE_LONG_DEG,
        SITE_HEIGHT_KM,
        cosm.frame("IAU Moon"),
        cosm.clone(),
    );

    // Get the Moon to Sun vector at that time
    let sun_in_moonj2k = cosm.celestial_state(
        Bodies::Sun.ephem_path(),
        landing_epoch,
        cosm.frame("Moon J2000"),
        LightTimeCalc::None,
    );

    // We want the orbital plane of XL-1 to be in the direction of the Sun
    // The negative of that vector allows us to have zero
    let h_hat = -sun_in_moonj2k.r_hat();

    let inc = h_hat[2].acos().to_degrees();
    let n = Vector3::new(0.0, 0.0, 1.0).cross(&h_hat);
    let raan = n[0].acos().to_degrees() + 30.0;
    // Let's define the AoP to be the longitude of the landing site
    let aop = 180.0;

    // DRM
    // let orbit = Orbit::cartesian(
    //     -77.19951157,
    //     637.0714101,
    //     -1639.05629855,
    //     0.72861834,
    //     1.47489041,
    //     0.48984635,
    //     landing_epoch,
    //     cosm.frame("Moon J2000"),
    // );

    let orbit = Orbit::keplerian_apsis_altitude(
        250.0,
        12.0,
        inc,
        raan,
        aop,
        0.0,
        landing_epoch,
        cosm.frame("Moon J2000"),
    );

    // Propagate this back in time by 15 orbits, and then forward in time
    // Currently, a traj can only be generated when propagating forward.
    let prop = Propagator::default(OrbitalDynamics::point_masses(
        &[Bodies::Sun, Bodies::Earth],
        cosm.clone(),
    ));

    let start_orbit = prop.with(orbit).for_duration(-15 * orbit.period()).unwrap();

    // And forward with traj
    let (final_orbit, traj) = prop
        .with(start_orbit)
        .for_duration_with_traj(15 * orbit.period())
        .unwrap();

    println!("{:x}\t{}", final_orbit, final_orbit.period());

    // Check the eclipse info
    let eloc = EclipseLocator {
        light_source: cosm.frame("Sun J2000"),
        shadow_bodies: vec![cosm.frame("Earth J2000"), cosm.frame("Moon J2000")],
        cosm: cosm.clone(),
    };

    let last_revolution_epoch = final_orbit.epoch() - final_orbit.period();
    println!("Hayworth: {}", ls.to_orbit(last_revolution_epoch));
    let mut total_checks = 0;
    let mut total_eclipses = 0;
    let mut highest_elevation = 0.0_f64;
    let mut highest_epoch = landing_epoch;
    for state in traj.every(0.2 * TimeUnit::Minute) {
        total_checks += 1;
        if eloc.compute(&state) != EclipseState::Visibilis {
            total_eclipses += 1;
        }

        if state.epoch() > last_revolution_epoch {
            // Check the elevation with the landing site
            let (elevation, _, _) = ls.elevation_of(&state);
            if elevation >= highest_elevation {
                highest_elevation = elevation;
                highest_epoch = state.epoch();
            }
        }
    }

    println!(
        "Percentage eclipse: {}",
        total_eclipses * 100 / total_checks
    );

    let (landing_elevation, _, _) = ls.elevation_of(&traj.at(landing_epoch).unwrap());
    println!("Max elevation during last orbit: {} degrees at {}\nLanding epoch planned for {}, elevation then is {}", highest_elevation, highest_epoch, landing_epoch, landing_elevation);
}
