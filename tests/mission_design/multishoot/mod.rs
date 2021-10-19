extern crate nyx_space as nyx;

use nyx::dynamics::guidance::Thruster;
use nyx::md::ui::*;
use nyx::opti::multishoot::*;
use std::str::FromStr;

#[test]
fn landing_demo() {
    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    const SITE_LAT_DEG: f64 = -86.798;
    const SITE_LONG_DEG: f64 = -21.150;
    const SITE_HEIGHT_KM: f64 = 0.4;
    let cosm = Cosm::de438_gmat();
    let moonj2k = cosm.frame("Luna");

    /* *** */
    /* Define the landing epoch, landing site, and initial orbit. */
    /* *** */
    // Epoch of passage at zenith
    let e = Epoch::from_str("2023-11-25T03:02:14.904144287").unwrap();
    // Landing site
    let ls = Orbit::from_geodesic(
        SITE_LAT_DEG,
        SITE_LONG_DEG,
        SITE_HEIGHT_KM,
        e,
        cosm.frame("IAU Moon"),
    );

    let start_orbit = Orbit::keplerian(
        5873.100000,
        0.661491,
        69.327288,
        176.514124,
        283.206000,
        144.089731,
        e,
        moonj2k,
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
    println!("Start: {}", xl1);

    /* *** */
    /* Run the differential corrector for the initial guess of the velocity vector. */
    /* *** */
    // Convert the landing site into the same frame as the spacecraft and use that as targeting values
    let ls_luna = cosm.frame_chg(&ls, moonj2k);
    println!("LANDING SITE: {}", ls_luna);

    let prop = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()));

    let pdi_start = prop.with(xl1).for_duration(-17 * TimeUnit::Minute).unwrap();

    // And run the multiple shooting algorithm

    let mut opti = MultipleShooting::equidistant_nodes(pdi_start, ls, 3, &prop).unwrap();
    opti.solve(CostFunction::MinimumEnergy).unwrap();
}
