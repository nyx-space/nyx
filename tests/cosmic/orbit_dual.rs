extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use nyx::cosmic::Orbit;
use nyx::time::Epoch;

use crate::propagation::GMAT_EARTH_GM;
use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;

#[fixture]
fn almanac() -> Almanac {
    use crate::test_almanac;
    test_almanac()
}

#[rstest]
fn orbit_dual_test(almanac: Almanac) {
    use nyx::cosmic::OrbitDual;
    use nyx::md::StateParameter;

    let eme2k = almanac
        .frame_info(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let dt = Epoch::from_mjd_tai(21_545.0);
    let cart = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 1.0, dt, eme2k,
    );

    println!("{cart:x}");

    let cart_dual = OrbitDual::from(cart);
    println!("{}", eme2k.mu_km3_s2().unwrap());
    println!("|v| = {}", cart_dual.vmag_km_s());
    println!(
        "|v|^2 = {}",
        cart_dual.vmag_km_s().dual * cart_dual.vmag_km_s().dual
    );
    println!("x = {}", cart_dual.x);
    println!("y = {}", cart_dual.y);
    println!("z = {}", cart_dual.z);
    println!("\\xi = {}\n\n", cart_dual.energy_km2_s2().unwrap());

    // Now print the table
    let params = vec![
        StateParameter::AoL,
        StateParameter::AoP,
        StateParameter::Apoapsis,
        StateParameter::C3,
        StateParameter::Declination,
        StateParameter::EccentricAnomaly,
        StateParameter::Eccentricity,
        StateParameter::Energy,
        StateParameter::FlightPathAngle,
        StateParameter::Height,
        StateParameter::Latitude,
        StateParameter::Longitude,
        StateParameter::Hmag,
        StateParameter::HX,
        StateParameter::HY,
        StateParameter::HZ,
        StateParameter::Inclination,
        StateParameter::MeanAnomaly,
        StateParameter::Periapsis,
        StateParameter::RightAscension,
        StateParameter::RAAN,
        StateParameter::Rmag,
        StateParameter::SemiParameter,
        StateParameter::SMA,
        StateParameter::TrueLongitude,
        StateParameter::Vmag,
        StateParameter::X,
        StateParameter::Y,
        StateParameter::Z,
        StateParameter::VX,
        StateParameter::VY,
        StateParameter::VZ,
    ];
    for param in params {
        let dual = cart_dual.partial_for(param).unwrap();
        println!(
            "{:?} & {} & {} & {} & {} & {} & {} \\\\",
            param,
            if dual.wtr_x().abs() > 1e-2 {
                "large"
            } else if dual.wtr_x().abs() > 0.0 {
                "small"
            } else {
                "none"
            },
            if dual.wtr_y().abs() > 1e-2 {
                "large"
            } else if dual.wtr_y().abs() > 0.0 {
                "small"
            } else {
                "none"
            },
            if dual.wtr_z().abs() > 1e-2 {
                "large"
            } else if dual.wtr_z().abs() > 0.0 {
                "small"
            } else {
                "none"
            },
            if dual.wtr_vx().abs() > 1e-2 {
                "large"
            } else if dual.wtr_vx().abs() > 0.0 {
                "small"
            } else {
                "none"
            },
            if dual.wtr_vy().abs() > 1e-2 {
                "large"
            } else if dual.wtr_vy().abs() > 0.0 {
                "small"
            } else {
                "none"
            },
            if dual.wtr_vz().abs() > 1e-2 {
                "large"
            } else if dual.wtr_vz().abs() > 0.0 {
                "small"
            } else {
                "none"
            },
            // dual,
        );
    }
}
