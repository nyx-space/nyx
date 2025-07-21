extern crate nyx_space as nyx;

use anise::constants::frames::MOON_J2000;
use nyx::cosmic::{Orbit, Spacecraft};
use nyx::dynamics::{Drag, OrbitalDynamics, SolarPressure, SpacecraftDynamics};
use nyx::linalg::Vector6;
use nyx::propagators::Propagator;
use nyx::time::{Epoch, Unit};
use nyx::utils::rss_orbit_vec_errors;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use nyx_space::md::prelude::{Interpolatable, State};
use rstest::*;
use std::sync::Arc;

use crate::propagation::GMAT_EARTH_GM;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn srp_earth_full_vis(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24 * Unit::Day;

    // Define the dynamics

    let srp = SolarPressure::default(eme2k, almanac.clone()).unwrap();

    let dry_mass = 300.0;

    // Create a spacecraft with SRP model
    let sc_dyn = SpacecraftDynamics::from_model(OrbitalDynamics::two_body(), srp);
    println!("{orbit:x}");

    let sc = Spacecraft::from_srp_defaults(orbit, dry_mass, 16.0);

    let setup = Propagator::default(sc_dyn);
    let mut prop = setup.with(sc, almanac);
    let final_state = prop.for_duration(prop_time).unwrap();

    println!("{final_state}");

    // GMAT result
    let rslt = Vector6::new(
        -10269.72958057943,
        -22135.59895717367,
        0.008121155161511498,
        3.666192351652537,
        -1.699972409878607,
        -8.640_729_256_514_233e-7,
    );

    let (err_r, err_v) = rss_orbit_vec_errors(&final_state.orbit.to_cartesian_pos_vel(), &rslt);
    println!(
        "Error accumulated in full sunlight over {} : {:.6} m \t{:.6} m/s",
        prop_time,
        err_r * 1e3,
        err_v * 1e3
    );
    assert!(err_r < 5e-4, "position error too large for SRP");
    assert!(err_v < 9e-8, "velocity error too large for SRP");
}

#[rstest]
fn srp_earth_leo(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(7_000.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24 * Unit::Day;

    // Define the dynamics

    let srp = SolarPressure::default(eme2k, almanac.clone()).unwrap();

    let dry_mass = 300.0;

    // Create a spacecraft with SRP model
    let sc_dyn = SpacecraftDynamics::from_model(OrbitalDynamics::two_body(), srp);
    println!("{orbit:x}");

    let sc = Spacecraft::from_srp_defaults(orbit, dry_mass, 16.0);

    let setup = Propagator::default(sc_dyn);
    let mut prop = setup.with(sc, almanac);
    let final_state = prop.for_duration(prop_time).unwrap();

    println!("{final_state}");

    // GMAT result
    let rslt = Vector6::new(
        791.0288295053131,
        -6955.409986553813,
        -0.02614433515330551,
        7.497359631253262,
        0.8535219376877066,
        9.281_283_498_115_046e-5,
    );

    let (err_r, err_v) = rss_orbit_vec_errors(&final_state.orbit.to_cartesian_pos_vel(), &rslt);
    println!(
        "Error accumulated in circular equatorial LEO (with penumbras) over {} : {:.6} m \t{:.6} m/s",
        prop_time,
        err_r * 1e3,
        err_v * 1e3
    );
    assert!(err_r < 7e-3, "position error too large for SRP");
    assert!(err_v < 8e-6, "velocity error too large for SRP");
}

#[rstest]
fn srp_earth_meo_ecc_inc(almanac: Arc<Almanac>) {
    use std::env::var as envvar;

    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(14_000.0, 0.5, 20.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24 * Unit::Day;

    // Define the dynamics

    let srp = SolarPressure::default(eme2k, almanac.clone()).unwrap();

    let dry_mass = 300.0;

    // Create a spacecraft with SRP model
    let sc_dyn = SpacecraftDynamics::from_model(OrbitalDynamics::two_body(), srp);
    println!("{orbit:x}");

    let sc = Spacecraft::from_srp_defaults(orbit, dry_mass, 16.0);

    let setup = Propagator::default(sc_dyn.clone());
    let mut prop = setup.with(sc, almanac.clone());
    let final_state = prop.for_duration(prop_time).unwrap();

    println!("{final_state}");

    // GMAT result
    let rslt = Vector6::new(
        -10536.43092609598,
        -11023.53096010533,
        -4012.155888856515,
        4.584027977980024,
        -0.9729922811514549,
        -0.3541244448596867,
    );

    let (err_r, err_v) = rss_orbit_vec_errors(&final_state.orbit.to_cartesian_pos_vel(), &rslt);
    println!(
        "Error accumulated in ecc+inc MEO (with penumbras) over {} : {:.6} m \t{:.6} m/s",
        prop_time,
        err_r * 1e3,
        err_v * 1e3
    );
    assert!(err_r < 2e-3, "position error too large for SRP");
    assert!(err_v < 1e-6, "velocity error too large for SRP");

    // Test that we can specify an integration frame in the options and that the result is the same.
    // Note that we also test here that we're setting the GM and shape of the integration frame as defined
    // and not just grabbing that data from the almanac.
    let orbit = almanac.transform_to(orbit, MOON_J2000, None).unwrap();
    println!("{orbit:x}");

    let sc_moon = Spacecraft::from_srp_defaults(orbit, dry_mass, 16.0);

    let mut setup_moon = Propagator::default(sc_dyn);
    setup_moon.opts.integration_frame = Some(eme2k);
    let mut prop = setup_moon.with(sc_moon, almanac.clone());
    let final_moon_state = prop.for_duration(prop_time).unwrap();
    assert!(
        final_moon_state.frame().ephem_origin_match(MOON_J2000),
        "expected a result in the Moon frame"
    );
    println!("{final_moon_state}");

    let final_earth_orbit = almanac
        .transform_to(final_moon_state.orbit, EARTH_J2000, None)
        .unwrap();

    let (fw_err_r, fw_err_v) =
        rss_orbit_vec_errors(&final_earth_orbit.to_cartesian_pos_vel(), &rslt);
    println!(
        "Error accumulated in ecc+inc MEO (with penumbras) over {} after integration frame swap: {:.6} m \t{:.6} m/s",
        prop_time,
        fw_err_r * 1e3,
        fw_err_v * 1e3
    );
    assert!(dbg!((fw_err_r - err_r).abs() / err_r) < 1e-9);
    assert!(dbg!((fw_err_v - err_v).abs() / err_v) < 1e-12);

    // Compare the case with the hyperdual EOMs (computation uses another part of the code)
    let mut prop = setup.with(sc.with_stm(), almanac);
    let final_state_dual = prop.for_duration(prop_time).unwrap();

    let (err_r, err_v) = rss_orbit_vec_errors(
        &final_state.orbit.to_cartesian_pos_vel(),
        &final_state_dual.orbit.to_cartesian_pos_vel(),
    );
    println!(
        "Error between reals and duals accumulated over {} : {:.3e} m \t{:.3e} m/s",
        prop_time,
        err_r * 1e3,
        err_v * 1e3
    );
    // This should be zero ... but I'm guessing that a successive set of rounding leads to the small accumulation we see
    // So we're allowing for 20 micrometers of difference over 24 days, or less than 1 picometer per second of integration time
    match envvar("CI") {
        Ok(_) => {
            // We're running on Gitlab. It seems to have more rounding error than my computer...
            assert!(
                err_r < 1e-3,
                "Error between reals and duals too large for SRP for CI"
            );
            assert!(
                err_v < 1e-6,
                "Error between reals and duals too large for SRP for CI"
            );
        }
        Err(_) => {
            // Running on a better machine, allow less error
            assert!(
                err_r < 2e-8,
                "Error between reals and duals too large for SRP"
            );
            assert!(
                err_v < 1e-11,
                "Error between reals and duals too large for SRP"
            );
        }
    }
}

#[rstest]
fn exp_drag_earth(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24 * Unit::Day;

    // Define the dynamics

    let srp = SolarPressure::default(eme2k, almanac.clone()).unwrap();
    let drag = Drag::earth_exp(almanac.clone()).unwrap();

    let dry_mass = 300.0;

    // Build a spacecraft with SRP and Drag enabled.
    let sc_dyn = SpacecraftDynamics::from_models(OrbitalDynamics::two_body(), vec![srp, drag]);
    println!("{orbit:x}");

    let sc = Spacecraft::from_srp_defaults(orbit, dry_mass, 1.0).with_drag(1.0, 2.0);

    let setup = Propagator::default_dp78(sc_dyn);
    let mut prop = setup.with(sc, almanac);
    prop.for_duration(prop_time).unwrap();

    let final_state = prop.state;
    println!("{final_state}");
    println!("{}", final_state.orbit);
}

#[rstest]
fn std_atm_drag_earth(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24 * Unit::Day;

    // Define the dynamics

    let srp = SolarPressure::default(eme2k, almanac.clone()).unwrap();
    let drag = Drag::std_atm1976(almanac.clone()).unwrap();

    let dry_mass = 300.0;

    // Build a spacecraft with SRP and Drag enabled.
    let sc_dyn = SpacecraftDynamics::from_models(OrbitalDynamics::two_body(), vec![srp, drag]);
    println!("{orbit:x}");

    let sc = Spacecraft::from_srp_defaults(orbit, dry_mass, 1.0).with_drag(1.0, 2.0);

    let setup = Propagator::default_dp78(sc_dyn);
    let mut prop = setup.with(sc, almanac);
    prop.for_duration(prop_time).unwrap();

    let final_state = prop.state;
    println!("{final_state}");
    println!("{}", final_state.orbit);

    /*
    Test: compared with exponential drag model, and the final states are similar:

    exp_drag_earth
    [Earth J2000] 2000-01-25T00:00:00 TAI   sma = 24394.167595 km   ecc = 0.000019  inc = 0.000001 deg      raan = 299.937993 deg   aop = 264.180317 deg    ta = 42.152440 deg      300 kg
    [Earth J2000] 2000-01-25T00:00:00 TAI   position = [-9816.442834, -22331.499423, -0.000477] km  velocity = [3.700562, -1.626744, 0.000000] km/s

    std_atm_drag_earth
    [Earth J2000] 2000-01-25T00:00:00 TAI   sma = 24396.000574 km   ecc = 0.000020  inc = 0.000001 deg      raan = 299.400845 deg   aop = 264.478503 deg    ta = 41.265314 deg      300 kg
    [Earth J2000] 2000-01-25T00:00:00 TAI   position = [-10254.183112, -22135.911958, -0.000484] km velocity = [3.667742, -1.699095, 0.000000] km/s

    */
}

#[rstest]
fn std_atm_drag_earth_low(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::try_keplerian_altitude(300.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k).unwrap();

    let prop_time = 24 * Unit::Day;

    // Define the dynamics

    let srp = SolarPressure::default(eme2k, almanac.clone()).unwrap();
    let drag = Drag::std_atm1976(almanac.clone()).unwrap();

    let dry_mass = 300.0;

    // Build a spacecraft with SRP and Drag enabled.
    let sc_dyn = SpacecraftDynamics::from_models(OrbitalDynamics::two_body(), vec![srp, drag]);
    println!("{orbit:x}");

    let sc = Spacecraft::from_srp_defaults(orbit, dry_mass, 1.0).with_drag(1.0, 2.0);

    let setup = Propagator::default_dp78(sc_dyn);
    let mut prop = setup.with(sc, almanac);
    prop.for_duration(prop_time).unwrap();

    let final_state = prop.state;
    println!("{final_state}");
    println!("{}", final_state.orbit);

    /*
    Test: compared with exponential drag model, and the final states are similar:

    exp_drag_earth
    [Earth J2000] 2000-01-25T00:00:00 TAI   sma = 24394.167595 km   ecc = 0.000019  inc = 0.000001 deg      raan = 299.937993 deg   aop = 264.180317 deg    ta = 42.152440 deg      300 kg
    [Earth J2000] 2000-01-25T00:00:00 TAI   position = [-9816.442834, -22331.499423, -0.000477] km  velocity = [3.700562, -1.626744, 0.000000] km/s

    std_atm_drag_earth
    [Earth J2000] 2000-01-25T00:00:00 TAI   sma = 24396.000574 km   ecc = 0.000020  inc = 0.000001 deg      raan = 299.400845 deg   aop = 264.478503 deg    ta = 41.265314 deg      300 kg
    [Earth J2000] 2000-01-25T00:00:00 TAI   position = [-10254.183112, -22135.911958, -0.000484] km velocity = [3.667742, -1.699095, 0.000000] km/s

    */
}
