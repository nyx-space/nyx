extern crate nyx_space as nyx;

use anise::constants::frames::{IAU_EARTH_FRAME, MOON_J2000};
use anise::f64_eq_tol;
use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use nyx::cosmic::{Orbit, Spacecraft};
use nyx::dynamics::{
    AtmDensity, Drag, GravityField, OrbitalDynamics, SolarPressure, SpacecraftDynamics,
};
use nyx::io::gravity::*;
use nyx::io::space_weather::{SpaceWeatherData, StaticSpaceWeather};
use nyx::linalg::Vector6;
use nyx::md::prelude::{Interpolatable, State};
use nyx::propagators::Propagator;
use nyx::time::{Epoch, Unit};
use nyx::utils::rss_orbit_vec_errors;
use nyx_space::cosmic::{DragData, Mass};
use nyx_space::od::Dynamics;
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
        .frame_info(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24 * Unit::Day;

    // Define the dynamics

    let srp = SolarPressure::default_flux(eme2k, &almanac).unwrap();

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
        .frame_info(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(7_000.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24 * Unit::Day;

    // Define the dynamics

    let srp = SolarPressure::default_flux(eme2k, &almanac).unwrap();

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
        .frame_info(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(14_000.0, 0.5, 20.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24 * Unit::Day;

    // Define the dynamics

    let srp = SolarPressure::default_flux(eme2k, &almanac).unwrap();

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
        .frame_info(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24 * Unit::Day;

    // Define the dynamics

    let srp = SolarPressure::default_flux(eme2k, &almanac).unwrap();
    let drag = Drag::earth_exp(&almanac).unwrap();

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
        .frame_info(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::keplerian(24396.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k);

    let prop_time = 24 * Unit::Day;

    // Define the dynamics

    let srp = SolarPressure::default_flux(eme2k, &almanac).unwrap();
    let drag = Drag::std_atm1976(&almanac).unwrap();

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
        .frame_info(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let dt = Epoch::from_gregorian_tai_at_midnight(2000, 1, 1);

    let orbit = Orbit::try_keplerian_altitude(300.0, 0.0, 0.0, 0.0, 0.0, 0.0, dt, eme2k).unwrap();

    let prop_time = 24 * Unit::Day;

    // Define the dynamics

    let srp = SolarPressure::default_flux(eme2k, &almanac).unwrap();
    let drag = Drag::std_atm1976(&almanac).unwrap();

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
fn test_prop_nrlmsise00_from_weather(almanac: Arc<Almanac>) {
    use std::path::PathBuf;

    let manifest_dir: PathBuf = [env!("CARGO_MANIFEST_DIR"), "../data/01_planetary"]
        .iter()
        .collect();

    let weather = SpaceWeatherData::from_csv_file(
        manifest_dir.join("SpaceWeather-2021-01-01_2026-09-06.csv.gz"),
        StaticSpaceWeather::SolarAverage(),
    )
    .unwrap();

    // Test that we can serialize the space weather BTreeMap
    let weather_toml =
        toml::to_string(&weather).expect("should be able to serialize the weather as TOML");
    let weather_rtn: SpaceWeatherData = toml::from_str(&weather_toml).unwrap();

    assert_eq!(weather_rtn, weather);

    let eme2k = almanac
        .frame_info(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let epoch = Epoch::from_gregorian_utc(2024, 3, 20, 12, 0, 0, 0);

    let orbit = Orbit::keplerian(6778.137, 0.001, 51.6, 0.0, 0.0, 0.0, epoch, eme2k);
    let spacecraft = Spacecraft::from_drag_defaults(orbit, 100.0, 1.0);

    let dynamics = SpacecraftDynamics::from_models(
        OrbitalDynamics::two_body(),
        vec![Arc::new(Drag {
            density: AtmDensity::NRLMSISE00 { weather },
            frame: IAU_EARTH_FRAME,
            estimate: false,
        })],
    );

    let accel_km_s2_with_msise00 = dynamics
        .eom(0.0, &spacecraft.to_vector(), &spacecraft, &almanac)
        .unwrap();

    let accel_km_s2_with_msise00 = accel_km_s2_with_msise00.fixed_rows::<3>(3);

    let state_with_msise00 = Propagator::default(dynamics)
        .with(spacecraft, almanac.clone())
        .for_duration(1.0 * Unit::Hour)
        .unwrap();

    // Build another propagator without drag and ensure drag was correctly applied.
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());

    let accel_km_s2_without_msise00 = dynamics
        .eom(0.0, &spacecraft.to_vector(), &spacecraft, &almanac)
        .unwrap();
    let accel_km_s2_without_msise00 = accel_km_s2_without_msise00.fixed_rows::<3>(3);

    let state_without_msise00 = Propagator::default(dynamics)
        .with(spacecraft, almanac.clone())
        .for_duration(1.0 * Unit::Hour)
        .unwrap();

    assert_ne!(state_without_msise00, state_with_msise00);

    let diff_accel_m_s2 = (accel_km_s2_without_msise00 - accel_km_s2_with_msise00) * 1e3;

    println!("Accel diff in m/s^2: {diff_accel_m_s2:.3e}");

    f64_eq_tol!(
        diff_accel_m_s2.norm(),
        5e-6,
        1e-6,
        "expected an acceleration difference"
    );
}

#[rstest]
fn val_ioastro_nrlmsise00(almanac: Arc<Almanac>) {
    let weather = SpaceWeatherData::from_static_weather(StaticSpaceWeather::Custom {
        f107: 150.0,
        ap: 18.6,
        kp: 3.0,
    });

    let eme2k = almanac
        .frame_info(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);
    let iau_earth = almanac.frame_info(IAU_EARTH_FRAME).unwrap();
    // .with_mu_km3_s2(GMAT_EARTH_GM);

    let epoch = Epoch::from_gregorian_utc_hms(2024, 2, 29, 1, 2, 3);

    let orbit = Orbit::new(
        -6210.6003575395861844,
        -445.9010437211380236,
        2353.4286399045840881,
        -0.6027299684384100,
        -7.2371397713442924,
        -2.7084864602730194,
        epoch,
        eme2k,
    );

    println!("Initial: {orbit:x}");

    let spacecraft = Spacecraft::builder()
        .orbit(orbit)
        .mass(Mass::from_dry_mass(1000.0))
        .drag(DragData {
            area_m2: 20.0,
            coeff_drag: 2.2,
        })
        .build();

    // Define the dynamics
    let earth_sph_harm = GravityFieldData::from_j2(EARTH_J2, iau_earth);
    let j2_mdl = GravityField::new(earth_sph_harm);

    let dynamics = SpacecraftDynamics::from_models(
        OrbitalDynamics::from_model(j2_mdl),
        vec![Arc::new(Drag {
            density: AtmDensity::NRLMSISE00 { weather },
            frame: iau_earth,
            estimate: false,
        })],
    );

    let final_state = Propagator::default(dynamics)
        .with(spacecraft, almanac.clone())
        .for_duration(1.0 * Unit::Day)
        .unwrap()
        .orbit;

    let expected_state = Orbit::new(
        -5904.4748605974446036,
        1330.9542946601800395,
        2687.0267716822277180,
        -2.6688074464373877,
        -6.9600802661844208,
        -2.2412448735527049,
        Epoch::from_gregorian_utc_hms(2024, 3, 1, 1, 2, 3),
        eme2k,
    );

    let ric_error = expected_state.ric_difference(&final_state).unwrap();

    println!("=== Cartesian ===\nGot:  {final_state}");
    println!("Want: {expected_state}");

    println!("=== Keplerian ===\nGot:  {final_state:x}");
    println!("Want: {expected_state:x}");

    println!(
        "RIC pos error (km) = {:.6}\n{:.6}",
        ric_error.rmag_km(),
        ric_error.radius_km
    );
    println!(
        "RIC vel error (km/s) = {:.6}\n{:.6}",
        ric_error.vmag_km_s(),
        ric_error.velocity_km_s
    );

    // NOTE There is a small in-track error that seems to be due to the Ap value.
    // Tuning that value from a guess of 15.0 to 18.6 reduces that in track error
    // from 14 km to < 1 km. While the values don't match as closesly as other
    // validation cases, I still consider this model validated because a single
    // parameter change in the Space Weather is sufficient to meet the test case.

    assert!(dbg!(ric_error.rmag_km()) < 0.72);
    assert!(dbg!(ric_error.vmag_km_s()) < 2e-3);
}

#[rstest]
fn regression_harris_drag(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_info(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let iau_earth = almanac
        .frame_info(IAU_EARTH_FRAME)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let epoch = Epoch::from_gregorian_utc_hms(2024, 2, 29, 1, 2, 3);

    let orbit = Orbit::new(
        -6210.6003575395861844,
        -445.9010437211380236,
        2353.4286399045840881,
        -0.6027299684384100,
        -7.2371397713442924,
        -2.7084864602730194,
        epoch,
        eme2k,
    );

    println!("Initial: {orbit:x}");

    let spacecraft = Spacecraft::builder()
        .orbit(orbit)
        .mass(Mass::from_dry_mass(1000.0))
        .drag(DragData {
            area_m2: 20.0,
            coeff_drag: 2.2,
        })
        .build();

    // Define the dynamics
    let earth_sph_harm = GravityFieldData::from_j2(EARTH_J2, iau_earth);
    let j2_mdl = GravityField::new(earth_sph_harm);

    let dynamics = SpacecraftDynamics::from_models(
        OrbitalDynamics::from_model(j2_mdl),
        vec![Arc::new(Drag {
            density: AtmDensity::HarrisPriester { n_parameter: 2 },
            frame: iau_earth,
            estimate: false,
        })],
    );

    let final_state = Propagator::default(dynamics)
        .with(spacecraft, almanac.clone())
        .for_duration(1.0 * Unit::Day)
        .unwrap()
        .orbit;

    let expected_state = Orbit::new(
        -5866.260366,
        1429.893370,
        2719.496468,
        -2.783483,
        -6.933258,
        -2.187250,
        Epoch::from_gregorian_utc_hms(2024, 3, 1, 1, 2, 3),
        eme2k,
    );

    let ric_error = expected_state.ric_difference(&final_state).unwrap();

    println!("=== Cartesian ===\nGot:  {final_state}");
    println!("Want: {expected_state}");

    println!("=== Keplerian ===\nGot:  {final_state:x}");
    println!("Want: {expected_state:x}");

    println!(
        "RIC pos error (km) = {:.6}\n{:.6}",
        ric_error.rmag_km(),
        ric_error.radius_km
    );
    println!(
        "RIC vel error (km/s) = {:.6}\n{:.6}",
        ric_error.vmag_km_s(),
        ric_error.velocity_km_s
    );

    assert!(dbg!(ric_error.rmag_km()) < 2e-5);
    assert!(dbg!(ric_error.vmag_km_s()) < 2e-6);

    // For checking purposes, add the diff to NRLMSISE00.

    let expected_state = Orbit::new(
        -5904.4748605974446036,
        1330.9542946601800395,
        2687.0267716822277180,
        -2.6688074464373877,
        -6.9600802661844208,
        -2.2412448735527049,
        Epoch::from_gregorian_utc_hms(2024, 3, 1, 1, 2, 3),
        eme2k,
    );
    let ric_error = expected_state.ric_difference(&final_state).unwrap();
    println!("*** NRLMSISE-00 COMPARISON ***");
    println!("=== Cartesian ===\nGot:  {final_state}");
    println!("Want: {expected_state}");

    println!("=== Keplerian ===\nGot:  {final_state:x}");
    println!("Want: {expected_state:x}");

    println!(
        "RIC pos error (km) = {:.6}\n{:.6}",
        ric_error.rmag_km(),
        ric_error.radius_km
    );
    println!(
        "RIC vel error (km/s) = {:.6}\n{:.6}",
        ric_error.vmag_km_s(),
        ric_error.velocity_km_s
    );
}
