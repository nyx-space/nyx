extern crate nalgebra as na;
extern crate nyx_space as nyx;

use anise::constants::celestial_objects::{EARTH, JUPITER_BARYCENTER, MOON, SUN};
use anise::constants::frames::IAU_EARTH_FRAME;
use hifitime::MJD_J2000;
use nalgebra::{Const, OMatrix};
use nyx::cosmic::{assert_orbit_eq_or_abs, Orbit};
use nyx::dynamics::{Dynamics, OrbitalDynamics, PointMasses, SpacecraftDynamics};
use nyx::linalg::Vector6;
use nyx::time::{Epoch, Unit};
use nyx::utils::{rss_orbit_errors, rss_orbit_vec_errors};
use nyx::State;
use nyx::{propagators::*, Spacecraft};

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[fixture]
fn almanac_gmat() -> Arc<Almanac> {
    use crate::test_almanac_gmat_arcd;
    test_almanac_gmat_arcd()
}

#[allow(clippy::identity_op)]
#[rstest]
fn energy_conservation(almanac: Arc<Almanac>) {
    let prop_time = 1 * Unit::Day;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let dt = Epoch::from_gregorian_utc_hms(2022, 2, 15, 17, 30, 37);
    let start_state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );

    let rslt = Orbit::cartesian(
        -5_971.194_376_797_643,
        3_945.517_912_574_178_4,
        2_864.620_957_744_429_2,
        0.049_083_101_605_507_95,
        -4.185_084_125_817_658,
        5.848_947_462_472_877,
        dt + prop_time,
        eme2k,
    );

    let rk89_final = Propagator::new(
        SpacecraftDynamics::new(OrbitalDynamics::two_body()),
        IntegratorMethod::RungeKutta89,
        IntegratorOptions::default(),
    )
    .with(Spacecraft::from(start_state), almanac.clone())
    .for_duration(prop_time)
    .unwrap();

    let rk89_energy_bleed =
        rk89_final.orbit.energy_km2_s2().unwrap() - start_state.energy_km2_s2().unwrap();

    println!(
        "[RK89] ==> energy_conservation absolute errors with RK89 val state\tenergy bleed = {rk89_energy_bleed:e}"
    );
    let delta = rk89_final.orbit.to_cartesian_pos_vel() - rslt.to_cartesian_pos_vel();
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    let dp78_final = Propagator::new(
        SpacecraftDynamics::new(OrbitalDynamics::two_body()),
        IntegratorMethod::DormandPrince78,
        IntegratorOptions::default(),
    )
    .with(start_state.into(), almanac)
    .for_duration(prop_time)
    .unwrap();

    let dp78_energy_bleed =
        dp78_final.orbit.energy_km2_s2().unwrap() - start_state.energy_km2_s2().unwrap();

    println!(
        "[DP78] ==> energy_conservation absolute errors with RK89 val state\tenergy bleed = {dp78_energy_bleed:e}"
    );
    let delta = dp78_final.orbit.to_cartesian_pos_vel() - rslt.to_cartesian_pos_vel();
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();
}

#[allow(clippy::identity_op)]
#[rstest]
fn val_two_body_dynamics(almanac: Arc<Almanac>) {
    let prop_time = 1 * Unit::Day;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let dt = Epoch::from_mjd_tai(MJD_J2000);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );

    let rslt = Orbit::cartesian(
        -5971.194375461378,
        3945.517831291771,
        2864.6210708007134,
        0.04908320163379219,
        -4.1850841921806206,
        5.848947414864886,
        dt + prop_time,
        eme2k,
    );

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());
    let setup = Propagator::rk89(dynamics, IntegratorOptions::default());
    let mut prop = setup.with(state.into(), almanac);
    prop.for_duration(prop_time).unwrap();

    let details = prop.details;

    println!("==> val_two_body_dynamics absolute errors");
    let delta = prop.state.orbit.to_cartesian_pos_vel() - rslt.to_cartesian_pos_vel();
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!("\nFinal step info {details}");

    assert_orbit_eq_or_abs(&prop.state.orbit, &rslt, 2e-9, "two body prop failed");

    // And now do the backprop by re-initializing a propagator to ensure correct step size
    prop.for_duration(-prop_time).unwrap();
    let (err_r, err_v) = rss_orbit_vec_errors(
        &prop.state.orbit.to_cartesian_pos_vel(),
        &state.to_cartesian_pos_vel(),
    );
    println!("RTN:  {}\nINIT: {}\n{:x}", prop.state, state, state);
    dbg!(err_r);
    assert!(
        err_r < 1e-5,
        "two body back prop failed to return to the initial state in position"
    );
    assert!(
        err_v < 1e-8,
        "two body back prop failed to return to the initial state in velocity"
    );
    assert_eq!(prop.state.epoch(), dt);
    // Forward propagation again to confirm that we can do repeated calls
    prop.for_duration(prop_time).unwrap();
    assert_eq!(prop.state.epoch(), dt + prop_time);
    let (err_r, err_v) = rss_orbit_vec_errors(
        &prop.state.orbit.to_cartesian_pos_vel(),
        &rslt.to_cartesian_pos_vel(),
    );
    assert!(
        err_r < 1e-5,
        "two body back+fwd prop failed to return to the initial state in position"
    );
    assert!(
        err_v < 1e-8,
        "two body back+fwd prop failed to return to the initial state in velocity"
    );
}

#[allow(clippy::identity_op)]
#[rstest]
fn val_halo_earth_moon_dynamics(almanac_gmat: Arc<Almanac>) {
    let almanac = almanac_gmat;
    /*
    We validate against GMAT after switching the GMAT script to use de438s.bsp. We are using GMAT's default GM values.
    The state in `rslt` is exactly the GMAT output.
    */
    let prop_time = 1 * Unit::Day;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let halo_rcvr = Orbit::cartesian(
        333_321.004_516,
        -76_134.198_887,
        -20_873.831_939,
        0.257_153_712,
        0.930_284_066,
        0.346_177,
        start_time,
        eme2k,
    );

    // GMAT data
    let rslt = Orbit::cartesian(
        345_395.216_758_754_4,
        5_967.890_264_751_025,
        7_350.734_617_702_599,
        0.022_370_754_768_832_33,
        0.957_450_818_399_485_1,
        0.303_172_019_604_272_5,
        start_time + prop_time,
        eme2k,
    );

    let bodies = vec![MOON];
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));

    let setup = Propagator::rk89(
        dynamics,
        IntegratorOptions::with_fixed_step(10 * Unit::Second),
    );
    let mut prop = setup.with(halo_rcvr.into(), almanac);
    prop.for_duration(prop_time).unwrap();
    let (err_r, err_v) = rss_orbit_errors(&prop.state.orbit, &rslt);

    println!("==> val_halo_earth_moon_dynamics absolute errors");
    let delta = prop.state.orbit.to_cartesian_pos_vel() - rslt.to_cartesian_pos_vel();
    for i in 0..3 {
        print!("{:.0e} m\t", delta[i].abs() * 1e3);
    }
    for i in 3..6 {
        print!("{:.0e} m/s\t", delta[i].abs() * 1e3);
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} m\tvel = {:.5e} m/s\ninit\t{}\nfinal\t{}",
        err_r * 1e3,
        err_v * 1e3,
        halo_rcvr,
        prop.state
    );

    assert!(err_r < 5e-5, "multi body failed in position: {err_r:.5e}");
    assert!(err_v < 1e-9, "multi body failed in velocity: {err_v:.5e}");
}

#[allow(clippy::identity_op)]
#[rstest]
fn val_halo_earth_moon_dynamics_adaptive(almanac_gmat: Arc<Almanac>) {
    let almanac = almanac_gmat;
    /*
    We validate against GMAT after switching the GMAT script to use de438s.bsp. We are using GMAT's default GM values.
    The state in `rslt` is exactly the GMAT output.
    */
    let prop_time = 1 * Unit::Day;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2002, 2, 7);

    let halo_rcvr = Orbit::cartesian(
        333_321.004_516,
        -76_134.198_887,
        -20_873.831_939,
        0.257_153_712,
        0.930_284_066,
        0.346_177,
        start_time,
        eme2k,
    );

    let rslt = Vector6::new(
        343_016.028_193_306_2,
        6_118.870_782_679_712,
        9_463.253_311_291_08,
        -0.033_885_504_418_292_03,
        0.961_942_577_960_542_2,
        0.351_738_121_709_363_5,
    );

    let bodies = vec![MOON];
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));

    let setup = Propagator::rk89(dynamics, IntegratorOptions::default());
    let mut prop = setup.with(halo_rcvr.into(), almanac);
    prop.for_duration(prop_time).unwrap();
    let (err_r, err_v) = rss_orbit_vec_errors(&prop.state.orbit.to_cartesian_pos_vel(), &rslt);

    println!("==> val_halo_earth_moon_dynamics_adaptive absolute errors");
    let delta = prop.state.orbit.to_cartesian_pos_vel() - rslt;
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} m\tvel = {:.5e} m/s\ninit\t{}\nfinal\t{}",
        err_r * 1e3,
        err_v * 1e3,
        halo_rcvr,
        prop.state
    );

    assert!(err_r < 1e-6, "multi body failed in position: {err_r:.5e}");
    assert!(err_v < 1e-11, "multi body failed in velocity: {err_v:.5e}");
}

#[allow(clippy::identity_op)]
#[rstest]
fn val_llo_earth_moon_dynamics_adaptive(almanac_gmat: Arc<Almanac>) {
    let almanac = almanac_gmat;
    /*
    We validate against GMAT after switching the GMAT script to use de438s.bsp. We are using GMAT's default GM values.
    The state in `rslt` is exactly the GMAT output.
    */
    let prop_time = 1 * Unit::Day;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2002, 2, 7);

    let llo_xmtr = Orbit::cartesian(
        3.919_869_89e5,
        -7.493_039_70e4,
        -7.022_605_11e4,
        -6.802_604_18e-1,
        1.992_053_61,
        4.369_389_94e-1,
        start_time,
        eme2k,
    );

    // GMAT data
    let rslt = Vector6::new(
        322_883.886_835_433_2,
        97_580.280_858_158,
        -30_871.085_807_431_58,
        -0.934_039_629_727_003_5,
        1.980_106_615_205_608,
        0.472_630_895_504_854_4,
    );

    let bodies = vec![MOON];
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));

    let setup = Propagator::rk89(dynamics, IntegratorOptions::default());
    let mut prop = setup.with(llo_xmtr.into(), almanac);
    prop.for_duration(prop_time).unwrap();
    let (err_r, err_v) = rss_orbit_vec_errors(&prop.state.orbit.to_cartesian_pos_vel(), &rslt);

    println!("==> val_llo_earth_moon_dynamics_adaptive absolute errors");
    let delta = prop.state.orbit.to_cartesian_pos_vel() - rslt;
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} m\tvel = {:.5e} m/s\ninit\t{}\nfinal\t{}",
        err_r * 1e3,
        err_v * 1e3,
        llo_xmtr,
        prop.state
    );

    assert!(err_r < 1e-5, "multi body failed in position: {err_r:.5e}");
    assert!(err_v < 1e-8, "multi body failed in velocity: {err_v:.5e}");
}

#[allow(clippy::identity_op)]
#[rstest]
fn val_halo_multi_body_dynamics(almanac_gmat: Arc<Almanac>) {
    let almanac = almanac_gmat;
    /*
    We validate against GMAT after switching the GMAT script to use de438s.bsp. We are using GMAT's default GM values.
    The state in `rslt` is exactly the GMAT output.
    */
    let prop_time = 1 * Unit::Day;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let halo_rcvr = Orbit::cartesian(
        333_321.004_516,
        -76_134.198_887,
        -20_873.831_939,
        0.257_153_712,
        0.930_284_066,
        0.346_177,
        start_time,
        eme2k,
    );

    // GMAT data
    let rslt = Vector6::new(
        345_350.664_306_402_7,
        5_930.672_402_473_843,
        7_333.283_870_811_47,
        0.021_298_196_465_430_16,
        0.956_678_964_966_812_2,
        0.302_817_582_487_008_6,
    );

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));

    let setup = Propagator::rk89(
        dynamics,
        IntegratorOptions::with_fixed_step(10 * Unit::Second),
    );
    let mut prop = setup.with(halo_rcvr.into(), almanac);
    prop.for_duration(prop_time).unwrap();
    let (err_r, err_v) = rss_orbit_vec_errors(&prop.state.orbit.to_cartesian_pos_vel(), &rslt);

    println!("==> val_halo_multi_body_dynamics absolute errors");
    let delta = prop.state.orbit.to_cartesian_pos_vel() - rslt;
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} m\tvel = {:.5e} m/s\ninit\t{}\nfinal\t{}",
        err_r * 1e3,
        err_v * 1e3,
        halo_rcvr,
        prop.state
    );

    assert!(err_r < 5e-5, "multi body failed in position: {err_r:.5e}");
    assert!(err_v < 1e-9, "multi body failed in velocity: {err_v:.5e}");
}

#[allow(clippy::identity_op)]
#[rstest]
fn val_halo_multi_body_dynamics_adaptive(almanac_gmat: Arc<Almanac>) {
    let almanac = almanac_gmat;
    /*
    We validate against GMAT after switching the GMAT script to use de438s.bsp. We are using GMAT's default GM values.
    The state in `rslt` is exactly the GMAT output.
    */

    let prop_time = 1 * Unit::Day;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2002, 2, 7);

    let halo_rcvr = Orbit::cartesian(
        333_321.004_516,
        -76_134.198_887,
        -20_873.831_939,
        0.257_153_712,
        0.930_284_066,
        0.346_177,
        start_time,
        eme2k,
    );

    // GMAT data
    let rslt = Vector6::new(
        343_063.315_079_726_9,
        6_045.912_866_799_058,
        9_430.044_002_816_507,
        -0.032_841_040_500_475_27,
        0.960_272_613_530_677_2,
        0.350_981_431_322_089_4,
    );

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));

    let setup = Propagator::default(dynamics);
    let mut prop = setup.with(halo_rcvr.into(), almanac);
    prop.for_duration(prop_time).unwrap();
    let (err_r, err_v) = rss_orbit_vec_errors(&prop.state.orbit.to_cartesian_pos_vel(), &rslt);

    println!("==> val_halo_multi_body_dynamics_adaptive absolute errors");
    let delta = prop.state.orbit.to_cartesian_pos_vel() - rslt;
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} m\tvel = {:.5e} m/s\ninit\t{}\nfinal\t{}",
        err_r * 1e3,
        err_v * 1e3,
        halo_rcvr,
        prop.state
    );

    assert!(err_r < 1e-6, "multi body failed in position: {err_r:.5e}");
    assert!(err_v < 1e-11, "multi body failed in velocity: {err_v:.5e}");
}

#[allow(clippy::identity_op)]
#[rstest]
fn val_llo_multi_body_dynamics_adaptive(almanac_gmat: Arc<Almanac>) {
    let almanac = almanac_gmat;
    /*
    We validate against GMAT after switching the GMAT script to use de438s.bsp. We are using GMAT's default GM values.
    The state in `rslt` is exactly the GMAT output.
    */

    let prop_time = 1 * Unit::Day;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2002, 2, 7);

    let llo_xmtr = Orbit::cartesian(
        3.919_869_89e5,
        -7.493_039_70e4,
        -7.022_605_11e4,
        -6.802_604_18e-1,
        1.992_053_61,
        4.369_389_94e-1,
        start_time,
        eme2k,
    );

    // GMAT data
    let rslt = Vector6::new(
        322_931.851_760_741_2,
        97_497.699_738_811_13,
        -30_899.323_820_367_2,
        -0.933_095_202_143_736_8,
        1.978_291_140_770_421,
        0.472_036_197_968_369_3,
    );

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));

    let setup = Propagator::default(dynamics);
    let mut prop = setup.with(llo_xmtr.into(), almanac);
    prop.for_duration(prop_time).unwrap();
    let (err_r, err_v) = rss_orbit_vec_errors(&prop.state.orbit.to_cartesian_pos_vel(), &rslt);

    println!("==> val_llo_multi_body_dynamics_adaptive absolute errors");
    let delta = prop.state.orbit.to_cartesian_pos_vel() - rslt;
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} m\tvel = {:.5e} m/s\ninit\t{}\nfinal\t{}",
        err_r * 1e3,
        err_v * 1e3,
        llo_xmtr,
        prop.state
    );

    assert!(err_r < 2e-6, "multi body failed in position: {err_r:.5e}");
    assert!(err_v < 1e-9, "multi body failed in velocity: {err_v:.5e}");
}

#[allow(clippy::identity_op)]
#[rstest]
fn val_leo_multi_body_dynamics_adaptive_wo_moon(almanac_gmat: Arc<Almanac>) {
    let almanac = almanac_gmat;
    /*
    We validate against GMAT after switching the GMAT script to use de438s.bsp. We are using GMAT's default GM values.
    The state in `rslt` is exactly the GMAT output.
    */

    let prop_time = 1 * Unit::Day;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let leo = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_time, eme2k,
    );

    // GMAT data
    let rslt = Vector6::new(
        -5_971.190_141_842_914,
        3_945.572_972_028_369,
        2_864.554_642_502_679,
        0.049_014_376_371_383_95,
        -4.185_051_832_316_421,
        5.848_971_837_743_221,
    );

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));

    let setup = Propagator::default(dynamics);
    let mut prop = setup.with(leo.into(), almanac);
    prop.for_duration(prop_time).unwrap();
    let (err_r, err_v) = rss_orbit_vec_errors(&prop.state.orbit.to_cartesian_pos_vel(), &rslt);

    println!("==> val_leo_multi_body_dynamics_adaptive_wo_moon absolute errors");
    let delta = prop.state.orbit.to_cartesian_pos_vel() - rslt;
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} m\tvel = {:.5e} m/s\ninit\t{}\nfinal\t{}",
        err_r * 1e3,
        err_v * 1e3,
        leo,
        prop.state
    );

    assert!(err_r < 5e-7, "multi body failed in position: {err_r:.5e}");
    assert!(err_v < 5e-10, "multi body failed in velocity: {err_v:.5e}");
}

#[allow(clippy::identity_op)]
#[rstest]
fn val_leo_multi_body_dynamics_adaptive(almanac_gmat: Arc<Almanac>) {
    let almanac = almanac_gmat;
    /*
    We validate against GMAT after switching the GMAT script to use de438s.bsp. We are using GMAT's default GM values.
    The state in `rslt` is exactly the GMAT output.
    */

    let prop_time = 1 * Unit::Day;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let leo = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, start_time, eme2k,
    );

    // GMAT data
    let rslt = Vector6::new(
        -5_971.190_491_039_24,
        3_945.529_211_711_111,
        2_864.613_171_213_388,
        0.049_086_325_111_121_92,
        -4.185_065_854_096_239,
        5.848_960_991_136_447,
    );

    let bodies = vec![SUN, JUPITER_BARYCENTER];
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));

    let setup = Propagator::default(dynamics);
    let mut prop = setup.with(leo.into(), almanac);
    prop.for_duration(prop_time).unwrap();
    let (err_r, err_v) = rss_orbit_vec_errors(&prop.state.orbit.to_cartesian_pos_vel(), &rslt);

    println!("==> val_leo_multi_body_dynamics_adaptive absolute errors");
    let delta = prop.state.orbit.to_cartesian_pos_vel() - rslt;
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    println!(
        "RSS errors:\tpos = {:.5e} m\tvel = {:.5e} m/s\ninit\t{}\nfinal\t{}",
        err_r * 1e3,
        err_v * 1e3,
        leo,
        prop.state
    );

    assert!(err_r < 3e-6, "multi body failed in position: {err_r:.5e}");
    assert!(err_v < 3e-9, "multi body failed in velocity: {err_v:.5e}");
}

#[allow(clippy::identity_op)]
#[rstest]
fn two_body_dual(almanac: Arc<Almanac>) {
    // This is a duplicate of the differentials test in hyperdual.

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let init = Orbit::cartesian(
        -9_042.862_233_600_335,
        18_536.333_069_123_244,
        6_999.957_069_486_411_5,
        -3.288_789_003_770_57,
        -2.226_285_193_102_822,
        1.646_738_380_722_676_5,
        Epoch::from_mjd_tai(21_546.0),
        eme2k,
    );

    let expected_fx = Vector6::new(
        -3.288_789_003_770_57,
        -2.226_285_193_102_822,
        1.646_738_380_722_676_5,
        0.000_348_875_166_711_715_13,
        -0.000_715_134_890_110_951_6,
        -0.000_270_059_537_180_366_4,
    );

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());

    let init_sc = Spacecraft::from(init).with_stm();
    let fx_real = dynamics
        .eom(0.0, &init_sc.to_vector(), &init_sc, almanac.clone())
        .unwrap();

    let fx_orbit_real = fx_real.fixed_rows::<6>(0).to_owned();

    let (fx, grad) = dynamics
        .dual_eom(0.0, &Spacecraft::from(init).with_stm(), almanac.clone())
        .unwrap();

    let fx = fx.fixed_rows::<6>(0).to_owned();

    println!("{fx_orbit_real}\n{fx}\n{expected_fx}");

    assert!(
        (fx_orbit_real - fx).norm() < 1e-16,
        "Dual and real computations differ"
    );

    assert!(
        (fx - expected_fx).norm() < 1e-16,
        "f(x) computation is incorrect {:e}",
        (fx - expected_fx).norm()
    );

    let mut expected = OMatrix::<f64, Const<9>, Const<9>>::zeros();

    expected[(0, 3)] = 1.0;
    expected[(1, 4)] = 1.0;
    expected[(2, 5)] = 1.0;
    expected[(3, 0)] = -0.000_000_018_628_398_391_083_86;
    expected[(4, 0)] = -0.000_000_040_897_747_124_379_53;
    expected[(5, 0)] = -0.000_000_015_444_396_313_003_294;
    expected[(3, 1)] = -0.000_000_040_897_747_124_379_53;
    expected[(4, 1)] = 0.000_000_045_253_271_058_430_05;
    expected[(5, 1)] = 0.000_000_031_658_391_636_846_51;
    expected[(3, 2)] = -0.000_000_015_444_396_313_003_294;
    expected[(4, 2)] = 0.000_000_031_658_391_636_846_51;
    expected[(5, 2)] = -0.000_000_026_624_872_667_346_21;

    assert!(
        (grad - expected).norm() < 1e-16,
        "gradient computation is incorrect {:e}",
        (grad - expected).norm()
    );

    // [Earth J2000] 1917-11-14T00:00:00 UTC   sma = 22000.000344 km   ecc = 0.010000  inc = 30.000000 deg     raan = 80.000000 deg    aop = 40.000000 deg     ta = 0.000000 deg
    // Quite non-linear near periapsis
    let prop_time = 2 * Unit::Minute;
    let step_size = 10 * Unit::Second;

    let setup = Propagator::rk89(dynamics, IntegratorOptions::with_fixed_step(step_size));
    let mut prop = setup.with(init_sc, almanac.clone());
    let final_state = prop.for_duration(prop_time).unwrap();

    // Check that the STM is correct by back propagating by the previous step, and multiplying by the STM.
    let stm_k_to_0 = final_state.stm.unwrap();

    let prev_state = setup
        .with(Spacecraft::from(init).with_stm(), almanac)
        .for_duration(prop_time - step_size)
        .unwrap();
    let stm_km1_to_0 = prev_state.stm.unwrap();

    let stm_k_to_km1 = stm_k_to_0 * stm_km1_to_0.try_inverse().unwrap();

    // And check the difference
    let stm_err = stm_k_to_km1 * prev_state.to_vector().fixed_rows::<9>(0)
        - final_state.to_vector().fixed_rows::<9>(0);
    let radius_err = stm_err.fixed_rows::<3>(0).into_owned();
    let velocity_err = stm_err.fixed_rows::<3>(3).into_owned();

    assert!(dbg!(radius_err.norm()) < 1e-1);
    assert!(dbg!(velocity_err.norm()) < 1e-1);
}

#[allow(clippy::identity_op)]
#[rstest]
fn multi_body_dynamics_dual(almanac: Arc<Almanac>) {
    // After trial and error, it seems that the linearization breaks after 45 minutes for this example.
    // Specifically, inverting the previous STM and multiplying it with the next STM to compute the STM
    // of this one step will round too much and cause an error greater than one meter.
    let prop_time = 45 * Unit::Minute;
    let step_size = 10 * Unit::Second;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let halo_rcvr = Orbit::cartesian(
        333_321.004_516,
        -76_134.198_887,
        -20_873.831_939,
        0.257_153_712,
        0.930_284_066,
        0.346_177,
        start_time,
        eme2k,
    );

    // let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    // let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());

    let setup = Propagator::rk89(dynamics, IntegratorOptions::with_fixed_step(step_size));
    let final_state = setup
        .with(halo_rcvr.into(), almanac.clone())
        .for_duration(prop_time)
        .unwrap();
    let mut prop = setup.with(Spacecraft::from(halo_rcvr).with_stm(), almanac.clone());
    let final_state_dual = prop.for_duration(prop_time).unwrap();
    println!("Final STM {}", final_state_dual.stm().unwrap());

    // Test that reset_stm() and a single step will lead to the correct STM diagonals
    prop.state.reset_stm();
    let post_reset = prop.for_duration(step_size).unwrap();
    println!("{}", post_reset.stm().unwrap());

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
    // This should be zero!
    assert!(
        err_r < f64::EPSILON,
        "position error too large for multibody gravity"
    );
    assert!(
        err_v < f64::EPSILON,
        "velocity error too large for multibody gravity"
    );

    // Check that the STM is correct by back propagating by the previous step, and multiplying by the STM.
    let stm_k_to_0 = final_state_dual.stm.unwrap();

    let prev_state = setup
        .with(Spacecraft::from(halo_rcvr).with_stm(), almanac)
        .for_duration(prop_time - step_size)
        .unwrap();
    let stm_km1_to_0 = prev_state.stm.unwrap();

    let stm_k_to_km1 = stm_k_to_0 * stm_km1_to_0.try_inverse().unwrap();

    // And check the difference
    let stm_err = stm_k_to_km1 * prev_state.to_vector().fixed_rows::<9>(0)
        - final_state.to_vector().fixed_rows::<9>(0);
    let radius_stm_delta = stm_err.fixed_rows::<3>(0).into_owned();
    let velocity_stm_delta = stm_err.fixed_rows::<3>(3).into_owned();

    assert!(dbg!(radius_stm_delta.norm()) < 1e-2);
    assert!(dbg!(velocity_stm_delta.norm()) < 1e-3);
}

#[allow(clippy::identity_op)]
#[rstest]
fn val_earth_sph_harmonics_j2(almanac: Arc<Almanac>) {
    // NOTE: GMAT and Monte are within about 0.1 meters of difference in position. Hence, we're checking we're in the same bracket.
    use nyx::dynamics::Harmonics;
    use nyx::io::gravity::*;

    let monte_earth_gm = 3.986_004_328_969_392e5;
    let monte_earth_j2 = -0.000_484_169_325_971;

    let eme2k = almanac
        .frame_info(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(monte_earth_gm);

    let iau_earth = almanac
        .frame_info(IAU_EARTH_FRAME)
        .unwrap()
        .with_mu_km3_s2(monte_earth_gm);

    let earth_sph_harm = HarmonicsMem::from_j2(monte_earth_j2);
    let harmonics = Harmonics::from_stor(iau_earth, earth_sph_harm);

    let dt = Epoch::from_mjd_tai(MJD_J2000);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );
    // GMAT validation case
    // -5751.473991327555          4721.214035135832           2045.947768608806           -0.7977402618596746          -3.656451983636495           6.139637921620765
    /*
    let rslt_gmat = Vector6::new(
        -5_751.473_991_327_555,
        4721.214035135832,
        2045.947768608806,
        -0.7977402618596746,
        -3.656451983636495,
        6.139637921620765,
    );*/

    // Monte validation case
    // Orbit (km, km/sec)
    // 'Earth' -> 'test' in 'EME2000' at '02-JAN-2000 12:00:00.0000 TAI'
    // Pos: -5.751472565170783e+03  4.721183256208691e+03  2.046020865167045e+03
    // Vel: -7.976895830677169e-01 -3.656498994998706e+00  6.139616747276084e+00
    let rslt_monte = Vector6::new(
        -5.751_472_565_170_783e3,
        4.721_183_256_208_691e3,
        2.046_020_865_167_045e3,
        -7.976_895_830_677_169e-1,
        -3.656_498_994_998_706,
        6.139_616_747_276_084,
    );

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::from_model(harmonics));

    let prop_state = Propagator::rk89(dynamics, IntegratorOptions::default())
        .with(state.into(), almanac)
        .for_duration(1 * Unit::Day)
        .unwrap();

    println!("{prop_state}");

    println!("==> val_earth_sph_harmonics_j2 absolute errors (MONTE)");
    let delta = prop_state.orbit.to_cartesian_pos_vel() - rslt_monte;
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    let (err_r, err_v) =
        rss_orbit_vec_errors(&prop_state.orbit.to_cartesian_pos_vel(), &rslt_monte);

    assert!(err_r < 1e-1, "J2 failed in position: {err_r:.5e}");
    assert!(err_v < 1e-4, "J2 failed in velocity: {err_v:.5e}");
}

#[allow(clippy::identity_op)]
#[rstest]
fn val_earth_sph_harmonics_12x12(almanac_gmat: Arc<Almanac>) {
    let almanac = almanac_gmat;
    extern crate pretty_env_logger;
    let _ = pretty_env_logger::try_init();
    use nyx::dynamics::sph_harmonics::Harmonics;
    use nyx::io::gravity::*;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();
    let iau_earth = almanac.frame_info(IAU_EARTH_FRAME).unwrap();

    let earth_sph_harm =
        HarmonicsMem::from_cof("data/01_planetary/JGM3.cof.gz", 12, 12, true).unwrap();
    let harmonics = Harmonics::from_stor(iau_earth, earth_sph_harm);

    let dt = Epoch::from_mjd_tai(MJD_J2000);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );

    // GMAT validation case
    let rslt_gmat = Vector6::new(
        -5_751.935_197_673_059,
        4_719.330_857_046_409,
        2_048.776_230_999_391,
        -0.795_315_465_634_082_6,
        -3.658_346_256_468_031,
        6.138_852_391_455_04,
    );

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::from_model(harmonics));

    let setup = Propagator::rk89(dynamics.clone(), IntegratorOptions::with_tolerance(1e-9));
    let prop_time = 1 * Unit::Day;
    let final_state = setup
        .with(state.into(), almanac.clone())
        .for_duration(prop_time)
        .unwrap();

    println!("{final_state}");

    println!("==> val_earth_sph_harmonics_12x12 absolute errors");
    let delta = final_state.orbit.to_cartesian_pos_vel() - rslt_gmat;
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    let (err_r, err_v) =
        rss_orbit_vec_errors(&final_state.orbit.to_cartesian_pos_vel(), &rslt_gmat);

    assert!(err_r < 1e-1, "12x12 failed in position: {err_r:.5e}");
    assert!(err_v < 1e-4, "12x12 failed in velocity: {err_v:.5e}");

    // We set up a new propagator with a fixed step. Without the fixed step, the error control
    // on the STM leads to a difference of 1.04 meters in this one day propagation.
    let setup = Propagator::rk89(dynamics, IntegratorOptions::with_fixed_step_s(30.0));
    let prop_time = 6 * Unit::Hour;
    let final_state = setup
        .with(state.into(), almanac.clone())
        .for_duration(prop_time)
        .unwrap();
    // Compare the case with the hyperdual EOMs (computation uses another part of the code)
    let mut prop = setup.with(Spacecraft::from(state).with_stm(), almanac);
    let final_state_dual = prop.for_duration(prop_time).unwrap();

    let (err_r, err_v) = rss_orbit_vec_errors(
        &final_state.orbit.to_cartesian_pos_vel(),
        &final_state_dual.orbit.to_cartesian_pos_vel(),
    );
    println!(
        "Error between reals and duals accumulated over {} : {:.6} m \t{:.6} m/s",
        prop_time,
        err_r * 1e3,
        err_v * 1e3
    );
    // This should be zero!
    assert!(
        err_r < f64::EPSILON,
        "position error too large for 12x12 gravity"
    );
    assert!(
        err_v < f64::EPSILON,
        "velocity error too large for 12x12 gravity"
    );
}

#[allow(clippy::identity_op)]
#[rstest]
fn val_earth_sph_harmonics_70x70(almanac_gmat: Arc<Almanac>) {
    let almanac = almanac_gmat;
    extern crate pretty_env_logger;
    let _ = pretty_env_logger::try_init();
    use nyx::dynamics::Harmonics;
    use nyx::io::gravity::*;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();
    let iau_earth = almanac.frame_info(IAU_EARTH_FRAME).unwrap();

    let earth_sph_harm =
        HarmonicsMem::from_cof("data/01_planetary/JGM3.cof.gz", 70, 70, true).unwrap();
    let harmonics = Harmonics::from_stor(iau_earth, earth_sph_harm);

    let dt = Epoch::from_mjd_tai(MJD_J2000);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );

    // GMAT validation case
    let rslt_gmat = Vector6::new(
        -5_751.924_618_076_704,
        4_719.386_612_440_923,
        2_048.696_011_823_441,
        -0.795_383_404_365_819_8,
        -3.658_301_183_319_466,
        6.138_865_498_487_843,
    );

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::from_model(harmonics));

    let prop_rslt = Propagator::default(dynamics)
        .with(state.into(), almanac)
        .for_duration(1 * Unit::Day)
        .unwrap();

    println!("{prop_rslt}");

    println!("==> val_earth_sph_harmonics_70x70 absolute errors");
    let delta = prop_rslt.orbit.to_cartesian_pos_vel() - rslt_gmat;
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    let (err_r, err_v) = rss_orbit_vec_errors(&prop_rslt.orbit.to_cartesian_pos_vel(), &rslt_gmat);

    assert!(dbg!(err_r) < 0.2, "70x70 failed in position: {err_r:.5e}");
    assert!(dbg!(err_v) < 1e-3, "70x70 failed in velocity: {err_v:.5e}");
}

#[allow(clippy::identity_op)]
#[rstest]
fn val_earth_sph_harmonics_70x70_partials(almanac_gmat: Arc<Almanac>) {
    let almanac = almanac_gmat;
    extern crate pretty_env_logger;
    let _ = pretty_env_logger::try_init();
    use nyx::dynamics::Harmonics;
    use nyx::io::gravity::*;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();
    let iau_earth = almanac.frame_info(IAU_EARTH_FRAME).unwrap();

    let earth_sph_harm =
        HarmonicsMem::from_cof("data/01_planetary/JGM3.cof.gz", 70, 70, true).unwrap();
    let harmonics = Harmonics::from_stor(iau_earth, earth_sph_harm);

    let dt = Epoch::from_mjd_tai(MJD_J2000);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );

    // GMAT validation case
    let rslt_gmat = Vector6::new(
        -5_751.924_618_076_704,
        4_719.386_612_440_923,
        2_048.696_011_823_441,
        -0.795_383_404_365_819_8,
        -3.658_301_183_319_466,
        6.138_865_498_487_843,
    );

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::from_model(harmonics));

    let prop_rslt = Propagator::default(dynamics)
        .with(state.into(), almanac)
        .for_duration(1 * Unit::Day)
        .unwrap();

    println!("{prop_rslt}");

    println!("==> val_earth_sph_harmonics_70x70_partials absolute errors");
    let delta = prop_rslt.orbit.to_cartesian_pos_vel() - rslt_gmat;
    for i in 0..6 {
        print!("{:.0e}\t", delta[i].abs());
    }
    println!();

    let (err_r, err_v) = rss_orbit_vec_errors(&prop_rslt.orbit.to_cartesian_pos_vel(), &rslt_gmat);

    assert!(dbg!(err_r) < 0.2, "12x12 failed in position: {err_r:.5e}");
    assert!(dbg!(err_v) < 1e-3, "12x12 failed in velocity: {err_v:.5e}");
}

#[rstest]
fn hf_prop(almanac: Arc<Almanac>) {
    // Tests a high fidelity propagation over several days for performance analysis.

    extern crate pretty_env_logger;
    let _ = pretty_env_logger::try_init();
    use nyx::dynamics::sph_harmonics::Harmonics;
    use nyx::io::gravity::*;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();
    let iau_earth = almanac.frame_info(IAU_EARTH_FRAME).unwrap();

    let earth_sph_harm =
        HarmonicsMem::from_cof("data/01_planetary/JGM3.cof.gz", 21, 21, true).unwrap();
    let harmonics = Harmonics::from_stor(iau_earth, earth_sph_harm);

    let dt = Epoch::from_mjd_tai(MJD_J2000);
    let state = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );

    let bodies = vec![MOON, SUN, JUPITER_BARYCENTER];
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::new(vec![
        PointMasses::new(bodies),
        harmonics,
    ]));

    let setup = Propagator::rk89(dynamics, IntegratorOptions::with_tolerance(1e-9));
    let rslt = setup
        .with(state.into(), almanac)
        .for_duration(30.0 * Unit::Day)
        .unwrap();

    println!("{rslt}\n{rslt:x}");
}

#[rstest]
fn val_cislunar_dynamics(almanac_gmat: Arc<Almanac>) {
    let almanac = almanac_gmat;
    let prop_time = 36 * Unit::Hour;

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();

    // "2022 NOV 27 05:55:49"
    let dt = Epoch::from_gregorian_utc_hms(2022, 11, 27, 5, 55, 49);
    let state = Orbit::cartesian(
        -7.529_485_277_404_609e2,
        5.624_035_455_855_085e3,
        3.278_632_833_875_311e3,
        -7.683_161_946_015_461,
        -0.860_670_301_418_699_3,
        -0.085_614_035_370_280_35,
        dt,
        eme2k,
    );

    let rslt = Orbit::cartesian(
        -6151.500843512164,
        1833.914559118256,
        1384.255278442759,
        -2.604989451776925,
        -6.353097736110432,
        -3.519264757529829,
        dt + prop_time,
        eme2k,
    );

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(vec![EARTH, SUN, MOON]));
    let setup = Propagator::new(
        dynamics,
        IntegratorMethod::RungeKutta4,
        IntegratorOptions::with_fixed_step_s(0.5),
    );
    let mut prop = setup.with(state.into(), almanac);
    prop.for_duration(prop_time).unwrap();

    println!("==> val_cislunar_dynamics absolute errors");
    let delta = prop.state.orbit.to_cartesian_pos_vel() - rslt.to_cartesian_pos_vel();
    for i in 0..3 {
        print!("{:.0e} m\t", delta[i].abs() * 1e3);
    }
    for i in 3..6 {
        print!("{:.0e} m/s\t", delta[i].abs() * 1e3);
    }
    println!();

    let (err_r, err_v) = rss_orbit_errors(&prop.state.orbit, &rslt);

    println!(
        "RSS errors:\tpos = {:.5e} m\tvel = {:.5e} m/s\ninit\t{}\nfinal\t{}",
        err_r * 1e3,
        err_v * 1e3,
        state,
        prop.state
    );

    assert!(
        err_r < 3e-6,
        "val_cislunar_dynamics failed in position: {err_r:.5e}"
    );
    assert!(
        err_v < 3e-9,
        "val_cislunar_dynamics failed in velocity: {err_v:.5e}"
    );
}
