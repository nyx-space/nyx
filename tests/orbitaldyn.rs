extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::na::{Vector6, U3};
fn rss_state_errors(prop_err: &Vector6<f64>, cur_state: &Vector6<f64>) -> (f64, f64) {
    let err_radius = (prop_err.fixed_rows::<U3>(0) - cur_state.fixed_rows::<U3>(0)).norm();

    let err_velocity = (prop_err.fixed_rows::<U3>(3) - cur_state.fixed_rows::<U3>(3)).norm();

    (err_radius, err_velocity)
}

#[test]
fn two_body_custom() {
    extern crate nalgebra as na;
    use self::na::Vector6;
    use hifitime::julian::ModifiedJulian;
    use nyx::celestia::{bodies, Cosm, Geoid, State};
    use nyx::dynamics::celestial::CelestialDynamics;
    use nyx::propagators::error_ctrl::RSSStepPV;
    use nyx::propagators::*;

    let prop_time = 24.0 * 3_600.0;

    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(bodies::EARTH_BARYCENTER).unwrap();

    let dt = ModifiedJulian::j2000();
    let mut state = State::<Geoid>::from_cartesian(-2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, earth_geoid);
    state.frame.gm = 398_600.441_5;

    let rslt = Vector6::new(
        -5_971.194_191_684_024,
        3_945.506_653_624_737_3,
        2_864.636_617_867_270_6,
        0.049_096_957_141_044_464,
        -4.185_093_318_149_689,
        5.848_940_867_979_176,
    );

    let mut dynamics = CelestialDynamics::two_body(state);

    let mut prop = Propagator::new::<RK89>(&mut dynamics, &PropOpts::<RSSStepPV>::default());
    prop.until_time_elapsed(prop_time);
    assert_eq!(prop.state(), rslt, "two body prop failed");
    // And now do the backprop
    prop.until_time_elapsed(-prop_time);
    let (err_r, err_v) = rss_state_errors(&prop.state(), &state.to_cartesian_vec());
    assert!(
        err_r < 1e-5,
        "two body back prop failed to return to the initial state in position"
    );
    assert!(
        err_v < 1e-8,
        "two body back prop failed to return to the initial state in velocity"
    );
}

#[test]
fn two_body_dynamics() {
    extern crate nalgebra as na;
    use self::na::Vector6;
    use hifitime::julian::ModifiedJulian;
    use nyx::celestia::{bodies, Cosm, Geoid, State};
    use nyx::dynamics::celestial::CelestialDynamics;
    use nyx::propagators::error_ctrl::RSSStepPV;
    use nyx::propagators::*;
    use std::f64::EPSILON;

    let prop_time = 24.0 * 3_600.0;

    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(bodies::EARTH_BARYCENTER).unwrap();

    let dt = ModifiedJulian::j2000();
    let state = State::<Geoid>::from_cartesian(-2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, earth_geoid);

    let rslt = Vector6::new(
        -5_971.194_376_797_642,
        3_945.517_912_574_167_4,
        2_864.620_957_744_445,
        0.049_083_101_605_521_72,
        -4.185_084_125_817_668,
        5.848_947_462_472_871,
    );

    let mut dynamics = CelestialDynamics::two_body(state);

    let mut prop = Propagator::new::<RK89>(&mut dynamics, &PropOpts::<RSSStepPV>::default());
    prop.until_time_elapsed(prop_time);
    assert_eq!(prop.state(), rslt, "two body prop failed");
    assert!((prop.dynamics.state.dt_as_modified_julian().days - dt.days - 1.0).abs() <= EPSILON);
    // And now do the backprop
    prop.until_time_elapsed(-prop_time);
    let (err_r, err_v) = rss_state_errors(&prop.state(), &state.to_cartesian_vec());
    assert!(
        err_r < 1e-6,
        "two body back prop failed to return to the initial state in position"
    );
    assert!(
        err_v < 1e-9,
        "two body back prop failed to return to the initial state in velocity"
    );
}

#[test]
fn three_body_dynamics() {
    extern crate nalgebra as na;
    use self::na::Vector6;
    use hifitime::datetime::Datetime;
    use hifitime::julian::ModifiedJulian;
    use hifitime::TimeSystem;
    use nyx::celestia::{bodies, Cosm, Geoid, State};
    use nyx::dynamics::celestial::CelestialDynamics;
    use nyx::propagators::error_ctrl::RSSStepPV;
    use nyx::propagators::*;

    let prop_time = 24.0 * 3_600.0;

    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(bodies::EARTH).unwrap();

    let start_time = ModifiedJulian::from_instant(Datetime::at_midnight(2020, 1, 1).unwrap().into_instant());

    let halo_rcvr = State::<Geoid>::from_cartesian(
        333_321.004_516,
        -76_134.198_887,
        -20_873.831_939,
        0.257_153_712,
        0.930_284_066,
        0.346_177,
        start_time,
        earth_geoid,
    );

    // Monte data:
    let rslt = Vector6::new(
        345_350.661_525_660,
        5_930.672_633_019_7,
        7_333.285_591_307,
        2.129_812_933e-2,
        9.566_789_680e-1,
        3.028_176_198e-1,
    );
    // GMAT data (uses different GMs)
    // 28850.5                   345350.6640304797           5930.672047088849           7333.28377928682            0.02129818943860685          0.9566789568516441           0.302817581134027

    // HAD:   left: `Matrix { data: [343206.9934741252, 6027.320454567925,  9430.188737830082, -0.02934707799392118,  0.9600169353090493, 0.35103917505038895] }`,
    // GOT:   left: `Matrix { data: [343162.7343812757, 5989.6805335097215, 9412.246006634898, -0.030412263432568358, 0.959229681167265,  0.3506632777757227] }`,
    // NOW:   left: `Matrix { data: [343007.25430964294, 6011.974861696354, 9417.141146684777, -0.034068311795565645, 0.9595987074407781, 0.35072154117333154] }`,\
    // WANT: right: `Matrix { data: [345350.66152566,   5930.6726330197,    7333.285591307,     0.02129812933,        0.956678968,        0.3028176198] }`: two body prop failed', tests/orbitaldyn.rs:158:5

    // Without third bodies:
    //  ┌                      ┐
    //  │    343051.5163061142 │
    //  │    6049.607891125282 │
    //  │    9435.080870474838 │
    //  │ -0.03300296559877033 │
    //  │   0.9603856580047353 │
    //  │  0.35109730571850106 │
    //  └                      ┘

    let bodies = vec![bodies::EARTH_MOON, bodies::SUN, bodies::JUPITER_BARYCENTER];
    let mut dynamics = CelestialDynamics::new(halo_rcvr, bodies, &cosm);
    // let mut dynamics = CelestialDynamics::new(halo_rcvr, vec![bodies::EARTH_MOON, bodies::JUPITER_BARYCENTER], &cosm);

    let mut prop = Propagator::new::<RK89>(&mut dynamics, &PropOpts::with_adaptive_step(1.0, 2700.0, 1e-13, RSSStepPV {}));
    prop.until_time_elapsed(prop_time);
    println!("{}", prop.state());
    let (err_r, err_v) = rss_state_errors(&prop.state(), &rslt);
    assert!(err_r < 1e-3, format!("multi body failed in position: {:.5e}", err_r));
    assert!(err_v < 1e-6, format!("multi body failed in velocity: {:.5e}", err_v));
}

#[test]
fn two_body_dual() {
    // This is a duplicate of the differentials test in hyperdual.
    extern crate nalgebra as na;
    use self::na::{Matrix6, Vector6};
    use hifitime::julian::ModifiedJulian;
    use nyx::celestia::{Cosm, Geoid, State};
    use nyx::dynamics::celestial::TwoBodyWithDualStm;
    use nyx::od::AutoDiffDynamics;

    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(3).unwrap();

    let init = State::<Geoid>::from_cartesian(
        -9_042.862_233_600_335,
        18_536.333_069_123_244,
        6_999.957_069_486_411_5,
        -3.288_789_003_770_57,
        -2.226_285_193_102_822,
        1.646_738_380_722_676_5,
        ModifiedJulian { days: 21546.0 },
        earth_geoid,
    );

    let expected_fx = Vector6::new(
        -3.288_789_003_770_57,
        -2.226_285_193_102_822,
        1.646_738_380_722_676_5,
        0.000_348_875_166_673_120_14,
        -0.000_715_134_890_031_838_4,
        -0.000_270_059_537_150_490_5,
    );

    let dynamics = TwoBodyWithDualStm::from_state(&init);
    let (fx, grad) = dynamics.compute(0.0, &init.to_cartesian_vec());

    assert!(
        (fx - expected_fx).norm() < 1e-16,
        "f(x) computation is incorrect {:e}",
        (fx - expected_fx).norm()
    );

    let mut expected = Matrix6::zeros();

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

    assert_eq!(dynamics.to_state(), init);
}
