extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::na::{Vector6, U3};
fn backprop_rss_state_errors(prop_err: &Vector6<f64>, cur_state: &Vector6<f64>) -> (f64, f64) {
    let err_radius = (prop_err.fixed_rows::<U3>(0) - cur_state.fixed_rows::<U3>(0)).norm();

    let err_velocity = (prop_err.fixed_rows::<U3>(3) - cur_state.fixed_rows::<U3>(3)).norm();

    (err_radius, err_velocity)
}

#[test]
fn two_body_parametrized() {
    extern crate nalgebra as na;
    use self::na::Vector6;
    use nyx::celestia::Cosm;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::propagators::error_ctrl::RSSStepPV;
    use nyx::propagators::*;

    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(3).unwrap();

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 60.0;

    let init = Vector6::new(-2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0);

    let rslt = Vector6::from_row_slice(&[
        -5_971.194_376_784_884,
        3_945.517_912_191_541,
        2_864.620_958_267_658_4,
        0.049_083_102_073_914_83,
        -4.185_084_126_130_087_5,
        5.848_947_462_252_259_5,
    ]);

    let mut dynamics = TwoBody::from_state_vec(init, earth_geoid);
    let mut prop = Propagator::new::<RK89>(
        &mut dynamics,
        &PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSStepPV {}),
    );
    prop.until_time_elapsed(prop_time);
    assert_eq!(prop.state(), rslt, "two body prop failed");
    // And now do the backprop
    prop.until_time_elapsed(-prop_time);
    let (err_r, err_v) = backprop_rss_state_errors(&prop.state(), &init);
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
fn two_body_custom() {
    extern crate nalgebra as na;
    use self::na::Vector6;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::propagators::error_ctrl::RSSStepPV;
    use nyx::propagators::*;

    let prop_time = 24.0 * 3_600.0;

    let init = Vector6::new(-2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0);

    let rslt = Vector6::new(
        -5_971.194_191_684_024,
        3_945.506_653_624_737_3,
        2_864.636_617_867_270_6,
        0.049_096_957_141_044_464,
        -4.185_093_318_149_689,
        5.848_940_867_979_176,
    );

    let mut dynamics = TwoBody::from_state_vec_with_gm(init, 398_600.441_5);
    let mut prop = Propagator::new::<RK89>(&mut dynamics, &PropOpts::<RSSStepPV>::default());
    prop.until_time_elapsed(prop_time);
    assert_eq!(prop.state(), rslt, "two body prop failed");
    // And now do the backprop
    prop.until_time_elapsed(-prop_time);
    let (err_r, err_v) = backprop_rss_state_errors(&prop.state(), &init);
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
fn two_body_state_parametrized() {
    extern crate nalgebra as na;
    use hifitime::julian::ModifiedJulian;
    use hifitime::SECONDS_PER_DAY;
    use nyx::celestia::{Cosm, Geoid, State};
    use nyx::dynamics::celestial::TwoBody;
    use nyx::propagators::error_ctrl::RSSStepPV;
    use nyx::propagators::{PropOpts, Propagator, RK89};

    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(3).unwrap();

    let dt = ModifiedJulian { days: 21545.0 };
    let initial_state =
        State::<Geoid>::from_cartesian(-2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, earth_geoid);

    println!("Initial state:\n{0}\n{0:o}\n", initial_state);

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 60.0;

    let rslt = State::<Geoid>::from_cartesian(
        -5_971.194_376_784_884,
        3_945.517_912_191_541,
        2_864.620_958_267_658_4,
        0.049_083_102_073_914_83,
        -4.185_084_126_130_087_5,
        5.848_947_462_252_259_5,
        ModifiedJulian { days: 21546.0 },
        earth_geoid,
    );

    let mut dynamics = TwoBody::from_state_vec(initial_state.to_cartesian_vec(), earth_geoid);
    let mut prop = Propagator::new::<RK89>(
        &mut dynamics,
        &PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSStepPV {}),
    );
    let (final_t, final_state0) = prop.until_time_elapsed(prop_time);

    let final_dt = ModifiedJulian {
        days: dt.days + final_t / SECONDS_PER_DAY,
    };
    let final_state = State::from_cartesian_vec(&prop.state(), final_dt, earth_geoid);
    assert_eq!(final_state, rslt, "two body prop failed");
    assert_eq!(prop.state(), final_state0, "until_time_elapsed returns the wrong value");

    println!("Final state:\n{0}\n{0:o}", final_state);

    // And now do the backprop
    prop.until_time_elapsed(-prop_time);
    let (err_r, err_v) = backprop_rss_state_errors(&prop.state(), &initial_state.to_cartesian_vec());
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
fn two_body_dual() {
    // This is a duplicate of the differentials test in dual_num.
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
}
