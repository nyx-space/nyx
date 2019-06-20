extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::na::{Vector6, U3};
fn backprop_rss_state_errors(prop_err: &Vector6<f64>, cur_state: &Vector6<f64>) -> (f64, f64) {
    let err_radius = (&prop_err.fixed_rows::<U3>(0).into_owned() - &cur_state.fixed_rows::<U3>(0).into_owned()).norm();

    let err_velocity = (&prop_err.fixed_rows::<U3>(3).into_owned() - &cur_state.fixed_rows::<U3>(3).into_owned()).norm();

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

    let init = Vector6::new(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0);

    let rslt = Vector6::from_row_slice(&[
        -5971.1941916712285,
        3945.5066532419537,
        2864.636618390466,
        0.04909695760948815,
        -4.1850933184621315,
        5.848940867758592,
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

    let init = Vector6::new(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0);

    let rslt = Vector6::new(
        -5971.194191684024,
        3945.5066536247373,
        2864.6366178672706,
        0.049096957141044464,
        -4.185093318149689,
        5.848940867979176,
    );

    let mut dynamics = TwoBody::from_state_vec_with_gm(init, 398600.4415);
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
    let initial_state = State::<Geoid>::from_cartesian(
        -2436.45,
        -2436.45,
        6891.037,
        5.088611,
        -5.088611,
        0.0,
        dt,
        earth_geoid.clone(),
    );

    println!("Initial state:\n{0}\n{0:o}\n", initial_state);

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 60.0;

    let rslt = State::<Geoid>::from_cartesian(
        -5971.1941916712285,
        3945.5066532419537,
        2864.636618390466,
        0.04909695760948815,
        -4.1850933184621315,
        5.848940867758592,
        ModifiedJulian { days: 21546.0 },
        earth_geoid.clone(),
    );

    let mut dynamics = TwoBody::from_state_vec(initial_state.to_cartesian_vec(), earth_geoid.clone());
    let mut prop = Propagator::new::<RK89>(
        &mut dynamics,
        &PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSStepPV {}),
    );
    let (final_t, final_state0) = prop.until_time_elapsed(prop_time);

    let final_dt = ModifiedJulian {
        days: dt.days + final_t / SECONDS_PER_DAY,
    };
    let final_state = State::from_cartesian_vec(&prop.state(), final_dt, earth_geoid.clone());
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
        -9042.862233600335,
        18536.333069123244,
        6999.9570694864115,
        -3.28878900377057,
        -2.226285193102822,
        1.6467383807226765,
        ModifiedJulian { days: 21546.0 },
        earth_geoid,
    );

    let expected_fx = Vector6::new(
        -3.28878900377057,
        -2.226285193102822,
        1.6467383807226765,
        0.0003488751720191492,
        -0.0007151349009902908,
        -0.00027005954128877916,
    );

    let dynamics = TwoBodyWithDualStm::from_state(&init);
    let (fx, grad) = dynamics.compute(0.0, &init.to_cartesian_vec());

    assert!(
        (fx - expected_fx).norm() < 1e-16,
        "f(x) computation is incorrect ({})",
        (fx - expected_fx).norm()
    );

    let mut expected = Matrix6::zeros();

    expected[(0, 3)] = 1.0;
    expected[(1, 4)] = 1.0;
    expected[(2, 5)] = 1.0;
    expected[(3, 0)] = -0.000000018628398676538285;
    expected[(4, 0)] = -0.00000004089774775108092;
    expected[(5, 0)] = -0.0000000154443965496673;
    expected[(3, 1)] = -0.00000004089774775108092;
    expected[(4, 1)] = 0.000000045253271751873843;
    expected[(5, 1)] = 0.00000003165839212196757;
    expected[(3, 2)] = -0.0000000154443965496673;
    expected[(4, 2)] = 0.00000003165839212196757;
    expected[(5, 2)] = -0.000000026624873075335538;

    assert!(
        (grad - expected).norm() < 1e-16,
        "gradient computation is incorrect {}",
        (grad - expected).norm()
    );
}
