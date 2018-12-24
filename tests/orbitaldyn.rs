extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

#[test]
fn two_body_parametrized() {
    extern crate nalgebra as na;
    use self::na::Vector6;
    use nyx::celestia::EARTH;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::dynamics::Dynamics;
    use nyx::propagators::*;

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

    let mut prop = Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy));
    let mut dyn = TwoBody::from_state_vec::<EARTH>(init);
    prop.until_time_elapsed(prop_time, &mut dyn, error_ctrl::rss_step_pos_vel);
    assert_eq!(dyn.state(), rslt, "two body prop failed");
    // And now do the backprop
    prop.until_time_elapsed(-prop_time, &mut dyn, error_ctrl::rss_step_pos_vel);
    let delta = (dyn.state() - init).norm();
    assert!(delta < 1e-5, "two body back prop failed to return to the initial state");
}

#[test]
fn two_body_custom() {
    extern crate nalgebra as na;
    use self::na::Vector6;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::dynamics::Dynamics;
    use nyx::propagators::*;

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 60.0;

    let rslt = Vector6::new(
        -5971.1941916712285,
        3945.5066532419537,
        2864.636618390466,
        0.04909695760948815,
        -4.1850933184621315,
        5.848940867758592,
    );

    let mut prop = Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy));
    let mut dyn = TwoBody::from_state_vec_with_gm(
        Vector6::new(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0),
        398600.4415,
    );
    prop.until_time_elapsed(prop_time, &mut dyn, error_ctrl::rss_step_pos_vel);
    assert_eq!(dyn.state(), rslt, "two body prop failed");
}

#[test]
fn two_body_state_parametrized() {
    extern crate nalgebra as na;
    use hifitime::julian::ModifiedJulian;
    use hifitime::SECONDS_PER_DAY;
    use nyx::celestia::{State, EARTH, ECI};
    use nyx::dynamics::celestial::TwoBody;
    use nyx::dynamics::Dynamics;
    use nyx::propagators::{error_ctrl, Options, Propagator, RK89};

    let dt = ModifiedJulian { days: 21545.0 };
    let initial_state = State::from_cartesian_eci(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0, dt);

    println!("Initial state:\n{0}\n{0:o}\n", initial_state);

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 60.0;

    let rslt = State::from_cartesian_eci(
        -5971.1941916712285,
        3945.5066532419537,
        2864.636618390466,
        0.04909695760948815,
        -4.1850933184621315,
        5.848940867758592,
        ModifiedJulian { days: 21546.0 },
    );

    let mut prop = Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy));
    let mut dyn = TwoBody::from_state_vec::<EARTH>(initial_state.to_cartesian_vec());
    let (final_t, _) = prop.until_time_elapsed(prop_time, &mut dyn, error_ctrl::rss_step_pos_vel);

    let final_dt = ModifiedJulian {
        days: dt.days + final_t / SECONDS_PER_DAY,
    };
    let final_state = State::from_cartesian_vec::<EARTH, ModifiedJulian>(&dyn.state(), final_dt, ECI {});
    assert_eq!(final_state, rslt, "two body prop failed",);

    println!("Final state:\n{0}\n{0:o}", final_state);
}

#[test]
fn two_body_dual() {
    // This is a duplicate of the differentials test in dual_num.
    extern crate nalgebra as na;
    use self::na::{Matrix6, Vector6};
    use hifitime::julian::ModifiedJulian;
    use nyx::celestia::{State, EARTH, ECI};
    use nyx::dynamics::celestial::TwoBodyWithDualStm;
    use nyx::od::AutoDiffDynamics;

    let init = State::from_cartesian_eci(
        -9042.862233600335,
        18536.333069123244,
        6999.9570694864115,
        -3.28878900377057,
        -2.226285193102822,
        1.6467383807226765,
        ModifiedJulian { days: 21546.0 },
    );

    let expected_fx = Vector6::new(
        -3.28878900377057,
        -2.226285193102822,
        1.6467383807226765,
        0.0003488751720191492,
        -0.0007151349009902908,
        -0.00027005954128877916,
    );

    let dyn = TwoBodyWithDualStm::from_state::<EARTH, ECI>(init);
    let (fx, grad) = dyn.compute(0.0, &init.to_cartesian_vec());

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
