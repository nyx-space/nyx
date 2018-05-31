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

    let rslt = Vector6::from_row_slice(&[
        -5971.1941916712285,
        3945.5066532419537,
        2864.636618390466,
        0.04909695760948815,
        -4.1850933184621315,
        5.848940867758592,
    ]);

    let mut prop = Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy));
    let mut dyn = TwoBody::from_state_vec::<EARTH>(&Vector6::new(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0));
    let (final_t, final_state) = prop.prop_for(prop_time, dyn, error_ctrl::rss_step_pos_vel);
    dyn.set_state(final_t, &final_state);
    assert_eq!(dyn.state(), rslt, "two body prop failed");
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

    let rslt = Vector6::from_row_slice(&[
        -5971.1941916712285,
        3945.5066532419537,
        2864.636618390466,
        0.04909695760948815,
        -4.1850933184621315,
        5.848940867758592,
    ]);

    let mut prop = Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy));
    let mut dyn = TwoBody::from_state_vec_with_gm(
        &Vector6::new(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0),
        398600.4415,
    );
    let (final_t, final_state) = prop.prop_for(prop_time, dyn, error_ctrl::rss_step_pos_vel);
    dyn.set_state(final_t, &final_state);
    assert_eq!(dyn.state(), rslt, "two body prop failed");
}

#[test]
fn two_body_state_parametrized() {
    extern crate nalgebra as na;
    use hifitime::julian::ModifiedJulian;
    use hifitime::SECONDS_PER_DAY;
    use nyx::celestia::{State, EARTH};
    use nyx::dynamics::celestial::TwoBody;
    use nyx::dynamics::Dynamics;
    use nyx::propagators::{error_ctrl, Options, Propagator, RK89};

    let dt = ModifiedJulian { days: 21545.0 };
    let initial_state = State::from_cartesian::<EARTH>(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0, dt);

    println!("Initial state:\n{0}\n{0:o}\n", initial_state);

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 60.0;

    let rslt = State::from_cartesian::<EARTH>(
        -5971.1941916712285,
        3945.5066532419537,
        2864.636618390466,
        0.04909695760948815,
        -4.1850933184621315,
        5.848940867758592,
        ModifiedJulian { days: 21546.0 },
    );

    let mut prop = Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy));
    let mut dyn = TwoBody::from_state_vec::<EARTH>(&initial_state.to_cartesian_vec());
    let (final_t, final_state_vec) = prop.prop_for(prop_time, dyn, error_ctrl::rss_step_pos_vel);
    dyn.set_state(final_t, &final_state_vec);

    let final_dt = ModifiedJulian {
        days: dt.days + final_t / SECONDS_PER_DAY,
    };
    let final_state = State::from_cartesian_vec::<EARTH>(&dyn.state(), final_dt);
    assert_eq!(final_state, rslt, "two body prop failed",);

    println!("Final state:\n{0}\n{0:o}", final_state);
}
