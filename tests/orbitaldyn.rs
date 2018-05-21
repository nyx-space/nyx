extern crate nalgebra as na;
extern crate nyx_space as nyx;
use std::f64;

#[test]
fn two_body_parametrized() {
    extern crate nalgebra as na;
    use nyx::propagators::*;
    use nyx::dynamics::Dynamics;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::celestia::EARTH;
    use self::na::Vector6;

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

    let mut dyn = TwoBody::from_state_vec::<EARTH>(&Vector6::new(
        -2436.45,
        -2436.45,
        6891.037,
        5.088611,
        -5.088611,
        0.0,
    ));
    loop {
        let (t, state) = prop.derive(
            dyn.time(),
            &dyn.state(),
            |t_: f64, state_: &Vector6<f64>| dyn.eom(t_, state_),
            error_ctrl::rss_step_pos_vel,
        );

        if t < prop_time {
            // We haven't passed the time based stopping condition.
            dyn.set_state(t, &state);
        } else {
            let prev_details = prop.latest_details().clone();
            let overshot = t - prop_time;
            if overshot > 0.0 {
                println!("overshot by {} seconds", overshot);
                prop.set_fixed_step(prev_details.step - overshot);
                // Take one final step
                let (t, state) = prop.derive(
                    dyn.time(),
                    &dyn.state(),
                    |t_: f64, state_: &Vector6<f64>| dyn.eom(t_, state_),
                    error_ctrl::rss_step_pos_vel,
                );
                dyn.set_state(t, &state);
            } else {
                dyn.set_state(t, &state);
            }

            if prev_details.error > accuracy {
                assert!(
                    prev_details.step - min_step < f64::EPSILON,
                    "step size should be at its minimum because error is higher than tolerance: {:?}",
                    prev_details
                );
            }

            assert_eq!(dyn.state(), rslt, "two body prop failed",);
            break;
        }
    }
}

#[test]
fn two_body_custom() {
    extern crate nalgebra as na;
    use nyx::propagators::*;
    use nyx::dynamics::Dynamics;
    use nyx::dynamics::celestial::TwoBody;
    use self::na::Vector6;

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
    loop {
        let (t, state) = prop.derive(
            dyn.time(),
            &dyn.state(),
            |t_: f64, state_: &Vector6<f64>| dyn.eom(t_, state_),
            error_ctrl::rss_step_pos_vel,
        );

        if t < prop_time {
            // We haven't passed the time based stopping condition.
            dyn.set_state(t, &state);
        } else {
            let prev_details = prop.latest_details().clone();
            let overshot = t - prop_time;
            if overshot > 0.0 {
                println!("overshot by {} seconds", overshot);
                prop.set_fixed_step(prev_details.step - overshot);
                // Take one final step
                let (t, state) = prop.derive(
                    dyn.time(),
                    &dyn.state(),
                    |t_: f64, state_: &Vector6<f64>| dyn.eom(t_, state_),
                    error_ctrl::rss_step_pos_vel,
                );
                dyn.set_state(t, &state);
            } else {
                dyn.set_state(t, &state);
            }

            if prev_details.error > accuracy {
                assert!(
                    prev_details.step - min_step < f64::EPSILON,
                    "step size should be at its minimum because error is higher than tolerance: {:?}",
                    prev_details
                );
            }

            assert_eq!(dyn.state(), rslt, "two body prop failed",);
            break;
        }
    }
}

#[test]
fn two_body_state_parametrized() {
    extern crate nalgebra as na;
    use nyx::propagators::{error_ctrl, Options, Propagator, RK89};
    use nyx::dynamics::Dynamics;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::celestia::{State, EARTH};
    use self::na::Vector6;

    let initial_state = State::from_cartesian::<EARTH>(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0);

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
    );

    let mut prop = Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy));

    let mut dyn = TwoBody::from_state_vec::<EARTH>(&initial_state.to_cartesian_vec());
    let final_state: State;

    loop {
        let (t, state) = prop.derive(
            dyn.time(),
            &dyn.state(),
            |t_: f64, state_: &Vector6<f64>| dyn.eom(t_, state_),
            error_ctrl::rss_step_pos_vel,
        );

        if t < prop_time {
            // We haven't passed the time based stopping condition.
            dyn.set_state(t, &state);
        } else {
            let prev_details = prop.latest_details().clone();
            let overshot = t - prop_time;
            if overshot > 0.0 {
                println!("overshot by {} seconds", overshot);
                prop.set_fixed_step(prev_details.step - overshot);
                // Take one final step
                let (t, state) = prop.derive(
                    dyn.time(),
                    &dyn.state(),
                    |t_: f64, state_: &Vector6<f64>| dyn.eom(t_, state_),
                    error_ctrl::rss_step_pos_vel,
                );
                dyn.set_state(t, &state);
            } else {
                dyn.set_state(t, &state);
            }

            if prev_details.error > accuracy {
                assert!(
                    prev_details.step - min_step < f64::EPSILON,
                    "step size should be at its minimum because error is higher than tolerance: {:?}",
                    prev_details
                );
            }

            final_state = State::from_cartesian_vec::<EARTH>(&dyn.state());
            assert_eq!(final_state, rslt, "two body prop failed",);
            break;
        }
    }
    println!("Final state:\n{0}\n{0:o}", final_state);
}
