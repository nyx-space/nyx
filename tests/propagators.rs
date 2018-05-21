extern crate nalgebra as na;
extern crate nyx_space as nyx;
use std::f64;
use self::na::{U1, U3, Vector6};

fn two_body_dynamics(_t: f64, state: &Vector6<f64>) -> Vector6<f64> {
    let radius = state.fixed_slice::<U3, U1>(0, 0);
    let velocity = state.fixed_slice::<U3, U1>(3, 0);
    let body_acceleration = (-398_600.4415 / radius.norm().powi(3)) * radius;
    Vector6::from_iterator(velocity.iter().chain(body_acceleration.iter()).cloned())
}

#[test]
fn regress_leo_day_adaptive() {
    // Regression test for propagators not available in GMAT.
    extern crate nalgebra as na;
    use nyx::propagators::*;
    use self::na::Vector6;

    let prop_time = 24.0 * 3600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 30.0;

    // NOTE: In this test we only use the propagators which also exist in GMAT.
    // Refer to `regress_leo_day_adaptive` for the additional propagators.
    let mut all_props = vec![
        Propagator::new::<CashKarp45>(&Options::with_adaptive_step(min_step, max_step, accuracy)),
        Propagator::new::<Fehlberg45>(&Options::with_adaptive_step(min_step, max_step, accuracy)),
    ];

    let all_it_cnt = vec![5178, 6817]; // NOTE: This is a decent estimate of which propagators to use depending if speed is important.

    let all_rslts = vec![
        Vector6::from_row_slice(&[
            -5971.194190251459,
            3945.5066108335227,
            2864.636676354708,
            0.049097009502150256,
            -4.185093353076805,
            5.848940843324859,
        ]),
        Vector6::from_row_slice(&[
            -5971.1941909881425,
            3945.506632794103,
            2864.6366463408517,
            0.04909698263114359,
            -4.185093335151727,
            5.848940855975827,
        ]),
    ];

    for (p_id, prop) in all_props.iter_mut().enumerate() {
        let mut init_state = Vector6::from_row_slice(&[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);
        let mut cur_t = 0.0;
        let mut iterations = 0;
        loop {
            let (t, state) = prop.derive(
                cur_t,
                &init_state,
                two_body_dynamics,
                error_ctrl::rss_state_pos_vel,
            );
            iterations += 1;
            if t < prop_time {
                // We haven't passed the time based stopping condition.
                cur_t = t;
                init_state = state;
            } else {
                // NOTE: The refined propagation time for here is different from what GMAT does.
                // In fact, GMAT uses the [secant method](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/command/Propagate.cpp#L5450)
                // even for time computation. This means there are several iterations over the derivatives, and a non-exact iteration time.
                // GMAT still uses a fixed step propagation though.
                // At this point, we've passed the condition, so let's switch to a fixed step of _exactly_ the
                // previous time step minus the amount by which we overshot. This allows us to propagate in time for
                // _exactly_ the time we want to propagate for.
                let prev_details = prop.latest_details().clone();
                let overshot = t - prop_time;
                prop.set_fixed_step(prev_details.step - overshot);
                // Take one final step
                let (t, state) = prop.derive(
                    cur_t,
                    &init_state,
                    two_body_dynamics,
                    error_ctrl::rss_state_pos_vel,
                );

                assert!(
                    (t - prop_time).abs() < 1e-12,
                    "propagated for {} instead of {}",
                    t,
                    prop_time
                );

                // Let's check that, prior to the refined step, we either hit the accuracy wanted,
                // or we are using the minimum step size.
                if prev_details.error > accuracy {
                    assert!(
                        prev_details.step - min_step < f64::EPSILON,
                        "step size should be at its minimum because error is higher than tolerance (p_id = {}): {:?}",
                        p_id,
                        prev_details
                    );
                }

                assert_eq!(
                    state,
                    all_rslts[p_id],
                    "leo prop failed for p_id = {}",
                    p_id
                );

                assert_eq!(
                    iterations,
                    all_it_cnt[p_id],
                    "wrong number of iterations (p_id = {})",
                    p_id
                );
                break;
            }
        }
    }
}

#[test]
fn gmat_val_leo_day_adaptive() {
    extern crate nalgebra as na;
    use nyx::propagators::*;
    use self::na::Vector6;

    let prop_time = 24.0 * 3600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 30.0;

    // NOTE: In this test we only use the propagators which also exist in GMAT.
    // Refer to `regress_leo_day_adaptive` for the additional propagators.
    let mut all_props = vec![
        Propagator::new::<Verner56>(&Options::with_adaptive_step(min_step, max_step, accuracy)),
        // Propagator::new::<CashKarp45>(&Options::with_adaptive_step(min_step, max_step, 1e-12)),
        // Propagator::new::<Fehlberg45>(&Options::with_adaptive_step(min_step, max_step, 1e-12)),
        Propagator::new::<Dormand45>(&Options::with_adaptive_step(min_step, max_step, accuracy)),
        Propagator::new::<Dormand78>(&Options::with_adaptive_step(min_step, max_step, accuracy)),
        Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy)),
    ];
    let all_rslts = vec![
        Vector6::from_row_slice(&[
            -5971.194191669716,
            3945.506653225098,
            2864.6366184126337,
            0.04909695762959668,
            -4.185093318476122,
            5.848940867749713,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191673688,
            3945.506653306904,
            2864.636618302827,
            0.04909695753042703,
            -4.185093318409565,
            5.848940867795624,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191668567,
            3945.5066531626767,
            2864.636618498951,
            0.04909695770740798,
            -4.185093318527218,
            5.848940867713008,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191668705,
            3945.5066532187348,
            2864.636618421759,
            0.04909695763815305,
            -4.185093318482125,
            5.84894086774621,
        ]),
    ];
    let all_it_cnt = vec![3347, 864_000, 864_000, 864_000, 864_000, 2_880, 2_880];

    for (p_id, prop) in all_props.iter_mut().enumerate() {
        let mut init_state = Vector6::from_row_slice(&[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);
        let mut cur_t = 0.0;
        let mut iterations = 0;
        loop {
            let (t, state) = prop.derive(
                cur_t,
                &init_state,
                two_body_dynamics,
                error_ctrl::rss_state_pos_vel,
            );
            iterations += 1;
            if t < prop_time {
                // We haven't passed the time based stopping condition.
                cur_t = t;
                init_state = state;
            } else {
                // NOTE: The refined propagation time for here is different from what GMAT does.
                // In fact, GMAT uses the [secant method](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/command/Propagate.cpp#L5450)
                // even for time computation. This means there are several iterations over the derivatives, and a non-exact iteration time.
                // GMAT still uses a fixed step propagation though.
                // At this point, we've passed the condition, so let's switch to a fixed step of _exactly_ the
                // previous time step minus the amount by which we overshot. This allows us to propagate in time for
                // _exactly_ the time we want to propagate for.
                let prev_details = prop.latest_details().clone();
                let overshot = t - prop_time;
                prop.set_fixed_step(prev_details.step - overshot);
                // Take one final step
                let (t, state) = prop.derive(
                    cur_t,
                    &init_state,
                    two_body_dynamics,
                    error_ctrl::rss_state_pos_vel,
                );

                assert!(
                    (t - prop_time).abs() < 1e-12,
                    "propagated for {} instead of {}",
                    t,
                    prop_time
                );

                // Let's check that, prior to the refined step, we either hit the accuracy wanted,
                // or we are using the minimum step size.
                if prev_details.error > accuracy {
                    assert!(
                        prev_details.step - min_step < f64::EPSILON,
                        "step size should be at its minimum because error is higher than tolerance (p_id = {}): {:?}",
                        p_id,
                        prev_details
                    );
                }

                assert_eq!(
                    state,
                    all_rslts[p_id],
                    "leo prop failed for p_id = {}",
                    p_id
                );

                assert_eq!(
                    iterations,
                    all_it_cnt[p_id],
                    "wrong number of iterations (p_id = {})",
                    p_id
                );
                break;
            }
        }
    }
}

#[test]
fn gmat_val_leo_day_fixed() {
    extern crate nalgebra as na;
    use nyx::propagators::*;
    use self::na::Vector6;
    let mut all_props = vec![
        Propagator::new::<RK4Fixed>(&Options::with_fixed_step(1.0)),
        Propagator::new::<Verner56>(&Options::with_fixed_step(10.0)),
        Propagator::new::<Dormand45>(&Options::with_fixed_step(10.0)),
        Propagator::new::<Dormand78>(&Options::with_fixed_step(10.0)),
        Propagator::new::<RK89>(&Options::with_fixed_step(10.0)),
    ];
    let all_rslts = vec![
        Vector6::from_row_slice(&[
            -5971.194191670768,
            3945.506653227154,
            2864.6366184109706,
            0.04909695762764177,
            -4.18509331847428,
            5.8489408677500965,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191670203,
            3945.5066532190967,
            2864.636618421618,
            0.04909695763733907,
            -4.185093318480867,
            5.848940867745654,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191699656,
            3945.50665408017,
            2864.63661724545,
            0.04909695658406228,
            -4.185093317777894,
            5.848940868241106,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191670044,
            3945.5066532117953,
            2864.636618431374,
            0.049096957645996114,
            -4.185093318486724,
            5.848940867741533,
        ]),
        Vector6::from_row_slice(&[
            -5971.19419167081,
            3945.5066532332503,
            2864.6366184022418,
            0.049096957620019005,
            -4.185093318469214,
            5.848940867753748,
        ]),
    ];

    // let mut p_id: usize = 0; // We're using this as a propagation index in order to avoid modifying borrowed content
    for (p_id, prop) in all_props.iter_mut().enumerate() {
        let mut init_state = Vector6::from_row_slice(&[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);
        let mut cur_t = 0.0;
        loop {
            let (t, state) = prop.derive(
                cur_t,
                &init_state,
                two_body_dynamics,
                error_ctrl::rss_state_pos_vel,
            );
            cur_t = t;
            init_state = state;
            if cur_t >= 3600.0 * 24.0 {
                let details = prop.latest_details();
                if details.error > 1e-2 {
                    assert!(
                        details.step - 1e-1 < f64::EPSILON,
                        "step size should be at its minimum because error is higher than tolerance (p_id = {}): {:?}",
                        p_id,
                        details
                    );
                }
                println!("p_id={} => {:?}", p_id, prop.latest_details());

                assert_eq!(
                    state,
                    all_rslts[p_id],
                    "leo fixed prop failed for p_id = {}",
                    p_id
                );

                break;
            }
        }
    }
}
