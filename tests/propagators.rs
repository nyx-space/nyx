extern crate nalgebra as na;
extern crate nyx_space as nyx;
use self::na::{Vector6, U1, U3};
use std::f64;

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
    use self::na::Vector6;
    use nyx::propagators::*;

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 30.0;

    // NOTE: In this test we only use the propagators which also exist in GMAT.
    // Refer to `regress_leo_day_adaptive` for the additional propagators.
    let mut all_props = vec![
        Propagator::new::<RK2Fixed>(&PropOpts::with_fixed_step(1.0)),
        Propagator::new::<CashKarp45>(&PropOpts::with_adaptive_step(min_step, max_step, accuracy)),
        Propagator::new::<Fehlberg45>(&PropOpts::with_adaptive_step(min_step, max_step, accuracy)),
    ];

    let all_it_cnt = vec![86_400, 5_178, 5_935]; // NOTE: This is a decent estimate of which propagators to use depending if speed is important.

    let all_rslts = vec![
        Vector6::from_row_slice(&[
            -5971.198524908157,
            3945.7755093263054,
            2864.262542023422,
            0.04876621287986919,
            -4.184873956982518,
            5.849098380963502,
        ]),
        Vector6::from_row_slice(&[
            -5971.194190251459,
            3945.5066108335227,
            2864.636676354708,
            0.049097009502150256,
            -4.185093353076805,
            5.848940843324859,
        ]),
        Vector6::from_row_slice(&[
            -5971.194190305514,
            3945.5066123572196,
            2864.6366742764963,
            0.04909700763937082,
            -4.185093351832398,
            5.848940844199317,
        ]),
    ];

    for (p_id, prop) in all_props.iter_mut().enumerate() {
        let mut init_state = Vector6::from_row_slice(&[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);
        let mut cur_t = 0.0;
        let mut iterations = 0;
        loop {
            let (t, state) = prop.derive(cur_t, &init_state, two_body_dynamics, error_ctrl::rss_state_pos_vel);
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
                let (t, state) = prop.derive(cur_t, &init_state, two_body_dynamics, error_ctrl::rss_state_pos_vel);

                assert!((t - prop_time).abs() < 1e-12, "propagated for {} instead of {}", t, prop_time);

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

                assert_eq!(state, all_rslts[p_id], "leo prop failed for p_id = {}", p_id);

                assert_eq!(iterations, all_it_cnt[p_id], "wrong number of iterations (p_id = {})", p_id);
                break;
            }
        }
    }
}

#[test]
fn gmat_val_leo_day_adaptive() {
    extern crate nalgebra as na;
    use self::na::Vector6;
    use nyx::propagators::*;

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 30.0;

    // NOTE: In this test we only use the propagators which also exist in GMAT.
    // Refer to `regress_leo_day_adaptive` for the additional propagators.
    let mut all_props = vec![
        Propagator::new::<Dormand45>(&PropOpts::with_adaptive_step(min_step, max_step, accuracy)),
        Propagator::new::<Verner56>(&PropOpts::with_adaptive_step(min_step, max_step, accuracy)),
        Propagator::new::<Dormand78>(&PropOpts::with_adaptive_step(min_step, max_step, accuracy)),
        Propagator::new::<RK89>(&PropOpts::with_adaptive_step(min_step, max_step, accuracy)),
    ];

    let all_it_cnt = vec![5_412, 3_130, 2_880, 2_880]; // This number of iterations does not include the final refined fixed step.

    let all_rslts = vec![
        Vector6::from_row_slice(&[
            -5971.194191972316,
            3945.5066620394823,
            2864.636606375189,
            0.049096946846225495,
            -4.185093311278742,
            5.848940872821119,
        ]),
        Vector6::from_row_slice(&[
            -5971.19419167894,
            3945.5066538720525,
            2864.6366175103476,
            0.049096956828390714,
            -4.1850933179466505,
            5.848940868134205,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191670392,
            3945.506653218658,
            2864.63661842225,
            0.049096957637897856,
            -4.185093318481106,
            5.8489408677453,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191670676,
            3945.506653225158,
            2864.6366184134445,
            0.04909695762999346,
            -4.185093318475795,
            5.848940867748944,
        ]),
    ];

    for (p_id, prop) in all_props.iter_mut().enumerate() {
        let mut init_state = Vector6::from_row_slice(&[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);
        let mut cur_t = 0.0;
        let mut iterations = 0;
        loop {
            let (t, state) = prop.derive(cur_t, &init_state, two_body_dynamics, error_ctrl::rss_state_pos_vel);
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
                let (t, state) = prop.derive(cur_t, &init_state, two_body_dynamics, error_ctrl::rss_state_pos_vel);

                assert!((t - prop_time).abs() < 1e-12, "propagated for {} instead of {}", t, prop_time);

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

                assert_eq!(state, all_rslts[p_id], "leo prop failed for p_id = {}", p_id);

                assert_eq!(iterations, all_it_cnt[p_id], "wrong number of iterations (p_id = {})", p_id);
                break;
            }
        }
    }
}

#[test]
fn gmat_val_leo_day_fixed() {
    extern crate nalgebra as na;
    use self::na::Vector6;
    use nyx::propagators::*;
    let mut all_props = vec![
        Propagator::new::<RK4Fixed>(&PropOpts::with_fixed_step(1.0)),
        Propagator::new::<Verner56>(&PropOpts::with_fixed_step(10.0)),
        Propagator::new::<Dormand45>(&PropOpts::with_fixed_step(10.0)),
        Propagator::new::<Dormand78>(&PropOpts::with_fixed_step(10.0)),
        Propagator::new::<RK89>(&PropOpts::with_fixed_step(10.0)),
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
            let (t, state) = prop.derive(cur_t, &init_state, two_body_dynamics, error_ctrl::rss_state_pos_vel);
            cur_t = t;
            init_state = state;
            if cur_t >= 3_600.0 * 24.0 {
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

                assert_eq!(state, all_rslts[p_id], "leo fixed prop failed for p_id = {}", p_id);

                break;
            }
        }
    }
}
