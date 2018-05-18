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
fn leo_day_prop() {
    extern crate nalgebra as na;
    use nyx::propagators::*;
    use self::na::Vector6;
    let all_props = vec![
        Propagator::new::<RK4Fixed>(&Options::with_fixed_step(1.0)),
        Propagator::new::<CashKarp45>(&Options::with_adaptive_step(0.1, 30.0, 1e-12)),
        Propagator::new::<Fehlberg45>(&Options::with_adaptive_step(0.1, 30.0, 1e-12)),
        Propagator::new::<Dormand45>(&Options::with_adaptive_step(0.1, 30.0, 1e-12)),
        Propagator::new::<Verner56>(&Options::with_adaptive_step(0.1, 30.0, 1e-12)),
        Propagator::new::<Dormand78>(&Options::with_adaptive_step(0.1, 30.0, 1e-12)),
        Propagator::new::<RK89>(&Options::with_adaptive_step(0.1, 30.0, 1e-12)),
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
            -5971.194191665485,
            3945.506653102251,
            2864.6366185816864,
            0.049096957780097136,
            -4.1850933185772945,
            5.8489408676776895,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191669591,
            3945.5066531527896,
            2864.6366185138795,
            0.0490969577185288,
            -4.185093318534312,
            5.848940867706425,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191667916,
            3945.506653191027,
            2864.636618460324,
            0.04909695767239903,
            -4.185093318504468,
            5.848940867730048,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191668705,
            3945.5066532187348,
            2864.636618421759,
            0.04909695763815305,
            -4.185093318482125,
            5.84894086774621,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191674274,
            3945.506653344657,
            2864.6366182499773,
            0.04909695748166678,
            -4.185093318377827,
            5.848940867819183,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191672271,
            3945.5066531862226,
            2864.636618471764,
            0.04909695768011727,
            -4.1850933185073345,
            5.848940867722855,
        ]),
    ];
    let all_it_cnt = vec![86_400, 864_000, 864_000, 864_000, 864_000, 864_000, 864_000];

    let mut p_id: usize = 0; // We're using this as a propagation index in order to avoid modifying borrowed content
    for mut prop in all_props {
        let mut init_state =
            Vector6::from_row_slice(&[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);
        let mut cur_t = 0.0;
        let mut iterations = 0;
        loop {
            let (t, state) = prop.derive(cur_t, &init_state, two_body_dynamics);
            iterations += 1;
            cur_t = t;
            init_state = state;
            if cur_t >= 3600.0 * 24.0 {
                if p_id > 0 {
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
                }
                assert_eq!(
                    state,
                    all_rslts[p_id],
                    "geo prop failed for p_id = {}",
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
        p_id += 1;
    }
}
