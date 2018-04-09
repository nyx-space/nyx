extern crate nalgebra as na;
extern crate nyx;
use std::f64;
use self::na::{U1, U3, Vector6};

fn two_body_dynamics(_t: f64, state: &Vector6<f64>) -> Vector6<f64> {
    let radius = state.fixed_slice::<U3, U1>(0, 0);
    let velocity = state.fixed_slice::<U3, U1>(3, 0);
    let body_acceleration = (-398_600.441500000015366822 / radius.norm().powi(3)) * radius;
    Vector6::from_iterator(velocity.iter().chain(body_acceleration.iter()).cloned())
}

#[test]
fn geo_day_prop() {
    extern crate nyx;
    extern crate nalgebra as na;
    use nyx::propagators::*;
    use self::na::Vector6;
    let all_props = vec![
        Propagator::new::<RK4Fixed>(&Options::with_fixed_step(1.0)),
        Propagator::new::<CashKarp45>(&Options::with_adaptive_step(0.1, 30.0, 1e-12)),
        Propagator::new::<Fehlberg45>(&Options::with_adaptive_step(0.1, 30.0, 1e-12)),
        Propagator::new::<Dormand45>(&Options::with_adaptive_step(0.1, 30.0, 1e-12)),
        Propagator::new::<Fehlberg56>(&Options::with_adaptive_step(0.1, 30.0, 1e-12)),
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
            -5971.1763060616095,
            3944.0101437173926,
            2866.727622126912,
            0.05095441538128432,
            -4.1863204152348725,
            5.848049439457298,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191670006,
            3945.5066531852003,
            2864.6366184678586,
            0.04909695767907051,
            -4.185093318508363,
            5.848940867725646,
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
    let all_it_cnt = vec![
        86_400, 864_000, 864_000, 864_000, 63_590, 864_000, 2_880, 2_880
    ];

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
                    let details = prop.clone().latest_details();
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
