extern crate nalgebra as na;
extern crate nyx;
use std::f64;
use self::na::{U1, U3, Vector6};

fn two_body_dynamics(_t: f64, state: &Vector6<f64>) -> Vector6<f64> {
    let radius = state.fixed_slice::<U3, U1>(0, 0);
    let velocity = state.fixed_slice::<U3, U1>(3, 0);
    let body_acceleration = (-398_600.4 / radius.norm().powi(3)) * radius;
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
        Propagator::new::<CashKarp45>(&Options::with_adaptive_step(0.1, 30.0, 1e-2)),
        Propagator::new::<Fehlberg45>(&Options::with_adaptive_step(0.01, 30.0, 1e-2)),
        Propagator::new::<Dormand45>(&Options::with_adaptive_step(0.1, 30.0, 1e-2)),
        Propagator::new::<Fehlberg56>(&Options::with_adaptive_step(0.1, 30.0, 1e-2)),
        Propagator::new::<Verner56>(&Options::with_adaptive_step(0.1, 30.0, 1e-2)),
        Propagator::new::<Dormand78>(&Options::with_adaptive_step(0.1, 30.0, 1e-2)),
        Propagator::new::<RK89>(&Options::with_adaptive_step(0.1, 30.0, 1e-2)),
    ];
    let all_rslts = vec![
        Vector6::from_row_slice(&[
            -5971.195448672869,
            3945.5831501786265,
            2864.530217443299,
            0.04900281804671473,
            -4.185030861894159,
            5.848985672431022,
        ]),
        Vector6::from_row_slice(&[
            -5971.195422075937,
            3945.5823545380795,
            2864.5313049901433,
            0.0490037916139182,
            -4.185031511302435,
            5.84898521402028,
        ]),
        Vector6::from_row_slice(&[
            -5971.19539794548,
            3945.581628280032,
            2864.5322979076436,
            0.04900468039431707,
            -4.185032104074844,
            5.848984795418991,
        ]),
        Vector6::from_row_slice(&[
            -5971.1954486729965,
            3945.5831501881476,
            2864.5302174302246,
            0.04900281803620559,
            -4.185030861887087,
            5.848985672436005,
        ]),
        Vector6::from_row_slice(&[
            -5971.195448672514,
            3945.5831501109246,
            2864.5302175392503,
            0.04900281813249359,
            -4.185030861949536,
            5.84898567238942,
        ]),
        Vector6::from_row_slice(&[
            -5971.195460543446,
            3945.5835028716965,
            2864.529735467152,
            0.04900238652205916,
            -4.185030573980658,
            5.848985875519705,
        ]),
        Vector6::from_row_slice(&[
            -5971.195448678458,
            3945.5831503373893,
            2864.530217225496,
            0.04900281785017413,
            -4.185030861762282,
            5.848985672523384,
        ]),
        Vector6::from_row_slice(&[
            -5971.195448675211,
            3945.583150208691,
            2864.530217404535,
            0.049002818010257174,
            -4.185030861868303,
            5.848985672447304,
        ]),
    ];
    let all_it_cnt = vec![86_400, 2880, 2880, 864_000, 864_000, 2880, 864_000, 864_000];

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
