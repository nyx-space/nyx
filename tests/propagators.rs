extern crate nalgebra as na;
extern crate nyx;
use self::na::DVector;

fn two_body_dynamics(_t: f64, state: &DVector<f64>) -> DVector<f64> {
    let radius = state.slice((0, 0), (3, 1)); // TODO: Change to compile time slice
    let velocity = state.slice((3, 0), (3, 1));
    let body_acceleration = (-398600.4 / radius.norm().powi(3)) * radius;
    DVector::from_iterator(6, velocity.iter().chain(body_acceleration.iter()).cloned())
}

#[test]
fn geo_day_prop() {
    extern crate nyx;
    extern crate nalgebra as na;
    use nyx::propagators::*;
    use self::na::DVector;
    let all_props = vec![
        Propagator::new::<RK4Fixed>(Options::with_fixed_step(1.0)),
        Propagator::new::<RKF54>(Options::with_adaptive_step(0.1, 30.0, 1e-2)),
        Propagator::new::<RKCK54>(Options::with_adaptive_step(0.1, 30.0, 1e-2)),
    ];
    let all_rslts = vec![
        DVector::from_row_slice(
            6,
            &[
                -5971.195448672869,
                3945.5831501786265,
                2864.530217443299,
                0.04900281804671473,
                -4.185030861894159,
                5.848985672431022,
            ],
        ),
        DVector::from_row_slice(
            6,
            &[
                -5971.195448672315,
                3945.583150181142,
                2864.530217437198,
                0.0490028180426569,
                -4.185030861891767,
                5.848985672433508,
            ],
        ),
        DVector::from_row_slice(
            6,
            &[
                -5971.195448670072,
                3945.5831501339185,
                2864.5302175025718,
                0.0490028181000867,
                -4.185030861930267,
                5.848985672407197,
            ],
        ),
    ];
    let all_it_cnt = vec![86400, 864000, 864000];
    let mut p_id: usize = 0; // We're using this as a propagation index in order to avoid modifying borrowed content
    for mut prop in all_props {
        let mut init_state =
            DVector::from_row_slice(6, &[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);
        let mut cur_t = 0.0;
        let mut iterations = 0;
        loop {
            let (t, state) = prop.derive(cur_t, init_state.clone(), two_body_dynamics);
            iterations += 1;
            cur_t = t;
            init_state = state.clone();
            if p_id > 0 {
                // Check that the error is less than the max error.
                let details = prop.clone().latest_details();
                assert!(
                    details.error < 1e-1,
                    "error larger than expected (p_id = {})",
                    p_id
                );
                assert_eq!(
                    details.step,
                    1e-1,
                    "step size should be at its minimum (p_id = {})",
                    p_id
                );
            }
            if cur_t >= 3600.0 * 24.0 {
                assert_eq!(
                    iterations,
                    all_it_cnt[p_id],
                    "wrong number of iterations (p_id = {})",
                    p_id
                );
                assert_eq!(
                    state,
                    all_rslts[p_id],
                    "geo prop failed for p_id = {}",
                    p_id
                );
                break;
            }
        }
        p_id += 1;
    }
}
