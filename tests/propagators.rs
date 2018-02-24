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
    let opts = Options::with_fixed_step(1.0);
    let mut prop = Propagator::new::<RK4Fixed>(opts);
    let mut init_state =
        DVector::from_row_slice(6, &[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);
    let mut cur_t = 0.0;
    loop {
        let (t, state) = prop.derive(cur_t, init_state.clone(), two_body_dynamics);
        cur_t = t;
        init_state = state.clone();
        if cur_t >= 3600.0 * 24.0 {
            let exp_state = DVector::from_row_slice(
                6,
                &[
                    -5971.195448672869,
                    3945.5831501786265,
                    2864.530217443299,
                    0.04900281804671473,
                    -4.185030861894159,
                    5.848985672431022,
                ],
            );
            assert_eq!(state, exp_state, "geo prop failed");
            break;
        }
    }
}
