extern crate nalgebra as na;
extern crate nyx_space as nyx;
use std::f64;

#[cfg(feature = "unvalidated")]
#[test]
fn basic_drag() {
    extern crate nalgebra as na;
    use self::na::{U6, VectorN};
    use nyx::celestia::{State, EARTH};
    use nyx::dynamics::celestial::TwoBody;
    use nyx::dynamics::drag::BasicDrag;
    use nyx::dynamics::Dynamics;
    use nyx::propagators::{error_ctrl, Options, Propagator, RK89};

    #[derive(Clone)]
    pub struct SimpleDrag {
        pub twobody: TwoBody,
        pub drag: BasicDrag,
    }

    impl Dynamics for SimpleDrag {
        type StateSize = U6;
        fn time(&self) -> f64 {
            // Both dynamical models have the same time because they share the propagator.
            self.twobody.time()
        }

        fn state(&self) -> VectorN<f64, Self::StateSize> {
            self.twobody.state() // Harmonics do not have a state of their own
        }

        fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
            self.twobody.set_state(new_t, new_state);
        }

        fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
            self.twobody.eom(_t, state) + self.drag.eom(_t, state)
        }
    }

    let initial_state = State::from_cartesian::<EARTH>(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0);

    println!("Initial state:\n{0}\n{0:o}\n", initial_state);

    let rslt = State::from_cartesian::<EARTH>(
        -5971.1941916712285,
        3945.5066532419537,
        2864.636618390466,
        0.04909695760948815,
        -4.1850933184621315,
        5.848940867758592,
    );

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 60.0;

    let mut prop = Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy));

    let mut dyn = SimpleDrag {
        twobody: TwoBody::from_state_vec::<EARTH>(&initial_state.to_cartesian_vec()),
        drag: BasicDrag {
            rho: 7e-15, // XXX: This is a dummy value which needs to be updated to something more realistic
            cd: 2.2,
            area: 15.0,
            mass: 850.0,
        },
    };

    let final_state: State;
    loop {
        let (t, state) = prop.derive(
            dyn.time(),
            &dyn.state(),
            |t_: f64, state_: &VectorN<f64, U6>| dyn.eom(t_, state_),
            error_ctrl::rss_step_pos_vel,
        );
        if t < prop_time {
            // We haven't passed the time based stopping condition.
            dyn.set_state(t, &state);
        } else {
            let prev_details = prop.latest_details().clone();
            let overshot = t - prop_time;
            if overshot > 0.0 {
                prop.set_fixed_step(prev_details.step - overshot);
                // Take one final step
                let (t, state) = prop.derive(
                    dyn.time(),
                    &dyn.state(),
                    |t_: f64, state_: &VectorN<f64, U6>| dyn.eom(t_, state_),
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
            let rss_pos =
                ((rslt.x - final_state.x).powi(2) + (rslt.y - final_state.y).powi(2) + (rslt.z - final_state.z).powi(2)).sqrt();

            assert_eq!(
                rslt, final_state,
                "Drag: RSS = {} km\nexpected: {:o}\ncomputed: {:o}",
                rss_pos, rslt, final_state
            );
            break;
        }
    }
    println!("Final state:\n{0}\n{0:o}", final_state);
}
