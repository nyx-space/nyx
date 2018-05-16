extern crate nalgebra as na;
extern crate nyx_space as nyx;
use std::f64;

#[test]
fn two_body_j2_state_parametrized() {
    /* NOTE: We only test the J2 paramaters here for the JGM3 models. */
    extern crate nalgebra as na;
    use nyx::propagators::{Dormand45, Options, Propagator};
    use nyx::dynamics::Dynamics;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::dynamics::gravity::Harmonics;
    use nyx::io::gravity::MemoryBackend;
    use nyx::celestia::{State, EARTH};
    use self::na::{U6, VectorN};

    // TODO: Provide a cleaner way to wrap these, probably by implementing the std::ops::Add.
    // Or at the very least provide some template scenario which already has the dynamics included.
    #[derive(Clone)]
    pub struct J2Dyn {
        pub twobody: TwoBody,
        pub harmonics: Harmonics<MemoryBackend>,
    }

    impl Dynamics for J2Dyn {
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

        fn eom(
            &self,
            _t: f64,
            state: &VectorN<f64, Self::StateSize>,
        ) -> VectorN<f64, Self::StateSize> {
            self.twobody.eom(_t, state) + self.harmonics.eom(_t, state)
        }
    }

    let initial_state =
        State::from_cartesian::<EARTH>(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0);

    println!("Initial state:\n{0}\n{0:o}\n", initial_state);

    let rslt = State::from_cartesian::<EARTH>(
        -5751.473991026703,
        4721.214036506542,
        2045.947766308757,
        -0.797740263736765,
        -3.656451982107864,
        6.139637922282735,
    );

    let mut prop = Propagator::new::<Dormand45>(&Options::with_adaptive_step(0.1, 30.0, 1e-12));

    let mut dyn = J2Dyn {
        twobody: TwoBody::from_state_vec::<EARTH>(&initial_state.to_cartesian_vec()),
        harmonics: Harmonics::from_stor::<EARTH>(MemoryBackend::j2_jgm3()),
    };
    let final_state: State;
    loop {
        let (t, state) = prop.derive(
            dyn.time(),
            &dyn.state(),
            |t_: f64, state_: &VectorN<f64, U6>| dyn.eom(t_, state_),
        );
        dyn.set_state(t, &state);
        if dyn.time() >= 3600.0 * 24.0 {
            let details = prop.latest_details();
            if details.error > 1e-2 {
                assert!(
                        details.step - 1e-1 < f64::EPSILON,
                        "step size should be at its minimum because error is higher than tolerance: {:?}",
                        details
                    );
            }
            final_state = State::from_cartesian_vec::<EARTH>(&dyn.state());
            assert_eq!(
                final_state,
                rslt,
                "two body prop failed\nexpected: {0:o}\ncomputed: {1:o}",
                rslt,
                final_state
            );
            break;
        }
    }
    println!("Final state:\n{0}\n{0:o}", final_state);
}
