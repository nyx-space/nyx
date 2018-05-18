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
        count: u64,
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
            self.count += 1;
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

    // One second
    /*let rslt = State::from_cartesian::<EARTH>(
        -2431.360331498824,
        -2441.537552025363,
        6891.034000290618,
        5.090725265477195,
        -5.086492314158094,
        -0.00599941833051744,
    );*/

    // One minute
    let rslt = State::from_cartesian::<EARTH>(
        -2127.483782786232,
        -2737.798886616517,
        6880.240862171737,
        5.207578930784956,
        -4.953733998024182,
        -0.3597773824106297,
    );

    // Ten thousand seconds
    /*let rslt = State::from_cartesian::<EARTH>(
        3075.869682748079,
        1776.782315795178,
        -6870.826570868569,
        -4.744831035692495,
        5.320693700850615,
        -0.7495049625159994,
    );*/
    /*
    // One day:
    let rslt = State::from_cartesian::<EARTH>(
        -5751.47399154785,
        4721.214034143153,
        2045.947770285066,
        -0.7977402605029352,
        -3.656451984758527,
        6.139637921125921,
    );*/
    // 10 days:
    /*
    let rslt = State::from_cartesian::<EARTH>(
        -3453.82913576055,
        -613.1058873700827,
        6859.70170958988,
        2.418815931304214,
        -6.749934134121439,
        0.6126057059401954,
    );*/
    // 100 days:
    /*
    let rslt = State::from_cartesian::<EARTH>(
        1616.297423706993,
        -6402.348787298117,
        3984.68419222919,
        2.861916092325618,
        4.003867304278716,
        5.246167325614458,
    );*/

    let mut prop = Propagator::new::<Dormand45>(&Options::with_adaptive_step(0.1, 30.0, 1e-12));

    let mut dyn = J2Dyn {
        count: 0,
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
        if dyn.time() >= 60.0 {
            println!("{:?}", 21545.00000039794 + dyn.time() / 86400.0);
            let details = prop.latest_details();
            if details.error > 1e-2 {
                assert!(
                        details.step - 1e-1 < f64::EPSILON,
                        "step size should be at its minimum because error is higher than tolerance: {:?}",
                        details
                    );
            }
            println!("{} ==> {:?}", dyn.count, details);
            final_state = State::from_cartesian_vec::<EARTH>(&dyn.state());
            let rss_pos = ((rslt.x - final_state.x).powi(2) + (rslt.y - final_state.y).powi(2)
                + (rslt.z - final_state.z).powi(2))
                .sqrt();
            assert_eq!(
                rslt, final_state,
                "J2 JGM3 prop failed: RSS = {} km\nexpected: {:o}\ncomputed: {:o}",
                rss_pos, rslt, final_state
            );
            break;
        }
    }
    println!("Final state:\n{0}\n{0:o}", final_state);
}
