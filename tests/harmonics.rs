extern crate nalgebra as na;
extern crate nyx_space as nyx;
use std::f64;

#[test]
fn gmat_val_harmonics_basicj2() {
    // WARNING: This test shows a significant difference between GMAT and `nyx`. Refer to the VALIDATION.md file for a discussion.
    extern crate nalgebra as na;
    use self::na::{U6, VectorN};
    use nyx::celestia::{State, EARTH};
    use nyx::dynamics::Dynamics;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::dynamics::gravity::BasicJ2;
    use nyx::propagators::{error_ctrl, Options, Propagator, RK89};

    #[derive(Clone)]
    pub struct J2Dyn {
        count: u64,
        pub twobody: TwoBody,
        pub j2: BasicJ2,
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

        fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
            self.twobody.eom(_t, state) + self.j2.eom(_t, state)
        }
    }

    let initial_state = State::from_cartesian::<EARTH>(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0);

    println!("Initial state:\n{0}\n{0:o}\n", initial_state);

    let rslt = State::from_cartesian::<EARTH>(
        -5751.473991327453,
        4721.214035136977,
        2045.947768607068,
        -0.7977402618610976,
        -3.656451983635231,
        6.139637921621191,
    );

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 60.0;

    let mut prop = Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy));

    let mut dyn = J2Dyn {
        count: 0,
        twobody: TwoBody::from_state_vec::<EARTH>(&initial_state.to_cartesian_vec()),
        j2: BasicJ2::around_earth(),
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
                println!("overshot by {} seconds", overshot);
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

            println!("{} ==> {:?}", dyn.time(), prev_details);
            final_state = State::from_cartesian_vec::<EARTH>(&dyn.state());
            let rss_pos =
                ((rslt.x - final_state.x).powi(2) + (rslt.y - final_state.y).powi(2) + (rslt.z - final_state.z).powi(2)).sqrt();

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

#[cfg(feature = "broken-harmonics")]
#[test]
fn gmat_val_harmonics_j2jgm3() {
    // WARNING: This test shows a significant difference between GMAT and `nyx`. Refer to the VALIDATION.md file for a discussion.
    extern crate nalgebra as na;
    use self::na::{U6, VectorN};
    use nyx::celestia::{State, EARTH};
    use nyx::dynamics::Dynamics;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::dynamics::gravity::Harmonics;
    use nyx::io::gravity::MemoryBackend;
    use nyx::propagators::{error_ctrl, Options, Propagator, RK89};

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

        fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
            self.twobody.eom(_t, state) + self.harmonics.eom(_t, state)
        }
    }

    let initial_state = State::from_cartesian::<EARTH>(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0);

    println!("Initial state:\n{0}\n{0:o}\n", initial_state);

    let rslt = State::from_cartesian::<EARTH>(
        -5751.473991328457,
        4721.214035131441,
        2045.947768616195,
        -0.7977402618536559,
        -3.656451983641397,
        6.139637921618609,
    );

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 60.0;

    let mut prop = Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy));

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
            error_ctrl::rss_step_pos_vel,
        );
        if t < prop_time {
            // We haven't passed the time based stopping condition.
            dyn.set_state(t, &state);
        } else {
            let prev_details = prop.latest_details().clone();
            let overshot = t - prop_time;
            if overshot > 0.0 {
                println!("overshot by {} seconds", overshot);
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

            println!("{} ==> {:?}", dyn.time(), prev_details);
            final_state = State::from_cartesian_vec::<EARTH>(&dyn.state());
            let rss_pos =
                ((rslt.x - final_state.x).powi(2) + (rslt.y - final_state.y).powi(2) + (rslt.z - final_state.z).powi(2)).sqrt();

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

#[cfg(feature = "broken-harmonics")]
#[test]
fn gmat_val_harmonics_21x21() {
    /* NOTE: We only test the J2 paramaters here for the JGM3 models. */
    extern crate nalgebra as na;
    use self::na::{U6, VectorN};
    use nyx::celestia::{State, EARTH};
    use nyx::dynamics::Dynamics;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::dynamics::gravity::Harmonics;
    use nyx::io::gravity::MemoryBackend;
    use nyx::propagators::{error_ctrl, Options, Propagator, RK89};

    let filepath = "./data/JGM3.cof.gz";

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

        fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
            self.twobody.eom(_t, state) + self.harmonics.eom(_t, state)
        }
    }

    let initial_state = State::from_cartesian::<EARTH>(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0);

    println!("Initial state:\n{0}\n{0:o}\n", initial_state);

    let rslt = State::from_cartesian::<EARTH>(
        -5751.473991328457,
        4721.214035131441,
        2045.947768616195,
        -0.7977402618536559,
        -3.656451983641397,
        6.139637921618609,
    );

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 60.0;

    let mut prop = Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy));

    let mut dyn = J2Dyn {
        count: 0,
        twobody: TwoBody::from_state_vec::<EARTH>(&initial_state.to_cartesian_vec()),
        harmonics: Harmonics::from_stor::<EARTH>(MemoryBackend::from_cof(filepath.to_owned(), 21, 21, true)),
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
                println!("overshot by {} seconds", overshot);
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

            println!("{} ==> {:?}", dyn.count, prev_details);
            final_state = State::from_cartesian_vec::<EARTH>(&dyn.state());
            let rss_pos =
                ((rslt.x - final_state.x).powi(2) + (rslt.y - final_state.y).powi(2) + (rslt.z - final_state.z).powi(2)).sqrt();

            assert_eq!(
                rslt, final_state,
                "21x21 prop failed: RSS = {} km\nexpected: {:o}\ncomputed: {:o}",
                rss_pos, rslt, final_state
            );
            break;
        }
    }
    println!("Final state:\n{0}\n{0:o}", final_state);
}

#[cfg(feature = "broken-harmonics")]
#[test]
fn gmat_val_harmonics_70x70() {
    extern crate nalgebra as na;
    use self::na::{U6, VectorN};
    use nyx::celestia::{State, EARTH};
    use nyx::dynamics::Dynamics;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::dynamics::gravity::Harmonics;
    use nyx::io::gravity::MemoryBackend;
    use nyx::propagators::{error_ctrl, Options, Propagator, RK89};

    let filepath = "./data/JGM3.cof.gz";

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

        fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
            self.twobody.eom(_t, state) + self.harmonics.eom(_t, state)
        }
    }

    let initial_state = State::from_cartesian::<EARTH>(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0);

    println!("Initial state:\n{0}\n{0:o}\n", initial_state);

    let rslt = State::from_cartesian::<EARTH>(
        -5751.924618075773,
        4719.386612445624,
        2048.696011815599,
        -0.795383404372165,
        -3.658301183314282,
        6.138865498490013,
    );

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 60.0;

    let mut prop = Propagator::new::<RK89>(&Options::with_adaptive_step(min_step, max_step, accuracy));

    let mut dyn = J2Dyn {
        count: 0,
        twobody: TwoBody::from_state_vec::<EARTH>(&initial_state.to_cartesian_vec()),
        harmonics: Harmonics::from_stor::<EARTH>(MemoryBackend::from_cof(filepath.to_owned(), 70, 70, true)),
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
                println!("overshot by {} seconds", overshot);
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

            println!("{} ==> {:?}", dyn.count, prev_details);
            final_state = State::from_cartesian_vec::<EARTH>(&dyn.state());
            let rss_pos =
                ((rslt.x - final_state.x).powi(2) + (rslt.y - final_state.y).powi(2) + (rslt.z - final_state.z).powi(2)).sqrt();

            assert_eq!(
                rslt, final_state,
                "70x70 prop failed: RSS = {} km\nexpected: {:o}\ncomputed: {:o}",
                rss_pos, rslt, final_state
            );
            break;
        }
    }
    println!("Final state:\n{0}\n{0:o}", final_state);
}
