extern crate nalgebra as na;
extern crate nyx_space as nyx;

#[test]
fn combining_dynamics() {
    // Testing that we can combine dynamics such that they use the same propagator.
    // Warning: this is arguably a bad example: attitude dynamics very significantly
    // faster than orbital mechanics. Hence you really should use different propagators
    // for the attitude and orbital position and velocity. However, this demonstrates
    // that combining several provided dynamics from `nyx::dynamics` works and is quite easy.
    use nyx::dynamics::Dynamics;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::celestia::EARTH;
    use self::na::{Matrix3, U9, Vector3, Vector6, VectorN};
    use nyx::dynamics::momentum::AngularMom;
    use nyx::propagators::{CashKarp45, Options, Propagator};

    // In the following struct, we only store the dynamics because this is only a proof of concept.
    // An engineer could add more useful information to this struct, such as a short cut to the position
    // or an attitude.
    #[derive(Copy, Clone)]
    pub struct PosVelAttMom {
        pub twobody: TwoBody,
        pub momentum: AngularMom,
        full_state: VectorN<f64, U9>, // TODO: Add a quaternion of a given orientation.
    }

    impl PosVelAttMom {
        pub fn init_state(&mut self) {
            let twobody_state = self.twobody.state();
            let momentum_state = self.momentum.state();
            // We're channing both states ti create a combined state.
            // The most important part here is make sure that the `state` and `set_state` handle the state in the same order.
            self.full_state = <VectorN<f64, U9>>::from_iterator(
                twobody_state.iter().chain(momentum_state.iter()).cloned(),
            );
        }
    }

    impl Dynamics for PosVelAttMom {
        type StateSize = U9;
        fn time(&self) -> f64 {
            // Both dynamical models have the same time because they share the propagator.
            self.twobody.time()
        }

        fn state(&self) -> &VectorN<f64, Self::StateSize> {
            &self.full_state
        }

        fn set_state(&mut self, new_t: f64, new_state: &VectorN<f64, Self::StateSize>) {
            // XXX: Reconstructing the Vector6 from scratch because for some reason it isn't the correct type... Must figure this out
            /*
                error[E0308]: mismatched types
      --> tests/lib.rs:45:43
       |
    45 |             self.twobody.set_state(new_t, pos_vel);
       |                                           ^^^^^^^ expected f64, found &f64
       |
       = note: expected type `&na::Matrix<f64, na::U6, na::U1, na::MatrixArray<f64, na::U6, na::U1>>`
              found type `&na::Matrix<&f64, na::U6, na::U1, na::MatrixArray<&f64, na::U6, na::U1>>`
              */
            self.full_state = *new_state;
            let mut pos_vel_vals = [0.0; 6];
            let mut mom_vals = [0.0; 3];
            for (i, val) in new_state.iter().enumerate() {
                if i < 6 {
                    pos_vel_vals[i] = *val;
                } else {
                    mom_vals[i - 6] = *val;
                }
            }
            self.twobody
                .set_state(new_t, &Vector6::from_row_slice(&pos_vel_vals));
            self.momentum
                .set_state(new_t, &Vector3::from_row_slice(&mom_vals));
        }

        fn eom(
            &self,
            _t: f64,
            state: &VectorN<f64, Self::StateSize>,
        ) -> VectorN<f64, Self::StateSize> {
            let mut pos_vel_vals = [0.0; 6];
            let mut mom_vals = [0.0; 3];
            for (i, val) in state.iter().enumerate() {
                if i < 6 {
                    pos_vel_vals[i] = *val;
                } else {
                    mom_vals[i - 6] = *val;
                }
            }
            let dpos_vel_dt = self.twobody
                .eom(_t, &Vector6::from_row_slice(&pos_vel_vals));
            let domega_dt = self.momentum.eom(_t, &Vector3::from_row_slice(&mom_vals));
            <VectorN<f64, U9>>::from_iterator(dpos_vel_dt.iter().chain(domega_dt.iter()).cloned())
        }
    }

    // Let's initialize our combined dynamics.

    let dyn_twobody = TwoBody::around(
        &Vector6::from_row_slice(&[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]),
        &EARTH,
    );

    let omega = Vector3::new(0.1, 0.4, -0.2);
    let tensor = Matrix3::new(10.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 2.0);
    let dyn_mom = AngularMom::from_tensor_matrix(&tensor, &omega);

    let mut full_model = PosVelAttMom {
        twobody: dyn_twobody,
        momentum: dyn_mom,
        full_state: <VectorN<f64, U9>>::zeros(),
    };
    full_model.init_state();

    let init_momentum = full_model.momentum.momentum().norm();
    let mom_tolerance = 1e-8;

    // And now let's define the propagator and propagate for a short amount of time.
    let mut prop = Propagator::new::<CashKarp45>(&Options::with_adaptive_step(0.01, 30.0, 1e-12));
    // let mut prop = Propagator::new::<RK4Fixed>(&Options::with_fixed_step(0.1));

    // And propagate
    loop {
        let (t, state) = prop.derive(
            full_model.time(),
            full_model.state(),
            |t_: f64, state_: &VectorN<f64, U9>| full_model.eom(t_, state_),
        );
        full_model.set_state(t, &state);
        if full_model.time() >= 3600.0 {
            println!("{:?}", prop.latest_details());
            println!("{}", full_model.state());
            let delta_mom =
                ((full_model.momentum.momentum().norm() - init_momentum) / init_momentum).abs();
            if delta_mom > mom_tolerance {
                panic!(
                    "angular momentum prop failed: momentum changed by {:e} (> {:e})",
                    delta_mom, mom_tolerance
                );
            }
            break;
        }
    }
}
