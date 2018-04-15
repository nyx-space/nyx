extern crate nalgebra as na;
extern crate nyx_space as nyx;

#[test]
fn const_mom() {
    use nyx::propagators::{Options, Propagator, RK4Fixed};
    use nyx::dynamics::Dynamics;
    use nyx::dynamics::momentum::AngularMom;
    use self::na::{Matrix3, Vector3};

    let mut prop = Propagator::new::<RK4Fixed>(&Options::with_fixed_step(0.1));

    let omega = Vector3::new(0.1, 0.4, -0.2);
    let tensor = Matrix3::new(10.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 2.0);
    let tolerance = 1e-8;

    let mut dyn = AngularMom::from_tensor_matrix(&tensor, &omega);
    let init_momentum = dyn.momentum().norm();
    loop {
        let (t, state) = prop.derive(
            dyn.time(),
            &dyn.state(),
            |t_: f64, state_: &Vector3<f64>| dyn.eom(t_, state_),
        );
        dyn.set_state(t, &state);
        if dyn.time() >= 5.0 {
            println!("{:?}", prop.latest_details());
            let delta_mom = ((dyn.momentum().norm() - init_momentum) / init_momentum).abs();
            if delta_mom > tolerance {
                panic!(
                    "angular momentum prop failed: momentum changed by {:e} (> {:e})",
                    delta_mom, tolerance
                );
            }
            break;
        }
    }
}

#[test]
#[should_panic]
fn non_diag_tensor() {
    use self::na::{Matrix3, Vector3};
    use nyx::dynamics::momentum::AngularMom;

    let omega = Vector3::new(0.1, 0.4, -0.2);
    let tensor = Matrix3::new(10.0, 1.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 2.0);

    AngularMom::from_tensor_matrix(&tensor, &omega);
}
