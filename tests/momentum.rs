extern crate nalgebra as na;
extern crate nyx_space as nyx;

#[test]
fn const_mom() {
    use self::na::{Matrix3, Vector3};
    use nyx::dynamics::momentum::AngularMom;
    use nyx::dynamics::Dynamics;
    use nyx::propagators::error_ctrl::LargestStep;
    use nyx::propagators::{CashKarp45, PropOpts, Propagator};

    let omega = Vector3::new(0.1, 0.4, -0.2);
    let tensor = Matrix3::new(10.0, 0.0, 0.0, 0.0, 5.0, 0.0, 0.0, 0.0, 2.0);
    let tolerance = 1e-8;

    let mut dynamics = AngularMom::from_tensor_matrix(&tensor, &omega);
    let init_momentum = dynamics.momentum().norm();

    let mut prop = Propagator::new::<CashKarp45>(&mut dynamics, &PropOpts::with_adaptive_step(0.1, 5.0, 1e-8, LargestStep {}));

    prop.until_time_elapsed(5.0);

    println!("{:?}", prop.latest_details());
    let delta_mom = ((prop.dynamics.momentum().norm() - init_momentum) / init_momentum).abs();
    if delta_mom > tolerance {
        panic!(
            "angular momentum prop failed: momentum changed by {:e} (> {:e})",
            delta_mom, tolerance
        );
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
