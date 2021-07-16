extern crate nyx_space as nyx;
use nyx::cosmic::{Cosm, Orbit};
use nyx::dimensions::Matrix6;
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::propagators::*;
use nyx::time::{Epoch, TimeUnit};

#[test]
fn linear_stm_step_traj() {
    // From an initial state in a linear-like regime, check that the STM computed in a single step is a decent approximation of the time-varying dynamics

    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let init = Orbit::keplerian(8000.0, 1e-5, 10.0, 5.0, 25.0, 0.0, epoch, eme2k);

    println!("init.to_cartesian_vec() = {}", init.to_cartesian_vec());

    let prop = Propagator::new::<RK4Fixed>(
        OrbitalDynamics::two_body(),
        PropOpts::with_fixed_step(10 * TimeUnit::Second),
    );

    let ten_steps = prop
        .with(init.with_stm())
        .for_duration(100 * TimeUnit::Second)
        .unwrap();

    let nominal = ten_steps.to_cartesian_vec();

    // let stm_err_ten_steps =
    //     ten_steps.stm() * init.to_cartesian_vec() - ten_steps.to_cartesian_vec();

    let stm_err_ten_steps = ten_steps.stm() * init.to_cartesian_vec() - nominal;

    println!(
        "ten_steps.to_cartesian_vec() = {}",
        ten_steps.to_cartesian_vec()
    );
    println!("ten_steps.stm() = {}", ten_steps.stm());

    println!("stm_err_ten_steps = {}", stm_err_ten_steps);

    // Compare with finite differencing
    let mut stm_fd = Matrix6::<f64>::zeros();
    let pert = 0.0001;

    for i in 0..6 {
        let mut this_init = init.with_stm();
        match i {
            0 => this_init.x += pert,
            1 => this_init.y += pert,
            2 => this_init.z += pert,
            3 => this_init.vx += pert,
            4 => this_init.vy += pert,
            5 => this_init.vz += pert,
            _ => unreachable!(),
        }

        let these_ten_steps = prop
            .with(this_init)
            .for_duration(100 * TimeUnit::Second)
            .unwrap();

        let jac_val = (these_ten_steps.to_cartesian_vec() - nominal) / pert;

        for j in 0..6 {
            stm_fd[(j, i)] = jac_val[j];
        }
    }

    println!("stm_fd = {}", stm_fd);

    let stm_fd_err_ten_steps = stm_fd * init.to_cartesian_vec() - nominal;

    println!("stm_fd_err_ten_steps = {}", stm_fd_err_ten_steps);

    println!(
        "Error between hyperdual and finite diff: {}",
        stm_err_ten_steps - stm_fd_err_ten_steps
    );
}
