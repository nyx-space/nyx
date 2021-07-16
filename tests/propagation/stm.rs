extern crate nyx_space as nyx;
use nyx::cosmic::{Bodies, Cosm, Orbit};
use nyx::dimensions::{Matrix2, Matrix6};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::propagators::*;
use nyx::time::{Epoch, TimeUnit};
use nyx::State;

// These tests compare the computation of the state transition matrix between the finite differencing methoid (common) and hyperdual numbers.
// Conclusion: hyperdual numbers lead to less error than finite differencing.

#[test]
fn stm_fixed_step() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let prop = Propagator::new::<RK4Fixed>(
        OrbitalDynamics::two_body(),
        PropOpts::with_fixed_step(10 * TimeUnit::Second),
    );

    // First test is in mostly linear regime (low eccentricity)
    // The second is quite non linear because the eccentricity is 0.2 and we start at periapsis
    let eccs = vec![1e-5, 0.2];

    for ecc in eccs {
        let init = Orbit::keplerian(8000.0, ecc, 10.0, 5.0, 25.0, 0.0, epoch, eme2k);

        let ten_steps = prop
            .with(init.with_stm())
            .for_duration(100 * TimeUnit::Second)
            .unwrap();

        let nominal = ten_steps.to_cartesian_vec();

        let stm_err_ten_steps = ten_steps.stm() * init.to_cartesian_vec() - nominal;

        println!("HD ten_steps.stm() = {}", ten_steps.stm());

        println!("HD stm_err_ten_steps = {}", stm_err_ten_steps);

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

        println!("FD stm_fd = {}", stm_fd);

        let stm_fd_err_ten_steps = stm_fd * init.to_cartesian_vec() - nominal;

        println!("FD stm_fd_err_ten_steps = {}", stm_fd_err_ten_steps);

        let delta = stm_err_ten_steps - stm_fd_err_ten_steps;

        println!("Error between hyperdual and finite diff: {}", delta);

        // Check that each component is less than 1 km or 1 km/s of difference, OR that the hyperdual computation is better than the finite differencing
        for i in 0..6 {
            if delta[i].abs() > 1.0 {
                assert!(delta[i] < 0.0, "FD more precise than HD?! {}", delta);
            }
        }
    }
}

#[test]
fn stm_variable_step() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let prop = Propagator::default(OrbitalDynamics::two_body());

    let eccs = vec![1e-5, 0.2];

    // First test is in mostly linear regime (low eccentricity)
    // The second is quite non linear because the eccentricity is 0.2 and we start at periapsis
    for ecc in eccs {
        let init = Orbit::keplerian(8000.0, ecc, 10.0, 5.0, 25.0, 0.0, epoch, eme2k);

        let ten_steps = prop
            .with(init.with_stm())
            .for_duration(100 * TimeUnit::Second)
            .unwrap();

        let nominal = ten_steps.to_cartesian_vec();

        let stm_err_ten_steps = ten_steps.stm() * init.to_cartesian_vec() - nominal;

        println!("HD ten_steps.stm() = {}", ten_steps.stm());

        println!("HD stm_err_ten_steps = {}", stm_err_ten_steps);

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

        println!("FD stm_fd = {}", stm_fd);

        let stm_fd_err_ten_steps = stm_fd * init.to_cartesian_vec() - nominal;

        println!("FD stm_fd_err_ten_steps = {}", stm_fd_err_ten_steps);

        let delta = stm_err_ten_steps - stm_fd_err_ten_steps;

        println!("Error between hyperdual and finite diff: {}", delta);

        // Check that each component is less than 1 km or 1 km/s of difference, OR that the hyperdual computation is better than the finite differencing
        for i in 0..6 {
            if delta[i].abs() > 1.0 {
                assert!(delta[i] < 0.0, "FD more precise than HD?! {}", delta);
            }
        }
    }
}

#[test]
fn stm_between_steps() {
    // Check that \Phi(t_2, t_1) = \Phi(t_2, t_0) * \Phi^{-1}(t_1, t_0)

    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let prop = Propagator::default(OrbitalDynamics::two_body());

    let eccs = vec![1e-5, 0.2];

    // First test is in mostly linear regime (low eccentricity)
    // The second is quite non linear because the eccentricity is 0.2 and we start at periapsis
    for ecc in eccs {
        let init = Orbit::keplerian(8000.0, ecc, 10.0, 5.0, 25.0, 0.0, epoch, eme2k).with_stm();

        let t100 = prop
            .with(init)
            .for_duration(100 * TimeUnit::Second)
            .unwrap();

        let phi_t100_t0 = t100.stm();

        let t50 = prop.with(init).for_duration(50 * TimeUnit::Second).unwrap();

        let phi_t50_t0 = t50.stm();

        let t50_to_t100 = prop.with(t50).for_duration(50 * TimeUnit::Second).unwrap();

        let phi_t100_t50 = t50_to_t100.stm();

        let phi_t100_t50_prime = phi_t100_t0 * phi_t50_t0.try_inverse().unwrap();

        println!("{}{}", phi_t50_t0, phi_t100_t50);
        let delta = phi_t100_t50_prime - phi_t100_t50;

        // NOTE: We don't check the state partials wrt velocity because the upper right is identity and therefore depends on the step size.
        let delta_pos = delta.fixed_columns::<3>(0);
        println!("{}", delta_pos.norm());
        assert!(
            delta_pos.norm() < 1.0,
            "state partials wrt position are too high"
        )
    }
}

#[test]
fn stm_hifi_variable_step() {
    // Using higher fidelity dynamics for STM testing
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let prop = Propagator::default(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun],
        cosm.clone(),
    ));

    let eccs = vec![1e-5, 0.2];

    // First test is in mostly linear regime (low eccentricity)
    // The second is quite non linear because the eccentricity is 0.2 and we start at periapsis
    for ecc in eccs {
        let init = Orbit::keplerian(8000.0, ecc, 10.0, 5.0, 25.0, 0.0, epoch, eme2k);

        let ten_steps = prop
            .with(init.with_stm())
            .for_duration(100 * TimeUnit::Second)
            .unwrap();

        let nominal = ten_steps.to_cartesian_vec();

        let stm_err_ten_steps = ten_steps.stm() * init.to_cartesian_vec() - nominal;

        println!("HD ten_steps.stm() = {}", ten_steps.stm());

        println!("HD stm_err_ten_steps = {}", stm_err_ten_steps);

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

        println!("FD stm_fd = {}", stm_fd);

        let stm_fd_err_ten_steps = stm_fd * init.to_cartesian_vec() - nominal;

        println!("FD stm_fd_err_ten_steps = {}", stm_fd_err_ten_steps);

        let delta = stm_err_ten_steps - stm_fd_err_ten_steps;

        println!("Error between hyperdual and finite diff: {}", delta);

        // Check that each component is less than 1 km or 1 km/s of difference, OR that the hyperdual computation is better than the finite differencing
        for i in 0..6 {
            if delta[i].abs() > 1.0 {
                assert!(delta[i] < 0.0, "FD more precise than HD?! {}", delta);
            }
        }
    }
}

#[test]
fn orbit_set_unset() {
    let m0 = Matrix2::new(1.0, 2.0, 3.0, 4.0);

    let mut m1 = Matrix2::zeros();
    m1.copy_from_slice(m0.as_slice());

    assert_eq!(m1, m0);

    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let init = Orbit::keplerian(8000.0, 0.5, 10.0, 5.0, 25.0, 0.0, epoch, eme2k).with_stm();

    let init_vec = init.as_vector().unwrap();

    let mut init2 = init;
    init2.set(epoch, &init_vec).unwrap();

    assert_eq!(init, init2);
}
