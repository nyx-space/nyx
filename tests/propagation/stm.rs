extern crate nyx_space as nyx;
use nyx::cosmic::{Bodies, Cosm, Orbit, Spacecraft};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::linalg::{Const, Matrix6, OVector};
use nyx::propagators::*;
use nyx::time::{Epoch, Unit};
use nyx::State;
use nyx_space::md::ui::SpacecraftDynamics;

// These tests compare the computation of the state transition matrix between the finite differencing methoid (common) and hyperdual numbers.
// Conclusion: hyperdual numbers lead to less error than finite differencing.

#[test]
fn stm_fixed_step() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let prop = Propagator::new::<RK4Fixed>(
        OrbitalDynamics::two_body(),
        PropOpts::with_fixed_step(10 * Unit::Second),
    );

    // First test is in mostly linear regime (low eccentricity)
    // The second is quite non linear because the eccentricity is 0.2 and we start at periapsis
    let eccs = vec![1e-5, 0.2];

    for ecc in eccs {
        let init = Orbit::keplerian(8000.0, ecc, 10.0, 5.0, 25.0, 0.0, epoch, eme2k);

        let ten_steps = prop
            .with(init.with_stm())
            .for_duration(100 * Unit::Second)
            .unwrap();

        let nominal = ten_steps.to_cartesian_vec();

        let stm_err_ten_steps = ten_steps.stm().unwrap() * init.to_cartesian_vec() - nominal;

        println!("HD ten_steps.stm() = {}", ten_steps.stm().unwrap());

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
                .for_duration(100 * Unit::Second)
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

    let prop = Propagator::default_dp78(OrbitalDynamics::two_body());

    let eccs = vec![1e-5, 0.2];

    // First test is in mostly linear regime (low eccentricity)
    // The second is quite non linear because the eccentricity is 0.2 and we start at periapsis
    for ecc in eccs {
        let init = Orbit::keplerian(8000.0, ecc, 10.0, 5.0, 25.0, 0.0, epoch, eme2k);

        let ten_steps = prop
            .with(init.with_stm())
            .for_duration(100 * Unit::Second)
            .unwrap();

        let nominal = ten_steps.to_cartesian_vec();

        let stm_err_ten_steps = ten_steps.stm().unwrap() * init.to_cartesian_vec() - nominal;

        println!("HD ten_steps.stm() = {}", ten_steps.stm().unwrap());

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
                .for_duration(100 * Unit::Second)
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

    let prop = Propagator::default_dp78(OrbitalDynamics::two_body());

    let eccs = vec![1e-5, 0.2];

    // First test is in mostly linear regime (low eccentricity)
    // The second is quite non linear because the eccentricity is 0.2 and we start at periapsis
    for ecc in eccs {
        let init = Orbit::keplerian(8000.0, ecc, 10.0, 5.0, 25.0, 0.0, epoch, eme2k).with_stm();

        let t100 = prop.with(init).for_duration(100 * Unit::Second).unwrap();

        let phi_t100_t0 = t100.stm().unwrap();

        let t50 = prop.with(init).for_duration(50 * Unit::Second).unwrap();

        let phi_t50_t0 = t50.stm().unwrap();

        let t50_to_t100 = prop.with(t50).for_duration(50 * Unit::Second).unwrap();

        let phi_t100_t50 = t50_to_t100.stm().unwrap();

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

    let prop = Propagator::default_dp78(OrbitalDynamics::point_masses(
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
            .for_duration(100 * Unit::Second)
            .unwrap();

        let nominal = ten_steps.to_cartesian_vec();

        let stm_err_ten_steps = ten_steps.stm().unwrap() * init.to_cartesian_vec() - nominal;

        println!("HD ten_steps.stm() = {}", ten_steps.stm().unwrap());

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
                .for_duration(100 * Unit::Second)
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
fn orbit_set_unset_static() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let mut init = Orbit::keplerian(8000.0, 0.5, 10.0, 5.0, 25.0, 0.0, epoch, eme2k).with_stm();

    // Change STM
    let stm_data = (0..36).map(|x| x as f64).collect::<Vec<f64>>();
    init.stm.as_mut().unwrap().copy_from_slice(&stm_data);

    let init_vec = init.as_vector().unwrap();

    let mut init2 = init;
    init2.set(epoch, &init_vec).unwrap();

    assert_eq!(init, init2);
}

#[test]
fn orbit_set_unset() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let init = Orbit::keplerian(8000.0, 0.5, 10.0, 5.0, 25.0, 0.0, epoch, eme2k).with_stm();

    let prop = Propagator::default_dp78(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun],
        cosm.clone(),
    ));

    let orbit = prop.with(init).for_duration(2 * Unit::Hour).unwrap();

    let vec = orbit.as_vector().unwrap();

    let mut init2 = orbit;
    init2.set(orbit.epoch(), &vec).unwrap();

    assert_eq!(orbit, init2);
}

#[test]
fn sc_set_unset_static() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let init = Orbit::keplerian(8000.0, 0.5, 10.0, 5.0, 25.0, 0.0, epoch, eme2k).with_stm();
    let mut init_sc = Spacecraft::from_srp_defaults(init, 100.0, 1.0);

    // Change the full vector
    let data = (0..90).map(|x| x as f64).collect::<Vec<f64>>();
    init_sc
        .set(
            init.epoch(),
            &OVector::<f64, Const<90>>::from_column_slice(&data),
        )
        .unwrap();

    let init_vec = init_sc.as_vector().unwrap();

    let mut init2 = init_sc;
    init2.set(epoch, &init_vec).unwrap();

    assert_eq!(init_sc, init2);
}

#[test]
fn sc_and_orbit_stm_chk() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let init_orbit = Orbit::keplerian(8000.0, 0.5, 10.0, 5.0, 25.0, 0.0, epoch, eme2k).with_stm();
    let init_sc = Spacecraft::from_srp_defaults(init_orbit, 0.0, 0.0);

    let prop_orbit = Propagator::default_dp78(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun],
        cosm.clone(),
    ));

    let prop_sc = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        &[Bodies::Luna, Bodies::Sun],
        cosm.clone(),
    )));

    let final_orbit = prop_orbit
        .with(init_orbit)
        .for_duration(2 * Unit::Hour)
        .unwrap();
    let final_sc = prop_sc.with(init_sc).for_duration(2 * Unit::Hour).unwrap();

    assert_eq!(
        final_orbit, final_sc.orbit,
        "Identical dynamics for Spacecraft and Orbit lead to different states"
    );

    let sc_stm = final_sc.stm().unwrap();
    let sc_orbit_stm = sc_stm.fixed_slice::<6, 6>(0, 0);
    assert_eq!(
        sc_orbit_stm,
        final_orbit.stm().unwrap(),
        "Identical dynamics for Spacecraft and Orbit lead to different STM"
    );
}
