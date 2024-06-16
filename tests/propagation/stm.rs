extern crate nyx_space as nyx;
use std::sync::Arc;

use anise::constants::celestial_objects::{MOON, SUN};
use nyx::cosmic::{Orbit, Spacecraft};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::linalg::{Const, Matrix6, OVector};
use nyx::propagators::*;
use nyx::time::{Epoch, Unit};
use nyx::State;
use nyx_space::md::prelude::SpacecraftDynamics;

// These tests compare the computation of the state transition matrix between the finite differencing methoid (common) and hyperdual numbers.
// Conclusion: hyperdual numbers lead to less error than finite differencing.

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;

use crate::propagation::GMAT_EARTH_GM;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn stm_fixed_step(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let prop = Propagator::new::<RK4Fixed>(
        SpacecraftDynamics::new(OrbitalDynamics::two_body()),
        PropOpts::with_fixed_step(10 * Unit::Second),
    );

    // First test is in mostly linear regime (low eccentricity)
    // The second is quite non linear because the eccentricity is 0.2 and we start at periapsis
    let eccs = vec![1e-5, 0.2];

    for ecc in eccs {
        let init = Spacecraft::from(Orbit::keplerian(
            8000.0, ecc, 10.0, 5.0, 25.0, 0.0, epoch, eme2k,
        ));

        let ten_steps = prop
            .with(init.with_stm(), almanac.clone())
            .for_duration(100 * Unit::Second)
            .unwrap();

        let nominal = ten_steps.orbit.to_cartesian_pos_vel();

        let stm_err_ten_steps = ten_steps.stm().unwrap().fixed_resize::<6, 6>(0.0)
            * init.orbit.to_cartesian_pos_vel()
            - nominal;

        println!("HD ten_steps.stm() = {}", ten_steps.stm().unwrap());

        println!("HD stm_err_ten_steps = {}", stm_err_ten_steps);

        // Compare with finite differencing
        let mut stm_fd = Matrix6::<f64>::zeros();
        let pert = 0.0001;

        for i in 0..6 {
            let mut this_init = init.with_stm();
            match i {
                0 => this_init.orbit.radius_km.x += pert,
                1 => this_init.orbit.radius_km.y += pert,
                2 => this_init.orbit.radius_km.z += pert,
                3 => this_init.orbit.velocity_km_s.x += pert,
                4 => this_init.orbit.velocity_km_s.y += pert,
                5 => this_init.orbit.velocity_km_s.z += pert,
                _ => unreachable!(),
            }

            let these_ten_steps = prop
                .with(this_init, almanac.clone())
                .for_duration(100 * Unit::Second)
                .unwrap();

            let jac_val = (these_ten_steps.orbit.to_cartesian_pos_vel() - nominal) / pert;

            for j in 0..6 {
                stm_fd[(j, i)] = jac_val[j];
            }
        }

        println!("FD stm_fd = {}", stm_fd);

        let stm_fd_err_ten_steps = stm_fd * init.orbit.to_cartesian_pos_vel() - nominal;

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

#[rstest]
fn stm_variable_step(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let prop = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::two_body()));

    let eccs = vec![1e-5, 0.2];

    // First test is in mostly linear regime (low eccentricity)
    // The second is quite non linear because the eccentricity is 0.2 and we start at periapsis
    for ecc in eccs {
        let init = Spacecraft::from(Orbit::keplerian(
            8000.0, ecc, 10.0, 5.0, 25.0, 0.0, epoch, eme2k,
        ));

        let ten_steps = prop
            .with(init.with_stm(), almanac.clone())
            .for_duration(100 * Unit::Second)
            .unwrap();

        let nominal = ten_steps.orbit.to_cartesian_pos_vel();

        let stm_err_ten_steps = ten_steps.stm().unwrap().fixed_resize::<6, 6>(0.0)
            * init.orbit.to_cartesian_pos_vel()
            - nominal;

        println!("HD ten_steps.stm() = {}", ten_steps.stm().unwrap());

        println!("HD stm_err_ten_steps = {}", stm_err_ten_steps);

        // Compare with finite differencing
        let mut stm_fd = Matrix6::<f64>::zeros();
        let pert = 0.0001;

        for i in 0..6 {
            let mut this_init = init.with_stm();
            match i {
                0 => this_init.orbit.radius_km.x += pert,
                1 => this_init.orbit.radius_km.y += pert,
                2 => this_init.orbit.radius_km.z += pert,
                3 => this_init.orbit.velocity_km_s.x += pert,
                4 => this_init.orbit.velocity_km_s.y += pert,
                5 => this_init.orbit.velocity_km_s.z += pert,
                _ => unreachable!(),
            }

            let these_ten_steps = prop
                .with(this_init, almanac.clone())
                .for_duration(100 * Unit::Second)
                .unwrap();

            let jac_val = (these_ten_steps.orbit.to_cartesian_pos_vel() - nominal) / pert;

            for j in 0..6 {
                stm_fd[(j, i)] = jac_val[j];
            }
        }

        println!("FD stm_fd = {}", stm_fd);

        let stm_fd_err_ten_steps = stm_fd * init.orbit.to_cartesian_pos_vel() - nominal;

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

#[rstest]
fn stm_between_steps(almanac: Arc<Almanac>) {
    // Check that \Phi(t_2, t_1) = \Phi(t_2, t_0) * \Phi^{-1}(t_1, t_0)

    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let prop = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::two_body()));

    let eccs = vec![1e-5, 0.2];

    // First test is in mostly linear regime (low eccentricity)
    // The second is quite non linear because the eccentricity is 0.2 and we start at periapsis
    for ecc in eccs {
        let init = Spacecraft::from(Orbit::keplerian(
            8000.0, ecc, 10.0, 5.0, 25.0, 0.0, epoch, eme2k,
        ))
        .with_stm();

        let t100 = prop
            .with(init, almanac.clone())
            .for_duration(100 * Unit::Second)
            .unwrap();

        let phi_t100_t0 = t100.stm().unwrap();

        let t50 = prop
            .with(init, almanac.clone())
            .for_duration(50 * Unit::Second)
            .unwrap();

        let phi_t50_t0 = t50.stm().unwrap();

        let t50_to_t100 = prop
            .with(t50, almanac.clone())
            .for_duration(50 * Unit::Second)
            .unwrap();

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

#[rstest]
fn stm_hifi_variable_step(almanac: Arc<Almanac>) {
    // Using higher fidelity dynamics for STM testing

    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let prop = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        vec![MOON, SUN],
    )));

    let eccs = vec![1e-5, 0.2];

    // First test is in mostly linear regime (low eccentricity)
    // The second is quite non linear because the eccentricity is 0.2 and we start at periapsis
    for ecc in eccs {
        let init = Spacecraft::from(Orbit::keplerian(
            8000.0, ecc, 10.0, 5.0, 25.0, 0.0, epoch, eme2k,
        ));

        let ten_steps = prop
            .with(init.with_stm(), almanac.clone())
            .for_duration(100 * Unit::Second)
            .unwrap();

        let nominal = ten_steps.orbit.to_cartesian_pos_vel();

        let stm_err_ten_steps = ten_steps.stm().unwrap().fixed_resize::<6, 6>(0.0)
            * init.orbit.to_cartesian_pos_vel()
            - nominal;

        println!("HD ten_steps.stm() = {}", ten_steps.stm().unwrap());

        println!("HD stm_err_ten_steps = {}", stm_err_ten_steps);

        // Compare with finite differencing
        let mut stm_fd = Matrix6::<f64>::zeros();
        let pert = 0.0001;

        for i in 0..6 {
            let mut this_init = init.with_stm();
            match i {
                0 => this_init.orbit.radius_km.x += pert,
                1 => this_init.orbit.radius_km.y += pert,
                2 => this_init.orbit.radius_km.z += pert,
                3 => this_init.orbit.velocity_km_s.x += pert,
                4 => this_init.orbit.velocity_km_s.y += pert,
                5 => this_init.orbit.velocity_km_s.z += pert,
                _ => unreachable!(),
            }

            let these_ten_steps = prop
                .with(this_init, almanac.clone())
                .for_duration(100 * Unit::Second)
                .unwrap();

            let jac_val = (these_ten_steps.orbit.to_cartesian_pos_vel() - nominal) / pert;

            for j in 0..6 {
                stm_fd[(j, i)] = jac_val[j];
            }
        }

        println!("FD stm_fd = {}", stm_fd);

        let stm_fd_err_ten_steps = stm_fd * init.orbit.to_cartesian_pos_vel() - nominal;

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

#[rstest]
fn orbit_set_unset_static(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let mut init = Spacecraft::from(Orbit::keplerian(
        8000.0, 0.5, 10.0, 5.0, 25.0, 0.0, epoch, eme2k,
    ))
    .with_stm();

    // Change STM
    let stm_data = (0..81).map(|x| x as f64).collect::<Vec<f64>>();
    init.stm.as_mut().unwrap().copy_from_slice(&stm_data);

    let init_vec = init.to_vector();

    let mut init2 = init;
    init2.set(epoch, &init_vec);

    assert_eq!(init, init2);
}

#[rstest]
fn orbit_set_unset(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let init = Spacecraft::from(Orbit::keplerian(
        8000.0, 0.5, 10.0, 5.0, 25.0, 0.0, epoch, eme2k,
    ))
    .with_stm();

    let prop = Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        vec![MOON, SUN],
    )));

    let orbit = prop
        .with(init, almanac.clone())
        .for_duration(2 * Unit::Hour)
        .unwrap();

    let vec = orbit.to_vector();

    let mut init2 = orbit;
    init2.set(orbit.epoch(), &vec);

    assert_eq!(orbit, init2);
}

#[rstest]
fn sc_set_unset_static(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);
    let epoch = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    let init = Orbit::keplerian(8000.0, 0.5, 10.0, 5.0, 25.0, 0.0, epoch, eme2k);
    let mut init_sc = Spacecraft::from_srp_defaults(init, 100.0, 1.0).with_stm();

    // Change the full vector
    let data = (0..90).map(|x| x as f64).collect::<Vec<f64>>();
    init_sc.set(
        init.epoch,
        &OVector::<f64, Const<90>>::from_column_slice(&data),
    );

    let init_vec = init_sc.to_vector();

    let mut init2 = init_sc;
    init2.set(epoch, &init_vec);

    assert_eq!(init_sc, init2);
}
