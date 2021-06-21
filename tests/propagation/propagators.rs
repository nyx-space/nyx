extern crate nyx_space as nyx;
use nyx::cosmic::{assert_orbit_eq_or_abs, assert_orbit_eq_or_rel, Cosm, Orbit};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::propagators::error_ctrl::RSSCartesianState;
use nyx::propagators::*;
use nyx::time::{Epoch, TimeUnit, J2000_OFFSET};
use nyx::utils::rss_orbit_errors;

#[allow(clippy::identity_op)]
#[test]
fn regress_leo_day_adaptive() {
    // Regression test for propagators not available in GMAT.
    let cosm = Cosm::de438();
    let eme2k = cosm.frame("EME2000");

    let prop_time = 1 * TimeUnit::Day;
    let accuracy = 1e-12;
    let min_step = 0.1 * TimeUnit::Second;
    let max_step = 30.0 * TimeUnit::Second;
    let dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let init = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );
    let final_dt = dt + prop_time;

    let all_rslts = vec![
        Orbit::cartesian(
            -5_971.198_709_133_600_5,
            3_945.786_767_659_806_6,
            2_864.246_881_515_823,
            0.048_752_357_390_149_66,
            -4.184_864_764_063_978,
            5.849_104_974_563_176_5,
            final_dt,
            eme2k,
        ),
        Orbit::cartesian(
            -5_971.194_375_364_978,
            3_945.517_869_775_919_7,
            2_864.621_016_241_924,
            0.049_083_153_975_562_65,
            -4.185_084_160_750_815,
            5.848_947_437_814_39,
            final_dt,
            eme2k,
        ),
        Orbit::cartesian(
            -5_971.194_375_418_999,
            3_945.517_871_298_253_3,
            2_864.621_014_165_613_4,
            0.049_083_152_114_520_266,
            -4.185_084_159_507_545,
            5.848_947_438_688_043,
            final_dt,
            eme2k,
        ),
    ];

    let dynamics = OrbitalDynamics::two_body();

    {
        let setup = Propagator::new::<RK2Fixed>(
            dynamics.clone(),
            PropOpts::with_fixed_step(1.0 * TimeUnit::Second),
        );
        let mut prop = setup.with(init);
        prop.for_duration(prop_time).unwrap();
        assert_eq!(
            prop.state.to_cartesian_vec(),
            all_rslts[0].to_cartesian_vec(),
            "two body prop failed"
        );
        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert_eq!(
                prev_details.step, min_step,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
    }

    {
        let setup = Propagator::new::<CashKarp45>(
            dynamics.clone(),
            PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSCartesianState {}),
        );
        let mut prop = setup.with(init);
        prop.for_duration(prop_time).unwrap();
        assert_orbit_eq_or_abs(&prop.state, &all_rslts[1], 1e-7, "two body prop failed");
        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert_eq!(
                prev_details.step, min_step,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
    }

    {
        let setup = Propagator::new::<Fehlberg45>(
            dynamics,
            PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSCartesianState {}),
        );
        let mut prop = setup.with(init);
        prop.for_duration(prop_time).unwrap();
        assert_orbit_eq_or_rel(&prop.state, &all_rslts[2], 1e-7, "two body prop failed");
        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert_eq!(
                prev_details.step, min_step,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
    }
}

#[allow(clippy::identity_op)]
#[test]
fn gmat_val_leo_day_adaptive() {
    // NOTE: In this test we only use the propagators which also exist in GMAT.
    // Refer to `regress_leo_day_adaptive` for the additional propagators.

    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    let prop_time = 1 * TimeUnit::Day;
    let accuracy = 1e-12;
    let min_step = 0.1 * TimeUnit::Second;
    let max_step = 30.0 * TimeUnit::Second;
    let dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let init = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );
    let final_dt = dt + prop_time;

    let all_rslts = vec![
        Orbit::cartesian(
            -5_971.194_191_972_314,
            3_945.506_662_039_457,
            2_864.636_606_375_225_7,
            0.049_096_946_846_257_56,
            -4.185_093_311_278_763,
            5.848_940_872_821_106,
            final_dt,
            eme2k,
        ),
        Orbit::cartesian(
            -5_971.194_191_678_94,
            3_945.506_653_872_037_5,
            2_864.636_617_510_367,
            0.049_096_956_828_408_46,
            -4.185_093_317_946_663,
            5.848_940_868_134_195_4,
            final_dt,
            eme2k,
        ),
        Orbit::cartesian(
            -5_971.194_191_670_392,
            3_945.506_653_218_658,
            2_864.636_618_422_25,
            0.049_096_957_637_897_856,
            -4.185_093_318_481_106,
            5.848_940_867_745_3,
            final_dt,
            eme2k,
        ),
        Orbit::cartesian(
            -5_971.194_191_670_676,
            3_945.506_653_225_158,
            2_864.636_618_413_444_5,
            0.049_096_957_629_993_46,
            -4.185_093_318_475_795,
            5.848_940_867_748_944,
            final_dt,
            eme2k,
        ),
    ];

    let dynamics = OrbitalDynamics::two_body();

    {
        let setup = Propagator::new::<Dormand45>(
            dynamics.clone(),
            PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSCartesianState {}),
        );
        let mut prop = setup.with(init);
        prop.for_duration(prop_time).unwrap();
        assert_orbit_eq_or_abs(&prop.state, &all_rslts[0], 1e-8, "two body prop failed");

        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert_eq!(
                prev_details.step, min_step,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
        assert_orbit_eq_or_abs(
            &prop.state,
            &all_rslts[0],
            1e-8,
            "first forward two body prop failed",
        );

        println!("==> Dormand45 adaptive");
        let delta = prop.state.to_cartesian_vec() - all_rslts[0].to_cartesian_vec();
        for i in 0..6 {
            print!("{:.0e}\t", delta[i].abs());
        }
        println!();

        prop.for_duration(-prop_time).unwrap();
        prop.for_duration(prop_time).unwrap();
        prop.for_duration(prop_time).unwrap();
        prop.for_duration(-prop_time).unwrap();
        let (err_r, err_v) = rss_orbit_errors(&prop.state, &all_rslts[0]);
        assert!(
            err_r < 1e-5,
            "two body 2*(fwd+back) prop failed to return to the initial state in position"
        );
        assert!(
            err_v < 1e-8,
            "two body 2*(fwd+back) prop failed to return to the initial state in velocity"
        );
    }

    {
        let setup = Propagator::new::<Verner56>(
            dynamics.clone(),
            PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSCartesianState {}),
        );
        let mut prop = setup.with(init);
        prop.for_duration(prop_time).unwrap();
        assert_orbit_eq_or_abs(&prop.state, &all_rslts[1], 1e-8, "two body prop failed");
        println!("==> Verner56 adaptive");
        let delta = prop.state.to_cartesian_vec() - all_rslts[1].to_cartesian_vec();
        for i in 0..6 {
            print!("{:.0e}\t", delta[i].abs());
        }
        println!();
        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert_eq!(
                prev_details.step, min_step,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
    }

    {
        let setup = Propagator::new::<Dormand78>(
            dynamics.clone(),
            PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSCartesianState {}),
        );
        let mut prop = setup.with(init);
        prop.for_duration(prop_time).unwrap();
        assert_eq!(
            prop.state.to_cartesian_vec(),
            all_rslts[2].to_cartesian_vec(),
            "two body prop failed"
        );
        println!("==> Dormand78 adaptive");
        let delta = prop.state.to_cartesian_vec() - all_rslts[2].to_cartesian_vec();
        for i in 0..6 {
            print!("{:.0e}\t", delta[i].abs());
        }
        println!();
        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert_eq!(
                prev_details.step, min_step,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
    }

    {
        let setup = Propagator::new::<RK89>(
            dynamics,
            PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSCartesianState {}),
        );
        let mut prop = setup.with(init);
        prop.for_duration(prop_time).unwrap();
        assert_eq!(
            prop.state.to_cartesian_vec(),
            all_rslts[3].to_cartesian_vec(),
            "two body prop failed"
        );
        println!("==> RK89 adaptive");
        let delta = prop.state.to_cartesian_vec() - all_rslts[3].to_cartesian_vec();
        for i in 0..6 {
            print!("{:.0e}\t", delta[i].abs());
        }
        println!();
        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert_eq!(
                prev_details.step, min_step,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
    }
}

#[allow(clippy::identity_op)]
#[test]
fn gmat_val_leo_day_fixed() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    let prop_time = 1 * TimeUnit::Day;
    let dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let init = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );
    let final_dt = dt + prop_time;

    let all_rslts = vec![
        Orbit::cartesian(
            -5_971.194_191_670_768,
            3_945.506_653_227_154,
            2_864.636_618_410_970_6,
            0.049_096_957_627_641_77,
            -4.185_093_318_474_28,
            5.848_940_867_750_096_5,
            final_dt,
            eme2k,
        ),
        Orbit::cartesian(
            -5_971.194_191_670_203,
            3_945.506_653_219_096_7,
            2_864.636_618_421_618,
            0.049_096_957_637_339_07,
            -4.185_093_318_480_867,
            5.848_940_867_745_654,
            final_dt,
            eme2k,
        ),
        Orbit::cartesian(
            -5_971.194_191_699_656,
            3_945.506_654_080_17,
            2_864.636_617_245_45,
            0.049_096_956_584_062_28,
            -4.185_093_317_777_894,
            5.848_940_868_241_106,
            final_dt,
            eme2k,
        ),
        Orbit::cartesian(
            -5_971.194_191_670_044,
            3_945.506_653_211_795_3,
            2_864.636_618_431_374,
            0.049_096_957_645_996_114,
            -4.185_093_318_486_724,
            5.848_940_867_741_533,
            final_dt,
            eme2k,
        ),
        Orbit::cartesian(
            -5_971.194_191_670_81,
            3_945.506_653_233_250_3,
            2_864.636_618_402_241_8,
            0.049_096_957_620_019_005,
            -4.185_093_318_469_214,
            5.848_940_867_753_748,
            final_dt,
            eme2k,
        ),
    ];

    let dynamics = OrbitalDynamics::two_body();

    {
        let setup = Propagator::new::<RK4Fixed>(
            dynamics.clone(),
            PropOpts::with_fixed_step(1.0 * TimeUnit::Second),
        );
        let mut prop = setup.with(init);
        prop.for_duration(prop_time).unwrap();
        assert_eq!(
            prop.state, all_rslts[0],
            "first forward two body prop failed"
        );
        prop.for_duration(-prop_time).unwrap();
        prop.for_duration(prop_time).unwrap();
        prop.for_duration(prop_time).unwrap();
        prop.for_duration(-prop_time).unwrap();
        let (err_r, err_v) = rss_orbit_errors(&prop.state, &all_rslts[0]);
        assert!(
            err_r < 1e-5,
            "two body 2*(fwd+back) prop failed to return to the initial state in position"
        );
        assert!(
            err_v < 1e-8,
            "two body 2*(fwd+back) prop failed to return to the initial state in velocity"
        );
    }

    {
        let setup = Propagator::new::<Verner56>(
            dynamics.clone(),
            PropOpts::with_fixed_step(10.0 * TimeUnit::Second),
        );
        let mut prop = setup.with(init);
        prop.for_duration(prop_time).unwrap();
        assert_orbit_eq_or_rel(&prop.state, &all_rslts[1], 1e-7, "two body prop failed");
        println!("==> Verner56");
        let delta = prop.state.to_cartesian_vec() - all_rslts[1].to_cartesian_vec();
        for i in 0..6 {
            print!("{:.0e}\t", delta[i].abs());
        }
        println!();
    }

    {
        let setup = Propagator::new::<Dormand45>(
            dynamics.clone(),
            PropOpts::with_fixed_step(10.0 * TimeUnit::Second),
        );
        let mut prop = setup.with(init);
        prop.for_duration(prop_time).unwrap();
        assert_eq!(prop.state, all_rslts[2], "two body prop failed");
        println!("==> Dormand45");
        let delta = prop.state.to_cartesian_vec() - all_rslts[2].to_cartesian_vec();
        for i in 0..6 {
            print!("{:.0e}\t", delta[i].abs());
        }
        println!();
    }

    {
        let setup = Propagator::new::<Dormand78>(
            dynamics.clone(),
            PropOpts::with_fixed_step(10.0 * TimeUnit::Second),
        );
        let mut prop = setup.with(init);
        prop.for_duration(prop_time).unwrap();
        assert_eq!(prop.state, all_rslts[3], "two body prop failed");
        println!("==> Dormand78");
        let delta = prop.state.to_cartesian_vec() - all_rslts[3].to_cartesian_vec();
        for i in 0..6 {
            print!("{:.0e}\t", delta[i].abs());
        }
        println!();
    }

    {
        let setup =
            Propagator::new::<RK89>(dynamics, PropOpts::with_fixed_step(10.0 * TimeUnit::Second));
        let mut prop = setup.with(init);
        prop.for_duration(prop_time).unwrap();
        assert_eq!(prop.state, all_rslts[4], "two body prop failed");
        println!("==> RK89");
        let delta = prop.state.to_cartesian_vec() - all_rslts[4].to_cartesian_vec();
        for i in 0..6 {
            print!("{:.0e}\t", delta[i].abs());
        }
        println!();
    }
}
