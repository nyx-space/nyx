extern crate nyx_space as nyx;
use std::sync::Arc;

use hifitime::JD_J2000;
use nyx::cosmic::{assert_orbit_eq_or_abs, Orbit};
use nyx::dynamics::orbital::OrbitalDynamics;
use nyx::dynamics::SpacecraftDynamics;
use nyx::propagators::error_ctrl::RSSCartesianState;
use nyx::time::{Epoch, Unit};
use nyx::utils::rss_orbit_errors;
use nyx::{propagators::*, Spacecraft};

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;

use crate::propagation::GMAT_EARTH_GM;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[allow(clippy::identity_op)]
#[rstest]
fn regress_leo_day_adaptive(almanac: Arc<Almanac>) {
    // Regression test for propagators not available in GMAT.
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let prop_time = 1 * Unit::Day;
    let accuracy = 1e-12;
    let min_step = 0.1 * Unit::Second;
    let max_step = 30.0 * Unit::Second;
    let dt = Epoch::from_mjd_tai(JD_J2000);
    let init = Spacecraft::from(Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    ));
    let final_dt = dt + prop_time;

    let all_rslts = vec![
        Orbit::cartesian(
            -5_971.198_524_908_157,
            3_945.775_509_326_305_4,
            2_864.262_542_023_422,
            0.048_766_212_879_869_19,
            -4.184_873_956_982_518,
            5.849_098_380_963_502,
            final_dt,
            eme2k,
        ),
        Orbit::cartesian(
            -5_971.194_190_197_366,
            3_945.506_606_221_459_6,
            2_864.636_682_800_498_4,
            0.049_097_015_227_526_38,
            -4.185_093_356_859_808,
            5.848_940_840_578_1,
            final_dt,
            eme2k,
        ),
        Orbit::cartesian(
            -5_971.194_190_305_766,
            3_945.506_612_356_549_3,
            2_864.636_674_277_756_4,
            0.049_097_007_640_393_29,
            -4.185_093_351_832_897,
            5.848_940_844_198_66,
            final_dt,
            eme2k,
        ),
    ];

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());

    {
        let setup = Propagator::new::<RK2Fixed>(
            dynamics.clone(),
            PropOpts::with_fixed_step(1.0 * Unit::Second),
        );
        let mut prop = setup.with(init, almanac.clone());
        prop.for_duration(prop_time).unwrap();
        assert_eq!(
            prop.state.orbit.to_cartesian_pos_vel(),
            all_rslts[0].to_cartesian_pos_vel(),
            "RK2Fixed two body prop failed"
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
        let mut prop = setup.with(init, almanac.clone());
        prop.for_duration(prop_time).unwrap();
        assert_orbit_eq_or_abs(
            &prop.state.orbit,
            &all_rslts[1],
            1e-7,
            "CashKarp45 two body prop failed",
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
        let setup = Propagator::new::<Fehlberg45>(
            dynamics,
            PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSCartesianState {}),
        );
        let mut prop = setup.with(init, almanac.clone());
        prop.for_duration(prop_time).unwrap();
        // TODO(ANISE): This was a rel check!
        assert_orbit_eq_or_abs(
            &prop.state.orbit,
            &all_rslts[2],
            1e-7,
            "Fehlberg45 two body prop failed",
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
}

#[allow(clippy::identity_op)]
#[rstest]
fn gmat_val_leo_day_adaptive(almanac: Arc<Almanac>) {
    // NOTE: In this test we only use the propagators which also exist in GMAT.
    // Refer to `regress_leo_day_adaptive` for the additional propagators.

    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let prop_time = 1 * Unit::Day;
    let accuracy = 1e-12;
    let min_step = 0.1 * Unit::Second;
    let max_step = 30.0 * Unit::Second;
    let dt = Epoch::from_mjd_tai(JD_J2000);
    let init = Spacecraft::from(Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    ));
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

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());

    {
        let setup = Propagator::new::<Dormand45>(
            dynamics.clone(),
            PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSCartesianState {}),
        );
        let mut prop = setup.with(init, almanac.clone());
        prop.for_duration(prop_time).unwrap();
        assert_orbit_eq_or_abs(
            &prop.state.orbit,
            &all_rslts[0],
            1e-8,
            "two body prop failed",
        );

        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert_eq!(
                prev_details.step, min_step,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
        assert_orbit_eq_or_abs(
            &prop.state.orbit,
            &all_rslts[0],
            1e-8,
            "==> Dormand45: two body prop failed",
        );

        println!("==> Dormand45 adaptive");
        let delta = prop.state.orbit.to_cartesian_pos_vel() - all_rslts[0].to_cartesian_pos_vel();
        for i in 0..6 {
            print!("{:.0e}\t", delta[i].abs());
        }
        println!();

        prop.for_duration(-prop_time).unwrap();
        prop.for_duration(prop_time).unwrap();
        prop.for_duration(prop_time).unwrap();
        prop.for_duration(-prop_time).unwrap();
        let (err_r, err_v) = rss_orbit_errors(&prop.state.orbit, &all_rslts[0]);
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
        let mut prop = setup.with(init, almanac.clone());
        prop.for_duration(prop_time).unwrap();
        assert_orbit_eq_or_abs(
            &prop.state.orbit,
            &all_rslts[1],
            1e-7,
            "==> Verner56: two body prop failed",
        );
        println!("==> Verner56 adaptive");
        let delta = prop.state.orbit.to_cartesian_pos_vel() - all_rslts[1].to_cartesian_pos_vel();
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
        let mut prop = setup.with(init, almanac.clone());
        prop.for_duration(prop_time).unwrap();
        assert_eq!(
            prop.state.orbit.to_cartesian_pos_vel(),
            all_rslts[2].to_cartesian_pos_vel(),
            "==> Dormand78: two body prop failed"
        );
        println!("==> Dormand78 adaptive");
        let delta = prop.state.orbit.to_cartesian_pos_vel() - all_rslts[2].to_cartesian_pos_vel();
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
        let mut prop = setup.with(init, almanac.clone());
        prop.for_duration(prop_time).unwrap();
        assert_eq!(
            prop.state.orbit.to_cartesian_pos_vel(),
            all_rslts[3].to_cartesian_pos_vel(),
            "==> RK89 adaptive: two body prop failed"
        );
        println!("==> RK89 adaptive");
        let delta = prop.state.orbit.to_cartesian_pos_vel() - all_rslts[3].to_cartesian_pos_vel();
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
#[rstest]
fn gmat_val_leo_day_fixed(almanac: Arc<Almanac>) {
    let eme2k = almanac
        .frame_from_uid(EARTH_J2000)
        .unwrap()
        .with_mu_km3_s2(GMAT_EARTH_GM);

    let prop_time = 1 * Unit::Day;
    let dt = Epoch::from_mjd_tai(JD_J2000);
    let init = Spacecraft::from(Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    ));
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

    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::two_body());

    {
        let setup = Propagator::new::<RK4Fixed>(
            dynamics.clone(),
            PropOpts::with_fixed_step(1.0 * Unit::Second),
        );
        let mut prop = setup.with(init, almanac.clone());
        prop.for_duration(prop_time).unwrap();
        assert_eq!(
            prop.state.orbit, all_rslts[0],
            "first forward two body prop failed"
        );
        prop.for_duration(-prop_time).unwrap();
        prop.for_duration(prop_time).unwrap();
        prop.for_duration(prop_time).unwrap();
        prop.for_duration(-prop_time).unwrap();
        let (err_r, err_v) = rss_orbit_errors(&prop.state.orbit, &all_rslts[0]);
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
            PropOpts::with_fixed_step(10.0 * Unit::Second),
        );
        let mut prop = setup.with(init, almanac.clone());
        prop.for_duration(prop_time).unwrap();
        // TODO(ANISE): This was a rel check!
        assert_orbit_eq_or_abs(
            &prop.state.orbit,
            &all_rslts[1],
            1e-7,
            "two body prop failed",
        );
        println!("==> Verner56");
        let delta = prop.state.orbit.to_cartesian_pos_vel() - all_rslts[1].to_cartesian_pos_vel();
        for i in 0..6 {
            print!("{:.0e}\t", delta[i].abs());
        }
        println!();
    }

    {
        let setup = Propagator::new::<Dormand45>(
            dynamics.clone(),
            PropOpts::with_fixed_step(10.0 * Unit::Second),
        );
        let mut prop = setup.with(init, almanac.clone());
        prop.for_duration(prop_time).unwrap();
        assert_eq!(prop.state.orbit, all_rslts[2], "two body prop failed");
        println!("==> Dormand45");
        let delta = prop.state.orbit.to_cartesian_pos_vel() - all_rslts[2].to_cartesian_pos_vel();
        for i in 0..6 {
            print!("{:.0e}\t", delta[i].abs());
        }
        println!();
    }

    {
        let setup = Propagator::new::<Dormand78>(
            dynamics.clone(),
            PropOpts::with_fixed_step(10.0 * Unit::Second),
        );
        let mut prop = setup.with(init, almanac.clone());
        prop.for_duration(prop_time).unwrap();
        assert_eq!(prop.state.orbit, all_rslts[3], "two body prop failed");
        println!("==> Dormand78");
        let delta = prop.state.orbit.to_cartesian_pos_vel() - all_rslts[3].to_cartesian_pos_vel();
        for i in 0..6 {
            print!("{:.0e}\t", delta[i].abs());
        }
        println!();
    }

    {
        let setup =
            Propagator::new::<RK89>(dynamics, PropOpts::with_fixed_step(10.0 * Unit::Second));
        let mut prop = setup.with(init, almanac.clone());
        prop.for_duration(prop_time).unwrap();
        assert_eq!(prop.state.orbit, all_rslts[4], "two body prop failed");
        println!("==> RK89");
        let delta = prop.state.orbit.to_cartesian_pos_vel() - all_rslts[4].to_cartesian_pos_vel();
        for i in 0..6 {
            print!("{:.0e}\t", delta[i].abs());
        }
        println!();
    }
}
