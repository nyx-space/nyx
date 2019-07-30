extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;
use std::f64;

#[test]
fn regress_leo_day_adaptive() {
    // Regression test for propagators not available in GMAT.
    use self::na::Vector6;
    use hifitime::{Epoch, J2000_OFFSET};
    use nyx::celestia::{Cosm, Geoid, State};
    use nyx::dynamics::celestial::CelestialDynamics;
    use nyx::propagators::error_ctrl::RSSStatePV;
    use nyx::propagators::*;
    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(3).unwrap();

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 30.0;
    let dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let init = State::<Geoid>::from_cartesian(-2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, earth_geoid);

    let all_rslts = vec![
        Vector6::from_row_slice(&[
            -5_971.198_709_133_600_5,
            3_945.786_767_659_806_6,
            2_864.246_881_515_823,
            0.048_752_357_390_149_66,
            -4.184_864_764_063_978,
            5.849_104_974_563_176_5,
        ]),
        Vector6::from_row_slice(&[
            -5_971.194_375_364_978,
            3_945.517_869_775_942,
            2_864.621_016_241_891,
            0.049_083_153_975_533_804,
            -4.185_084_160_750_795,
            5.848_947_437_814_404,
        ]),
        Vector6::from_row_slice(&[
            -5_971.194_375_418_999,
            3_945.517_871_298_25,
            2_864.621_014_165_619,
            0.049_083_152_114_524_85,
            -4.185_084_159_507_549,
            5.848_947_438_688_040_5,
        ]),
    ];

    {
        let mut dynamics = CelestialDynamics::two_body(init);
        let mut prop = Propagator::new::<RK2Fixed>(&mut dynamics, &PropOpts::with_fixed_step(1.0, RSSStatePV {}));
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[0], "two body prop failed");
        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert!(
                prev_details.step - min_step < f64::EPSILON,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
    }

    {
        let mut dynamics = CelestialDynamics::two_body(init);
        let mut prop = Propagator::new::<CashKarp45>(
            &mut dynamics,
            &PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSStatePV {}),
        );
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[1], "two body prop failed");
        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert!(
                prev_details.step - min_step < f64::EPSILON,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
    }

    {
        let mut dynamics = CelestialDynamics::two_body(init);
        let mut prop = Propagator::new::<Fehlberg45>(
            &mut dynamics,
            &PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSStatePV {}),
        );
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[2], "two body prop failed");
        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert!(
                prev_details.step - min_step < f64::EPSILON,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
    }
}

#[test]
fn gmat_val_leo_day_adaptive() {
    // NOTE: In this test we only use the propagators which also exist in GMAT.
    // Refer to `regress_leo_day_adaptive` for the additional propagators.

    use self::na::Vector6;
    use hifitime::{Epoch, J2000_OFFSET};
    use nyx::celestia::{Cosm, Geoid, State};
    use nyx::dynamics::celestial::CelestialDynamics;
    use nyx::propagators::error_ctrl::RSSStatePV;
    use nyx::propagators::*;

    let cosm = Cosm::from_xb("./de438s");
    let mut earth_geoid = cosm.geoid_from_id(3).unwrap();
    earth_geoid.gm = 398_600.441_5; // Using GMAT's value

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 30.0;
    let dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let init = State::<Geoid>::from_cartesian(-2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, earth_geoid);

    let all_rslts = vec![
        Vector6::from_row_slice(&[
            -5_971.194_191_972_316,
            3_945.506_662_039_482_3,
            2_864.636_606_375_189,
            0.049_096_946_846_225_495,
            -4.185_093_311_278_742,
            5.848_940_872_821_119,
        ]),
        Vector6::from_row_slice(&[
            -5_971.194_191_678_94,
            3_945.506_653_872_052_5,
            2_864.636_617_510_347_6,
            0.049_096_956_828_390_714,
            -4.185_093_317_946_650_5,
            5.848_940_868_134_205,
        ]),
        Vector6::from_row_slice(&[
            -5_971.194_191_670_392,
            3_945.506_653_218_658,
            2_864.636_618_422_25,
            0.049_096_957_637_897_856,
            -4.185_093_318_481_106,
            5.848_940_867_745_3,
        ]),
        Vector6::from_row_slice(&[
            -5_971.194_191_670_676,
            3_945.506_653_225_158,
            2_864.636_618_413_444_5,
            0.049_096_957_629_993_46,
            -4.185_093_318_475_795,
            5.848_940_867_748_944,
        ]),
    ];

    {
        let mut dynamics = CelestialDynamics::two_body(init);
        let mut prop = Propagator::new::<Dormand45>(
            &mut dynamics,
            &PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSStatePV {}),
        );
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[0], "two body prop failed");
        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert!(
                prev_details.step - min_step < f64::EPSILON,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
    }

    {
        let mut dynamics = CelestialDynamics::two_body(init);
        let mut prop = Propagator::new::<Verner56>(
            &mut dynamics,
            &PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSStatePV {}),
        );
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[1], "two body prop failed");
        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert!(
                prev_details.step - min_step < f64::EPSILON,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
    }

    {
        let mut dynamics = CelestialDynamics::two_body(init);
        let mut prop = Propagator::new::<Dormand78>(
            &mut dynamics,
            &PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSStatePV {}),
        );
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[2], "two body prop failed");
        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert!(
                prev_details.step - min_step < f64::EPSILON,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
    }

    {
        let mut dynamics = CelestialDynamics::two_body(init);
        let mut prop = Propagator::new::<RK89>(
            &mut dynamics,
            &PropOpts::with_adaptive_step(min_step, max_step, accuracy, RSSStatePV {}),
        );
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[3], "two body prop failed");
        let prev_details = prop.latest_details();
        if prev_details.error > accuracy {
            assert!(
                prev_details.step - min_step < f64::EPSILON,
                "step size should be at its minimum because error is higher than tolerance: {:?}",
                prev_details
            );
        }
    }
}

#[test]
fn gmat_val_leo_day_fixed() {
    use crate::na::Vector6;
    use hifitime::{Epoch, J2000_OFFSET};
    use nyx::celestia::{Cosm, Geoid, State};
    use nyx::dynamics::celestial::CelestialDynamics;
    use nyx::propagators::error_ctrl::RSSStatePV;
    use nyx::propagators::*;

    let cosm = Cosm::from_xb("./de438s");
    let mut earth_geoid = cosm.geoid_from_id(3).unwrap();
    earth_geoid.gm = 398_600.441_5; // Using GMAT's value

    let prop_time = 3_600.0 * 24.0;
    let dt = Epoch::from_mjd_tai(J2000_OFFSET);
    let init = State::<Geoid>::from_cartesian(-2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, earth_geoid);

    let all_rslts = vec![
        Vector6::from_row_slice(&[
            -5_971.194_191_670_768,
            3_945.506_653_227_154,
            2_864.636_618_410_970_6,
            0.049_096_957_627_641_77,
            -4.185_093_318_474_28,
            5.848_940_867_750_096_5,
        ]),
        Vector6::from_row_slice(&[
            -5_971.194_191_670_203,
            3_945.506_653_219_096_7,
            2_864.636_618_421_618,
            0.049_096_957_637_339_07,
            -4.185_093_318_480_867,
            5.848_940_867_745_654,
        ]),
        Vector6::from_row_slice(&[
            -5_971.194_191_699_656,
            3_945.506_654_080_17,
            2_864.636_617_245_45,
            0.049_096_956_584_062_28,
            -4.185_093_317_777_894,
            5.848_940_868_241_106,
        ]),
        Vector6::from_row_slice(&[
            -5_971.194_191_670_044,
            3_945.506_653_211_795_3,
            2_864.636_618_431_374,
            0.049_096_957_645_996_114,
            -4.185_093_318_486_724,
            5.848_940_867_741_533,
        ]),
        Vector6::from_row_slice(&[
            -5_971.194_191_670_81,
            3_945.506_653_233_250_3,
            2_864.636_618_402_241_8,
            0.049_096_957_620_019_005,
            -4.185_093_318_469_214,
            5.848_940_867_753_748,
        ]),
    ];

    {
        let mut dynamics = CelestialDynamics::two_body(init);
        let mut prop = Propagator::new::<RK4Fixed>(&mut dynamics, &PropOpts::with_fixed_step(1.0, RSSStatePV {}));
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[0], "two body prop failed");
    }

    {
        let mut dynamics = CelestialDynamics::two_body(init);
        let mut prop = Propagator::new::<Verner56>(&mut dynamics, &PropOpts::with_fixed_step(10.0, RSSStatePV {}));
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[1], "two body prop failed");
    }

    {
        let mut dynamics = CelestialDynamics::two_body(init);
        let mut prop = Propagator::new::<Dormand45>(&mut dynamics, &PropOpts::with_fixed_step(10.0, RSSStatePV {}));
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[2], "two body prop failed");
    }

    {
        let mut dynamics = CelestialDynamics::two_body(init);
        let mut prop = Propagator::new::<Dormand78>(&mut dynamics, &PropOpts::with_fixed_step(10.0, RSSStatePV {}));
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[3], "two body prop failed");
    }

    {
        let mut dynamics = CelestialDynamics::two_body(init);
        let mut prop = Propagator::new::<RK89>(&mut dynamics, &PropOpts::with_fixed_step(10.0, RSSStatePV {}));
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[4], "two body prop failed");
    }
}
