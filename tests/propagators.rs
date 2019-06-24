extern crate nalgebra as na;
extern crate nyx_space as nyx;
use std::f64;

#[test]
fn regress_leo_day_adaptive() {
    // Regression test for propagators not available in GMAT.
    use self::na::Vector6;
    use nyx::celestia::Cosm;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::propagators::error_ctrl::RSSStatePV;
    use nyx::propagators::*;
    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(3).unwrap();

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 30.0;
    let init = Vector6::from_row_slice(&[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);

    let all_rslts = vec![
        Vector6::from_row_slice(&[
            -5971.198524908157,
            3945.7755093263054,
            2864.262542023422,
            0.04876621287986919,
            -4.184873956982518,
            5.849098380963502,
        ]),
        Vector6::from_row_slice(&[
            -5971.194190251459,
            3945.5066108335227,
            2864.636676354708,
            0.049097009502150256,
            -4.185093353076805,
            5.848940843324859,
        ]),
        Vector6::from_row_slice(&[
            -5971.194190305514,
            3945.5066123572196,
            2864.6366742764963,
            0.04909700763937082,
            -4.185093351832398,
            5.848940844199317,
        ]),
    ];

    {
        let mut dynamics = TwoBody::from_state_vec(init.clone(), earth_geoid.clone());
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
        let mut dynamics = TwoBody::from_state_vec(init.clone(), earth_geoid.clone());
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
        let mut dynamics = TwoBody::from_state_vec(init.clone(), earth_geoid);
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
    use nyx::celestia::Cosm;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::propagators::error_ctrl::RSSStatePV;
    use nyx::propagators::*;

    let cosm = Cosm::from_xb("./de438s");
    let earth_geoid = cosm.geoid_from_id(3).unwrap();

    let prop_time = 24.0 * 3_600.0;
    let accuracy = 1e-12;
    let min_step = 0.1;
    let max_step = 30.0;
    let init = Vector6::from_row_slice(&[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);

    let all_rslts = vec![
        Vector6::from_row_slice(&[
            -5971.194191972316,
            3945.5066620394823,
            2864.636606375189,
            0.049096946846225495,
            -4.185093311278742,
            5.848940872821119,
        ]),
        Vector6::from_row_slice(&[
            -5971.19419167894,
            3945.5066538720525,
            2864.6366175103476,
            0.049096956828390714,
            -4.1850933179466505,
            5.848940868134205,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191670392,
            3945.506653218658,
            2864.63661842225,
            0.049096957637897856,
            -4.185093318481106,
            5.8489408677453,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191670676,
            3945.506653225158,
            2864.6366184134445,
            0.04909695762999346,
            -4.185093318475795,
            5.848940867748944,
        ]),
    ];

    {
        let mut dynamics = TwoBody::from_state_vec(init.clone(), earth_geoid.clone());
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
        let mut dynamics = TwoBody::from_state_vec(init.clone(), earth_geoid.clone());
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
        let mut dynamics = TwoBody::from_state_vec(init.clone(), earth_geoid.clone());
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
        let mut dynamics = TwoBody::from_state_vec(init.clone(), earth_geoid);
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
    use nyx::celestia::Cosm;
    use nyx::dynamics::celestial::TwoBody;
    use nyx::propagators::error_ctrl::RSSStatePV;
    use nyx::propagators::*;

    let cosm = Cosm::from_xb("./de438");
    let earth_geoid = cosm.geoid_from_id(3).unwrap();

    let prop_time = 3_600.0 * 24.0;
    let init = Vector6::from_row_slice(&[-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0]);

    let all_rslts = vec![
        Vector6::from_row_slice(&[
            -5971.194191670768,
            3945.506653227154,
            2864.6366184109706,
            0.04909695762764177,
            -4.18509331847428,
            5.8489408677500965,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191670203,
            3945.5066532190967,
            2864.636618421618,
            0.04909695763733907,
            -4.185093318480867,
            5.848940867745654,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191699656,
            3945.50665408017,
            2864.63661724545,
            0.04909695658406228,
            -4.185093317777894,
            5.848940868241106,
        ]),
        Vector6::from_row_slice(&[
            -5971.194191670044,
            3945.5066532117953,
            2864.636618431374,
            0.049096957645996114,
            -4.185093318486724,
            5.848940867741533,
        ]),
        Vector6::from_row_slice(&[
            -5971.19419167081,
            3945.5066532332503,
            2864.6366184022418,
            0.049096957620019005,
            -4.185093318469214,
            5.848940867753748,
        ]),
    ];

    {
        let mut dynamics = TwoBody::from_state_vec(init.clone(), earth_geoid.clone());
        dbg!(dynamics.mu);
        let mut prop = Propagator::new::<RK4Fixed>(&mut dynamics, &PropOpts::with_fixed_step(1.0, RSSStatePV {}));
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[0], "two body prop failed");
    }

    {
        let mut dynamics = TwoBody::from_state_vec(init.clone(), earth_geoid.clone());
        let mut prop = Propagator::new::<Verner56>(&mut dynamics, &PropOpts::with_fixed_step(10.0, RSSStatePV {}));
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[1], "two body prop failed");
    }

    {
        let mut dynamics = TwoBody::from_state_vec(init.clone(), earth_geoid.clone());
        let mut prop = Propagator::new::<Dormand45>(&mut dynamics, &PropOpts::with_fixed_step(10.0, RSSStatePV {}));
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[2], "two body prop failed");
    }

    {
        let mut dynamics = TwoBody::from_state_vec(init.clone(), earth_geoid.clone());
        let mut prop = Propagator::new::<Dormand78>(&mut dynamics, &PropOpts::with_fixed_step(10.0, RSSStatePV {}));
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[3], "two body prop failed");
    }

    {
        let mut dynamics = TwoBody::from_state_vec(init.clone(), earth_geoid);
        let mut prop = Propagator::new::<RK89>(&mut dynamics, &PropOpts::with_fixed_step(10.0, RSSStatePV {}));
        prop.until_time_elapsed(prop_time);
        assert_eq!(prop.state(), all_rslts[4], "two body prop failed");
    }
}
