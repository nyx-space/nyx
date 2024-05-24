extern crate nyx_space as nyx;

use anise::constants::celestial_objects::{EARTH, MOON, SUN};
use anise::constants::frames::IAU_EARTH_FRAME;
use nyx::cosmic::Orbit;
use nyx::od::noise::GaussMarkov;
use nyx::od::prelude::*;
use nyx::time::Epoch;
use nyx::{dynamics::OrbitalDynamics, propagators::Propagator};
use rand::SeedableRng;
use rand_pcg::Pcg64Mcg;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn nil_measurement(almanac: Arc<Almanac>) {
    use hifitime::J2000_OFFSET;
    // Let's create a station and make it estimate the range and range rate of something which is strictly in the same spot.

    let lat = -7.906_635_7;
    let long = 345.5975;
    let height = 0.0;
    let epoch = Epoch::from_mjd_tai(J2000_OFFSET);

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let mut station = GroundStation {
        name: "nil".to_string(),
        latitude_deg: lat,
        longitude_deg: long,
        height_km: height,
        frame: eme2k,
        elevation_mask_deg: 0.0,
        timestamp_noise_s: None,
        range_noise_km: Some(GaussMarkov::ZERO),
        doppler_noise_km_s: Some(GaussMarkov::ZERO),
        integration_time: None,
        light_time_correction: false,
    };

    let at_station = Orbit::from_geodesic(lat, long, height, epoch, eme2k);

    let (_, traj) = Propagator::default(OrbitalDynamics::two_body())
        .with(at_station, almanac.clone())
        .for_duration_with_traj(1.seconds())
        .unwrap();

    assert!(station
        .measure(epoch, &traj, None, almanac)
        .unwrap()
        .is_none());
}

/// Tests that the measurements generated from a topocentric frame are correct.
/// We're propagating a spacecraft in a cislunar trajectory as that will use all of the logic
/// for frame transformation that simpler cases might not use.
/// GMAT script: Cislunar_Measurement_Generation.script
#[allow(clippy::identity_op)]
#[rstest]
fn val_measurements_topo(almanac: Arc<Almanac>) {
    use self::nyx::cosmic::Orbit;
    use self::nyx::md::prelude::*;
    use self::nyx::od::prelude::*;
    use self::nyx::propagators::RK4Fixed;
    use std::str::FromStr;

    let cislunar1 = Orbit::cartesian(
        -6252.59501113,
        1728.23921802,
        1054.21399354,
        -3.86295539,
        -8.85806596,
        -5.08576325,
        Epoch::from_str("2023-11-16T06:36:30.232000 UTC").unwrap(),
        almanac.frame_from_uid(EARTH_J2000).unwrap(),
    );

    // In GMAT's uncommon MJD notation, this epoch corresponds to 29912.78296296296.
    let cislunar2 = Orbit::cartesian(
        4391.84282386,
        -8819.24914059,
        -5415.11431877,
        7.92817749977,
        -1.78800739052,
        -1.69330836191,
        Epoch::from_str("2022-11-29T06:47:28.0 TAI").unwrap(),
        almanac.frame_from_uid(EARTH_J2000).unwrap(),
    );

    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();

    let elevation_mask = 7.0; // in degrees
    let mut dss65_madrid = GroundStation::dss65_madrid(
        elevation_mask,
        GaussMarkov::ZERO,
        GaussMarkov::ZERO,
        iau_earth,
    );

    // Generate the measurements

    // Define the propagator information.
    let prop_time = 12 * Unit::Hour;
    let step_size = 10.0 * Unit::Second;
    let opts = PropOpts::with_fixed_step(step_size);

    let setup =
        Propagator::new::<RK4Fixed>(OrbitalDynamics::point_masses(&[EARTH, MOON, SUN]), opts);

    struct GmatMsrData {
        offset: Duration,
        range: f64,
        range_rate: f64,
    }

    let (_, traj1) = setup
        .with(cislunar1, almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Now iterate the trajectory to generate the measurements.
    let traj1_val_data = vec![
        GmatMsrData {
            offset: 0.29097222222117125 * Unit::Day,
            range: 9.145_755_787_575_61e4,
            range_rate: 2.199_227_723_432_48e0,
        },
        GmatMsrData {
            offset: 0.3368055555547471 * Unit::Day,
            range: 9.996_505_560_799_869e4,
            range_rate: 2.105_490_397_794_733e0,
        },
        GmatMsrData {
            offset: 0.37777777777591837 * Unit::Day,
            range: 1.073_229_118_411_670_2e5,
            range_rate: 2.056_308_226_930_496e0,
        },
        GmatMsrData {
            offset: 0.4187500000007276 * Unit::Day,
            range: 1.145_516_751_191_464_7e5,
            range_rate: 2.031_146_181_775_705_7e0,
        },
        GmatMsrData {
            offset: 0.4874999999992724 * Unit::Day,
            range: 1.265_739_190_638_930_7e5,
            range_rate: 2.021_375_530_901_736_7e0,
        },
    ];

    let mut rng = Pcg64Mcg::from_entropy();
    let mut traj1_msr_cnt = 0;
    for state in traj1.every(1 * Unit::Minute) {
        if dss65_madrid
            .measure(state.epoch(), &traj1, Some(&mut rng), almanac.clone())
            .unwrap()
            .is_some()
        {
            traj1_msr_cnt += 1;
        }
    }
    println!("Generated {} measurements for cislunar1", traj1_msr_cnt);
    // GMAT generates 302 measurements but Nyx generates 303. GMAT generates the same measurement twice at the start of the file... like wut?
    assert_eq!(
        traj1_msr_cnt, 303,
        "incorrect number of measurements generated"
    );
    // Now, let's check specific arbitrarily selected observations
    for truth in &traj1_val_data {
        let now = cislunar1.epoch() + truth.offset;
        let state = traj1.at(now).unwrap();
        // Will panic if the measurement is not visible
        let meas = dss65_madrid
            .measure(state.epoch(), &traj1, Some(&mut rng), almanac.clone())
            .unwrap()
            .unwrap();

        let obs = meas.observation();
        println!(
            "range difference {:e}\t range rate difference: {:e}",
            (obs[0] - truth.range).abs(),
            (obs[1] - truth.range_rate).abs()
        );
    }

    // Second cislunar test
    let traj2_val_data = vec![
        GmatMsrData {
            offset: 0.32777777778028394 * Unit::Day,
            range: 1.020_601_774_210_878_8e5,
            range_rate: 1.956_752_045_319_600_3e0,
        },
        GmatMsrData {
            offset: 0.37222222222408163 * Unit::Day,
            range: 1.093_894_902_936_570_1e5,
            range_rate: 1.867_718_050_780_170_7e0,
        },
        GmatMsrData {
            offset: 0.41319444444889086 * Unit::Day,
            range: 1.159_072_016_126_479_3e5,
            range_rate: 1.819_777_023_286_441_9e0,
        },
        GmatMsrData {
            offset: 0.4541666666700621 * Unit::Day,
            range: 1.223_057_077_408_475e5,
            range_rate: 1.799_383_353_751_318_2e0,
        },
        GmatMsrData {
            offset: 0.4993055555605679 * Unit::Day,
            range: 1.293_208_210_899_899_3e5,
            range_rate: 1.801_787_541_374_800_8e0,
        },
    ];
    let mut traj2_msr_cnt = 0;
    let (_, traj2) = setup
        .with(cislunar2)
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Now iterate the trajectory to count the measurements.
    for state in traj2.every(1 * Unit::Minute) {
        if dss65_madrid
            .measure(state.epoch(), &traj2, Some(&mut rng), almanac.clone())
            .unwrap()
            .is_some()
        {
            traj2_msr_cnt += 1;
        }
    }

    println!("Generated {} measurements for cislunar2", traj2_msr_cnt);
    assert_eq!(
        traj2_msr_cnt, 249,
        "incorrect number of measurements generated"
    );
    // Now, let's check specific arbitrarily selected observations
    for truth in &traj2_val_data {
        let now = cislunar2.epoch() + truth.offset;
        let state = traj2.at(now).unwrap();
        // Will panic if the measurement is not visible
        let meas = dss65_madrid
            .measure(state.epoch(), &traj2, Some(&mut rng), almanac.clone())
            .unwrap()
            .unwrap();
        let obs = meas.observation();
        println!(
            "range difference {:e}\t range rate difference: {:e}",
            (obs[0] - truth.range).abs(),
            (obs[1] - truth.range_rate).abs()
        );
    }
}
