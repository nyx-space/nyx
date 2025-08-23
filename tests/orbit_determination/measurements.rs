extern crate nyx_space as nyx;

use anise::constants::celestial_objects::{EARTH, MOON, SUN};
use anise::constants::frames::IAU_EARTH_FRAME;
use anise::constants::usual_planetary_constants::MEAN_EARTH_ANGULAR_VELOCITY_DEG_S;
use indexmap::{IndexMap, IndexSet};
use nalgebra::{Const, U2};
use nyx::cosmic::Orbit;
use nyx::dynamics::SpacecraftDynamics;
use nyx::md::prelude::*;
use nyx::od::prelude::*;
use nyx::time::Epoch;
use nyx::{dynamics::OrbitalDynamics, propagators::Propagator};
use nyx_space::propagators::IntegratorMethod;
use rand::SeedableRng;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::str::FromStr;

use rand_pcg::Pcg64Mcg;

use anise::{constants::frames::EARTH_J2000, prelude::Almanac};
use rstest::*;
use sensitivity::TrackerSensitivity;
use std::sync::Arc;

#[fixture]
fn almanac() -> Arc<Almanac> {
    use crate::test_almanac_arcd;
    test_almanac_arcd()
}

#[rstest]
fn nil_measurement(almanac: Arc<Almanac>) {
    use hifitime::JD_J2000;
    // Let's create a station and make it estimate the range and range rate of something which is strictly in the same spot.

    let lat = -7.906_635_7;
    let long = 345.5975;
    let height = 0.0;
    let epoch = Epoch::from_mjd_tai(JD_J2000);

    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();

    let mut stochastics = IndexMap::new();
    stochastics.insert(MeasurementType::Range, StochasticNoise::MIN);
    stochastics.insert(MeasurementType::Doppler, StochasticNoise::MIN);

    let mut station = GroundStation {
        name: "nil".to_string(),
        latitude_deg: lat,
        longitude_deg: long,
        height_km: height,
        frame: eme2k,
        elevation_mask_deg: 0.0,
        timestamp_noise_s: None,
        stochastic_noises: Some(stochastics),
        integration_time: None,
        light_time_correction: false,
        ..Default::default()
    };

    let at_station = Orbit::try_latlongalt(
        lat,
        long,
        height,
        MEAN_EARTH_ANGULAR_VELOCITY_DEG_S,
        epoch,
        eme2k,
    )
    .unwrap();

    let (_, traj) = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::two_body()))
        .with(at_station.into(), almanac.clone())
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
    let mut msr_types = IndexSet::new();
    msr_types.insert(MeasurementType::Range);
    msr_types.insert(MeasurementType::Doppler);

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
        StochasticNoise::MIN,
        StochasticNoise::MIN,
        iau_earth,
    );

    // Generate the measurements

    // Define the propagator information.
    let prop_time = 12 * Unit::Hour;
    let step_size = 10.0 * Unit::Second;
    let opts = IntegratorOptions::with_fixed_step(step_size);

    let setup = Propagator::new(
        SpacecraftDynamics::new(OrbitalDynamics::point_masses(vec![EARTH, MOON, SUN])),
        IntegratorMethod::RungeKutta4,
        opts,
    );

    struct GmatMsrData {
        offset: Duration,
        range: f64,
        range_rate: f64,
    }

    let (_, traj1) = setup
        .with(cislunar1.into(), almanac.clone())
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

    let mut rng = Pcg64Mcg::from_os_rng();
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
    println!("Generated {traj1_msr_cnt} measurements for cislunar1");
    // GMAT generates 302 measurements but Nyx generates 303. GMAT generates the same measurement twice at the start of the file... like wut?
    assert_eq!(
        traj1_msr_cnt, 303,
        "incorrect number of measurements generated"
    );
    // Now, let's check specific arbitrarily selected observations
    for truth in &traj1_val_data {
        let epoch = cislunar1.epoch + truth.offset;
        let state = traj1.at(epoch).unwrap();
        // Will panic if the measurement is not visible
        let meas = dss65_madrid
            .measure(state.epoch(), &traj1, Some(&mut rng), almanac.clone())
            .unwrap()
            .unwrap();

        let obs = meas.observation::<U2>(&msr_types);
        println!("@{epoch}");
        println!(
            "COMPUTED range = {:.3} km\trange-rate = {:.3} km/s",
            obs[0], obs[1]
        );
        println!(
            "EXPECTED range = {:.3} km\trange-rate = {:.3} km/s",
            truth.range, truth.range_rate
        );
        println!(
            "ERROR    range = {:.3} km\trange-rate = {:.3} km/s",
            obs[0] - truth.range,
            obs[1] - truth.range_rate
        );
        assert!(
            (obs[1] - truth.range_rate).abs() < 1e-3,
            "range rate error greater than 1 m/s"
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
        .with(cislunar2.into(), almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    // Now iterate the trajectory to count the measurements.
    for state in traj2.every(1 * Unit::Minute) {
        if let Some(msr) = dss65_madrid
            .measure(state.epoch(), &traj2, Some(&mut rng), almanac.clone())
            .unwrap()
        {
            if traj2_msr_cnt == 0 {
                println!("{msr}");
            }
            traj2_msr_cnt += 1;
        }
    }

    println!("Generated {traj2_msr_cnt} measurements for cislunar2");
    assert_eq!(
        traj2_msr_cnt, 249,
        "incorrect number of measurements generated"
    );
    // Now, let's check specific arbitrarily selected observations
    for truth in &traj2_val_data {
        let epoch = cislunar2.epoch + truth.offset;
        let state = traj2.at(epoch).unwrap();
        // Will panic if the measurement is not visible
        let meas = dss65_madrid
            .measure(state.epoch(), &traj2, Some(&mut rng), almanac.clone())
            .unwrap()
            .unwrap();
        let obs = meas.observation::<U2>(&msr_types);

        println!(
            "COMPUTED range = {:.3} km\trange-rate = {:.3} km/s",
            obs[0], obs[1]
        );
        println!(
            "EXPECTED range = {:.3} km\trange-rate = {:.3} km/s",
            truth.range, truth.range_rate
        );
        println!(
            "ERROR    range = {:.3} km\trange-rate = {:.3} km/s",
            obs[0] - truth.range,
            obs[1] - truth.range_rate
        );
        assert!(
            (obs[1] - truth.range_rate).abs() < 1e-3,
            "range rate error greater than 1 m/s"
        );
    }
}

/// Verifies that the sensitivity matrix is reasonably well.
#[allow(clippy::identity_op)]
#[rstest]
fn verif_sensitivity_mat(almanac: Arc<Almanac>) {
    let cislunar1 = Orbit::cartesian(
        58643.769540,
        -61696.435624,
        -36178.745722,
        2.148654,
        -1.202489,
        -0.714016,
        Epoch::from_str("2022-11-16T13:35:31.0 UTC").unwrap(),
        almanac.frame_from_uid(EARTH_J2000).unwrap(),
    );

    let cislunar_sc: Spacecraft = cislunar1.into();

    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();

    let mut dss65_madrid =
        GroundStation::dss65_madrid(0.0, StochasticNoise::ZERO, StochasticNoise::ZERO, iau_earth)
            .with_msr_type(MeasurementType::Azimuth, StochasticNoise::ZERO)
            .with_msr_type(MeasurementType::Elevation, StochasticNoise::ZERO);

    let mut cislunar_sc_pert = cislunar_sc;
    cislunar_sc_pert.orbit.radius_km.x += 1.0;
    cislunar_sc_pert.orbit.radius_km.y -= 1.0;
    cislunar_sc_pert.orbit.radius_km.z += 1.0;
    cislunar_sc_pert.orbit.velocity_km_s.x += 1.0e-3;
    cislunar_sc_pert.orbit.velocity_km_s.y -= 1.0e-3;
    cislunar_sc_pert.orbit.velocity_km_s.z += 1.0e-3;

    let truth_meas = dss65_madrid
        .measure_instantaneous(cislunar_sc, None, almanac.clone())
        .expect("successful measurement")
        .expect("a measurement");

    let pert_meas = dss65_madrid
        .measure_instantaneous(cislunar_sc_pert, None, almanac.clone())
        .expect("successful measurement")
        .expect("a measurement");

    for msr_type in [
        MeasurementType::Range,
        MeasurementType::Doppler,
        MeasurementType::Elevation,
        MeasurementType::Azimuth,
    ] {
        let mut msr_types = IndexSet::new();
        msr_types.insert(msr_type);

        let truth_obs = truth_meas.observation::<Const<1>>(&msr_types);

        let pert_obs = pert_meas.observation::<Const<1>>(&msr_types);

        // Given this observation, feed it to the sensitivity matrix, and we should find the original state.

        let h_tilde = dss65_madrid
            .h_tilde::<Const<1>>(&truth_meas, &msr_types, &cislunar_sc, almanac.clone())
            .expect("sensitivity should not fail");

        let delta_state = cislunar_sc.to_vector().fixed_rows::<9>(0)
            - cislunar_sc_pert.to_vector().fixed_rows::<9>(0);

        let delta_obs = h_tilde * delta_state;
        let computed_obs = truth_obs - delta_obs;

        let sensitivity_error = (pert_obs - computed_obs)[0];
        println!(
            "{msr_type:?} error = {sensitivity_error:.3e} {}",
            msr_type.unit()
        );

        assert!(sensitivity_error.abs() < 1e-3);
    }
}

/// Validate the white noise modeling on an arbitrary trajectory for all measurement types
/// is within the covariance and the sigma bounds as configured by the stochastic process.
#[allow(clippy::identity_op)]
#[rstest]
fn val_measurement_noise(almanac: Arc<Almanac>) {
    use arrow::datatypes::FieldRef;
    use arrow::datatypes::{DataType, Field, Schema};
    use parquet::arrow::ArrowWriter;
    use serde_arrow::schema::{SchemaLike, TracingOptions};

    // Build an example trajectory.
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    let epoch = Epoch::from_gregorian_tai_at_midnight(2025, 8, 22);

    let orbit =
        Orbit::keplerian_apsis_radii(99_000.0, 180.0, 32.0, 128.0, 256.0, 0.0, epoch, eme2k);

    let sc: Spacecraft = orbit.into();

    println!("{sc}");

    let (_, traj) = Propagator::default(SpacecraftDynamics::new(OrbitalDynamics::point_masses(
        vec![MOON, SUN],
    )))
    .with(sc, almanac.clone())
    .for_duration_with_traj(orbit.period().unwrap() * 2)
    .unwrap();

    println!("{traj}");

    // Set up the noise modeling.

    let iau_earth = almanac.frame_from_uid(IAU_EARTH_FRAME).unwrap();

    let range_wn = WhiteNoise::constant_white_noise(2.0e-3);
    let doppler_wn = WhiteNoise::constant_white_noise(3e-6);
    let angle_wn = WhiteNoise::constant_white_noise(2e-2);

    let noisy_ground_station = GroundStation::dss13_goldstone(
        10.0,
        StochasticNoise {
            white_noise: Some(range_wn),
            bias: None,
        },
        StochasticNoise {
            white_noise: Some(doppler_wn),
            bias: None,
        },
        iau_earth,
    )
    .with_msr_type(
        MeasurementType::Azimuth,
        StochasticNoise {
            white_noise: Some(angle_wn),
            bias: None,
        },
    )
    .with_msr_type(
        MeasurementType::Elevation,
        StochasticNoise {
            white_noise: Some(angle_wn),
            bias: None,
        },
    );

    let mut noisy_devices = BTreeMap::new();
    noisy_devices.insert("Station".to_string(), noisy_ground_station);

    let mut configs = BTreeMap::new();
    configs.insert(
        "Station".to_string(),
        TrkConfig::from_sample_rate(Unit::Minute * 1),
    );

    let mut noisy_trk_sim =
        TrackingArcSim::new(noisy_devices, traj.clone(), configs.clone()).unwrap();
    noisy_trk_sim.build_schedule(almanac.clone()).unwrap();
    let noisy_trk_data = noisy_trk_sim
        .generate_measurements(almanac.clone())
        .unwrap();

    let perfect_ground_station = GroundStation::dss13_goldstone(
        10.0,
        StochasticNoise::ZERO,
        StochasticNoise::ZERO,
        iau_earth,
    )
    .with_msr_type(MeasurementType::Azimuth, StochasticNoise::ZERO)
    .with_msr_type(MeasurementType::Elevation, StochasticNoise::ZERO);

    let mut noisy_devices = BTreeMap::new();
    noisy_devices.insert("Station".to_string(), perfect_ground_station);

    let mut perfect_trk_sim =
        TrackingArcSim::new(noisy_devices, traj.clone(), configs.clone()).unwrap();
    perfect_trk_sim.build_schedule(almanac.clone()).unwrap();
    let perfect_trk_data = perfect_trk_sim
        .generate_measurements(almanac.clone())
        .unwrap();

    assert_eq!(perfect_trk_data.len(), noisy_trk_data.len());

    #[derive(Serialize, Deserialize)]
    struct Noises {
        noise_err: f64,
        std_dev: f64,
        name: String,
        out_of_family: bool,
    }

    let mut records = vec![];

    // Compute the diff between both tracking data.
    for (msr_type, std_dev) in [
        (MeasurementType::Range, range_wn.sigma),
        (MeasurementType::Doppler, doppler_wn.sigma),
        (MeasurementType::Azimuth, angle_wn.sigma),
        (MeasurementType::Elevation, angle_wn.sigma),
    ] {
        let noisy_subset = noisy_trk_data.clone().filter_by_measurement_type(msr_type);
        let perfect_subset = perfect_trk_data
            .clone()
            .filter_by_measurement_type(msr_type);

        let mut out_of_family = 0;

        for (noisy_msr, perfect_msr) in noisy_subset
            .measurements
            .values()
            .zip(perfect_subset.measurements.values())
        {
            // First, test that the filter by measurement type works
            assert_eq!(noisy_msr.data.len(), 1);
            assert_eq!(perfect_msr.data.len(), 1);
            let noisy_val = noisy_msr.data[&msr_type];
            let perfect_val = perfect_msr.data[&msr_type];

            let noise_err = noisy_val - perfect_val;
            let is_oof = noise_err.abs() > 3.0 * std_dev;
            if is_oof {
                out_of_family += 1;
            }
            records.push(Noises {
                noise_err,
                std_dev,
                name: format!("{msr_type:?}"),
                out_of_family: is_oof,
            });
        }

        let prct_in_family = 100.0 - (out_of_family as f64) / (noisy_subset.len() as f64) * 100.0;

        // We expect 99.7% to be in-family for 3-sigma

        println!("percentage IN FAMILY for {msr_type:?} = {prct_in_family:.4} %");

        assert!((prct_in_family - 99.7).abs() < 0.3);
    }

    let hdrs = vec![
        Field::new("Noise err", DataType::Float64, false),
        Field::new("StdDev", DataType::Float64, false),
        Field::new("Name", DataType::LargeUtf8, false),
        Field::new("Out of family", DataType::Boolean, false),
    ];

    // Build the schema
    let schema = Arc::new(Schema::new(hdrs));
    // Determine Arrow schema
    let fields = Vec::<FieldRef>::from_type::<Noises>(TracingOptions::default()).unwrap();

    // Build the record batch
    let batch = serde_arrow::to_record_batch(&fields, &records).unwrap();

    let file = File::create("data/04_output/noise_val.parquet").unwrap();

    let mut writer = ArrowWriter::try_new(file, schema.clone(), None).unwrap();
    // let batch = RecordBatch::try_new(schema, record)
    //     .context(ArrowSnafu {
    //         action: "writing OD results (building batch record)",
    //     })
    //     .context(ODIOSnafu)?;

    writer.write(&batch).unwrap()
}
