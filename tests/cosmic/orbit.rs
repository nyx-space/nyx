extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use nyx::cosmic::{Cosm, Frame, Orbit};
use nyx::time::{Epoch, TimeUnit};

macro_rules! f64_eq {
    ($x:expr, $val:expr, $msg:expr) => {
        assert!(
            ($x - $val).abs() < 1e-10,
            "{}: {:.2e}",
            $msg,
            ($x - $val).abs()
        )
    };
}

#[test]
fn state_def_circ_inc() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_mjd_tai(21_545.0);
    let cart = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k,
    );
    let cart2 = Orbit::cartesian(
        -2436.45,
        -2436.45,
        6891.037,
        5.088_611,
        -5.088_611,
        0.0,
        Epoch::from_jde_tai(dt.as_jde_tai_days()),
        eme2k,
    );
    assert_eq!(
        cart, cart2,
        "different representations of the datetime are not considered equal"
    );
    f64_eq!(cart.x, -2436.45, "x");
    f64_eq!(cart.y, -2436.45, "y");
    f64_eq!(cart.z, 6891.037, "z");
    f64_eq!(cart.vx, 5.088_611, "vx");
    f64_eq!(cart.vy, -5.088_611, "vy");
    f64_eq!(cart.vz, 0.0, "vz");
    f64_eq!(cart.energy(), -25.842_247_282_849_137, "energy");
    assert_eq!(
        cart.period(),
        6_740.269_063_643_045 * TimeUnit::Second,
        "period"
    );
    f64_eq!(cart.hx(), 35_065.806_679_607_005, "HX");
    f64_eq!(cart.hy(), 35_065.806_679_607_005, "HY");
    f64_eq!(cart.hz(), 24_796.292_541_9, "HZ");
    f64_eq!(cart.sma(), 7_712.186_117_895_043, "sma");
    f64_eq!(cart.ecc(), 0.000_999_582_831_432_052_5, "ecc");
    f64_eq!(cart.inc(), 63.434_003_407_751_14, "inc");
    f64_eq!(cart.raan(), 135.0, "raan");
    f64_eq!(cart.aop(), 90.0, "aop");
    f64_eq!(cart.ta(), 0.0, "ta");
    f64_eq!(cart.tlong(), 225.0, "tlong");
    f64_eq!(cart.ea(), 0.0, "ea");
    f64_eq!(cart.ma(), 0.0, "ma");
    f64_eq!(cart.apoapsis(), 7_719.895_086_731_299, "apo");
    f64_eq!(cart.periapsis(), 7_704.477_149_058_786, "peri");
    f64_eq!(
        cart.semi_parameter(),
        7_712.178_412_142_147,
        "semi parameter"
    );

    let kep = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);
    f64_eq!(kep.x, 8_057.976_452_202_976, "x");
    f64_eq!(kep.y, -0.196_740_370_290_888_9, "y");
    f64_eq!(kep.z, 1_475.383_214_274_138, "z");
    f64_eq!(kep.vx, -0.166_470_488_584_076_31, "vx");
    f64_eq!(kep.vy, 6.913_868_638_275_646_5, "vy");
    f64_eq!(kep.vz, 0.910_157_981_443_279_1, "vz");
    f64_eq!(kep.sma(), 8_191.929_999_999_999, "sma");
    f64_eq!(kep.ecc(), 1.000_000_000_388_51e-06, "ecc");
    f64_eq!(kep.inc(), 12.849_999_999_999_987, "inc");
    f64_eq!(kep.raan(), 306.614, "raan");
    f64_eq!(kep.aop(), 314.189_999_994_618_1, "aop");
    f64_eq!(kep.ta(), 99.887_700_005_381_9, "ta");
    f64_eq!(kep.energy(), -24.328_848_116_377_95, "energy");
    assert_eq!(
        kep.period(),
        7_378.877_993_957_958 * TimeUnit::Second,
        "period"
    );
    f64_eq!(kep.hx(), -10_200.784_799_426_574, "HX");
    f64_eq!(kep.hy(), -7_579.639_346_783_497, "HY");
    f64_eq!(kep.hz(), 55_711.757_929_384_25, "HZ");
    f64_eq!(kep.tlong(), 0.691_700_000_000_082_6, "tlong");
    f64_eq!(kep.ea(), 99.887_643_560_656_85, "ea");
    f64_eq!(kep.ma(), 99.887_587_115_926_96, "ma");
    f64_eq!(kep.apoapsis(), 8_191.938_191_930_002, "apo");
    f64_eq!(kep.periapsis(), 8_191.921_808_069_997, "peri");
    f64_eq!(
        kep.semi_parameter(),
        8_191.929_999_991_808,
        "semi parameter"
    );

    let kep = Orbit::keplerian(8_191.93, 0.2, 12.85, 306.614, 314.19, -99.887_7, dt, eme2k);
    f64_eq!(kep.ta(), 260.1123, "ta");

    // Test that DCMs are valid
    let dcm = kep.dcm_from_traj_frame(Frame::VNC).unwrap();
    assert!(((dcm * dcm.transpose()).determinant() - 1.0).abs() < 1e-12);
    assert!(((dcm.transpose() * dcm).determinant() - 1.0).abs() < 1e-12);

    let dcm = kep.dcm_from_traj_frame(Frame::RCN).unwrap();
    assert!(((dcm * dcm.transpose()).determinant() - 1.0).abs() < 1e-12);
    assert!(((dcm.transpose() * dcm).determinant() - 1.0).abs() < 1e-12);

    let dcm = kep.dcm_from_traj_frame(Frame::RIC).unwrap();
    assert!(((dcm * dcm.transpose()).determinant() - 1.0).abs() < 1e-12);
    assert!(((dcm.transpose() * dcm).determinant() - 1.0).abs() < 1e-12);
}

#[test]
fn state_def_elliptical() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_mjd_tai(21_545.0);
    let cart = Orbit::cartesian(
        5_946.673_548_288_958,
        1_656.154_606_023_661,
        2_259.012_129_598_249,
        -3.098_683_050_943_824,
        4.579_534_132_135_011,
        6.246_541_551_539_432,
        dt,
        eme2k,
    );
    f64_eq!(cart.energy(), -25.842_247_282_849_144, "energy");
    assert_eq!(
        cart.period(),
        6_740.269_063_643_042_5 * TimeUnit::Second,
        "period"
    );
    f64_eq!(cart.hx(), 0.015_409_898_034_704_383, "HX");
    f64_eq!(cart.hy(), -44_146.106_010_690_01, "HY");
    f64_eq!(cart.hz(), 32_364.892_694_481_765, "HZ");
    f64_eq!(cart.sma(), 7_712.186_117_895_041, "sma");
    f64_eq!(cart.ecc(), 0.158_999_999_999_999_95, "ecc");
    f64_eq!(cart.inc(), 53.753_69, "inc");
    f64_eq!(cart.raan(), 1.998_632_864_211_17e-05, "raan");
    f64_eq!(cart.aop(), 359.787_880_000_004, "aop");
    f64_eq!(cart.ta(), 25.434_003_407_751_188, "ta");
    f64_eq!(cart.tlong(), 25.221_903_394_083_824, "tlong");
    f64_eq!(cart.ea(), 21.763_052_882_584_79, "ea");
    f64_eq!(cart.ma(), 18.385_336_330_516_39, "ma");
    f64_eq!(cart.apoapsis(), 8_938.423_710_640_353, "apo");
    f64_eq!(cart.periapsis(), 6_485.948_525_149_73, "peri");
    f64_eq!(
        cart.semi_parameter(),
        7_517.214_340_648_537,
        "semi parameter"
    );

    let kep = Orbit::keplerian(
        8_191.93, 0.024_5, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k,
    );
    f64_eq!(kep.x, 8_087.161_618_048_522_5, "x");
    f64_eq!(kep.y, -0.197_452_943_772_520_73, "y");
    f64_eq!(kep.z, 1_480.726_901_246_883, "z");
    f64_eq!(kep.vx, -0.000_168_592_186_843_952_16, "vx");
    f64_eq!(kep.vy, 6.886_845_792_370_852, "vy");
    f64_eq!(kep.vz, 0.936_931_260_302_891_8, "vz");
    f64_eq!(kep.sma(), 8_191.930_000_000_003, "sma");
    f64_eq!(kep.ecc(), 0.024_500_000_000_000_348, "ecc");
    f64_eq!(kep.inc(), 12.850_000_000_000_016, "inc");
    f64_eq!(kep.raan(), 306.614, "raan");
    f64_eq!(kep.aop(), 314.190_000_000_000_4, "aop");
    f64_eq!(kep.ta(), 99.887_699_999_999_58, "ta");
    f64_eq!(kep.energy(), -24.328_848_116_377_94, "energy");
    assert_eq!(
        kep.period(),
        7_378.877_993_957_964 * TimeUnit::Second,
        "period"
    );
    f64_eq!(kep.hx(), -10_197.722_829_337_885, "HX");
    f64_eq!(kep.hy(), -7_577.364_166_057_776, "HY");
    f64_eq!(kep.hz(), 55_695.034_928_191_49, "HZ");
    f64_eq!(kep.tlong(), 0.691_699_999_999_855_2, "tlong");
    f64_eq!(kep.ea(), 98.501_748_370_880_22, "ea");
    f64_eq!(kep.ma(), 97.113_427_049_323_43, "ma");
    f64_eq!(kep.apoapsis(), 8_392.632_285_000_007, "apo");
    f64_eq!(kep.periapsis(), 7_991.227_715_000_001, "peri");
    f64_eq!(
        kep.semi_parameter(),
        8_187.012_794_017_503,
        "semi parameter"
    );

    // Test that DCMs are valid
    let dcm = kep.dcm_from_traj_frame(Frame::VNC).unwrap();
    assert!(((dcm * dcm.transpose()).determinant() - 1.0).abs() < 1e-12);
    assert!(((dcm.transpose() * dcm).determinant() - 1.0).abs() < 1e-12);

    let dcm = kep.dcm_from_traj_frame(Frame::RCN).unwrap();
    assert!(((dcm * dcm.transpose()).determinant() - 1.0).abs() < 1e-12);
    assert!(((dcm.transpose() * dcm).determinant() - 1.0).abs() < 1e-12);

    let dcm = kep.dcm_from_traj_frame(Frame::RIC).unwrap();
    assert!(((dcm * dcm.transpose()).determinant() - 1.0).abs() < 1e-12);
    assert!(((dcm.transpose() * dcm).determinant() - 1.0).abs() < 1e-12);
}

#[test]
fn state_def_circ_eq() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_mjd_tai(21_545.0);
    let cart = Orbit::cartesian(
        -38_892.724_449_149_02,
        16_830.384_772_891_86,
        0.722_659_929_135_562_2,
        -1.218_008_333_846_6,
        -2.814_651_172_605_98,
        1.140_294_223_185_661e-5,
        dt,
        eme2k,
    );
    f64_eq!(cart.energy(), -4.702_902_670_552_006, "energy");
    assert_eq!(
        cart.period(),
        86_820.776_152_986_1 * TimeUnit::Second,
        "period"
    );
    f64_eq!(cart.hx(), 2.225_951_522_241_969_5, "HX");
    f64_eq!(cart.hy(), -0.436_714_326_090_944_6, "HY");
    f64_eq!(cart.hz(), 129_969.001_391_865_75, "HZ");
    f64_eq!(cart.sma(), 42_378.129_999_999_98, "sma");
    f64_eq!(cart.ecc(), 9.999_999_809_555_511e-9, "ecc");
    f64_eq!(cart.inc(), 0.001_000_000_401_564_538_6, "inc");
    f64_eq!(cart.raan(), 78.9, "raan");
    f64_eq!(cart.aop(), 65.399_999_847_186_78, "aop");
    f64_eq!(cart.ta(), 12.300_000_152_813_197, "ta");
    f64_eq!(cart.tlong(), 156.599_999_999_999_97, "tlong");
    f64_eq!(cart.ea(), 12.300_000_030_755_777, "ea");
    f64_eq!(cart.ma(), 12.299_999_908_698_359, "ma");
    f64_eq!(cart.apoapsis(), 42_378.130_423_781_27, "apo");
    f64_eq!(cart.periapsis(), 42_378.129_576_218_69, "peri");
    f64_eq!(
        cart.semi_parameter(),
        42_378.129_999_999_976,
        "semi parameter"
    );

    let kep = Orbit::keplerian(18191.098, 1e-6, 1e-6, 306.543, 314.32, 98.765, dt, eme2k);
    f64_eq!(kep.x, 18_190.717_357_886_37, "x");
    f64_eq!(kep.y, -118.107_162_539_218_69, "y");
    f64_eq!(kep.z, 0.000_253_845_647_633_053_35, "z");
    f64_eq!(kep.vx, 0.030_396_440_130_264_88, "vx");
    f64_eq!(kep.vy, 4.680_909_107_924_576, "vy");
    f64_eq!(kep.vz, 4.907_089_816_726_583e-8, "vz");
    f64_eq!(kep.sma(), 18_191.098_000_000_013, "sma");
    f64_eq!(kep.ecc(), 9.999_999_997_416_087e-7, "ecc");
    f64_eq!(kep.inc(), 1.207_418_269_725_733_3e-6, "inc");
    f64_eq!(kep.raan(), 306.543, "raan");
    f64_eq!(kep.aop(), 314.320_000_025_403_66, "aop");
    f64_eq!(kep.ta(), 98.764_999_974_596_28, "ta");
    f64_eq!(kep.energy(), -10.955_920_349_063_035, "energy");
    assert_eq!(
        kep.period(),
        24_417.396_242_570_256 * TimeUnit::Second,
        "period"
    );
    f64_eq!(kep.hx(), -0.001_194_024_028_558_358_7, "HX");
    f64_eq!(kep.hy(), -0.000_884_918_835_027_750_6, "HY");
    f64_eq!(kep.hz(), 85_152.684_597_507_06, "HZ");
    f64_eq!(kep.tlong(), 359.627_999_999_999_93, "tlong");
    f64_eq!(kep.ea(), 98.764_943_347_932_57, "ea");
    f64_eq!(kep.ma(), 98.764_886_721_264_56, "ma");
    f64_eq!(kep.apoapsis(), 18_191.116_191_098_008, "apo");
    f64_eq!(kep.periapsis(), 18_191.079_808_902_017, "peri");
    f64_eq!(
        kep.semi_parameter(),
        18_191.097_999_981_823,
        "semi parameter"
    );

    // Test that DCMs are valid
    let dcm = kep.dcm_from_traj_frame(Frame::VNC).unwrap();
    assert!(((dcm * dcm.transpose()).determinant() - 1.0).abs() < 1e-12);
    assert!(((dcm.transpose() * dcm).determinant() - 1.0).abs() < 1e-12);

    let dcm = kep.dcm_from_traj_frame(Frame::RCN).unwrap();
    assert!(((dcm * dcm.transpose()).determinant() - 1.0).abs() < 1e-12);
    assert!(((dcm.transpose() * dcm).determinant() - 1.0).abs() < 1e-12);

    let dcm = kep.dcm_from_traj_frame(Frame::RIC).unwrap();
    assert!(((dcm * dcm.transpose()).determinant() - 1.0).abs() < 1e-12);
    assert!(((dcm.transpose() * dcm).determinant() - 1.0).abs() < 1e-12);
}

#[test]
fn state_def_equatorial() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_mjd_tai(21_545.0);
    let cart = Orbit::cartesian(
        -7273.338970882,
        253.990592670,
        0.022164861,
        -0.258285289,
        -7.396322922,
        -0.000645451,
        dt,
        eme2k,
    );

    f64_eq!(cart.sma(), 7278.136188379306, "sma");
    f64_eq!(cart.ecc(), 4.99846643158263e-05, "ecc");
    f64_eq!(cart.inc(), 0.005000000478594339, "inc");
    f64_eq!(cart.raan(), 360.0, "raan");
    f64_eq!(cart.aop(), 177.9999736473912, "aop");
    f64_eq!(cart.ta(), 2.650826247094554e-05, "ta");
}

#[test]
fn state_def_reciprocity() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_mjd_tai(21_545.0);

    assert_eq!(
        Orbit::cartesian(
            -38_892.724_449_149_02,
            16_830.384_772_891_86,
            0.722_659_929_135_562_2,
            -1.218_008_333_846_6,
            -2.814_651_172_605_98,
            1.140_294_223_185_661e-5,
            dt,
            eme2k
        ),
        Orbit::keplerian(
            42_378.129_999_999_98,
            9.999_999_809_555_511e-9,
            0.001_000_000_401_564_538_6,
            78.9,
            65.399_999_847_186_78,
            12.300_000_152_813_197,
            dt,
            eme2k
        ),
        "circ_eq"
    );

    assert_eq!(
        Orbit::cartesian(
            5_946.673_548_288_958,
            1_656.154_606_023_661,
            2_259.012_129_598_249,
            -3.098_683_050_943_824,
            4.579_534_132_135_011,
            6.246_541_551_539_432,
            dt,
            eme2k
        ),
        Orbit::keplerian(
            7_712.186_117_895_041,
            0.158_999_999_999_999_95,
            53.75369,
            1.998_632_864_211_17e-5,
            359.787_880_000_004,
            25.434_003_407_751_188,
            dt,
            eme2k
        ),
        "elliptical"
    );

    assert_eq!(
        Orbit::cartesian(-2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, eme2k),
        Orbit::keplerian(
            7_712.186_117_895_043,
            0.000_999_582_831_432_052_5,
            63.434_003_407_751_14,
            135.0,
            90.0,
            0.0,
            dt,
            eme2k
        ),
        "circ_inc"
    );
}

#[test]
fn geodetic_vallado() {
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    dbg!(eme2k.semi_major_radius());
    let dt = Epoch::from_mjd_tai(51_545.0);
    // Test case from Vallado, 4th Ed., page 173, Example 3-3
    let ri = 6524.834;
    let ri_val = 6_524.833_999_999_999;
    let rj = 6862.875;
    let rj_val = 6_862.874_999_999_999;
    let rk = 6448.296;
    let lat = 34.352_495_139_917_26; // Valldo: 34.352496
    let long = 46.446_416_856_789_96; // Vallado 46.4464
    let height = 5_085.219_430_345_17; // Valldo: 5085.22
    let r = Orbit::from_position(ri, rj, rk, dt, eme2k);
    f64_eq!(r.geodetic_latitude(), lat, "latitude (φ)");
    f64_eq!(r.geodetic_longitude(), long, "longitude (λ)");
    f64_eq!(r.geodetic_height(), height, "height");
    let r = Orbit::from_geodesic(lat, long, height, dt, eme2k);
    f64_eq!(r.x, ri_val, "r_i");
    f64_eq!(r.y, rj_val, "r_j");
    f64_eq!(r.z, rk, "r_k");

    // Test case from Vallado, 4th Ed., page 173, Example 3-4
    let lat = -7.906_635_7;
    let lat_val = -7.906_635_699_999_994_5;
    let long = 345.5975;
    let height = 56.0e-3;
    let height_val = 0.056_000_000_000_494_765;
    let ri = 6_119.399_587_411_616;
    let rj = -1_571.479_380_333_195;
    let rk = -871.561_161_926_003_9;
    let r = Orbit::from_geodesic(lat, long, height, dt, eme2k);
    f64_eq!(r.x, ri, "r_i");
    f64_eq!(r.y, rj, "r_j");
    f64_eq!(r.z, rk, "r_k");
    let r = Orbit::from_position(ri, rj, rk, dt, eme2k);
    f64_eq!(r.geodetic_latitude(), lat_val, "latitude (φ)");
    f64_eq!(r.geodetic_longitude(), long, "longitude (λ)");
    f64_eq!(r.geodetic_height(), height_val, "height");
}

#[test]
fn with_init() {
    use nyx::utils::{between_0_360, between_pm_180};
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_tai_at_midnight(2021, 3, 4);
    let kep = Orbit::keplerian(
        8_191.93, 0.024_5, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k,
    );
    for sma_incr in 100..1000 {
        let new_sma = kep.sma() + f64::from(sma_incr);
        f64_eq!(kep.with_sma(new_sma).sma(), new_sma, "wrong sma");
    }
    for ecc_incr in 0..100 {
        let new_ecc = kep.ecc() + f64::from(ecc_incr) / 100.0;
        let new_state = kep.with_ecc(new_ecc);
        f64_eq!(
            new_state.ecc(),
            new_ecc,
            format!("wrong ecc: got {}\twanted {}", new_state.inc(), new_ecc)
        );
    }
    for angle_incr in 0..360 {
        let new_aop = between_0_360(kep.aop() + f64::from(angle_incr));
        let new_state = kep.add_aop(f64::from(angle_incr));
        f64_eq!(
            new_state.aop(),
            new_aop,
            format!("wrong aop: got {}\twanted {}", new_state.aop(), new_aop)
        );
    }
    for angle_incr in 0..360 {
        let new_raan = between_0_360(kep.raan() + f64::from(angle_incr));
        let new_state = kep.add_raan(f64::from(angle_incr));
        f64_eq!(
            new_state.raan(),
            new_raan,
            format!("wrong raan: got {}\twanted {}", new_state.raan(), new_raan)
        );
    }
    for angle_incr in 0..360 {
        let new_ta = between_0_360(kep.ta() + f64::from(angle_incr));
        let new_state = kep.with_ta(new_ta);
        f64_eq!(
            new_state.ta(),
            new_ta,
            format!("wrong ta: got {}\twanted {}", new_state.aop(), new_ta)
        );
    }
    for angle_incr in 0..360 {
        // NOTE: Inclination is bounded between 0 and 180, hence the slightly different logic here.
        let new_inc = between_pm_180(kep.inc() + f64::from(angle_incr)).abs();
        let new_state = kep.add_inc(f64::from(angle_incr));
        f64_eq!(
            new_state.inc(),
            new_inc,
            format!("wrong inc: got {}\twanted {}", new_state.inc(), new_inc)
        );
    }
    for apsis_delta in 100..1000 {
        let new_ra = kep.apoapsis() + f64::from(apsis_delta);
        let new_rp = kep.periapsis() - f64::from(apsis_delta);
        let new_orbit = kep.with_apoapsis_periapsis(new_ra, new_rp);
        f64_eq!(new_orbit.apoapsis(), new_ra, "wrong ra");
        f64_eq!(new_orbit.periapsis(), new_rp, "wrong rp");
    }
}

#[test]
fn orbit_dual_test() {
    use nyx::cosmic::OrbitDual;
    use nyx::md::StateParameter;
    let cosm = Cosm::de438_gmat();
    let eme2k = cosm.frame("EME2000");

    let dt = Epoch::from_mjd_tai(21_545.0);
    let cart = Orbit::cartesian(
        -2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 1.0, dt, eme2k,
    );

    println!("{:x}", cart);

    let cart_dual = OrbitDual::from(cart);
    println!("{}", eme2k.gm());
    println!("|v| = {}", cart_dual.vmag());
    println!("|v|^2 = {}", cart_dual.vmag().dual * cart_dual.vmag().dual);
    println!("x = {}", cart_dual.x);
    println!("y = {}", cart_dual.y);
    println!("z = {}", cart_dual.z);
    println!("\\xi = {}\n\n", cart_dual.energy());

    // Now print the table
    let params = vec![
        StateParameter::AoL,
        StateParameter::AoP,
        StateParameter::Apoapsis,
        StateParameter::C3,
        StateParameter::Declination,
        StateParameter::EccentricAnomaly,
        StateParameter::Eccentricity,
        StateParameter::Energy,
        StateParameter::FlightPathAngle,
        StateParameter::GeodeticHeight,
        StateParameter::GeodeticLatitude,
        StateParameter::GeodeticLongitude,
        StateParameter::Hmag,
        StateParameter::HX,
        StateParameter::HY,
        StateParameter::HZ,
        StateParameter::Inclination,
        StateParameter::MeanAnomaly,
        StateParameter::Periapsis,
        StateParameter::RightAscension,
        StateParameter::RAAN,
        StateParameter::Rmag,
        StateParameter::SemiParameter,
        StateParameter::SMA,
        StateParameter::TrueLongitude,
        StateParameter::Vmag,
        StateParameter::X,
        StateParameter::Y,
        StateParameter::Z,
        StateParameter::VX,
        StateParameter::VY,
        StateParameter::VZ,
    ];
    for param in &params {
        let dual = cart_dual.partial_for(param).unwrap();
        println!(
            "{:?} & {} & {} & {} & {} & {} & {} \\\\",
            param,
            if dual.wtr_x().abs() > 1e-2 {
                "large"
            } else if dual.wtr_x().abs() > 0.0 {
                "small"
            } else {
                "none"
            },
            if dual.wtr_y().abs() > 1e-2 {
                "large"
            } else if dual.wtr_y().abs() > 0.0 {
                "small"
            } else {
                "none"
            },
            if dual.wtr_z().abs() > 1e-2 {
                "large"
            } else if dual.wtr_z().abs() > 0.0 {
                "small"
            } else {
                "none"
            },
            if dual.wtr_vx().abs() > 1e-2 {
                "large"
            } else if dual.wtr_vx().abs() > 0.0 {
                "small"
            } else {
                "none"
            },
            if dual.wtr_vy().abs() > 1e-2 {
                "large"
            } else if dual.wtr_vy().abs() > 0.0 {
                "small"
            } else {
                "none"
            },
            if dual.wtr_vz().abs() > 1e-2 {
                "large"
            } else if dual.wtr_vz().abs() > 0.0 {
                "small"
            } else {
                "none"
            },
            // dual,
        );
    }
}
