extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use hifitime::Epoch;

macro_rules! f64_eq {
    ($x:expr, $val:expr, $msg:expr) => {
        assert!(($x - $val).abs() < 1e-10, $msg)
    };
}

#[test]
fn state_def_() {
    pel::init();
}

#[test]
fn state_def_circ_inc() {
    use nyx::celestia::{Cosm, Geoid, State};
    let cosm = Cosm::from_xb("./de438s");
    let mut earth_geoid = cosm.geoid_from_id(3).unwrap();
    // Let's use the GMAT GM value for which these tests we written.
    earth_geoid.gm = 398_600.441_5;
    let dt = Epoch::from_mjd_tai(21_545.0);
    let cart = State::<Geoid>::from_cartesian(-2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, earth_geoid);
    let cart2 = State::<Geoid>::from_cartesian(
        -2436.45,
        -2436.45,
        6891.037,
        5.088_611,
        -5.088_611,
        0.0,
        Epoch::from_jde_tai(dt.as_jde_tai_days()),
        earth_geoid,
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
    f64_eq!(cart.period(), 6_740.269_063_643_045, "period");
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
    f64_eq!(cart.semi_parameter(), 7_712.178_412_142_147, "semi parameter");

    let kep = State::<Geoid>::from_keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, earth_geoid);
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
    f64_eq!(kep.period(), 7_378.877_993_957_958, "period");
    f64_eq!(kep.hx(), -10_200.784_799_426_574, "HX");
    f64_eq!(kep.hy(), -7_579.639_346_783_497, "HY");
    f64_eq!(kep.hz(), 55_711.757_929_384_25, "HZ");
    f64_eq!(kep.tlong(), 0.691_700_000_000_082_6, "tlong");
    f64_eq!(kep.ea(), 99.887_643_560_656_85, "ea");
    f64_eq!(kep.ma(), 99.887_587_115_926_96, "ma");
    f64_eq!(kep.apoapsis(), 8_191.938_191_930_002, "apo");
    f64_eq!(kep.periapsis(), 8_191.921_808_069_997, "peri");
    f64_eq!(kep.semi_parameter(), 8_191.929_999_991_808, "semi parameter");
}

#[test]
fn state_def_elliptical() {
    use nyx::celestia::{Cosm, Geoid, State};
    let cosm = Cosm::from_xb("./de438s");
    let mut earth_geoid = cosm.geoid_from_id(3).unwrap();
    // Let's use the GMAT GM value for which these tests we written.
    earth_geoid.gm = 398_600.441_5;
    let dt = Epoch::from_mjd_tai(21_545.0);
    let cart = State::<Geoid>::from_cartesian(
        5_946.673_548_288_958,
        1_656.154_606_023_661,
        2_259.012_129_598_249,
        -3.098_683_050_943_824,
        4.579_534_132_135_011,
        6.246_541_551_539_432,
        dt,
        earth_geoid,
    );
    f64_eq!(cart.energy(), -25.842_247_282_849_144, "energy");
    f64_eq!(cart.period(), 6_740.269_063_643_042_5, "period");
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
    f64_eq!(cart.semi_parameter(), 7_517.214_340_648_537, "semi parameter");

    let kep = State::<Geoid>::from_keplerian(8_191.93, 0.024_5, 12.85, 306.614, 314.19, 99.887_7, dt, earth_geoid);
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
    f64_eq!(kep.period(), 7_378.877_993_957_964, "period");
    f64_eq!(kep.hx(), -10_197.722_829_337_885, "HX");
    f64_eq!(kep.hy(), -7_577.364_166_057_776, "HY");
    f64_eq!(kep.hz(), 55_695.034_928_191_49, "HZ");
    f64_eq!(kep.tlong(), 0.691_699_999_999_855_2, "tlong");
    f64_eq!(kep.ea(), 98.501_748_370_880_22, "ea");
    f64_eq!(kep.ma(), 97.113_427_049_323_43, "ma");
    f64_eq!(kep.apoapsis(), 8_392.632_285_000_007, "apo");
    f64_eq!(kep.periapsis(), 7_991.227_715_000_001, "peri");
    f64_eq!(kep.semi_parameter(), 8_187.012_794_017_503, "semi parameter");
}

#[test]
fn state_def_circ_eq() {
    use nyx::celestia::{Cosm, Geoid, State};
    let cosm = Cosm::from_xb("./de438s");
    let mut earth_geoid = cosm.geoid_from_id(3).unwrap();
    // Let's use the GMAT GM value for which these tests we written.
    earth_geoid.gm = 398_600.441_5;
    let dt = Epoch::from_mjd_tai(21_545.0);
    let cart = State::<Geoid>::from_cartesian(
        -38_892.724_449_149_02,
        16_830.384_772_891_86,
        0.722_659_929_135_562_2,
        -1.218_008_333_846_6,
        -2.814_651_172_605_98,
        1.140_294_223_185_661e-5,
        dt,
        earth_geoid,
    );
    f64_eq!(cart.energy(), -4.702_902_670_552_006, "energy");
    f64_eq!(cart.period(), 86_820.776_152_986_1, "period");
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
    f64_eq!(cart.semi_parameter(), 42_378.129_999_999_976, "semi parameter");

    let kep = State::<Geoid>::from_keplerian(18191.098, 1e-6, 1e-6, 306.543, 314.32, 98.765, dt, earth_geoid);
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
    f64_eq!(kep.period(), 24_417.396_242_570_256, "period");
    f64_eq!(kep.hx(), -0.001_194_024_028_558_358_7, "HX");
    f64_eq!(kep.hy(), -0.000_884_918_835_027_750_6, "HY");
    f64_eq!(kep.hz(), 85_152.684_597_507_06, "HZ");
    f64_eq!(kep.tlong(), 359.627_999_999_999_93, "tlong");
    f64_eq!(kep.ea(), 98.764_943_347_932_57, "ea");
    f64_eq!(kep.ma(), 98.764_886_721_264_56, "ma");
    f64_eq!(kep.apoapsis(), 18_191.116_191_098_008, "apo");
    f64_eq!(kep.periapsis(), 18_191.079_808_902_017, "peri");
    f64_eq!(kep.semi_parameter(), 18_191.097_999_981_823, "semi parameter");
}

#[test]
fn state_def_reciprocity() {
    use nyx::celestia::{Cosm, Geoid, State};
    let cosm = Cosm::from_xb("./de438s");
    let mut earth_geoid = cosm.geoid_from_id(3).unwrap();
    // Let's use the GMAT GM value for which these tests we written.
    earth_geoid.gm = 398_600.441_5;
    let dt = Epoch::from_mjd_tai(21_545.0);

    assert_eq!(
        State::<Geoid>::from_cartesian(
            -38_892.724_449_149_02,
            16_830.384_772_891_86,
            0.722_659_929_135_562_2,
            -1.218_008_333_846_6,
            -2.814_651_172_605_98,
            1.140_294_223_185_661e-5,
            dt,
            earth_geoid
        ),
        State::<Geoid>::from_keplerian(
            42_378.129_999_999_98,
            9.999_999_809_555_511e-9,
            0.001_000_000_401_564_538_6,
            78.9,
            65.399_999_847_186_78,
            12.300_000_152_813_197,
            dt,
            earth_geoid
        ),
        "circ_eq"
    );

    assert_eq!(
        State::<Geoid>::from_cartesian(
            5_946.673_548_288_958,
            1_656.154_606_023_661,
            2_259.012_129_598_249,
            -3.098_683_050_943_824,
            4.579_534_132_135_011,
            6.246_541_551_539_432,
            dt,
            earth_geoid
        ),
        State::<Geoid>::from_keplerian(
            7_712.186_117_895_041,
            0.158_999_999_999_999_95,
            53.75369,
            1.998_632_864_211_17e-5,
            359.787_880_000_004,
            25.434_003_407_751_188,
            dt,
            earth_geoid
        ),
        "elliptical"
    );

    assert_eq!(
        State::<Geoid>::from_cartesian(-2436.45, -2436.45, 6891.037, 5.088_611, -5.088_611, 0.0, dt, earth_geoid),
        State::<Geoid>::from_keplerian(
            7_712.186_117_895_043,
            0.000_999_582_831_432_052_5,
            63.434_003_407_751_14,
            135.0,
            90.0,
            0.0,
            dt,
            earth_geoid
        ),
        "circ_inc"
    );
}

#[test]
fn geodetic_vallado() {
    use nyx::celestia::{Cosm, Geoid, State};
    let cosm = Cosm::from_xb("./de438s");
    let mut earth_geoid = cosm.geoid_from_id(3).unwrap();
    // Let's use the GMAT GM value for which these tests we written.
    earth_geoid.gm = 398_600.441_5;
    let dt = Epoch::from_mjd_tai(51_545.0);
    // Test case from Vallado, 4th Ed., page 173, Example 3-3
    let ri = 6524.834;
    let ri_val = 6_524.833_999_999_999;
    let rj = 6862.875;
    let rj_val = 6_862.874_999_999_999;
    let rk = 6448.296;
    let lat = 34.352_495_150_861_564;
    let long = 46.446_416_856_789_96;
    let height = 5_085.218_731_091_624;
    let r = State::<Geoid>::from_position(ri, rj, rk, dt, earth_geoid);
    f64_eq!(r.geodetic_latitude(), lat, "latitude (φ)");
    f64_eq!(r.geodetic_longitude(), long, "longitude (λ)");
    f64_eq!(r.geodetic_height(), height, "height");
    let r = State::<Geoid>::from_geodesic(lat, long, height, dt, earth_geoid);
    f64_eq!(r.x, ri_val, "r_i");
    f64_eq!(r.y, rj_val, "r_j");
    f64_eq!(r.z, rk, "r_k");

    // Test case from Vallado, 4th Ed., page 173, Example 3-4
    let lat = -7.906_635_7;
    let lat_val = -7.906_635_699_999_994_5;
    let long = 345.5975;
    let height = 56.0e-3;
    let height_val = 0.056_000_000_000_494_765;
    let ri = 6_119.400_259_009_384;
    let rj = -1_571.479_552_801_429;
    let rk = -871.561_257_578_933_4;
    let r = State::<Geoid>::from_geodesic(lat, long, height, dt, earth_geoid);
    f64_eq!(r.x, ri, "r_i");
    f64_eq!(r.y, rj, "r_j");
    f64_eq!(r.z, rk, "r_k");
    let r = State::<Geoid>::from_position(ri, rj, rk, dt, earth_geoid);
    f64_eq!(r.geodetic_latitude(), lat_val, "latitude (φ)");
    f64_eq!(r.geodetic_longitude(), long, "longitude (λ)");
    f64_eq!(r.geodetic_height(), height_val, "height");
}
