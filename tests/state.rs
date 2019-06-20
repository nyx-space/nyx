extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use hifitime::julian::ModifiedJulian;
use std::f64::EPSILON;

macro_rules! f64_eq {
    ($x:expr, $val:expr, $msg:expr) => {
        assert!(($x - $val).abs() < EPSILON, $msg)
    };
}

#[test]
fn state_def_() {
    pel::init();
}

#[test]
fn state_def_circ_inc() {
    use hifitime::datetime::Datetime;
    use hifitime::TimeSystem;
    use nyx::celestia::State;
    let dt = ModifiedJulian { days: 21545.0 };
    let cart = State::from_cartesian_eci(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0, dt);
    let cart2 = State::from_cartesian_eci(
        -2436.45,
        -2436.45,
        6891.037,
        5.088611,
        -5.088611,
        0.0,
        Datetime::from_instant(dt.into_instant()),
    );
    assert_eq!(
        cart, cart2,
        "different representations of the datetime are not considered equal"
    );
    f64_eq!(cart.x, -2436.45, "x");
    f64_eq!(cart.y, -2436.45, "y");
    f64_eq!(cart.z, 6891.037, "z");
    f64_eq!(cart.vx, 5.088611, "vx");
    f64_eq!(cart.vy, -5.088611, "vy");
    f64_eq!(cart.vz, 0.0, "vz");
    f64_eq!(cart.energy(), -25.842247282849137, "energy");
    f64_eq!(cart.period(), 6740.269063643045, "period");
    f64_eq!(cart.hx(), 35065.806679607005, "HX");
    f64_eq!(cart.hy(), 35065.806679607005, "HY");
    f64_eq!(cart.hz(), 24796.2925419, "HZ");
    f64_eq!(cart.sma(), 7712.186117895043, "sma");
    f64_eq!(cart.ecc(), 0.0009995828314320525, "ecc");
    f64_eq!(cart.inc(), 63.43400340775114, "inc");
    f64_eq!(cart.raan(), 135.0, "raan");
    f64_eq!(cart.aop(), 90.0, "aop");
    f64_eq!(cart.ta(), 0.0, "ta");
    f64_eq!(cart.tlong(), 225.0, "tlong");
    f64_eq!(cart.ea(), 0.0, "ea");
    f64_eq!(cart.ma(), 0.0, "ma");
    f64_eq!(cart.apoapsis(), 7719.895086731299, "apo");
    f64_eq!(cart.periapsis(), 7704.477149058786, "peri");
    f64_eq!(cart.semi_parameter(), 7712.178412142147, "semi parameter");

    let kep = State::from_keplerian_eci(8191.93, 1e-6, 12.85, 306.614, 314.19, 99.8877, dt);
    f64_eq!(kep.x, 8057.976452202976, "x");
    f64_eq!(kep.y, -0.1967403702908889, "y");
    f64_eq!(kep.z, 1475.383214274138, "z");
    f64_eq!(kep.vx, -0.16647048858407631, "vx");
    f64_eq!(kep.vy, 6.9138686382756465, "vy");
    f64_eq!(kep.vz, 0.9101579814432791, "vz");
    f64_eq!(kep.sma(), 8191.929999999999, "sma");
    f64_eq!(kep.ecc(), 1.00000000038851e-06, "ecc");
    f64_eq!(kep.inc(), 12.849999999999987, "inc");
    f64_eq!(kep.raan(), 306.614, "raan");
    f64_eq!(kep.aop(), 314.1899999946181, "aop");
    f64_eq!(kep.ta(), 99.8877000053819, "ta");
    f64_eq!(kep.energy(), -24.32884811637795, "energy");
    f64_eq!(kep.period(), 7378.877993957958, "period");
    f64_eq!(kep.hx(), -10200.784799426574, "HX");
    f64_eq!(kep.hy(), -7579.639346783497, "HY");
    f64_eq!(kep.hz(), 55711.75792938425, "HZ");
    f64_eq!(kep.tlong(), 0.6917000000000826, "tlong");
    f64_eq!(kep.ea(), 99.88764356065685, "ea");
    f64_eq!(kep.ma(), 99.88758711592696, "ma");
    f64_eq!(kep.apoapsis(), 8191.938191930002, "apo");
    f64_eq!(kep.periapsis(), 8191.921808069997, "peri");
    f64_eq!(kep.semi_parameter(), 8191.929999991808, "semi parameter");
}

#[test]
fn state_def_elliptical() {
    use nyx::celestia::State;
    let dt = ModifiedJulian { days: 21545.0 };
    let cart = State::from_cartesian_eci(
        5946.673548288958,
        1656.154606023661,
        2259.012129598249,
        -3.098683050943824,
        4.579534132135011,
        6.246541551539432,
        dt,
    );
    f64_eq!(cart.energy(), -25.842247282849144, "energy");
    f64_eq!(cart.period(), 6740.2690636430425, "period");
    f64_eq!(cart.hx(), 0.015409898034704383, "HX");
    f64_eq!(cart.hy(), -44146.10601069001, "HY");
    f64_eq!(cart.hz(), 32364.892694481765, "HZ");
    f64_eq!(cart.sma(), 7712.186117895041, "sma");
    f64_eq!(cart.ecc(), 0.15899999999999995, "ecc");
    f64_eq!(cart.inc(), 53.75369, "inc");
    f64_eq!(cart.raan(), 1.99863286421117e-05, "raan");
    f64_eq!(cart.aop(), 359.787880000004, "aop");
    f64_eq!(cart.ta(), 25.434003407751188, "ta");
    f64_eq!(cart.tlong(), 25.221903394083824, "tlong");
    f64_eq!(cart.ea(), 21.76305288258479, "ea");
    f64_eq!(cart.ma(), 18.38533633051639, "ma");
    f64_eq!(cart.apoapsis(), 8938.423710640353, "apo");
    f64_eq!(cart.periapsis(), 6485.94852514973, "peri");
    f64_eq!(cart.semi_parameter(), 7517.214340648537, "semi parameter");

    let kep = State::from_keplerian_eci(8191.93, 0.0245, 12.85, 306.614, 314.19, 99.8877, dt);
    f64_eq!(kep.x, 8087.1616180485225, "x");
    f64_eq!(kep.y, -0.19745294377252073, "y");
    f64_eq!(kep.z, 1480.726901246883, "z");
    f64_eq!(kep.vx, -0.00016859218684395216, "vx");
    f64_eq!(kep.vy, 6.886845792370852, "vy");
    f64_eq!(kep.vz, 0.9369312603028918, "vz");
    f64_eq!(kep.sma(), 8191.930000000003, "sma");
    f64_eq!(kep.ecc(), 0.024500000000000348, "ecc");
    f64_eq!(kep.inc(), 12.850000000000016, "inc");
    f64_eq!(kep.raan(), 306.614, "raan");
    f64_eq!(kep.aop(), 314.1900000000004, "aop");
    f64_eq!(kep.ta(), 99.88769999999958, "ta");
    f64_eq!(kep.energy(), -24.32884811637794, "energy");
    f64_eq!(kep.period(), 7378.877993957964, "period");
    f64_eq!(kep.hx(), -10197.722829337885, "HX");
    f64_eq!(kep.hy(), -7577.364166057776, "HY");
    f64_eq!(kep.hz(), 55695.03492819149, "HZ");
    f64_eq!(kep.tlong(), 0.6916999999998552, "tlong");
    f64_eq!(kep.ea(), 98.50174837088022, "ea");
    f64_eq!(kep.ma(), 97.11342704932343, "ma");
    f64_eq!(kep.apoapsis(), 8392.632285000007, "apo");
    f64_eq!(kep.periapsis(), 7991.227715000001, "peri");
    f64_eq!(kep.semi_parameter(), 8187.012794017503, "semi parameter");
}

#[test]
fn state_def_circ_eq() {
    use nyx::celestia::State;
    let dt = ModifiedJulian { days: 21545.0 };
    let cart = State::from_cartesian_eci(
        -38892.72444914902,
        16830.38477289186,
        0.7226599291355622,
        -1.2180083338466,
        -2.81465117260598,
        1.140294223185661e-05,
        dt,
    );
    f64_eq!(cart.energy(), -4.702902670552006, "energy");
    f64_eq!(cart.period(), 86820.7761529861, "period");
    f64_eq!(cart.hx(), 2.2259515222419695, "HX");
    f64_eq!(cart.hy(), -0.4367143260909446, "HY");
    f64_eq!(cart.hz(), 129969.00139186575, "HZ");
    f64_eq!(cart.sma(), 42378.12999999998, "sma");
    f64_eq!(cart.ecc(), 9.999999809555511e-09, "ecc");
    f64_eq!(cart.inc(), 0.0010000004015645386, "inc");
    f64_eq!(cart.raan(), 78.90000000000001, "raan");
    f64_eq!(cart.aop(), 65.39999984718678, "aop");
    f64_eq!(cart.ta(), 12.300000152813197, "ta");
    f64_eq!(cart.tlong(), 156.59999999999997, "tlong");
    f64_eq!(cart.ea(), 12.300000030755777, "ea");
    f64_eq!(cart.ma(), 12.299999908698359, "ma");
    f64_eq!(cart.apoapsis(), 42378.13042378127, "apo");
    f64_eq!(cart.periapsis(), 42378.12957621869, "peri");
    f64_eq!(cart.semi_parameter(), 42378.129999999976, "semi parameter");

    let kep = State::from_keplerian_eci(18191.098, 1e-6, 1e-6, 306.543, 314.32, 98.765, dt);
    f64_eq!(kep.x, 18190.717357886369, "x");
    f64_eq!(kep.y, -118.10716253921869, "y");
    f64_eq!(kep.z, 0.00025384564763305335, "z");
    f64_eq!(kep.vx, 0.03039644013026488, "vx");
    f64_eq!(kep.vy, 4.680909107924576, "vy");
    f64_eq!(kep.vz, 4.907089816726583e-8, "vz");
    f64_eq!(kep.sma(), 18191.098000000013, "sma");
    f64_eq!(kep.ecc(), 9.999999997416087e-7, "ecc");
    f64_eq!(kep.inc(), 1.2074182697257333e-06, "inc");
    f64_eq!(kep.raan(), 306.543, "raan");
    f64_eq!(kep.aop(), 314.32000002540366, "aop");
    f64_eq!(kep.ta(), 98.76499997459628, "ta");
    f64_eq!(kep.energy(), -10.955920349063035, "energy");
    f64_eq!(kep.period(), 24417.396242570256, "period");
    f64_eq!(kep.hx(), -0.0011940240285583587, "HX");
    f64_eq!(kep.hy(), -0.0008849188350277506, "HY");
    f64_eq!(kep.hz(), 85152.68459750706, "HZ");
    f64_eq!(kep.tlong(), 359.62799999999993, "tlong");
    f64_eq!(kep.ea(), 98.76494334793257, "ea");
    f64_eq!(kep.ma(), 98.76488672126456, "ma");
    f64_eq!(kep.apoapsis(), 18191.116191098008, "apo");
    f64_eq!(kep.periapsis(), 18191.079808902017, "peri");
    f64_eq!(kep.semi_parameter(), 18191.097999981823, "semi parameter");
}

#[test]
fn state_def_reciprocity() {
    use nyx::celestia::State;
    let dt = ModifiedJulian { days: 21545.0 };

    assert_eq!(
        State::from_cartesian_eci(
            -38892.72444914902,
            16830.38477289186,
            0.7226599291355622,
            -1.2180083338466,
            -2.81465117260598,
            1.140294223185661e-05,
            dt,
        ),
        State::from_keplerian_eci(
            42378.12999999998,
            9.999999809555511e-09,
            0.0010000004015645386,
            78.90000000000001,
            65.39999984718678,
            12.300000152813197,
            dt,
        ),
        "circ_eq"
    );

    assert_eq!(
        State::from_cartesian_eci(
            5946.673548288958,
            1656.154606023661,
            2259.012129598249,
            -3.098683050943824,
            4.579534132135011,
            6.246541551539432,
            dt,
        ),
        State::from_keplerian_eci(
            7712.186117895041,
            0.15899999999999995,
            53.75369,
            1.99863286421117e-05,
            359.787880000004,
            25.434003407751188,
            dt,
        ),
        "elliptical"
    );

    assert_eq!(
        State::from_cartesian_eci(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0, dt),
        State::from_keplerian_eci(
            7712.186117895043,
            0.0009995828314320525,
            63.43400340775114,
            135.0,
            90.0,
            0.0,
            dt,
        ),
        "circ_inc"
    );
}

#[test]
fn geodetic_vallado() {
    use nyx::celestia::{Cosm, State};
    let cosm = Cosm::from_xb("./de438s");
    let dt = ModifiedJulian { days: 51545.0 };
    // Test case from Vallado, 4th Ed., page 173, Example 3-3
    let ri = 6524.834;
    let ri_val = 6524.833999999999;
    let rj = 6862.875;
    let rj_val = 6862.874999999999;
    let rk = 6448.296;
    let rk_val = 6448.295999999998;
    let lat = 34.352495150861564;
    let long = 46.44641685678996;
    let height = 5085.218731091624;
    let r = State::<ECEF>::from_position(ri, rj, rk, dt);
    f64_eq!(r.geodetic_latitude(), lat, "latitude (φ)");
    f64_eq!(r.geodetic_longitude(), long, "longitude (λ)");
    f64_eq!(r.geodetic_height(), height, "height");
    let r = State::<ECEF>::from_geodesic(lat, long, height, dt);
    f64_eq!(r.ri(), ri_val, "r_i");
    f64_eq!(r.rj(), rj_val, "r_j");
    f64_eq!(r.rk(), rk_val, "r_k");

    // Test case from Vallado, 4th Ed., page 173, Example 3-4
    let lat = -7.906_635_7;
    let lat_val = -7.906_635_699_999_994_5;
    let long = 345.5975;
    let height = 56.0e-3;
    let height_val = 0.056000000000494765;
    let ri = 6119.400259009384;
    let rj = -1571.4795528014297;
    let rk = -871.5612575789334;
    let r = State::<ECEF>::from_geodesic(lat, long, height, dt);
    f64_eq!(r.ri(), ri, "r_i");
    f64_eq!(r.rj(), rj, "r_j");
    f64_eq!(r.rk(), rk, "r_k");
    let r = State::<ECEF>::from_position(ri, rj, rk, dt);
    f64_eq!(r.geodetic_latitude(), lat_val, "latitude (φ)");
    f64_eq!(r.geodetic_longitude(), long, "longitude (λ)");
    f64_eq!(r.geodetic_height(), height_val, "height");
}
