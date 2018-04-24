extern crate nalgebra as na;
extern crate nyx_space as nyx;

#[test]
fn state_def_circ_inc() {
    use nyx::celestia::{State, EARTH};
    let leo =
        State::from_cartesian::<EARTH>(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0);
    assert_eq!(leo.x, -2436.45, "x");
    assert_eq!(leo.y, -2436.45, "y");
    assert_eq!(leo.z, 6891.037, "z");
    assert_eq!(leo.vx, 5.088611, "vx");
    assert_eq!(leo.vy, -5.088611, "vy");
    assert_eq!(leo.vz, 0.0, "vz");
    assert_eq!(leo.energy(), -25.842247282849137, "energy");
    assert_eq!(leo.period(), 6740.269063643045, "period");
    assert_eq!(leo.hx(), 35065.806679607005, "HX");
    assert_eq!(leo.hy(), 35065.806679607005, "HY");
    assert_eq!(leo.hz(), 24796.2925419, "HZ");
    assert_eq!(leo.sma(), 7712.186117895043, "sma");
    assert_eq!(leo.inc(), 63.43400340775114, "inc");
    assert_eq!(leo.ecc(), 0.0009995828314320525, "ecc");
    assert_eq!(leo.aop(), 90.0, "aop");
    assert_eq!(leo.raan(), 135.0, "raan");
    assert_eq!(leo.ta(), 0.0, "ta");
    assert_eq!(leo.tlong(), 225.0, "tlong");
}

#[test]
fn state_def_elliptical() {
    use nyx::celestia::{State, EARTH};
    let leo = State::from_cartesian::<EARTH>(
        5946.673548288958,
        1656.154606023661,
        2259.012129598249,
        -3.098683050943824,
        4.579534132135011,
        6.246541551539432,
    );
    assert_eq!(leo.energy(), -25.842247282849144, "energy");
    assert_eq!(leo.period(), 6740.2690636430425, "period");
    assert_eq!(leo.hx(), 0.015409898034704383, "HX");
    assert_eq!(leo.hy(), -44146.10601069001, "HY");
    assert_eq!(leo.hz(), 32364.892694481765, "HZ");
    assert_eq!(leo.sma(), 7712.186117895041, "sma");
    assert_eq!(leo.inc(), 53.75369, "inc");
    assert_eq!(leo.ecc(), 0.15899999999999995, "ecc");
    assert_eq!(leo.aop(), 359.787880000004, "aop");
    assert_eq!(leo.raan(), 1.99863286421117e-05, "raan");
    assert_eq!(leo.ta(), 25.434003407751188, "ta");
    assert_eq!(leo.tlong(), 25.221903394083824, "tlong");
}

#[test]
fn state_def_circ_eq() {
    use nyx::celestia::{State, EARTH};
    let geo = State::from_cartesian::<EARTH>(
        -38892.72444914902,
        16830.38477289186,
        0.7226599291355622,
        -1.2180083338466,
        -2.81465117260598,
        1.140294223185661e-05,
    );
    assert_eq!(geo.energy(), -4.702902670552006, "energy");
    assert_eq!(geo.period(), 86820.7761529861, "period");
    assert_eq!(geo.hx(), 2.2259515222419695, "HX");
    assert_eq!(geo.hy(), -0.4367143260909446, "HY");
    assert_eq!(geo.hz(), 129969.00139186575, "HZ");
    assert_eq!(geo.sma(), 42378.12999999998, "sma");
    assert_eq!(geo.inc(), 0.0010000004015645386, "inc");
    assert_eq!(geo.ecc(), 9.999999809555511e-09, "ecc");
    assert_eq!(geo.aop(), 65.39999984718678, "aop");
    assert_eq!(geo.raan(), 78.90000000000001, "raan");
    assert_eq!(geo.ta(), 12.300000152813197, "ta");
    assert_eq!(geo.tlong(), 156.59999999999997, "tlong");
}
