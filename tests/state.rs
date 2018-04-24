extern crate nalgebra as na;
extern crate nyx_space as nyx;

extern crate pretty_env_logger;
//
// #[test]
// fn state_def_geo() {
//     use nyx::celestia::{State, EARTH};
//     let geo =
//         State::from_cartesian::<EARTH>(6524.834, 6862.875, 6448.296, 4.901327, 5.533756, -1.976341);
//     // oT := NewOrbitFromOE(36127.343, 0.832853, 87.869126, 227.898260, 53.384931, 92.335157, Earth)
//     // assert_eq!(geo.aop(), 53.384931, "aop");
// }

#[test]
fn state_def_leo() {
    pretty_env_logger::init();
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
