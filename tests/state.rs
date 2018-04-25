extern crate nalgebra as na;
extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

#[test]
fn state_def_() {
    pel::init();
}

#[test]
fn state_def_circ_inc() {
    use nyx::celestia::{State, EARTH};
    let cart =
        State::from_cartesian::<EARTH>(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0);
    assert_eq!(cart.x, -2436.45, "x");
    assert_eq!(cart.y, -2436.45, "y");
    assert_eq!(cart.z, 6891.037, "z");
    assert_eq!(cart.vx, 5.088611, "vx");
    assert_eq!(cart.vy, -5.088611, "vy");
    assert_eq!(cart.vz, 0.0, "vz");
    assert_eq!(cart.energy(), -25.842247282849137, "energy");
    assert_eq!(cart.period(), 6740.269063643045, "period");
    assert_eq!(cart.hx(), 35065.806679607005, "HX");
    assert_eq!(cart.hy(), 35065.806679607005, "HY");
    assert_eq!(cart.hz(), 24796.2925419, "HZ");
    assert_eq!(cart.sma(), 7712.186117895043, "sma");
    assert_eq!(cart.ecc(), 0.0009995828314320525, "ecc");
    assert_eq!(cart.inc(), 63.43400340775114, "inc");
    assert_eq!(cart.raan(), 135.0, "raan");
    assert_eq!(cart.aop(), 90.0, "aop");
    assert_eq!(cart.ta(), 0.0, "ta");
    assert_eq!(cart.tlong(), 225.0, "tlong");
    assert_eq!(cart.ea(), 0.0, "ea");
    assert_eq!(cart.ma(), 0.0, "ma");
    assert_eq!(cart.apoapsis(), 7719.895086731299, "apo");
    assert_eq!(cart.periapsis(), 7704.477149058786, "peri");
    assert_eq!(cart.semi_parameter(), 7712.178412142147, "semi parameter");

    let kep = State::from_keplerian::<EARTH>(8191.93, 1e-6, 12.85, 306.614, 314.19, 99.8877);
    assert_eq!(kep.x, 8057.976452202976, "x");
    assert_eq!(kep.y, -0.1967403702908889, "y");
    assert_eq!(kep.z, 1475.383214274138, "z");
    assert_eq!(kep.vx, -0.16647048858407631, "vx");
    assert_eq!(kep.vy, 6.9138686382756465, "vy");
    assert_eq!(kep.vz, 0.9101579814432791, "vz");
    assert_eq!(kep.sma(), 8191.929999999999, "sma");
    assert_eq!(kep.ecc(), 1.00000000038851e-06, "ecc");
    assert_eq!(kep.inc(), 12.849999999999987, "inc");
    assert_eq!(kep.raan(), 306.614, "raan");
    assert_eq!(kep.aop(), 314.1899999946181, "aop");
    assert_eq!(kep.ta(), 99.8877000053819, "ta");
    assert_eq!(kep.energy(), -24.32884811637795, "energy");
    assert_eq!(kep.period(), 7378.877993957958, "period");
    assert_eq!(kep.hx(), -10200.784799426574, "HX");
    assert_eq!(kep.hy(), -7579.639346783497, "HY");
    assert_eq!(kep.hz(), 55711.75792938425, "HZ");
    assert_eq!(kep.tlong(), 0.6917000000000826, "tlong");
    assert_eq!(kep.ea(), 99.88764356065685, "ea");
    assert_eq!(kep.ma(), 99.88758711592696, "ma");
    assert_eq!(kep.apoapsis(), 8191.938191930002, "apo");
    assert_eq!(kep.periapsis(), 8191.921808069997, "peri");
    assert_eq!(kep.semi_parameter(), 8191.929999991808, "semi parameter");
}

#[test]
fn state_def_elliptical() {
    use nyx::celestia::{State, EARTH};
    let cart = State::from_cartesian::<EARTH>(
        5946.673548288958,
        1656.154606023661,
        2259.012129598249,
        -3.098683050943824,
        4.579534132135011,
        6.246541551539432,
    );
    assert_eq!(cart.energy(), -25.842247282849144, "energy");
    assert_eq!(cart.period(), 6740.2690636430425, "period");
    assert_eq!(cart.hx(), 0.015409898034704383, "HX");
    assert_eq!(cart.hy(), -44146.10601069001, "HY");
    assert_eq!(cart.hz(), 32364.892694481765, "HZ");
    assert_eq!(cart.sma(), 7712.186117895041, "sma");
    assert_eq!(cart.ecc(), 0.15899999999999995, "ecc");
    assert_eq!(cart.inc(), 53.75369, "inc");
    assert_eq!(cart.raan(), 1.99863286421117e-05, "raan");
    assert_eq!(cart.aop(), 359.787880000004, "aop");
    assert_eq!(cart.ta(), 25.434003407751188, "ta");
    assert_eq!(cart.tlong(), 25.221903394083824, "tlong");
    assert_eq!(cart.ea(), 21.76305288258479, "ea");
    assert_eq!(cart.ma(), 18.38533633051639, "ma");
    assert_eq!(cart.apoapsis(), 8938.423710640353, "apo");
    assert_eq!(cart.periapsis(), 6485.94852514973, "peri");
    assert_eq!(cart.semi_parameter(), 7517.214340648537, "semi parameter");

    let kep = State::from_keplerian::<EARTH>(8191.93, 0.0245, 12.85, 306.614, 314.19, 99.8877);
    assert_eq!(kep.x, 8087.1616180485225, "x");
    assert_eq!(kep.y, -0.19745294377252073, "y");
    assert_eq!(kep.z, 1480.726901246883, "z");
    assert_eq!(kep.vx, -0.00016859218684395216, "vx");
    assert_eq!(kep.vy, 6.886845792370852, "vy");
    assert_eq!(kep.vz, 0.9369312603028918, "vz");
    assert_eq!(kep.sma(), 8191.930000000003, "sma");
    assert_eq!(kep.ecc(), 0.024500000000000348, "ecc");
    assert_eq!(kep.inc(), 12.850000000000016, "inc");
    assert_eq!(kep.raan(), 306.614, "raan");
    assert_eq!(kep.aop(), 314.1900000000004, "aop");
    assert_eq!(kep.ta(), 99.88769999999958, "ta");
    assert_eq!(kep.energy(), -24.32884811637794, "energy");
    assert_eq!(kep.period(), 7378.877993957964, "period");
    assert_eq!(kep.hx(), -10197.722829337885, "HX");
    assert_eq!(kep.hy(), -7577.364166057776, "HY");
    assert_eq!(kep.hz(), 55695.03492819149, "HZ");
    assert_eq!(kep.tlong(), 0.6916999999998552, "tlong");
    assert_eq!(kep.ea(), 98.50174837088022, "ea");
    assert_eq!(kep.ma(), 97.11342704932343, "ma");
    assert_eq!(kep.apoapsis(), 8392.632285000007, "apo");
    assert_eq!(kep.periapsis(), 7991.227715000001, "peri");
    assert_eq!(kep.semi_parameter(), 8187.012794017503, "semi parameter");
}

#[test]
fn state_def_circ_eq() {
    use nyx::celestia::{State, EARTH};
    let cart = State::from_cartesian::<EARTH>(
        -38892.72444914902,
        16830.38477289186,
        0.7226599291355622,
        -1.2180083338466,
        -2.81465117260598,
        1.140294223185661e-05,
    );
    assert_eq!(cart.energy(), -4.702902670552006, "energy");
    assert_eq!(cart.period(), 86820.7761529861, "period");
    assert_eq!(cart.hx(), 2.2259515222419695, "HX");
    assert_eq!(cart.hy(), -0.4367143260909446, "HY");
    assert_eq!(cart.hz(), 129969.00139186575, "HZ");
    assert_eq!(cart.sma(), 42378.12999999998, "sma");
    assert_eq!(cart.ecc(), 9.999999809555511e-09, "ecc");
    assert_eq!(cart.inc(), 0.0010000004015645386, "inc");
    assert_eq!(cart.raan(), 78.90000000000001, "raan");
    assert_eq!(cart.aop(), 65.39999984718678, "aop");
    assert_eq!(cart.ta(), 12.300000152813197, "ta");
    assert_eq!(cart.tlong(), 156.59999999999997, "tlong");
    assert_eq!(cart.ea(), 12.300000030755777, "ea");
    assert_eq!(cart.ma(), 12.299999908698359, "ma");
    assert_eq!(cart.apoapsis(), 42378.13042378127, "apo");
    assert_eq!(cart.periapsis(), 42378.12957621869, "peri");
    assert_eq!(cart.semi_parameter(), 42378.129999999976, "semi parameter");

    let kep = State::from_keplerian::<EARTH>(18191.098, 1e-6, 1e-6, 306.543, 314.32, 98.765);
    assert_eq!(kep.x, 18190.717357886369, "x");
    assert_eq!(kep.y, -118.10716253921869, "y");
    assert_eq!(kep.z, 0.00025384564763305335, "z");
    assert_eq!(kep.vx, 0.03039644013026488, "vx");
    assert_eq!(kep.vy, 4.680909107924576, "vy");
    assert_eq!(kep.vz, 4.907089816726583e-8, "vz");
    assert_eq!(kep.sma(), 18191.098000000013, "sma");
    assert_eq!(kep.ecc(), 9.999999997416087e-7, "ecc");
    assert_eq!(kep.inc(), 1.2074182697257333e-06, "inc");
    assert_eq!(kep.raan(), 306.543, "raan");
    assert_eq!(kep.aop(), 314.32000002540366, "aop");
    assert_eq!(kep.ta(), 98.76499997459628, "ta");
    assert_eq!(kep.energy(), -10.955920349063035, "energy");
    assert_eq!(kep.period(), 24417.396242570256, "period");
    assert_eq!(kep.hx(), -0.0011940240285583587, "HX");
    assert_eq!(kep.hy(), -0.0008849188350277506, "HY");
    assert_eq!(kep.hz(), 85152.68459750706, "HZ");
    assert_eq!(kep.tlong(), 359.62799999999993, "tlong");
    assert_eq!(kep.ea(), 98.76494334793257, "ea");
    assert_eq!(kep.ma(), 98.76488672126456, "ma");
    assert_eq!(kep.apoapsis(), 18191.116191098008, "apo");
    assert_eq!(kep.periapsis(), 18191.079808902017, "peri");
    assert_eq!(kep.semi_parameter(), 18191.097999981823, "semi parameter");
}

#[test]
fn state_def_reciprocity() {
    use nyx::celestia::{State, EARTH};

    assert_eq!(
        State::from_cartesian::<EARTH>(
            -38892.72444914902,
            16830.38477289186,
            0.7226599291355622,
            -1.2180083338466,
            -2.81465117260598,
            1.140294223185661e-05,
        ),
        State::from_keplerian::<EARTH>(
            42378.12999999998,
            9.999999809555511e-09,
            0.0010000004015645386,
            78.90000000000001,
            65.39999984718678,
            12.300000152813197,
        ),
        "circ_eq"
    );

    assert_eq!(
        State::from_cartesian::<EARTH>(
            5946.673548288958,
            1656.154606023661,
            2259.012129598249,
            -3.098683050943824,
            4.579534132135011,
            6.246541551539432,
        ),
        State::from_keplerian::<EARTH>(
            7712.186117895041,
            0.15899999999999995,
            53.75369,
            1.99863286421117e-05,
            359.787880000004,
            25.434003407751188,
        ),
        "elliptical"
    );

    assert_eq!(
        State::from_cartesian::<EARTH>(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0),
        State::from_keplerian::<EARTH>(
            7712.186117895043,
            0.0009995828314320525,
            63.43400340775114,
            135.0,
            90.0,
            0.0,
        ),
        "circ_inc"
    );
}
