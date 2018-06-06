extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

#[test]
fn eci() {
    use hifitime::instant::{Era, Instant};
    use na::Matrix3;
    use nyx::celestia::{CoordinateFrame, ECI};

    assert_eq!(
        ECI::to_inertial(Instant::from_precise_seconds(1.0, Era::Present)),
        Matrix3::identity()
    );
}

#[test]
fn ecef_2_eci() {
    use hifitime::datetime::*;
    use nyx::celestia::State;
    let dt = Datetime::new(1992, 8, 20, 12, 14, 0, 0).expect("wut?");
    let eci = State::from_cartesian_eci(
        -38892.72444914902,
        16830.38477289186,
        0.7226599291355622,
        -1.2180083338466,
        -2.81465117260598,
        1.140294223185661e-05,
        dt,
    );
    assert_eq!(eci.in_ecef().in_eci(), eci, "reciprocity failed");
}

#[test]
fn ecef_theta_gmst() {
    // Vallado's example
    use hifitime::datetime::*;
    use nyx::celestia::ECEF;
    use std::f64::EPSILON;
    let dt = Datetime::new(1992, 8, 20, 12, 14, 0, 0).expect("wut?");
    assert!(
        (ECEF::gmst(dt.into_instant()) - 152.5787878104796).abs() < EPSILON,
        "wrong θ GMST computed"
    );
}
